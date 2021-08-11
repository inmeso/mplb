/**
 * Copyright 2019 United Kingdom Research and Innovation
 *
 * Authors: See AUTHORS
 *
 * Contact: [jianping.meng@stfc.ac.uk and/or jpmeng@gmail.com]
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice
 *    this list of conditions and the following disclaimer in the documentation
 *    and or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * ANDANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

/**
 * @brief The application is to demostrate how to calculate the Mie solution
 * of the Maxwell equation at a mesh.
 *
 * @author Sina Haeri (s.haeri@ed.ac.uk) and Jianping Meng
 */

#include <iostream>
#include <ostream>
#include "ops_seq_v2.h"
#include "block.h"
#include "field.h"
#include <string>
#include <vector>
//#include <map>
#include "MieConfig.h"
#include "scatterer.hpp"
#include "mie_theory_kernel.inc"

int main(int argc, const char** argv) {
    // OPS initialisation

    int blkNx[3];
    double blkLx[3];
    double minReal = std::numeric_limits<double>::min();
    ops_init(argc, argv, 4);
    double ct0, ct1, et0, et1;
    ops_timers(&ct0, &et0);
    
    // start a simulation by hard-coding
    BOOST_ASSERT_MSG(argc != 1,
        "Provide the name of the input configuration file <input name>.json");
        
    // start a new simulaton from a configuration file
    std::string configFileName(argv[1]);
    ReadConfiguration(configFileName);
     
    // Bessel Function Container
    int maxOrder = 5;
    complexd besselArgz = complexd(1.0, 0.5);

    bessels = new besselFunctions(maxOrder,besselArgz);

    partScatt = new scatterer(maxOrder, Config().partRadius, 
                Config().vaccumWaveLength,
                Config().partPermeability, Config().envPermeability,
                Config().envRefractiveIndex,
                Config().partRefractiveIndex);
                
    // Also define the variable in the GPU memory space.
    ops_decl_const("spaceDim", 1, "int", &(Config().spaceDim));

    // Set a stencil for numerical schemes
    // here we only need the local(current) grid point
    int currentNode[]{0, 0, 0};
    ops_stencil LOCALSTENCIL = ops_decl_stencil(3, 1, currentNode, "000");
   
    // input the mesh size here
    // we assume only one block in this simple case
    Block block(0, "Block", Config().blockSize);
    
    RealField CoordinateXYZ("CoordinateXYZ");
    CoordinateXYZ.SetDataDim(Config().spaceDim);
    CoordinateXYZ.CreateFieldFromScratch(block);

    RealField E("E");
    E.SetDataDim(Config().spaceDim);
    E.CreateFieldFromScratch(block);
    
    RealField H("H");
    H.SetDataDim(Config().spaceDim);
    H.CreateFieldFromScratch(block);
    
    ops_partition("Mie");
    ops_diagnostic_output();
    std::vector<int> iterRng;
    iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
    
    //Particle always at (0, 0, 0)
    blkLx[0] = Config().blockExtent[0];
    blkNx[0] = Config().blockSize[0];
    blkLx[1] = Config().blockExtent[1];
    blkNx[1] = Config().blockSize[1];
    blkLx[2] = Config().blockExtent[2];
    blkNx[2] = Config().blockSize[2];

    // set the CoordinateXYZ
    ops_par_loop(KerGridGen, "KerGridGen", block.Get(), 
                 Config().spaceDim, iterRng.data(),
                 ops_arg_dat(CoordinateXYZ[block.ID()], Config().spaceDim, LOCALSTENCIL,
                             "double", OPS_RW),
                 ops_arg_gbl(&blkNx[0], 1, "int", OPS_READ),
                 ops_arg_gbl(&blkLx[0], 1, "double", OPS_READ),
                 ops_arg_idx());

    // Convert to spherical
    ops_par_loop(KerCart2Sph, "KerCart2Sph", block.Get(), 
                 Config().spaceDim, iterRng.data(),
                 ops_arg_dat(CoordinateXYZ[block.ID()], Config().spaceDim, LOCALSTENCIL,
                             "double", OPS_RW),
                 ops_arg_gbl(&minReal, 1, "double", OPS_READ));

    // calculate the Mie solution
    ops_par_loop(
        KerCalculateMieSolution, "KerCalculateMieSolution", block.Get(),
        Config().spaceDim, iterRng.data(),
        ops_arg_dat(E[block.ID()], Config().spaceDim, LOCALSTENCIL, "double", OPS_RW),
        ops_arg_dat(H[block.ID()], Config().spaceDim, LOCALSTENCIL, "double", OPS_RW),
        ops_arg_dat(CoordinateXYZ[block.ID()], Config().spaceDim, LOCALSTENCIL, "double",
                    OPS_READ));

    ops_timers(&ct1, &et1);
    ops_printf("\nTotal Wall time %lf\n", et1 - et0);
    // Print OPS performance details to output stream
    ops_timing_output(std::cout);
    
    //back to cartesian before write
    ops_par_loop(KerSph2Cart, "KerSph2Cart", block.Get(), 
                 Config().spaceDim, iterRng.data(),
                 ops_arg_dat(CoordinateXYZ[block.ID()], Config().spaceDim, LOCALSTENCIL,
                             "double", OPS_RW));
                             
    CoordinateXYZ.WriteToHDF5(Config().caseName, 0);
    E.WriteToHDF5(Config().caseName, 0);
    H.WriteToHDF5(Config().caseName, 0);

    ops_exit();
    delete bessels;
}