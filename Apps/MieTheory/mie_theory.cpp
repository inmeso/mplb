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
#include <map>
#include "mie_theory_kernel.inc"

int main(int argc, const char** argv) {
    // OPS initialisation
    ops_init(argc, argv, 4);
    double ct0, ct1, et0, et1;
    ops_timers(&ct0, &et0);
    // set the space dimension
    const int spaceDim{3};
    // Also define the variable in the GPU memory space.
    ops_decl_const("spaceDim", 1, "int", &spaceDim);
    const std::string caseName{"Mie"};
    // Set a stencil for numerical schemes
    // here we only need the local(current) grid point
    int currentNode[]{0, 0, 0};
    ops_stencil LOCALSTENCIL = ops_decl_stencil(3, 1, currentNode, "000");
    // input the mesh size here
    // we assume only one block in this simple case
    std::vector<int> blockSize{101, 101, 101};
    Block block(0, "Block", blockSize);
    RealField coordinates("coordinates");
    coordinates.SetDataDim(spaceDim);
    coordinates.CreateFieldFromScratch(block);
    RealField E("E");
    E.SetDataDim(spaceDim);
    E.CreateFieldFromScratch(block);
    RealField H("H");
    H.SetDataDim(spaceDim);
    H.CreateFieldFromScratch(block);
    ops_partition("Mie");
    ops_diagnostic_output();
    std::vector<int> iterRng;
    iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
    // set the coordinates
    ops_par_loop(KerSetSphericalCoord, "KerSetSphericalCoord", block.Get(), spaceDim,
                 iterRng.data(),
                 ops_arg_dat(coordinates[block.ID()], spaceDim, LOCALSTENCIL,
                             "double", OPS_RW),
                 ops_arg_idx());

    // calculate the Mie solution
    ops_par_loop(
        KerCalculateMieSolution, "KerCalculateMieSolution", block.Get(),
        spaceDim, iterRng.data(),
        ops_arg_dat(E[block.ID()], spaceDim, LOCALSTENCIL, "double", OPS_RW),
        ops_arg_dat(H[block.ID()], spaceDim, LOCALSTENCIL, "double", OPS_RW),
        ops_arg_dat(coordinates[block.ID()], spaceDim, LOCALSTENCIL, "double",
                    OPS_READ));

    ops_timers(&ct1, &et1);
    ops_printf("\nTotal Wall time %lf\n", et1 - et0);
    // Print OPS performance details to output stream
    ops_timing_output(std::cout);
    coordinates.WriteToHDF5(caseName, 0);
    E.WriteToHDF5(caseName, 0);
    H.WriteToHDF5(caseName, 0);
    ops_exit();
}