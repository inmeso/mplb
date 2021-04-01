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

/** @brief Basic testing codes for the block and field module
 *  @author Jianping Meng
 **/

#include <iostream>
#include <ostream>
#define FIELD
#include "ops_seq_v2.h"
#ifdef FIELD
#include "block.h"
#include "field.h"
#include <string>
#include <vector>
#include <map>
#endif
#include "setvalue_kernel.inc"

int main(int argc, const char** argv) {
    // OPS initialisation
    ops_init(argc, argv, 3);
    double ct0, ct1, et0, et1;
    ops_timers(&ct0, &et0);
    const int nxyz{32};
#ifdef OPS_3D
    int currentNode[]{0, 0, 0};
    ops_stencil LOCALSTENCIL = ops_decl_stencil(3, 1, currentNode, "000");
#endif  // OPS_3D
#ifdef OPS_2D
    int currentNode[]{0, 0};
    ops_stencil LOCALSTENCIL = ops_decl_stencil(2, 1, currentNode, "00");
#endif  // OPS_2D

    // Test the ops native block and dat interface
    const int haloDepth{1};
#ifdef OPS_3D
    ops_block oblock = ops_decl_block(3, "otest");
    int d_p[]{haloDepth, haloDepth, haloDepth};
    int d_m[]{-haloDepth, -haloDepth, -haloDepth};
    // size of the dat
    int size[]{nxyz, nxyz, nxyz};
    int base[]{0, 0, 0};
#endif  // OPS_3D
#ifdef OPS_2D
    ops_block oblock = ops_decl_block(2, "otest");
    int d_p[]{haloDepth, haloDepth};
    int d_m[]{-haloDepth, -haloDepth};
    // size of the dat
    int size[]{nxyz, nxyz};
    int base[]{0, 0};
#endif  // OPS_2D

    double* temp = NULL;
    int* inttemp = NULL;
    ops_dat omar{
        ops_decl_dat(oblock, 1, size, base, d_m, d_p, inttemp, "int", "omar")};
    ops_dat ogeo{
        ops_decl_dat(oblock, 4, size, base, d_m, d_p, temp, "double", "ogeo")};
#ifdef FIELD
#ifdef OPS_3D
    std::vector<int> blockSize{nxyz, nxyz, nxyz};
#endif
#ifdef OPS_2D
    std::vector<int> blockSize{nxyz, nxyz};
#endif
    Block block(0, "test", blockSize);
    IntField mar("mar");
    RealField geo("geo");
    geo.SetDataDim(4);
    geo.CreateFieldFromScratch(block);
    mar.SetDataDim(1);
    mar.CreateFieldFromScratch(block);
#endif  // FIELD
    ops_partition("FieldBLockTest");
    ops_diagnostic_output();
#ifdef OPS_3D
    int range[]{0, nxyz, 0, nxyz, 0, nxyz};
    ops_par_loop(KerSetValuesSingle, "KerSetValuesSingle", oblock, 3, range,
                 ops_arg_dat(omar, 1, LOCALSTENCIL, "int", OPS_RW));
    ops_printf("Single-dimension test of ops native calls succeeded!\n");
#ifdef FIELD
    std::vector<int> iterRng;
    iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
    ops_par_loop(KerSetValuesSingle, "KerSetValuesSingle", block.Get(), 3,
                 iterRng.data(),
                 ops_arg_dat(mar[block.ID()], 1, LOCALSTENCIL, "int", OPS_RW));
    ops_printf("Single-dimension test of field and block calls succeeded!\n");
#endif  // FIELD
    ops_par_loop(KerSetValues, "KerSetValues", oblock, 3, range,
                 ops_arg_dat(ogeo, 4, LOCALSTENCIL, "double", OPS_RW));
    ops_printf("Multi-dimension test of ops native calls succeeded!\n");
#ifdef FIELD
    ops_par_loop(
        KerSetValues, "KerSetValues", block.Get(), 3, iterRng.data(),
        ops_arg_dat(geo[block.ID()], 4, LOCALSTENCIL, "double", OPS_RW));
    ops_printf("Multi-dimension test of field and block calls succeeded!\n");
#endif  // FIELD
#endif  // OPS_3D
#ifdef OPS_2D
    int range[]{0, nxyz, 0, nxyz};
    ops_par_loop(KerSetValuesSingle, "KerSetValuesSingle", oblock, 2, range,
                 ops_arg_dat(omar, 1, LOCALSTENCIL, "int", OPS_RW));
    ops_printf("Single-dimension test of ops native calls succeeded!\n");
#ifdef FIELD
    std::vector<int> iterRng;
    iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
    ops_par_loop(KerSetValuesSingle, "KerSetValuesSingle", block.Get(), 2,
                 iterRng.data(),
                 ops_arg_dat(mar[block.ID()], 1, LOCALSTENCIL, "int", OPS_RW));
    ops_printf("Single-dimension test of field and block calls succeeded!\n");
#endif  // FIELD
    ops_par_loop(KerSetValues, "KerSetValues", oblock, 2, range,
                 ops_arg_dat(ogeo, 4, LOCALSTENCIL, "double", OPS_RW));
    ops_printf("Multi-dimension test of ops native calls succeeded!\n");
#ifdef FIELD
    ops_par_loop(
        KerSetValues, "KerSetValues", block.Get(), 2, iterRng.data(),
        ops_arg_dat(geo[block.ID()], 4, LOCALSTENCIL, "double", OPS_RW));
    ops_printf("Multi-dimension test of field and block calls succeeded!\n");
#endif  // FIELD
#endif  // OPS_2D

    ops_timers(&ct1, &et1);
    ops_printf("\nTotal Wall time %lf\n", et1 - et0);
    // Print OPS performance details to output stream
    ops_timing_output(std::cout);
    ops_exit();
}