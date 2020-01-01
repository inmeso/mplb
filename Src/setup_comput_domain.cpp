
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

/*! @brief   Head files for importing geometry from HDF5 file
 * @author  Jianping Meng
 * @details Declaring kernel functions related to create computing domain
 */
#include "setup_comput_domain.h"
#include "setup_comput_domain_kernel.h"
std::string CASENAME;
int BLOCKNUM = 1;
int HALODEPTH = 0;
int SPACEDIM = 2;
ops_block* g_Block;
ops_dat* g_MacroVars;
ops_dat* g_CoordinateXYZ;
ops_dat* g_NodeType;
ops_dat* g_GeometryProperty;
int HaloRelationNum = 0;
ops_halo* HaloRelations;
ops_halo_group HaloGroups;
int* BlockIterRngWhole;
int* BlockIterRngJmin;
int* BlockIterRngJmax;
int* BlockIterRngImin;
int* BlockIterRngImax;
int* BlockIterRngKmax;
int* BlockIterRngKmin;
int* BlockIterRngBulk;
int* BLOCKSIZE;

const int* BlockSize(const int blockId) {
    return &BLOCKSIZE[blockId * SPACEDIM];
}

const int BlockNum() { return BLOCKNUM; }

void DestroyVariables() {
    delete[] g_Block;
    delete[] g_MacroVars;
    delete[] g_CoordinateXYZ;
    if (HaloRelationNum > 0) delete[] HaloRelations;
    delete[] g_NodeType;
    delete[] g_GeometryProperty;
    delete[] BlockIterRngWhole;
    delete[] BlockIterRngBulk;
    delete[] BlockIterRngImax;
    delete[] BlockIterRngImin;
    delete[] BlockIterRngJmax;
    delete[] BlockIterRngJmin;
    delete[] BLOCKSIZE;
}

void DefineVariables() {
    void* temp = nullptr;
    g_Block = new ops_block[BLOCKNUM];
    g_MacroVars = new ops_dat[BLOCKNUM];
    g_CoordinateXYZ = new ops_dat[BLOCKNUM];
    BlockIterRngWhole = new int[BLOCKNUM * 2 * SPACEDIM];
    BlockIterRngJmin = new int[BLOCKNUM * 2 * SPACEDIM];
    BlockIterRngJmax = new int[BLOCKNUM * 2 * SPACEDIM];
    BlockIterRngImax = new int[BLOCKNUM * 2 * SPACEDIM];
    BlockIterRngImin = new int[BLOCKNUM * 2 * SPACEDIM];
    if (3 == SPACEDIM) {
        BlockIterRngKmax = new int[BLOCKNUM * 2 * SPACEDIM];
        BlockIterRngKmin = new int[BLOCKNUM * 2 * SPACEDIM];
    }
    BlockIterRngBulk = new int[BLOCKNUM * 2 * SPACEDIM];
    int haloDepth = std::max(SchemeHaloNum(), BoundaryHaloNum());

#ifdef debug
    ops_printf("%s%i\n", "DefineVariable: haloDepth=", haloDepth);
#endif
    // max halo depths for the dat in the positive direction
    int* d_p = new int[SPACEDIM];
    // max halo depths for the dat in the negative direction
    int* d_m = new int[SPACEDIM];
    int* base = new int[SPACEDIM];
    for (int cordIdx = 0; cordIdx < SPACEDIM; cordIdx++) {
        d_p[cordIdx] = haloDepth;
        d_m[cordIdx] = -haloDepth;
        base[cordIdx] = 0;
    }
    // problem specific
    // if boundary fitting scheme
    // g_Metrics= new ops_dat[BLOCKNUM];
    // calculate the g_Metrics;
    // if cutting cell method
    g_NodeType = new ops_dat[BLOCKNUM];
    g_GeometryProperty = new ops_dat[BLOCKNUM];
    for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
        std::string label(std::to_string(blockIndex));
        std::string blockName("Block_" + label);
        // The name parameter is not properly typed in the definition of
        // ops_decl_block, so there is a minor warning here.
        g_Block[blockIndex] =
            ops_decl_block(SPACEDIM, (char*)blockName.c_str());
        int* size = new int[SPACEDIM];  // size of the dat
        for (int cordIdx = 0; cordIdx < SPACEDIM; cordIdx++) {
            size[cordIdx] = BlockSize(blockIndex)[cordIdx];
        }
        // we assume a problem is at least 2D
        BlockIterRngWhole[blockIndex * 2 * SPACEDIM] = 0;
        BlockIterRngWhole[blockIndex * 2 * SPACEDIM + 1] = size[0];
        BlockIterRngWhole[blockIndex * 2 * SPACEDIM + 2] = 0;
        BlockIterRngWhole[blockIndex * 2 * SPACEDIM + 3] = size[1];
        BlockIterRngBulk[blockIndex * 2 * SPACEDIM] = 1;
        BlockIterRngBulk[blockIndex * 2 * SPACEDIM + 1] = size[0] - 1;
        BlockIterRngBulk[blockIndex * 2 * SPACEDIM + 2] = 1;
        BlockIterRngBulk[blockIndex * 2 * SPACEDIM + 3] = size[1] - 1;
        BlockIterRngJmax[blockIndex * 2 * SPACEDIM] = 0;
        BlockIterRngJmax[blockIndex * 2 * SPACEDIM + 1] = size[0];
        BlockIterRngJmax[blockIndex * 2 * SPACEDIM + 2] = size[1] - 1;
        BlockIterRngJmax[blockIndex * 2 * SPACEDIM + 3] = size[1];
        BlockIterRngJmin[blockIndex * 2 * SPACEDIM] = 0;
        BlockIterRngJmin[blockIndex * 2 * SPACEDIM + 1] = size[0];
        BlockIterRngJmin[blockIndex * 2 * SPACEDIM + 2] = 0;
        BlockIterRngJmin[blockIndex * 2 * SPACEDIM + 3] = 1;
        BlockIterRngImax[blockIndex * 2 * SPACEDIM] = size[0] - 1;
        BlockIterRngImax[blockIndex * 2 * SPACEDIM + 1] = size[0];
        BlockIterRngImax[blockIndex * 2 * SPACEDIM + 2] = 0;
        BlockIterRngImax[blockIndex * 2 * SPACEDIM + 3] = size[1];
        BlockIterRngImin[blockIndex * 2 * SPACEDIM] = 0;
        BlockIterRngImin[blockIndex * 2 * SPACEDIM + 1] = 1;
        BlockIterRngImin[blockIndex * 2 * SPACEDIM + 2] = 0;
        BlockIterRngImin[blockIndex * 2 * SPACEDIM + 3] = size[1];
        if (3 == SPACEDIM) {
            BlockIterRngWhole[blockIndex * 2 * SPACEDIM + 4] = 0;
            BlockIterRngWhole[blockIndex * 2 * SPACEDIM + 5] = size[2];
            BlockIterRngBulk[blockIndex * 2 * SPACEDIM + 4] = 1;
            BlockIterRngBulk[blockIndex * 2 * SPACEDIM + 5] = size[2] - 1;
            BlockIterRngJmax[blockIndex * 2 * SPACEDIM + 4] = 0;
            BlockIterRngJmax[blockIndex * 2 * SPACEDIM + 5] = size[2];
            BlockIterRngJmin[blockIndex * 2 * SPACEDIM + 4] = 0;
            BlockIterRngJmin[blockIndex * 2 * SPACEDIM + 5] = size[2];
            BlockIterRngImax[blockIndex * 2 * SPACEDIM + 4] = 0;
            BlockIterRngImax[blockIndex * 2 * SPACEDIM + 5] = size[2];
            BlockIterRngImin[blockIndex * 2 * SPACEDIM + 4] = 0;
            BlockIterRngImin[blockIndex * 2 * SPACEDIM + 5] = size[2];
            BlockIterRngKmax[blockIndex * 2 * SPACEDIM] = 0;
            BlockIterRngKmax[blockIndex * 2 * SPACEDIM + 1] = size[0];
            BlockIterRngKmax[blockIndex * 2 * SPACEDIM + 2] = 0;
            BlockIterRngKmax[blockIndex * 2 * SPACEDIM + 3] = size[1];
            BlockIterRngKmax[blockIndex * 2 * SPACEDIM + 4] = size[2] - 1;
            BlockIterRngKmax[blockIndex * 2 * SPACEDIM + 5] = size[2];
            BlockIterRngKmin[blockIndex * 2 * SPACEDIM] = 0;
            BlockIterRngKmin[blockIndex * 2 * SPACEDIM + 1] = size[0];
            BlockIterRngKmin[blockIndex * 2 * SPACEDIM + 2] = 0;
            BlockIterRngKmin[blockIndex * 2 * SPACEDIM + 3] = size[1];
            BlockIterRngKmin[blockIndex * 2 * SPACEDIM + 4] = 0;
            BlockIterRngKmin[blockIndex * 2 * SPACEDIM + 5] = 1;
        }

        std::string dataName;
        dataName = "MacroVars_" + label;
        g_MacroVars[blockIndex] =
            ops_decl_dat(g_Block[blockIndex], NUMMACROVAR, size, base, d_m, d_p,
                         (Real*)temp, RealC, dataName.c_str());
        dataName = "Nodetype_" + label;
        // problem specific -- cut cell method
        g_NodeType[blockIndex] =
            ops_decl_dat(g_Block[blockIndex], 1, size, base, d_m, d_p,
                         (int*)temp, "int", dataName.c_str());
        dataName = "GeometryProperty_" + label;
        g_GeometryProperty[blockIndex] =
            ops_decl_dat(g_Block[blockIndex], 1, size, base, d_m, d_p,
                         (int*)temp, "int", dataName.c_str());
        dataName = "CoordinateXYZ_" + label;
        g_CoordinateXYZ[blockIndex] =
            ops_decl_dat(g_Block[blockIndex], SPACEDIM, size, base, d_m, d_p,
                         (Real*)temp, RealC, dataName.c_str());
        delete[] size;
    }
    delete[] d_p;
    delete[] d_m;
    delete[] base;
}

void WriteNodePropertyToHdf5() {
    for (int blockIndex = 0; blockIndex < BLOCKNUM; blockIndex++) {
        std::string blockName("Block");
        std::string label(std::to_string(blockIndex));
        blockName += label;
        std::string fileName = CASENAME + "_" + blockName + ".h5";
        ops_fetch_block_hdf5_file(g_Block[blockIndex], fileName.c_str());
        ops_fetch_dat_hdf5_file(g_GeometryProperty[blockIndex],
                                fileName.c_str());
        ops_fetch_dat_hdf5_file(g_NodeType[blockIndex], fileName.c_str());
    }
}

void WriteFlowfieldToHdf5() {
    for (int blockIndex = 0; blockIndex < BLOCKNUM; blockIndex++) {
        std::string blockName("Block");
        std::string label(std::to_string(blockIndex));
        blockName += label;
        std::string fileName = CASENAME + "_" + blockName + ".h5";
        ops_fetch_block_hdf5_file(g_Block[blockIndex], fileName.c_str());
        ops_fetch_dat_hdf5_file(g_MacroVars[blockIndex], fileName.c_str());
        ops_fetch_dat_hdf5_file(g_CoordinateXYZ[blockIndex], fileName.c_str());
    }
}

void SetupDomainNodeType(int blockIndex, VertexType* faceType,
                         VertexType* edgeType, VertexType* cornerType) {
    int nodeType = (int)VertexType::Fluid;
    int* iterRange = BlockIterRng(blockIndex, BlockIterRngBulk);
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        iterRange, ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
        ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_WRITE));

    // specify halo points
    nodeType = (int)VertexType::ImmersedSolid;
    iterRange = BlockIterRng(blockIndex, BlockIterRngJmin);
    int* haloIterRng = new int[2 * SPACEDIM];
    haloIterRng[0] = iterRange[0] - 1;
    haloIterRng[1] = iterRange[1] + 1;
    haloIterRng[2] = iterRange[2] - 1;
    haloIterRng[3] = iterRange[3] - 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        haloIterRng, ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
        ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_WRITE));
    iterRange = BlockIterRng(blockIndex, BlockIterRngJmax);
    haloIterRng[0] = iterRange[0] - 1;
    haloIterRng[1] = iterRange[1] + 1;
    haloIterRng[2] = iterRange[2] + 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        haloIterRng, ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
        ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_WRITE));
    iterRange = BlockIterRng(blockIndex, BlockIterRngImin);
    haloIterRng[0] = iterRange[0] - 1;
    haloIterRng[1] = iterRange[1] - 1;
    haloIterRng[2] = iterRange[2] - 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        haloIterRng, ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
        ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_WRITE));
    iterRange = BlockIterRng(blockIndex, BlockIterRngImax);
    haloIterRng[0] = iterRange[0] + 1;
    haloIterRng[1] = iterRange[1] + 1;
    haloIterRng[2] = iterRange[2] - 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        haloIterRng, ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
        ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_WRITE));
    if (3 == SPACEDIM) {
        iterRange = BlockIterRng(blockIndex, BlockIterRngKmin);
        haloIterRng[0] = iterRange[0] - 1;
        haloIterRng[1] = iterRange[1] + 1;
        haloIterRng[2] = iterRange[2] - 1;
        haloIterRng[3] = iterRange[3] + 1;
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] - 1;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, haloIterRng,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        iterRange = BlockIterRng(blockIndex, BlockIterRngKmax);
        haloIterRng[0] = iterRange[0] - 1;
        haloIterRng[1] = iterRange[1] + 1;
        haloIterRng[2] = iterRange[2] - 1;
        haloIterRng[3] = iterRange[3] + 1;
        haloIterRng[4] = iterRange[4] + 1;
        haloIterRng[5] = iterRange[5] + 1;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, haloIterRng,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
    }
    // specify faces for 2D cases, they are actually edges.
    // 0 VG_IP left, 1 VG_IM right
    // 2 VG_JP bottom, 3 VG_JM top
    // 4 VG_KP back, 5 VG_KM front
    nodeType = (int)faceType[0];
    iterRange = BlockIterRng(blockIndex, BlockIterRngImin);
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        iterRange, ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
        ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_WRITE));
    nodeType = (int)faceType[1];
    iterRange = BlockIterRng(blockIndex, BlockIterRngImax);
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        iterRange, ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
        ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_WRITE));
    nodeType = (int)faceType[2];
    iterRange = BlockIterRng(blockIndex, BlockIterRngJmin);
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        iterRange, ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
        ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_WRITE));
    nodeType = (int)faceType[3];
    iterRange = BlockIterRng(blockIndex, BlockIterRngJmax);
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        iterRange, ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
        ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_WRITE));

    if (3 == SPACEDIM) {
        nodeType = (int)faceType[4];
        iterRange = BlockIterRng(blockIndex, BlockIterRngKmin);
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iterRange,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        nodeType = (int)faceType[5];
        iterRange = BlockIterRng(blockIndex, BlockIterRngKmax);
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iterRange,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
    }

    const int nx = BlockSize(blockIndex)[0];
    const int ny = BlockSize(blockIndex)[1];
    // 2D Domain corner points 4 types
    // 0 IPJP leftBottom, 1 IPJM leftTop
    // 2 IMJP rightBottom, 3 IMJM rightTop
    if (2 == SPACEDIM) {
        int iminjmin[]{0, 1, 0, 1};
        nodeType = (int)cornerType[0];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmin,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int iminjmax[] = {0, 1, ny - 1, ny};
        nodeType = (int)cornerType[1];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmax,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int imaxjmax[] = {nx - 1, nx, ny - 1, ny};
        nodeType = (int)cornerType[3];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmax,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int imaxjmin[] = {nx - 1, nx, 0, 1};
        nodeType = (int)cornerType[2];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmin,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
    }

    if (3 == SPACEDIM) {
        const int nz = BlockSize(blockIndex)[2];
        // 3D Domain edges 12 types
        // 0 IPJP leftBottom, 1 IPJM leftTop, 2 IMJP rightBottom, 3 IMJM
        // rightTop 4 IPKP leftBack, 5  IPKM leftFront, 6 IMKP rightBack, 7 IMKM
        // rightFront 8 JPKP bottomBack, 9 JPKM bottomFront, 10 JMKP topBack, 11
        // JMKM topFront
        int iminjmin[]{0, 1, 0, 1, 0, nz};
        nodeType = (int)edgeType[0];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmin,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int iminjmax[]{0, 1, ny - 1, ny, 0, nz};
        nodeType = (int)edgeType[1];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmax,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int imaxjmin[]{nx - 1, nx, 0, 1, 0, nz};
        nodeType = (int)edgeType[2];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmin,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int imaxjmax[]{nx - 1, nx, ny - 1, ny, 0, nz};
        nodeType = (int)edgeType[3];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmax,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int iminkmin[]{0, 1, 0, ny, 0, 1};
        nodeType = (int)edgeType[4];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminkmin,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int iminkmax[]{0, 1, 0, ny, nz - 1, nz};
        nodeType = (int)edgeType[5];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminkmax,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int imaxkmin[]{nx - 1, nx, 0, ny, 0, 1};
        nodeType = (int)edgeType[6];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxkmin,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int imaxkmax[]{nx - 1, nx, 0, ny, nz - 1, nz};
        nodeType = (int)edgeType[7];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxkmax,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int jminkmin[]{0, nx, 0, 1, 0, 1};
        nodeType = (int)edgeType[8];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, jminkmin,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int jminkmax[]{0, nx, 0, 1, nz - 1, nz};
        nodeType = (int)edgeType[9];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, jminkmax,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int jmaxkmin[]{0, nx, ny - 1, ny, 0, 1};
        nodeType = (int)edgeType[10];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, jmaxkmin,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int jmaxkmax[]{0, nx, ny - 1, ny, nz - 1, nz};
        nodeType = (int)edgeType[11];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, jmaxkmax,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));

        // 3D domain corners 8 types
        // 0 IPJPKP LeftBottomBack 1 IPJPKM LeftBottomFront
        // 2 IPJMKP LeftTopBack 3 IPJMKM LeftTopFront
        // 4 IMJPKP RightBottomBack 5 IMJPKM RightBottomFront
        // 6 IMJMKP RightTopBack 7 IMJMKM RightTopFront
        int iminjminkmin[]{0, 1, 0, 1, 0, 1};
        nodeType = (int)cornerType[0];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjminkmin,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int iminjminkmax[]{0, 1, 0, 1, nz - 1, nz};
        nodeType = (int)cornerType[1];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjminkmax,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int iminjmaxkmin[]{0, 1, ny - 1, ny, 0, 1};
        nodeType = (int)cornerType[2];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmaxkmin,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int iminjmaxkmax[]{0, 1, ny-1, ny, nz - 1, nz};
        nodeType = (int)cornerType[3];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmaxkmax,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int imaxjminkmin[]{nx - 1, nx, 0, 1, 0, 1};
        nodeType = (int)cornerType[4];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjminkmin,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int imaxjminkmax[]{nx - 1, nx, 0, 1, nz - 1, nz};
        nodeType = (int)cornerType[5];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjminkmax,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int imaxjmaxkmin[]{nx - 1, nx, ny - 1, ny, 0, 1};
        nodeType = (int)cornerType[6];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmaxkmin,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
        int imaxjmaxkmax[]{nx - 1, nx, ny-1, ny, nz - 1, nz};
        nodeType = (int)cornerType[7];
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmaxkmax,
                     ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int",
                                 OPS_WRITE));
    }
}

void  SetBlockGeometryProperty(int blockIndex) {
    int geometryProperty = (int)VG_Fluid;
    int* iterRange = BlockIterRng(blockIndex, BlockIterRngBulk);
    ops_par_loop(KerAssignProperty, "KerAssignProperty", g_Block[blockIndex],
                 SPACEDIM, iterRange,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    // specify halo points
    geometryProperty = VG_ImmersedSolid;
    iterRange = BlockIterRng(blockIndex, BlockIterRngJmin);
    int* haloIterRng = new int[2 * SPACEDIM];
    haloIterRng[0] = iterRange[0] - 1;
    haloIterRng[1] = iterRange[1] + 1;
    haloIterRng[2] = iterRange[2] - 1;
    haloIterRng[3] = iterRange[3] - 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(KerAssignProperty, "KerAssignProperty", g_Block[blockIndex],
                 SPACEDIM, haloIterRng,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    iterRange = BlockIterRng(blockIndex, BlockIterRngJmax);
    haloIterRng[0] = iterRange[0] - 1;
    haloIterRng[1] = iterRange[1] + 1;
    haloIterRng[2] = iterRange[2] + 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(KerAssignProperty, "KerAssignProperty", g_Block[blockIndex],
                 SPACEDIM, haloIterRng,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    iterRange = BlockIterRng(blockIndex, BlockIterRngImin);
    haloIterRng[0] = iterRange[0] - 1;
    haloIterRng[1] = iterRange[1] - 1;
    haloIterRng[2] = iterRange[2] - 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(KerAssignProperty, "KerAssignProperty", g_Block[blockIndex],
                 SPACEDIM, haloIterRng,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));

    iterRange = BlockIterRng(blockIndex, BlockIterRngImax);
    haloIterRng[0] = iterRange[0] + 1;
    haloIterRng[1] = iterRange[1] + 1;
    haloIterRng[2] = iterRange[2] - 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(KerAssignProperty, "KerAssignProperty", g_Block[blockIndex],
                 SPACEDIM, haloIterRng,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    if (3 == SPACEDIM) {
        iterRange = BlockIterRng(blockIndex, BlockIterRngKmin);
        haloIterRng[0] = iterRange[0] - 1;
        haloIterRng[1] = iterRange[1] + 1;
        haloIterRng[2] = iterRange[2] - 1;
        haloIterRng[3] = iterRange[3] + 1;
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] - 1;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, haloIterRng,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        iterRange = BlockIterRng(blockIndex, BlockIterRngKmax);
        haloIterRng[0] = iterRange[0] - 1;
        haloIterRng[1] = iterRange[1] + 1;
        haloIterRng[2] = iterRange[2] - 1;
        haloIterRng[3] = iterRange[3] + 1;
        haloIterRng[4] = iterRange[4] + 1;
        haloIterRng[5] = iterRange[5] + 1;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, haloIterRng,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
    }

    // specify domain
    geometryProperty = VG_JP;
    iterRange = BlockIterRng(blockIndex, BlockIterRngJmin);
    ops_par_loop(KerAssignProperty, "KerAssignProperty", g_Block[blockIndex],
                 SPACEDIM, iterRange,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    geometryProperty = VG_JM;
    iterRange = BlockIterRng(blockIndex, BlockIterRngJmax);
    ops_par_loop(KerAssignProperty, "KerAssignProperty", g_Block[blockIndex],
                 SPACEDIM, iterRange,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    geometryProperty = VG_IP;
    iterRange = BlockIterRng(blockIndex, BlockIterRngImin);
    ops_par_loop(KerAssignProperty, "KerAssignProperty", g_Block[blockIndex],
                 SPACEDIM, iterRange,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    geometryProperty = VG_IM;
    iterRange = BlockIterRng(blockIndex, BlockIterRngImax);
    ops_par_loop(KerAssignProperty, "KerAssignProperty", g_Block[blockIndex],
                 SPACEDIM, iterRange,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    if (3 == SPACEDIM) {
        geometryProperty = VG_KP;
        iterRange = BlockIterRng(blockIndex, BlockIterRngKmin);
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iterRange,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        geometryProperty = VG_KM;
        iterRange = BlockIterRng(blockIndex, BlockIterRngKmax);
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iterRange,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
    }

    const int nx = BlockSize(blockIndex)[0];
    const int ny = BlockSize(blockIndex)[1];
    // 2D Domain corner points four types
    if (2 == SPACEDIM) {
        int iminjmin[]{0, 1, 0, 1};
        geometryProperty = VG_IPJP_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminjmax[] = {0, 1, ny - 1, ny};
        geometryProperty = VG_IPJM_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmax[] = {nx - 1, nx, ny - 1, ny};
        geometryProperty = VG_IMJM_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmin[] = {nx - 1, nx, 0, 1};
        geometryProperty = VG_IMJP_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
    }

    if (3 == SPACEDIM) {
        const int nz = BlockSize(blockIndex)[2];
        // 3D Domain edges 12 types
        int iminjmin[]{0, 1, 0, 1, 0, nz};
        geometryProperty = VG_IPJP_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminjmax[]{0, 1, ny - 1, ny, 0, nz};
        geometryProperty = VG_IPJM_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmax[]{nx - 1, nx, ny - 1, ny, 0, nz};
        geometryProperty = VG_IMJM_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmin[]{nx - 1, nx, 0, 1, 0, nz};
        geometryProperty = VG_IMJP_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));

        int iminkmin[]{0, 1, 0, ny, 0, 1};
        geometryProperty = VG_IPKP_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminkmax[]{0, 1, 0, ny, nz - 1, nz};
        geometryProperty = VG_IPKM_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxkmax[]{nx - 1, nx, 0, ny, nz - 1, nz};
        geometryProperty = VG_IMKM_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxkmin[]{nx - 1, nx, 0, ny, 0, 1};
        geometryProperty = VG_IMKP_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));

        int jminkmin[]{0, nx, 0, 1, 0, 1};
        geometryProperty = VG_JPKP_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, jminkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int jminkmax[]{0, nx, 0, 1, nz - 1, nz};
        geometryProperty = VG_JPKM_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, jminkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int jmaxkmax[]{0, nx, ny - 1, ny, nz - 1, nz};
        geometryProperty = VG_JMKM_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, jmaxkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int jmaxkmin[]{0, nx, ny - 1, ny, 0, 1};
        geometryProperty = VG_JMKP_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, jmaxkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));

        // 3D domain corners 8 types
        int iminjminkmin[]{0, 1, 0, 1, 0, 1};
        geometryProperty = VG_IPJPKP_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjminkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminjminkmax[]{0, 1, 0, 1, nz - 1, nz};
        geometryProperty = VG_IPJPKM_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjminkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminjmaxkmin[]{0, 1, ny - 1, ny, 0, 1};
        geometryProperty = VG_IPJMKP_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmaxkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminjmaxkmax[]{0, 1, ny-1, ny, nz - 1, nz};
        geometryProperty = VG_IPJMKM_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmaxkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjminkmin[]{nx - 1, nx, 0, 1, 0, 1};
        geometryProperty = VG_IMJPKP_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjminkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjminkmax[]{nx - 1, nx, 0, 1, nz - 1, nz};
        geometryProperty = VG_IMJPKM_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjminkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmaxkmin[]{nx - 1, nx, ny - 1, ny, 0, 1};
        geometryProperty = VG_IMJMKP_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmaxkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmaxkmax[]{nx - 1, nx, ny-1, ny, nz - 1, nz};
        geometryProperty = VG_IMJMKM_I;
        ops_par_loop(KerAssignProperty, "KerAssignProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmaxkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
    }
}
/*!
 * Example of manually setting up the embedded solid body
 * This case: flow around a circle cylinder
 * setting up both the geometry and note property
 * L Channel Length D Square Length H Channel height
 */
#ifdef OPS_2D
void SetupEmbeddedBodyFlowAroundCircle(const int ratioLD, const int front,
                                       const VertexType surface) {
    // specify the square
    const int nx = BLOCKSIZE[0];
    const int ny = BLOCKSIZE[1];
    // ratioLD=15 is the channel length: Circle diameter
    const int circleDiameter =
        (nx - 1) / ratioLD;  // nx must be the times of ratioLD
    const int circleCenter = (front + 1 / 2) * circleDiameter;
    int blockIndex{0};
    const Real diameter{1};
    Real circlePos[]{front + 0.5, 0.5 * (ny - 1) / circleDiameter};
    // mark all solid points inside the first circle to be ImmersedSolid
    int* bulkRng = BlockIterRng(blockIndex, BlockIterRngBulk);
    ops_par_loop(
        KerSetEmbeddedCircle, "KerSetEmbeddedCircle", g_Block[blockIndex],
        SPACEDIM, bulkRng, ops_arg_gbl(&diameter, 1, "double", OPS_READ),
        ops_arg_gbl(circlePos, SPACEDIM, "double", OPS_READ),
        ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM, LOCALSTENCIL,
                    "double", OPS_READ),
        ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_WRITE),
        ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL, "int",
                    OPS_WRITE));

    //circlePos[0] = front + 3.5;
    // mark all solid points inside the second circle to be ImmersedSolid
    // ops_par_loop(
    //     KerSetEmbeddedCircle, "KerSetEmbeddedCircle", g_Block[blockIndex],
    //     SPACEDIM, bulkRng, ops_arg_gbl(&diameter, 1, "double", OPS_READ),
    //     ops_arg_gbl(circlePos, SPACEDIM, "double", OPS_READ),
    //     ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM, LOCALSTENCIL,
    //                 "double", OPS_READ),
    //     ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_WRITE),
    //     ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL, "int",
    //                 OPS_WRITE));
    // wipe off some solid points that cannot be consideres
    // as a good surface point
    ops_par_loop(
        KerSweep, "KerSweep", g_Block[blockIndex], SPACEDIM, bulkRng,
        ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL, "int",
                    OPS_READ),
        ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_WRITE));
    // sync the Geometry property to reflect the impact of wiping off some
    // solid points
    ops_par_loop(
        KerSyncGeometryProperty, "KerSyncGeometryProperty", g_Block[blockIndex],
        SPACEDIM, bulkRng,
        ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_READ),
        ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL, "int",
                    OPS_RW));
    // set the correct  geometry property e.g., corner types
    // i.e., mark out the surface points
    ops_par_loop(KerSetEmbeddedBodyGeometry, "KerSetEmbeddedBodyGeometry",
                 g_Block[blockIndex], SPACEDIM, bulkRng,
                 ops_arg_dat(g_NodeType[blockIndex], 1, ONEPTLATTICESTENCIL,
                             "int", OPS_RW),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    // set the boundary type
    int nodeType{surface};
    ops_par_loop(
        KerSetEmbeddedBodyBoundary, "KerSetEmbeddedBodyBoundary",
        g_Block[blockIndex], SPACEDIM, bulkRng,
        ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
        ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL, "int",
                    OPS_READ),
        ops_arg_dat(g_NodeType[blockIndex], 1, LOCALSTENCIL, "int", OPS_RW));
}

/*!
 * Example of manually setting up the embedded solid body
 * This case: flow around a square
 * setting up both the geometry and note property
 * L Channel Length D Square Length H Channel height
 */
void SetupEmbeddedBodyFlowAroundSquare(const int ratioLD, const int front,
                                       const VertexType surface) {
    // specify the square
    const int nx = BLOCKSIZE[0];
    const int ny = BLOCKSIZE[1];
    // ratioLD=15 is the channel length: square length
    const int squareLen =
        (nx - 1) / ratioLD;  // nx must be the times of ratioLD
    const int squareXStart = front * squareLen;       // square Bulk
    const int squareXEnd = squareXStart + squareLen;  // square Bulk
    const int squareCenterY = (ny + 1) / 2 - 1;
    const int squareYStart = squareCenterY - squareLen / 2;
    const int squareYEnd = squareCenterY + squareLen / 2;
    int blockIndex{0};
    // GeometryProperty
    int geometryProperty{VG_JM};
    int squareJmin[] = {squareXStart, squareXEnd + 1, squareYStart,
                        squareYStart + 1};
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        squareJmin, ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
        ops_arg_dat(g_GeometryProperty[0], 1, LOCALSTENCIL, "int", OPS_WRITE));
    geometryProperty = VG_JP;
    int squareJmax[] = {squareXStart, squareXEnd + 1, squareYEnd,
                        squareYEnd + 1};
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        squareJmax, ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
        ops_arg_dat(g_GeometryProperty[0], 1, LOCALSTENCIL, "int", OPS_WRITE));
    geometryProperty = VG_IM;
    int squareImin[] = {squareXStart, squareXStart + 1, squareYStart,
                        squareYEnd + 1};
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        squareImin, ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
        ops_arg_dat(g_GeometryProperty[0], 1, LOCALSTENCIL, "int", OPS_WRITE));
    geometryProperty = VG_IP;
    int squareImax[] = {squareXEnd, squareXEnd + 1, squareYStart,
                        squareYEnd + 1};
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        squareImax, ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
        ops_arg_dat(g_GeometryProperty[0], 1, LOCALSTENCIL, "int", OPS_WRITE));
    int squareiminjmin[] = {squareXStart, squareXStart + 1, squareYStart,
                            squareYStart + 1};
    geometryProperty = VG_IMJM_O;
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        squareiminjmin, ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
        ops_arg_dat(g_GeometryProperty[0], 1, LOCALSTENCIL, "int", OPS_WRITE));
    int squareiminjmax[] = {squareXStart, squareXStart + 1, squareYEnd,
                            squareYEnd + 1};
    geometryProperty = VG_IMJP_O;
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        squareiminjmax, ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
        ops_arg_dat(g_GeometryProperty[0], 1, LOCALSTENCIL, "int", OPS_WRITE));
    int squareimaxjmax[] = {squareXEnd, squareXEnd + 1, squareYEnd,
                            squareYEnd + 1};
    geometryProperty = VG_IPJP_O;
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        squareimaxjmax, ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
        ops_arg_dat(g_GeometryProperty[0], 1, LOCALSTENCIL, "int", OPS_WRITE));
    int squareimaxjmin[] = {squareXEnd, squareXEnd + 1, squareYStart,
                            squareYStart + 1};
    geometryProperty = VG_IPJM_O;
    ops_par_loop(
        KerAssignProperty, "KerAssignProperty", g_Block[blockIndex], SPACEDIM,
        squareimaxjmin, ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
        ops_arg_dat(g_GeometryProperty[0], 1, LOCALSTENCIL, "int", OPS_WRITE));
    // nodeType
    int nodeType{VertexType::ImmersedSolid};
    int squareBulkRng[] = {squareXStart + 1, squareXEnd, squareYStart + 1,
                           squareYEnd};
    ops_par_loop(KerAssignProperty, "KerAssignProperty", g_Block[blockIndex],
                 SPACEDIM, squareBulkRng,
                 ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                 ops_arg_dat(g_NodeType[0], 1, LOCALSTENCIL, "int", OPS_WRITE));
    nodeType = surface;
    int squareLeftWallRng[] = {squareXStart, squareXStart + 1, squareYStart,
                               squareYEnd + 1};
    ops_par_loop(KerAssignProperty, "KerAssignProperty", g_Block[blockIndex],
                 SPACEDIM, squareLeftWallRng,
                 ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                 ops_arg_dat(g_NodeType[0], 1, LOCALSTENCIL, "int", OPS_WRITE));

    int squareRightWallRng[] = {squareXEnd, squareXEnd + 1, squareYStart,
                                squareYEnd + 1};
    ops_par_loop(KerAssignProperty, "KerAssignProperty", g_Block[blockIndex],
                 SPACEDIM, squareRightWallRng,
                 ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                 ops_arg_dat(g_NodeType[0], 1, LOCALSTENCIL, "int", OPS_WRITE));
    int squareTopWallRng[] = {squareXStart, squareXEnd + 1, squareYEnd,
                              squareYEnd + 1};
    ops_par_loop(KerAssignProperty, "KerAssignProperty", g_Block[blockIndex],
                 SPACEDIM, squareTopWallRng,
                 ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                 ops_arg_dat(g_NodeType[0], 1, LOCALSTENCIL, "int", OPS_WRITE));
    int squareBottomWallRng[] = {squareXStart, squareXEnd + 1, squareYStart,
                                 squareYStart + 1};
    ops_par_loop(KerAssignProperty, "KerAssignProperty", g_Block[blockIndex],
                 SPACEDIM, squareBottomWallRng,
                 ops_arg_gbl(&nodeType, 1, "int", OPS_READ),
                 ops_arg_dat(g_NodeType[0], 1, LOCALSTENCIL, "int", OPS_WRITE));
}
#endif  // OPS_2D
#ifdef OPS_3D

#endif  // OPS_3D
void InputCaseName() {
    std::cout << "Please input CASENAME:";
    std::cin >> CASENAME;
}

void InputBlcokDimensions(std::ofstream& output) {
    if (output.is_open()) {
        std::cout << "Please input the space dimensional:";
        std::cin >> SPACEDIM;
        std::cout << "Please input total number of blocks:";
        std::cin >> BLOCKNUM;
        output << "Totally " << BLOCKNUM << " of blocks" << std::endl;
        BLOCKSIZE = new int[BLOCKNUM * SPACEDIM];
        std::cout
            << "Please input the dimension of each block (grid point number):"
            << std::endl;
        for (int blockIndex = 0; blockIndex < BLOCKNUM; blockIndex++) {
            std::cout << "Block " << blockIndex << ":" << std::endl;
            output << "Block " << blockIndex << ":" << std::endl;
            for (int coordIndex = 0; coordIndex < SPACEDIM; coordIndex++) {
                std::cout << "At the coordinate " << coordIndex << ":";
                std::cin >> BLOCKSIZE[SPACEDIM * blockIndex + coordIndex];
                output << "At the coordinate " << coordIndex << ":"
                       << BLOCKSIZE[SPACEDIM * blockIndex + coordIndex]
                       << std::endl;
            }
        }
    }
}

void InputSegmentsAndConstructCoordinates(int blockIndex,
                                          Real* coordinates[SPACEDIM],
                                          std::ofstream& output) {
    if (output.is_open()) {
        int segmentNum[SPACEDIM];
        for (int coordIndex = 0; coordIndex < SPACEDIM; coordIndex++) {
            std::cout << "    Input total number of segments at Coordinate "
                      << coordIndex << ":" << std::endl;
            std::cout << "        At the coordinate " << coordIndex << ":";
            std::cin >> segmentNum[coordIndex];
            output << "        At the coordinate " << coordIndex << ":"
                   << segmentNum[coordIndex] << std::endl;
            int* cellNum = new int[segmentNum[coordIndex]];
            Real* endPos = new Real[segmentNum[coordIndex]];
            std::cout << "            Input the total number of cells and the "
                         "position of the final grid point at each segment:"
                      << std::endl;
            int totalCellNum{0};
            for (int segmentIndex = 0; segmentIndex < segmentNum[coordIndex];
                 segmentIndex++) {
                std::cout << "               Total cell number:";
                std::cin >> cellNum[segmentIndex];
                std::cout << "               Position of the final grid:";
                std::cin >> endPos[segmentIndex];
                totalCellNum += cellNum[segmentIndex];
            }
            if ((totalCellNum + 1) !=
                BlockSize(blockIndex)[SPACEDIM * blockIndex + coordIndex]) {
                std::cout << "               The segment setup is not "
                             "consistent with the "
                             "block size"
                          << std::endl;
                std::exit(1);
            }
            coordinates[coordIndex][0] = 0;
            Real startPos{0};
            int coordNodeIdx{1};
            for (int segmentIndex = 0; segmentIndex < segmentNum[coordIndex];
                 segmentIndex++) {
                Real step{(endPos[segmentIndex] - startPos) /
                          cellNum[segmentIndex]};
                for (int nodeIndex = 0; nodeIndex < cellNum[segmentIndex];
                     nodeIndex++) {
                    coordinates[coordIndex][coordNodeIdx] =
                        startPos + (nodeIndex + 1) * step;
                    coordNodeIdx++;
                }
                startPos = endPos[segmentIndex];
            }
            delete[] cellNum;
            delete[] endPos;
        }
    }
}

void SetInitialMacroVar(const int blockIdx) {
    int* iterRng = BlockIterRng(blockIdx, BlockIterRngWhole);
    ops_par_loop(KerSetInitialMacroVars, "KerSetInitialMacroVars",
                 g_Block[blockIdx], SPACEDIM, iterRng,
                 ops_arg_dat(g_CoordinateXYZ[blockIdx], SPACEDIM, LOCALSTENCIL,
                             "double", OPS_READ),
                 ops_arg_idx(),
                 ops_arg_dat(g_MacroVars[blockIdx], NUMMACROVAR, LOCALSTENCIL,
                             "double", OPS_RW));
}

void PrintCoordinates(const int blockIndex, Real* coordinates[SPACEDIM]) {
    std::cout << "****The coordinates***" << std::endl;
    for (int coordIndex = 0; coordIndex < SPACEDIM; coordIndex++) {
        std::cout << "Coordinate " << coordIndex << " :";
        for (int nodeIndex = 0;
             nodeIndex <
             BlockSize(blockIndex)[SPACEDIM * blockIndex + coordIndex];
             nodeIndex++) {
            std::cout << coordinates[coordIndex][nodeIndex] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "****The coordinates***" << std::endl;
}

void AssignCoordinates(int blockIndex, Real* coordinates[SPACEDIM]) {
#ifdef OPS_2D
    if (SPACEDIM == 2) {
        int* range = BlockIterRng(blockIndex, BlockIterRngWhole);
        ops_par_loop(KerSetCoordinates, "KerSetCoordinates",
                     g_Block[blockIndex], SPACEDIM, range,
                     ops_arg_gbl(coordinates[0], BlockSize(blockIndex)[0],
                                 "double", OPS_READ),
                     ops_arg_gbl(coordinates[1], BlockSize(blockIndex)[1],
                                 "double", OPS_READ),
                     ops_arg_idx(),
                     ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM,
                                 LOCALSTENCIL, "double", OPS_WRITE));
    }
#endif
#ifdef OPS_3D
    if (SPACEDIM == 3) {
        int* range = BlockIterRng(blockIndex, BlockIterRngWhole);
        ops_par_loop(KerSetCoordinates3D, "KerSetCoordinates3D",
                     g_Block[blockIndex], SPACEDIM, range,
                     ops_arg_gbl(coordinates[0], BlockSize(blockIndex)[0],
                                 "double", OPS_READ),
                     ops_arg_gbl(coordinates[1], BlockSize(blockIndex)[1],
                                 "double", OPS_READ),
                     ops_arg_gbl(coordinates[2], BlockSize(blockIndex)[2],
                                 "double", OPS_READ),
                     ops_arg_idx(),
                     ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM,
                                 LOCALSTENCIL, "double", OPS_WRITE));
    }
    std::cout << "Coordinates are assigned!" << std::endl;
#endif
}

void InputBoundaryType(const int blockIndex) {
    int tmpVt;
    // input face property
    VertexType faceType[4];
    std::cout << "Please input boundary at left:";
    std::cin >> tmpVt;
    faceType[0] = (VertexType)tmpVt;
    std::cout << "Please input boundary at right:";
    std::cin >> tmpVt;
    faceType[1] = (VertexType)tmpVt;
    std::cout << "Please input boundary at bottom:";
    std::cin >> tmpVt;
    faceType[2] = (VertexType)tmpVt;
    std::cout << "Please input boundary at top:";
    std::cin >> tmpVt;
    faceType[3] = (VertexType)tmpVt;

    // input corner property
    VertexType cornerType[4];
    std::cout << "Please input boundary at leftBottom:";
    std::cin >> tmpVt;
    cornerType[0] = (VertexType)tmpVt;
    std::cout << "Please input boundary at leftTop:";
    std::cin >> tmpVt;
    cornerType[1] = (VertexType)tmpVt;
    std::cout << "Please input boundary at rightBottom:";
    std::cin >> tmpVt;
    cornerType[2] = (VertexType)tmpVt;
    std::cout << "Please input boundary at rightTop:";
    std::cin >> tmpVt;
    cornerType[3] = (VertexType)tmpVt;

    VertexType* edgeType{nullptr};

    SetupDomainNodeType(blockIndex, faceType, edgeType, cornerType);
}

void InputBoundaryType3D(const int blockIndex) {
    int tmpVt;
    // input face property
    VertexType faceType[6];
    std::cout << "Please input boundary at left:";
    std::cin >> tmpVt;
    faceType[0] = (VertexType)tmpVt;
    std::cout << "Please input boundary at right:";
    std::cin >> tmpVt;
    faceType[1] = (VertexType)tmpVt;
    std::cout << "Please input boundary at bottom:";
    std::cin >> tmpVt;
    faceType[2] = (VertexType)tmpVt;
    std::cout << "Please input boundary at top:";
    std::cin >> tmpVt;
    faceType[3] = (VertexType)tmpVt;
    std::cout << "Please input boundary at front:";
    std::cin >> tmpVt;
    faceType[5] = (VertexType)tmpVt;
    std::cout << "Please input boundary at back:";
    std::cin >> tmpVt;
    faceType[4] = (VertexType)tmpVt;

    // input edge type
    VertexType edgeType[12];
    std::cout << "Please input boundary at leftBottom:";
    std::cin >> tmpVt;
    edgeType[0] = (VertexType)tmpVt;
    std::cout << "Please input boundary at LeftTop:";
    std::cin >> tmpVt;
    edgeType[1] = (VertexType)tmpVt;
    std::cout << "Please input boundary at rightBottom:";
    std::cin >> tmpVt;
    edgeType[2] = (VertexType)tmpVt;
    std::cout << "Please input boundary at rightTop:";
    std::cin >> tmpVt;
    edgeType[3] = (VertexType)tmpVt;
    std::cout << "Please input boundary at leftBack:";
    std::cin >> tmpVt;
    edgeType[4] = (VertexType)tmpVt;
    std::cout << "Please input boundary at leftFront:";
    std::cin >> tmpVt;
    edgeType[5] = (VertexType)tmpVt;
    std::cout << "Please input boundary at rightBack:";
    std::cin >> tmpVt;
    edgeType[6] = (VertexType)tmpVt;
    std::cout << "Please input boundary at rightFront:";
    std::cin >> tmpVt;
    edgeType[7] = (VertexType)tmpVt;
    std::cout << "Please input boundary at bottomBack:";
    std::cin >> tmpVt;
    edgeType[8] = (VertexType)tmpVt;
    std::cout << "Please input boundary at bottomFront:";
    std::cin >> tmpVt;
    edgeType[9] = (VertexType)tmpVt;
    std::cout << "Please input boundary at topBack:";
    std::cin >> tmpVt;
    edgeType[10] = (VertexType)tmpVt;
    std::cout << "Please input boundary at topFront:";
    std::cin >> tmpVt;
    edgeType[11] = (VertexType)tmpVt;

    // input corner property
    VertexType cornerType[8];
    std::cout << "Please input boundary at leftBottomBack:";
    std::cin >> tmpVt;
    cornerType[0] = (VertexType)tmpVt;
    std::cout << "Please input boundary at LeftBottomFront:";
    std::cin >> tmpVt;
    cornerType[1] = (VertexType)tmpVt;
    std::cout << "Please input boundary at LeftTopBack:";
    std::cin >> tmpVt;
    cornerType[2] = (VertexType)tmpVt;
    std::cout << "Please input boundary at LeftTopFront:";
    std::cin >> tmpVt;
    cornerType[3] = (VertexType)tmpVt;
    std::cout << "Please input boundary at RightBottomBack:";
    std::cin >> tmpVt;
    cornerType[4] = (VertexType)tmpVt;
    std::cout << "Please input boundary at RightBottomFront:";
    std::cin >> tmpVt;
    cornerType[5] = (VertexType)tmpVt;
    std::cout << "Please input boundary at RightTopBack:";
    std::cin >> tmpVt;
    cornerType[6] = (VertexType)tmpVt;
    std::cout << "Please input boundary at RightTopFront:";
    std::cin >> tmpVt;
    cornerType[7] = (VertexType)tmpVt;
    SetupDomainNodeType(blockIndex, faceType, edgeType, cornerType);
}

int main(int argc, char* argv[]) {
    ops_init(argc, argv, 5);
    std::vector<std::string> compoNames{"Fluid"};
    std::vector<int> compoid{0};
    std::vector<std::string> lattNames{"D3Q19"};
    DefineComponents(compoNames, compoid, lattNames);
    std::vector<VariableTypes> marcoVarTypes{Variable_Rho, Variable_U,
                                             Variable_V, Variable_W};
    std::vector<std::string> macroVarNames{"rho", "u", "v", "w"};
    std::vector<int> macroVarId{0, 1, 2, 3};
    std::vector<int> macroCompoId{0, 0, 0, 0};
    DefineMacroVars(marcoVarTypes, macroVarNames, macroVarId, macroCompoId);
    std::vector<CollisionType> equTypes{Equilibrium_BGKIsothermal2nd};
    std::vector<int> equCompoId{0};
    DefineCollision(equTypes, equCompoId);

    SetupBoundary();
    SetupScheme();
    InputCaseName();
    std::ofstream geometryList;
    geometryList.open(CASENAME + "_GeometryList.dat");
    geometryList << "CASENAME: " << CASENAME << std::endl;
    InputBlcokDimensions(geometryList);
    DefineVariables();
    ops_partition((char*)"LBMPreProcessor");

    std::cout << "Now constructing each block/domain..." << std::endl;
    geometryList << "Now constructing each block/domain..." << std::endl;
    for (int blockIndex = 0; blockIndex < BLOCKNUM; blockIndex++) {
        Real* coordinates[SPACEDIM];
        std::cout << "****Constructing coordinates for Block " << blockIndex
                  << "****" << std::endl;
        geometryList << "****Constructing coordinates for Block " << blockIndex
                     << "****" << std::endl;
        for (int coordIndex = 0; coordIndex < SPACEDIM; coordIndex++) {
            coordinates[coordIndex] = new Real[BlockSize(
                blockIndex)[SPACEDIM * blockIndex + coordIndex]];
        }
        InputSegmentsAndConstructCoordinates(blockIndex, coordinates,
                                             geometryList);
        AssignCoordinates(blockIndex, coordinates);
        std::cout << "****Constructing coordinates for Block " << blockIndex
                  << " finished****" << std::endl;
        geometryList << "****Constructing coordinates for Block " << blockIndex
                     << " finished****" << std::endl;
         SetBlockGeometryProperty(blockIndex);
        std::cout << "****Constructing boundary property for Block "
                  << blockIndex << "****" << std::endl;
        geometryList << "****Constructing boundary property for Block "
                     << blockIndex << "****" << std::endl;
        if (3 == SPACEDIM) {
            InputBoundaryType3D(blockIndex);
        }
        if (2 == SPACEDIM) {
            InputBoundaryType(blockIndex);
        }
        std::cout << "****Constructing boundary property for Block "
                  << blockIndex << "finished****" << std::endl;
        geometryList << "****Constructing boundary property for "
                        "Block "
                     << blockIndex << "finished****" << std::endl;
        for (int coordIndex = 0; coordIndex < SPACEDIM; coordIndex++) {
            if (coordinates[coordIndex] != nullptr)
                delete[] coordinates[coordIndex];
        }
    }
    geometryList.close();
    for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
        int* iterRng = BlockIterRng(blockIndex, BlockIterRngWhole);
        ops_par_loop(KerSetInitialMacroVars, "KerSetInitialMacroVars",
                     g_Block[blockIndex], SPACEDIM, iterRng,
                     ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM,
                                 LOCALSTENCIL, "double", OPS_READ),
                     ops_arg_idx(),
                     ops_arg_dat(g_MacroVars[blockIndex], NUMMACROVAR,
                                 LOCALSTENCIL, "double", OPS_RW));
    }
    //SetupEmbeddedBodyFlowAroundCircle(15, 4, Vertex_EQMDiffuseRefl);
    WriteFlowfieldToHdf5();
    WriteNodePropertyToHdf5();
    DestroyVariables();
    ops_timing_output(stdout);
    ops_exit();
}
