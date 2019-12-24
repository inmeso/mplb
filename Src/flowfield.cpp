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

/*! @brief   Implementing functions related to the flow field
 * @author  Jianping Meng
 * @details Implementing functions related to create the flow
 * field (allocate memory), set up the geometry and the boundary
 * property, and deallocate the memory.
 */

#include "flowfield.h"
std::string CASENAME;
bool TRANSIENT{false};
int BLOCKNUM{1};
/*!
 * SPACEDIM=2 for 2D 3 for three 3D
 */
int SPACEDIM{2};
/**
 * Layers of halos
 */
int HALODEPTH{0};
ops_block* g_Block{nullptr};
Field<Real> g_f{"f"};
Field<Real> g_fStage{"fStage"};
Field<Real> g_MacroVars{"MacroVars"};
Field<Real> g_MacroVarsCopy{"MacroVars_Copy"};
Real* g_ResidualError{nullptr};
ops_reduction* g_ResidualErrorHandle{nullptr};
Field<Real> g_MacroBodyforce{"MacroBodyForce"};
/**
 * DT: time step
 */
Real DT{1};
/**
 * TAUREF: the reference relaxation time
 * In appropriate non-dimensional system, it is the Knudsen number
 * It must be a constant during the run time
 */
Real* TAUREF{nullptr};
Field<Real> g_CoordinateXYZ{"CoordinateXYZ"};
std::vector<std::vector<std::vector<Real>>> COORDINATES;

Field<int> g_NodeType{"NodeType"};
Field<int> g_GeometryProperty{"GeometryProperty"};

/*!
 * Formal collection of halo relations required by the OPS library
 */
std::vector<ops_halo_group> HALOGROUPS;

int* BlockIterRngWhole{nullptr};
int* BlockIterRngJmin{nullptr};
int* BlockIterRngJmax{nullptr};
int* BlockIterRngImin{nullptr};
int* BlockIterRngImax{nullptr};
int* BlockIterRngKmax{nullptr};
int* BlockIterRngKmin{nullptr};
int* BlockIterRngBulk{nullptr};
/*!
 * The size of each block, i.e., each domain
 */
int* BLOCKSIZE{nullptr};

const int HaloPtNum() { return std::max(SchemeHaloNum(), BoundaryHaloNum()); }

void DefineCase(const std::string& caseName, const int spaceDim, const bool transient) {
    SetCaseName(caseName);
    SPACEDIM = spaceDim;
    TRANSIENT = transient;
}

bool IsTransient() { return TRANSIENT; }

/*!
 * Manually setting up all the variables necessary for the simulation
 * This function has not been updated for 3D problems.
 */
/*
void DefineVariables() {
    void* temp = NULL;
    g_Block = new ops_block[BLOCKNUM];
    g_f = new ops_dat[BLOCKNUM];
    g_Bodyforce = new ops_dat[BLOCKNUM];
    g_fStage = new ops_dat[BLOCKNUM];
    g_MacroVars = new ops_dat[BLOCKNUM];
    g_CoordinateXYZ = new ops_dat[BLOCKNUM];
    BlockIterRngWhole = new int[BLOCKNUM * 2 * SPACEDIM];
    BlockIterRngJmin = new int[BLOCKNUM * 2 * SPACEDIM];
    BlockIterRngJmax = new int[BLOCKNUM * 2 * SPACEDIM];
    BlockIterRngImax = new int[BLOCKNUM * 2 * SPACEDIM];
    BlockIterRngImin = new int[BLOCKNUM * 2 * SPACEDIM];
    BlockIterRngBulk = new int[BLOCKNUM * 2 * SPACEDIM];
    // if steady flow
    g_MacroVarsCopy = new ops_dat[BLOCKNUM];
    g_ResidualErrorHandle = new ops_reduction[MacroVarsNum()];
    g_ResidualError = new Real[2 * MacroVarsNum()];
    // end if steady flow
    int haloDepth = HaloDepth();
#ifdef debug
    ops_printf("%s%i\n", "DefineVariable: haloDepth=", haloDepth);
#endif
    // max halo depths for the dat in the positive direction
    int d_p[2] = {haloDepth, haloDepth};
    // max halo depths for the dat in the negative direction
    int d_m[2] = {-haloDepth, -haloDepth};
    int base[2] = {0, 0};
    // problem specific
    // if boundary fitting scheme
    // g_Metrics= new ops_dat[BLOCKNUM];
    // calculate the g_Metrics;
    // if cutting cell method
    g_NodeType = new ops_dat[BLOCKNUM];
    g_GeometryProperty = new ops_dat[BLOCKNUM];
    for (int blockIndex = 0; blockIndex < BLOCKNUM; blockIndex++) {
        std::string label(std::to_string(blockIndex));
        std::string blockName("Block_" + label);
        // The name parameter is not properly typed in the definition of
        // ops_decl_block, so there is a minor warning here.
        g_Block[blockIndex] =
            ops_decl_block(SPACEDIM, (char*)blockName.c_str());
        int size[2] = {BlockSize(blockIndex)[0],
                       BlockSize(blockIndex)[1]};  // size of the dat

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

        std::string dataName("f_");
        dataName += label;
        g_f[blockIndex] =
            ops_decl_dat(g_Block[blockIndex], NUMXI, size, base, d_m, d_p,
                         (Real*)temp, RealC, dataName.c_str());
        dataName = "fStage_" + label;
        g_fStage[blockIndex] =
            ops_decl_dat(g_Block[blockIndex], NUMXI, size, base, d_m, d_p,
                         (Real*)temp, RealC, dataName.c_str());
        dataName = "MacroBodyForce_" + label;
        g_Bodyforce[blockIndex] =
            ops_decl_dat(g_Block[blockIndex], NUMXI, size, base, d_m, d_p,
                         (Real*)temp, RealC, dataName.c_str());
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
            ops_decl_dat(g_Block[blockIndex], 2, size, base, d_m, d_p,
                         (Real*)temp, RealC, dataName.c_str());
        // if steady flow
        // in the future, we may consider to add an option for the "if"
        dataName = "MacroVars_Copy" + label;
        g_MacroVarsCopy[blockIndex] =
            ops_decl_dat(g_Block[blockIndex], NUMMACROVAR, size, base, d_m, d_p,
                         (Real*)temp, RealC, dataName.c_str());
        for (int localIdx = 0; localIdx < MacroVarsNum(); localIdx++) {
            g_ResidualErrorHandle[localIdx] = ops_decl_reduction_handle(
                // this is double
                sizeof(double), "double", MacroVarName()[localIdx].c_str());
        }
        // end if steady flow
    }
}
*/

// void AllocateMemory() {

//     int haloDepth{HaloPtNum()};
//     HALODEPTH = HaloPtNum();
// }

/*!
 * setting up all the variables necessary for the simulation from a HDF5 file
 * This function can be used for both 2D and 3D cases
 */
// void DefineVariablesFromHDF5() {
//     void* temp = NULL;
//     g_Block = new ops_block[BLOCKNUM];
//     g_f = new ops_dat[BLOCKNUM];
//     g_MacroBodyforce = new ops_dat[BLOCKNUM];
//     g_fStage = new ops_dat[BLOCKNUM];
//     g_MacroVars = new ops_dat[BLOCKNUM];
//     g_CoordinateXYZ = new ops_dat[BLOCKNUM];
//     BlockIterRngWhole = new int[BLOCKNUM * 2 * SPACEDIM];
//     BlockIterRngJmin = new int[BLOCKNUM * 2 * SPACEDIM];
//     BlockIterRngJmax = new int[BLOCKNUM * 2 * SPACEDIM];
//     BlockIterRngImax = new int[BLOCKNUM * 2 * SPACEDIM];
//     BlockIterRngImin = new int[BLOCKNUM * 2 * SPACEDIM];
//     if (3 == SPACEDIM) {
//         BlockIterRngKmax = new int[BLOCKNUM * 2 * SPACEDIM];
//         BlockIterRngKmin = new int[BLOCKNUM * 2 * SPACEDIM];
//     }
//     BlockIterRngBulk = new int[BLOCKNUM * 2 * SPACEDIM];
//     // if steady flow
//     // g_MacroVarsCopy = new ops_dat[BLOCKNUM];
//     g_ResidualErrorHandle = new ops_reduction[MacroVarsNum()];
//     g_ResidualError = new Real[2 * MacroVarsNum()];
//     // end if steady flow
//     int haloDepth = HaloDepth();
//     // max halo depths for the dat in the positive direction
//     int* d_p = new int[SPACEDIM];
//     // max halo depths for the dat in the negative direction
//     int* d_m = new int[SPACEDIM];
//     int* base = new int[SPACEDIM];
//     for (int cordIdx = 0; cordIdx < SPACEDIM; cordIdx++) {
//         d_p[cordIdx] = haloDepth;
//         d_m[cordIdx] = -haloDepth;
//         base[cordIdx] = 0;
//     }
//     // problem specific
//     // if boundary fitting scheme
//     // g_Metrics= new ops_dat[BLOCKNUM];
//     // calculate the g_Metrics;
//     // if cutting cell method
//     g_NodeType = new ops_dat[BLOCKNUM];
//     g_GeometryProperty = new ops_dat[BLOCKNUM];
//     for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
//         std::string label(std::to_string(blockIndex));
//         std::string blockName("Block_" + label);
//         std::string fileName = CASENAME + "_" + "Block" + label + ".h5";
//         // The name parameter is not properly typed in the definition of
//         // ops_decl_block, so there is a minor warning here.
//         g_Block[blockIndex] =
//             ops_decl_block(SPACEDIM, (char*)blockName.c_str());
//         int* size = new int[SPACEDIM];  // size of the dat
//         for (int cordIdx = 0; cordIdx < SPACEDIM; cordIdx++) {
//             size[cordIdx] = BlockSize(blockIndex)[cordIdx];
//         }
//         // we assume a problem is at least 2D
//         BlockIterRngWhole[blockIndex * 2 * SPACEDIM] = 0;
//         BlockIterRngWhole[blockIndex * 2 * SPACEDIM + 1] = size[0];
//         BlockIterRngWhole[blockIndex * 2 * SPACEDIM + 2] = 0;
//         BlockIterRngWhole[blockIndex * 2 * SPACEDIM + 3] = size[1];
//         BlockIterRngBulk[blockIndex * 2 * SPACEDIM] = 1;
//         BlockIterRngBulk[blockIndex * 2 * SPACEDIM + 1] = size[0] - 1;
//         BlockIterRngBulk[blockIndex * 2 * SPACEDIM + 2] = 1;
//         BlockIterRngBulk[blockIndex * 2 * SPACEDIM + 3] = size[1] - 1;
//         BlockIterRngJmax[blockIndex * 2 * SPACEDIM] = 0;
//         BlockIterRngJmax[blockIndex * 2 * SPACEDIM + 1] = size[0];
//         BlockIterRngJmax[blockIndex * 2 * SPACEDIM + 2] = size[1] - 1;
//         BlockIterRngJmax[blockIndex * 2 * SPACEDIM + 3] = size[1];
//         BlockIterRngJmin[blockIndex * 2 * SPACEDIM] = 0;
//         BlockIterRngJmin[blockIndex * 2 * SPACEDIM + 1] = size[0];
//         BlockIterRngJmin[blockIndex * 2 * SPACEDIM + 2] = 0;
//         BlockIterRngJmin[blockIndex * 2 * SPACEDIM + 3] = 1;
//         BlockIterRngImax[blockIndex * 2 * SPACEDIM] = size[0] - 1;
//         BlockIterRngImax[blockIndex * 2 * SPACEDIM + 1] = size[0];
//         BlockIterRngImax[blockIndex * 2 * SPACEDIM + 2] = 0;
//         BlockIterRngImax[blockIndex * 2 * SPACEDIM + 3] = size[1];
//         BlockIterRngImin[blockIndex * 2 * SPACEDIM] = 0;
//         BlockIterRngImin[blockIndex * 2 * SPACEDIM + 1] = 1;
//         BlockIterRngImin[blockIndex * 2 * SPACEDIM + 2] = 0;
//         BlockIterRngImin[blockIndex * 2 * SPACEDIM + 3] = size[1];
//         if (3 == SPACEDIM) {
//             BlockIterRngWhole[blockIndex * 2 * SPACEDIM + 4] = 0;
//             BlockIterRngWhole[blockIndex * 2 * SPACEDIM + 5] = size[2];
//             BlockIterRngBulk[blockIndex * 2 * SPACEDIM + 4] = 1;
//             BlockIterRngBulk[blockIndex * 2 * SPACEDIM + 5] = size[2] - 1;
//             BlockIterRngJmax[blockIndex * 2 * SPACEDIM + 4] = 0;
//             BlockIterRngJmax[blockIndex * 2 * SPACEDIM + 5] = size[2];
//             BlockIterRngJmin[blockIndex * 2 * SPACEDIM + 4] = 0;
//             BlockIterRngJmin[blockIndex * 2 * SPACEDIM + 5] = size[2];
//             BlockIterRngImax[blockIndex * 2 * SPACEDIM + 4] = 0;
//             BlockIterRngImax[blockIndex * 2 * SPACEDIM + 5] = size[2];
//             BlockIterRngImin[blockIndex * 2 * SPACEDIM + 4] = 0;
//             BlockIterRngImin[blockIndex * 2 * SPACEDIM + 5] = size[2];
//             BlockIterRngKmax[blockIndex * 2 * SPACEDIM] = 0;
//             BlockIterRngKmax[blockIndex * 2 * SPACEDIM + 1] = size[0];
//             BlockIterRngKmax[blockIndex * 2 * SPACEDIM + 2] = 0;
//             BlockIterRngKmax[blockIndex * 2 * SPACEDIM + 3] = size[1];
//             BlockIterRngKmax[blockIndex * 2 * SPACEDIM + 4] = size[2] - 1;
//             BlockIterRngKmax[blockIndex * 2 * SPACEDIM + 5] = size[2];
//             BlockIterRngKmin[blockIndex * 2 * SPACEDIM] = 0;
//             BlockIterRngKmin[blockIndex * 2 * SPACEDIM + 1] = size[0];
//             BlockIterRngKmin[blockIndex * 2 * SPACEDIM + 2] = 0;
//             BlockIterRngKmin[blockIndex * 2 * SPACEDIM + 3] = size[1];
//             BlockIterRngKmin[blockIndex * 2 * SPACEDIM + 4] = 0;
//             BlockIterRngKmin[blockIndex * 2 * SPACEDIM + 5] = 1;
//         }
//         std::string dataName("f_");
//         dataName += label;
//         g_f[blockIndex] =
//             ops_decl_dat(g_Block[blockIndex], NUMXI, size, base, d_m, d_p,
//                          (Real*)temp, RealC, dataName.c_str());
//         dataName = "fStage_" + label;
//         g_fStage[blockIndex] =
//             ops_decl_dat(g_Block[blockIndex], NUMXI, size, base, d_m, d_p,
//                          (Real*)temp, RealC, dataName.c_str());
//         dataName = "MacroBodyForce_" + label;
//         const int bodyForceSize{SPACEDIM * NUMCOMPONENTS};
//         g_MacroBodyforce[blockIndex] =
//             ops_decl_dat(g_Block[blockIndex], bodyForceSize, size, base, d_m,
//                          d_p, (Real*)temp, RealC, dataName.c_str());
//         dataName = "MacroVars_" + label;
//         g_MacroVars[blockIndex] =
//             ops_decl_dat_hdf5(g_Block[blockIndex], NUMMACROVAR, "double",
//                               dataName.c_str(), fileName.c_str());
//         dataName = "Nodetype_" + label;
//         // problem specific -- cut cell method
//         g_NodeType[blockIndex] = ops_decl_dat_hdf5(
//             g_Block[blockIndex], NUMCOMPONENTS, "int", dataName.c_str(), fileName.c_str());
//         dataName = "GeometryProperty_" + label;
//         g_GeometryProperty[blockIndex] = ops_decl_dat_hdf5(
//             g_Block[blockIndex], 1, "int", dataName.c_str(), fileName.c_str());
//         dataName = "CoordinateXYZ_" + label;
//         g_CoordinateXYZ[blockIndex] =
//             ops_decl_dat_hdf5(g_Block[blockIndex], SPACEDIM, RealC,
//                               dataName.c_str(), fileName.c_str());
//         // if steady flow
//         // in the future, we may consider to add an option for the "if"
//         // dataName = "MacroVars_Copy" + label;
//         // g_MacroVarsCopy[blockIndex] =
//         //     ops_decl_dat(g_Block[blockIndex], NUMMACROVAR, size, base, d_m, d_p,
//         //                  (Real*)temp, RealC, dataName.c_str());
//         for (int localIdx = 0; localIdx < MacroVarsNum(); localIdx++) {
//             g_ResidualErrorHandle[localIdx] = ops_decl_reduction_handle(
//                 // this is double
//                 sizeof(double), "double", MacroVarName()[localIdx].c_str());
//         }
//         // end if steady flow
//         delete[] size;
//     }
//     delete[] d_p;
//     delete[] d_m;
//     delete[] base;
// }
/*!
 * Manually define the halo relation between blocks.
 * When usingS
 */
void DefineHaloTransfer() {
    /*! @brief Defining the halo relationship
     *  @details Currently we need to manually define them,
     *  will be modified to read CGNF format in the future
     **/
    //     HaloRelationNum = 8;
    //     HaloRelations = new ops_halo[HaloRelationNum];
    //     int haloDepth = HaloDepth();
    //     int d_p[2] =
    //     { haloDepth, haloDepth }; //max halo depths for the dat in the
    //     positive
    //     direction
    //     int d_m[2] =
    //     { -haloDepth, -haloDepth }; //max halo depths for the dat in the
    //     negative direction
    //     int nx = BlockSize ( 0 ) [0];
    //     int ny = BlockSize ( 0 ) [0];
    //     int dir[] =
    //     { 1, 2 };
    //     {
    //         int halo_iter[] =
    //         { 1, ny[0] + d_p[1] - d_m[1] };
    //         int base_from[] =
    //         { 0, d_m[1] };
    //         int base_to[] =
    //         { nx, d_m[1] };
    //         HaloRelations[0] = ops_decl_halo ( g_f[0], g_f[0], halo_iter,
    //         base_from,
    //                                      base_to, dir, dir );
    //         base_from[0] = nx - 1; // need to be changed
    //         base_to[0] = d_m[0];
    //         HaloRelations[1] = ops_decl_halo ( g_f[0], g_f[0], halo_iter,
    //         base_from,
    //                                      base_to, dir, dir );
    //     }
    //     {
    //         int halo_iter[] =
    //         { nx + d_p[0] - d_m[0], 1 };
    //         int base_from[] =
    //         { d_m[0], 0 };
    //         int base_to[] =
    //         { d_m[0], ny };
    //         HaloRelations[2] = ops_decl_halo ( g_f[0], g_f[0], halo_iter,
    //         base_from,
    //                                      base_to, dir, dir );
    //         base_from[1] = ny - 1; //need to be changed
    //         base_to[1] = d_m[1];
    //         HaloRelations[3] = ops_decl_halo ( g_f[0], g_f[0], halo_iter,
    //         base_from,
    //                                      base_to, dir, dir );
    //     }
    //     // corner points
    //     {
    //         //
    //         int halo_iter[] =
    //         { 1, 1 };
    //         int base_from[] =
    //         { 0, 0 };
    //         int base_to[] =
    //         { nx, ny };
    //         HaloRelations[4] = ops_decl_halo ( g_f[0], g_f[0], halo_iter,
    //         base_from,
    //                                      base_to, dir, dir );
    //         base_from[0] = nx - 1; // need to be changed
    //         base_from[1] = ny - 1;
    //         base_to[0] = d_m[0];
    //         base_to[1] = d_m[1];
    //         HaloRelations[5] = ops_decl_halo ( g_f[0], g_f[0], halo_iter,
    //         base_from,
    //                                      base_to, dir, dir );
    //     }
    //     {
    //         //
    //         int halo_iter[] =
    //         { 1, 1 };
    //         int base_from[] =
    //         { 0, ny - 1 };
    //         int base_to[] =
    //         { nx, d_m[1] };
    //         HaloRelations[6] = ops_decl_halo ( g_f[0], g_f[0], halo_iter,
    //         base_from,
    //                                      base_to, dir, dir );
    //         base_from[0] = nx - 1; // need to be changed
    //         base_from[1] = 0;
    //         base_to[0] = d_m[0];
    //         base_to[1] = ny;
    //         HaloRelations[7] = ops_decl_halo ( g_f[0], g_f[0], halo_iter,
    //         base_from,
    //                                      base_to, dir, dir );
    //     }
    //     HaloGroups = ops_decl_halo_group ( HaloRelationNum,
    //     HaloRelations );

    // HaloRelationNum = 2;
    // HaloRelations = new ops_halo[HaloRelationNum];
    // int haloDepth = HaloDepth();
    // // max halo depths for the dat in the positive direction
    // int d_p[2] = {haloDepth, haloDepth};
    // // max halo depths for the dat in the negative direction
    // int d_m[2] = {-haloDepth, -haloDepth};
    // // The domain size in the Block 0
    // int nx = BlockSize(0)[0];
    // int ny = BlockSize(0)[1];
    // int dir[] = {1, 2};
    // {
    //     int halo_iter[] = {nx + d_p[0] - d_m[0], 1};
    //     int base_from[] = {d_m[0], 0};
    //     int base_to[] = {d_m[0], ny};
    //     HaloRelations[0] = ops_decl_halo(g_f[0], g_f[0], halo_iter, base_from,
    //                                      base_to, dir, dir);
    //     base_from[1] = ny - 1;  // need to be changed
    //     base_to[1] = d_m[1];
    //     HaloRelations[1] = ops_decl_halo(g_f[0], g_f[0], halo_iter, base_from,
    //                                      base_to, dir, dir);
    // }

    // HaloGroups = ops_decl_halo_group(HaloRelationNum, HaloRelations);
}
//This version is for the distribution function
//Other version to follow.
//TODO this function only suitable for single halo
#ifdef OPS_3D
void DefinePeriodicHaloPair3D(const std::vector<int>& haloPair) {
    DefinePeriodicHaloPair3D(haloPair, g_f[0], g_f.HaloDepth());
}

void DefinePeriodicHaloPair3D(const std::vector<int>& haloPair, ops_dat dat, const int haloDepth) {
    if (BlockNum() > 1) {
        ops_printf(
            "Periodic boundary conditions are only valid for single-block "
            "applications!");
        assert(BlockNum() > 1);
    }
    if (dat==nullptr){
        ops_printf(
            "Data must be allocated before defining halo relations!");
        assert(dat==nullptr);

    }
    // max halo depths for the dat in the positive direction
    int d_p[3] = {haloDepth, haloDepth, haloDepth};
    // max halo depths for the dat in the negative direction
    int d_m[3] = {-haloDepth, -haloDepth, -haloDepth};
    // The domain size in the Block 0
    int nx = BlockSize(0)[0];
    int ny = BlockSize(0)[1];
    int nz = BlockSize(0)[2];
    int dir[] = {1, 2, 3};
    for (auto& pair : haloPair) {
        if (0 == pair) {
            // left and right pair
            int halo_iter[] = {haloDepth, ny + d_p[0] - d_m[0],
                               nz + d_p[0] - d_m[0]};
            int base_from[] = {0, d_m[0], d_m[0]};
            int base_to[] = {nx, d_m[0], d_m[0]};
            ops_halo leftToRight = ops_decl_halo(dat, dat, halo_iter, base_from,
                                                 base_to, dir, dir);
            base_from[0] = nx + d_m[1];
            base_to[0] = d_m[1];
            ops_halo rightToLeft = ops_decl_halo(dat, dat, halo_iter, base_from,
                                                 base_to, dir, dir);
            ops_halo group[]{leftToRight, rightToLeft};
            HALOGROUPS.push_back(ops_decl_halo_group(2, group));
        }

        if (1 == pair) {
            // top and bottom pair
            int halo_iter[] = {nx + d_p[0] - d_m[0], haloDepth,
                               nz + d_p[0] - d_m[0]};
            int base_from[] = {d_m[0], 0, d_m[0]};
            int base_to[] = {d_m[0], ny, d_m[0]};
            ops_halo botToTop = ops_decl_halo(dat, dat, halo_iter, base_from,
                                              base_to, dir, dir);
            base_from[1] = ny + d_m[1];
            base_to[1] = d_m[1];
            ops_halo topToBot = ops_decl_halo(dat, dat, halo_iter, base_from,
                                              base_to, dir, dir);
            ops_halo group[]{botToTop, topToBot};
            HALOGROUPS.push_back(ops_decl_halo_group(2, group));
        }

        if (2 == pair) {
            // front and back pair
            int halo_iter[] = {nx + d_p[0] - d_m[0], ny + d_p[0] - d_m[0],
                               haloDepth};
            int base_from[] = {d_m[0], d_m[0], 0};
            int base_to[] = {d_m[0], d_m[0], nz};
            ops_halo backToFront = ops_decl_halo(dat, dat, halo_iter, base_from,
                                                 base_to, dir, dir);
            base_from[2] = nz + d_m[2];
            base_to[2] = d_m[2];
            ops_halo frontToBack = ops_decl_halo(dat, dat, halo_iter, base_from,
                                                 base_to, dir, dir);
            ops_halo group[]{backToFront, frontToBack};
            HALOGROUPS.push_back(ops_decl_halo_group(2, group));
        }
    }
}

#endif //OPS_3D
void Partition() {
    ops_partition((char*)"LBM Solver");
}

void DefineHaloTransfer3D() {
    // This is a hard coded version
    // could be used as an example for user-defined routines.

}
/*
 * We need a name to specify which file to input
 * To be decided: a single filename or an array of filenames
 */

void WriteFlowfieldToHdf5(const SizeType timeStep) {
    for (int blockIndex = 0; blockIndex < BLOCKNUM; blockIndex++) {
        std::string blockName("Block_");
        std::string label(std::to_string(blockIndex));
        std::string time(std::to_string(timeStep));
        blockName += (label + "_" + time);
        std::string fileName = CASENAME + "_" + blockName + ".h5";
        ops_fetch_block_hdf5_file(g_Block[blockIndex], fileName.c_str());
        ops_fetch_dat_hdf5_file(g_MacroVars[blockIndex], fileName.c_str());
        ops_fetch_dat_hdf5_file(g_CoordinateXYZ[blockIndex], fileName.c_str());
        ops_fetch_dat_hdf5_file(g_MacroBodyforce[blockIndex], fileName.c_str());
    }
}

void WriteDistributionsToHdf5(const SizeType timeStep) {
    for (int blockIndex = 0; blockIndex < BLOCKNUM; blockIndex++) {
        std::string blockName("Block_");
        std::string label(std::to_string(blockIndex));
        std::string time(std::to_string(timeStep));
        blockName += (label + "_" + time);
        std::string fileName = CASENAME + "_" + blockName + ".h5";
        ops_fetch_block_hdf5_file(g_Block[blockIndex], fileName.c_str());
        ops_fetch_dat_hdf5_file(g_f[blockIndex], fileName.c_str());
    }
}

void WriteNodePropertyToHdf5(const SizeType timeStep) {
    for (int blockIndex = 0; blockIndex < BLOCKNUM; blockIndex++) {
        std::string blockName("Block_");
        std::string label(std::to_string(blockIndex));
        std::string time(std::to_string(timeStep));
        blockName += (label + "_" + time);
        std::string fileName = CASENAME + "_" + blockName + ".h5";
        ops_fetch_block_hdf5_file(g_Block[blockIndex], fileName.c_str());
        ops_fetch_dat_hdf5_file(g_GeometryProperty[blockIndex],
                                fileName.c_str());
        ops_fetch_dat_hdf5_file(g_NodeType[blockIndex], fileName.c_str());
    }
}

void DefineHaloTransferFromHdf5() {}
/*
 * Importing geometry from an external HDF5 file
 */
// void SetupFlowfieldfromHdf5() {
//     CASENAME = "Cavity3D";  // Input parameter
//     SPACEDIM = 3;
//     BLOCKNUM = 1;  // Input parameter
//     BLOCKSIZE = new int[BLOCKNUM * SPACEDIM];
//     BLOCKSIZE[0] = 101;  // Input parameters
//     BLOCKSIZE[1] = 101;  // Input parameters
//     BLOCKSIZE[2] = 3;    // Input parameters

//     TAUREF = new Real[ComponentNum()];
//     TAUREF[0] = 0.001;  // Input parameters
//     // All above parameters should be written down by the python script
//     Real minDx{1. / 100};  // Input parameters at this moment
//     Real minDy{1. / 100};  // Input parameters at this moment
//     // DT = 0.01 * fmin(minDx, minDy) / MaximumSpeed();  // finite difference
//     // scheme
//     // DT = 0.0001414;
//     // DT = fmin(fmin(minDx, minDy) / MaximumSpeed(),
//     //              0.5 * TAUREF[0]);  // finite difference scheme
//     DT = minDx / SoundSpeed();  // stream-collision
//     HALODEPTH = HaloPtNum();
//     ops_printf("%s\n", "Starting to allocate...");
//     DefineVariablesFromHDF5();
//     DefineHaloTransfer3D();
//     // above calls must be before the ops_partition call
//     ops_partition((char*)"LBM");
// }

const std::string CaseName() { return CASENAME; }
void SetCaseName(const std::string caseName) { CASENAME = caseName; }
void setCaseName(const char* caseName) {
    std::string tmp(caseName);
    CASENAME = tmp;
}

const int BlockNum() { return BLOCKNUM; }
const int SpaceDim() { return SPACEDIM; }
const int HaloDepth() { return HALODEPTH; }

void SetHaloDepth(const int haloDepth) { HALODEPTH = haloDepth; }

void DestroyFlowfield() {
    FreeArrayMemory(g_Block);
    FreeArrayMemory(TAUREF);
    FreeArrayMemory(BlockIterRngWhole);
    FreeArrayMemory(BlockIterRngBulk);
    FreeArrayMemory(BlockIterRngImax);
    FreeArrayMemory(BlockIterRngImin);
    FreeArrayMemory(BlockIterRngJmax);
    FreeArrayMemory(BlockIterRngJmin);
    FreeArrayMemory(BLOCKSIZE);
    // if steady flow
    // FreeArrayMemory(g_MacroVarsCopy);
    FreeArrayMemory(g_ResidualErrorHandle);
    FreeArrayMemory(g_ResidualError);
    if (3 == SPACEDIM) {
        FreeArrayMemory(BlockIterRngKmax);
        FreeArrayMemory(BlockIterRngKmin);
    }
    // end if steady flow
    // delete[] halos;
}

int* IterRngWhole() { return BlockIterRngWhole; }
int* IterRngJmin() { return BlockIterRngJmin; }
int* IterRngJmax() { return BlockIterRngJmax; }
int* IterRngImin() { return BlockIterRngImin; }
int* IterRngImax() { return BlockIterRngImax; }
int* IterRngBulk() { return BlockIterRngBulk; }
int* IterRngKmax() { return BlockIterRngKmax; }
int* IterRngKmin() { return BlockIterRngKmin; }

const int* BlockSize(const int blockId) {
    return &BLOCKSIZE[blockId * SPACEDIM];
}

Real TotalMeshSize() {
    Real size = 1;
    for (int blockIdx = 0; blockIdx < BlockNum(); blockIdx++) {
        for (int cordIdx = 0; cordIdx < SPACEDIM; cordIdx++) {
            size *= BlockSize(blockIdx)[cordIdx];
        }
    }
    return size;
}

const Real TimeStep() { return DT; }
const Real* pTimeStep() { return &DT; }
void SetTimeStep(Real dt) { DT = dt; }
const Real* TauRef() { return TAUREF; }

void SetTauRef(const std::vector<Real> tauRef) {
    const int tauNum = SizeofTau();
    if (tauRef.size() == tauNum) {
        if (nullptr == TAUREF) {
            TAUREF = new Real[tauNum];
        }
        for (int idx = 0; idx < tauNum; idx++) {
            TAUREF[idx] = tauRef[idx];
            // ops_printf("\n tau is %f \n",TAUREF[idx]);
        }
    } else {
        ops_printf("Error! %i taus are required but there are %i!\n", tauNum,
                   tauRef.size());
        assert(tauRef.size() == tauNum);
    }
}

void SetBlockSize(const std::vector<SizeType> blockSize) {
    const int dim{SPACEDIM * BLOCKNUM};
    if (blockSize.size() == dim) {
        BLOCKSIZE = new int[BLOCKNUM * SPACEDIM];

        for (int blockIndex = 0; blockIndex < BLOCKNUM; blockIndex++) {
            for (int coordIndex = 0; coordIndex < SPACEDIM; coordIndex++) {
                BLOCKSIZE[SPACEDIM * blockIndex + coordIndex] =
                    blockSize[SPACEDIM * blockIndex + coordIndex];
            }
        }

    } else {
        ops_printf(
            "Error! %i numbers are required for specifying the size of %i blocks!\n",
            dim, BLOCKNUM);
            assert(blockSize.size() == dim);
    }
}

void SetBlockNum(const SizeType blockNum) {
    if (blockNum > 0) {
        BLOCKNUM = blockNum;
    } else {
        ops_printf("%s\n", "Error! There must be at least one block");
        assert(blockNum > 0);
    }
}

Real GetMaximumResidual(const SizeType checkPeriod) {
    Real maxResError{0};
    Real relResErrorMacroVar{0};
    for (int macroVarIdx = 0; macroVarIdx < MacroVarsNum(); macroVarIdx++) {
        relResErrorMacroVar = g_ResidualError[2 * macroVarIdx] /
                              g_ResidualError[2 * macroVarIdx + 1] /
                              (checkPeriod * TimeStep());

        if (maxResError <= relResErrorMacroVar) {
            maxResError = relResErrorMacroVar;
        }
    }
    return maxResError;
}

const std::vector<ops_halo_group>& HaloGroups() { return HALOGROUPS; }

void  SetBlockGeometryProperty(int blockIndex) {
    int geometryProperty = (int)VG_Fluid;
    int* iterRange = BlockIterRng(blockIndex, IterRngBulk());
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                 g_Block[blockIndex], SPACEDIM, iterRange,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    // specify halo points
    geometryProperty = VG_ImmersedSolid;
    iterRange = BlockIterRng(blockIndex, IterRngJmin());
    int* haloIterRng = new int[2 * SPACEDIM];
    haloIterRng[0] = iterRange[0] - 1;
    haloIterRng[1] = iterRange[1] + 1;
    haloIterRng[2] = iterRange[2] - 1;
    haloIterRng[3] = iterRange[3] - 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                 g_Block[blockIndex], SPACEDIM, haloIterRng,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    iterRange = BlockIterRng(blockIndex, IterRngJmax());
    haloIterRng[0] = iterRange[0] - 1;
    haloIterRng[1] = iterRange[1] + 1;
    haloIterRng[2] = iterRange[2] + 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                 g_Block[blockIndex], SPACEDIM, haloIterRng,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    iterRange = BlockIterRng(blockIndex, IterRngImin());
    haloIterRng[0] = iterRange[0] - 1;
    haloIterRng[1] = iterRange[1] - 1;
    haloIterRng[2] = iterRange[2] - 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                 g_Block[blockIndex], SPACEDIM, haloIterRng,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));

    iterRange = BlockIterRng(blockIndex, IterRngImax());
    haloIterRng[0] = iterRange[0] + 1;
    haloIterRng[1] = iterRange[1] + 1;
    haloIterRng[2] = iterRange[2] - 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                 g_Block[blockIndex], SPACEDIM, haloIterRng,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    if (3 == SPACEDIM) {
        iterRange = BlockIterRng(blockIndex, IterRngKmin());
        haloIterRng[0] = iterRange[0] - 1;
        haloIterRng[1] = iterRange[1] + 1;
        haloIterRng[2] = iterRange[2] - 1;
        haloIterRng[3] = iterRange[3] + 1;
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] - 1;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, haloIterRng,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        iterRange = BlockIterRng(blockIndex, IterRngKmax());
        haloIterRng[0] = iterRange[0] - 1;
        haloIterRng[1] = iterRange[1] + 1;
        haloIterRng[2] = iterRange[2] - 1;
        haloIterRng[3] = iterRange[3] + 1;
        haloIterRng[4] = iterRange[4] + 1;
        haloIterRng[5] = iterRange[5] + 1;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, haloIterRng,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
    }

    // specify domain
    geometryProperty = VG_JP;
    iterRange = BlockIterRng(blockIndex, IterRngJmin());
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                 g_Block[blockIndex], SPACEDIM, iterRange,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    geometryProperty = VG_JM;
    iterRange = BlockIterRng(blockIndex, IterRngJmax());
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                 g_Block[blockIndex], SPACEDIM, iterRange,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    geometryProperty = VG_IP;
    iterRange = BlockIterRng(blockIndex, IterRngImin());
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                 g_Block[blockIndex], SPACEDIM, iterRange,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    geometryProperty = VG_IM;
    iterRange = BlockIterRng(blockIndex, IterRngImax());
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                 g_Block[blockIndex], SPACEDIM, iterRange,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(g_GeometryProperty[blockIndex], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    if (3 == SPACEDIM) {
        geometryProperty = VG_KP;
        iterRange = BlockIterRng(blockIndex, IterRngKmin());
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, iterRange,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        geometryProperty = VG_KM;
        iterRange = BlockIterRng(blockIndex, IterRngKmax());
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
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
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminjmax[] = {0, 1, ny - 1, ny};
        geometryProperty = VG_IPJM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmax[] = {nx - 1, nx, ny - 1, ny};
        geometryProperty = VG_IMJM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmin[] = {nx - 1, nx, 0, 1};
        geometryProperty = VG_IMJP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
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
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminjmax[]{0, 1, ny - 1, ny, 0, nz};
        geometryProperty = VG_IPJM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmax[]{nx - 1, nx, ny - 1, ny, 0, nz};
        geometryProperty = VG_IMJM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmin[]{nx - 1, nx, 0, 1, 0, nz};
        geometryProperty = VG_IMJP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));

        int iminkmin[]{0, 1, 0, ny, 0, 1};
        geometryProperty = VG_IPKP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, iminkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminkmax[]{0, 1, 0, ny, nz - 1, nz};
        geometryProperty = VG_IPKM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, iminkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxkmax[]{nx - 1, nx, 0, ny, nz - 1, nz};
        geometryProperty = VG_IMKM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, imaxkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxkmin[]{nx - 1, nx, 0, ny, 0, 1};
        geometryProperty = VG_IMKP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, imaxkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));

        int jminkmin[]{0, nx, 0, 1, 0, 1};
        geometryProperty = VG_JPKP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, jminkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int jminkmax[]{0, nx, 0, 1, nz - 1, nz};
        geometryProperty = VG_JPKM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, jminkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int jmaxkmax[]{0, nx, ny - 1, ny, nz - 1, nz};
        geometryProperty = VG_JMKM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, jmaxkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int jmaxkmin[]{0, nx, ny - 1, ny, 0, 1};
        geometryProperty = VG_JMKP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, jmaxkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));

        // 3D domain corners 8 types
        int iminjminkmin[]{0, 1, 0, 1, 0, 1};
        geometryProperty = VG_IPJPKP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, iminjminkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminjminkmax[]{0, 1, 0, 1, nz - 1, nz};
        geometryProperty = VG_IPJPKM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, iminjminkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminjmaxkmin[]{0, 1, ny - 1, ny, 0, 1};
        geometryProperty = VG_IPJMKP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmaxkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminjmaxkmax[]{0, 1, ny - 1, ny, nz - 1, nz};
        geometryProperty = VG_IPJMKM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, iminjmaxkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjminkmin[]{nx - 1, nx, 0, 1, 0, 1};
        geometryProperty = VG_IMJPKP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjminkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjminkmax[]{nx - 1, nx, 0, 1, nz - 1, nz};
        geometryProperty = VG_IMJPKM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjminkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmaxkmin[]{nx - 1, nx, ny - 1, ny, 0, 1};
        geometryProperty = VG_IMJMKP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmaxkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmaxkmax[]{nx - 1, nx, ny - 1, ny, nz - 1, nz};
        geometryProperty = VG_IMJMKM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                     g_Block[blockIndex], SPACEDIM, imaxjmaxkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
    }
}

void SetBulkandHaloNodesType(int blockIndex, int compoId) {
    const int fluidType{(int)Vertex_Fluid};
    const int immersedSolidType{(int)Vertex_ImmersedSolid};
    const int boundaryType{(int)Vertex_Boundary};
    int* iterRange = BlockIterRng(blockIndex, IterRngBulk());
    ops_par_loop(KerSetNodeType, "KerSetNodeType", g_Block[blockIndex],
                 SPACEDIM, iterRange,
                 ops_arg_gbl(&fluidType, 1, "int", OPS_READ),
                 ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                             LOCALSTENCIL, "int", OPS_WRITE),
                 ops_arg_gbl(&compoId, 1, "int", OPS_READ));

    iterRange = BlockIterRng(blockIndex, IterRngJmin());
    // Specify general boundary type
    ops_par_loop(KerSetNodeType, "KerSetNodeType", g_Block[blockIndex],
                 SPACEDIM, iterRange,
                 ops_arg_gbl(&boundaryType, 1, "int", OPS_READ),
                 ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                             LOCALSTENCIL, "int", OPS_WRITE),
                 ops_arg_gbl(&compoId, 1, "int", OPS_READ));
    // Specify halo points
    int* haloIterRng = new int[2 * SPACEDIM];
    haloIterRng[0] = iterRange[0] - 1;
    haloIterRng[1] = iterRange[1] + 1;
    haloIterRng[2] = iterRange[2] - 1;
    haloIterRng[3] = iterRange[3] - 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(KerSetNodeType, "KerSetNodeType", g_Block[blockIndex],
                 SPACEDIM, haloIterRng,
                 ops_arg_gbl(&immersedSolidType, 1, "int", OPS_READ),
                 ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                             LOCALSTENCIL, "int", OPS_WRITE),
                 ops_arg_gbl(&compoId, 1, "int", OPS_READ));

    iterRange = BlockIterRng(blockIndex, IterRngJmax());
    // Specify general boundary type
    ops_par_loop(KerSetNodeType, "KerSetNodeType", g_Block[blockIndex],
                 SPACEDIM, iterRange,
                 ops_arg_gbl(&boundaryType, 1, "int", OPS_READ),
                 ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                             LOCALSTENCIL, "int", OPS_WRITE),
                 ops_arg_gbl(&compoId, 1, "int", OPS_READ));
    haloIterRng[0] = iterRange[0] - 1;
    haloIterRng[1] = iterRange[1] + 1;
    haloIterRng[2] = iterRange[2] + 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    // Specify halo points
    ops_par_loop(KerSetNodeType, "KerSetNodeType", g_Block[blockIndex],
                 SPACEDIM, haloIterRng,
                 ops_arg_gbl(&immersedSolidType, 1, "int", OPS_READ),
                 ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                             LOCALSTENCIL, "int", OPS_WRITE),
                 ops_arg_gbl(&compoId, 1, "int", OPS_READ));

    iterRange = BlockIterRng(blockIndex, IterRngImin());
    // Specify general boundary type
    ops_par_loop(KerSetNodeType, "KerSetNodeType", g_Block[blockIndex],
                 SPACEDIM, iterRange,
                 ops_arg_gbl(&boundaryType, 1, "int", OPS_READ),
                 ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                             LOCALSTENCIL, "int", OPS_WRITE),
                 ops_arg_gbl(&compoId, 1, "int", OPS_READ));
    haloIterRng[0] = iterRange[0] - 1;
    haloIterRng[1] = iterRange[1] - 1;
    haloIterRng[2] = iterRange[2] - 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    // Specify halo points
    ops_par_loop(KerSetNodeType, "KerSetNodeType", g_Block[blockIndex],
                 SPACEDIM, haloIterRng,
                 ops_arg_gbl(&immersedSolidType, 1, "int", OPS_READ),
                 ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                             LOCALSTENCIL, "int", OPS_WRITE),
                 ops_arg_gbl(&compoId, 1, "int", OPS_READ));

    iterRange = BlockIterRng(blockIndex, IterRngImax());
    // Specify general boundary type
    ops_par_loop(KerSetNodeType, "KerSetNodeType", g_Block[blockIndex],
                 SPACEDIM, iterRange,
                 ops_arg_gbl(&boundaryType, 1, "int", OPS_READ),
                 ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                             LOCALSTENCIL, "int", OPS_WRITE),
                 ops_arg_gbl(&compoId, 1, "int", OPS_READ));
    haloIterRng[0] = iterRange[0] + 1;
    haloIterRng[1] = iterRange[1] + 1;
    haloIterRng[2] = iterRange[2] - 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    // Specify halo points
    ops_par_loop(KerSetNodeType, "KerSetNodeType", g_Block[blockIndex],
                 SPACEDIM, haloIterRng,
                 ops_arg_gbl(&immersedSolidType, 1, "int", OPS_READ),
                 ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                             LOCALSTENCIL, "int", OPS_WRITE),
                 ops_arg_gbl(&compoId, 1, "int", OPS_READ));

    if (3 == SPACEDIM) {
        iterRange = BlockIterRng(blockIndex, IterRngKmin());
        // Specify general boundary type
        ops_par_loop(KerSetNodeType, "KerSetNodeType", g_Block[blockIndex],
                     SPACEDIM, iterRange,
                     ops_arg_gbl(&boundaryType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                                 LOCALSTENCIL, "int", OPS_WRITE),
                     ops_arg_gbl(&compoId, 1, "int", OPS_READ));
        haloIterRng[0] = iterRange[0] - 1;
        haloIterRng[1] = iterRange[1] + 1;
        haloIterRng[2] = iterRange[2] - 1;
        haloIterRng[3] = iterRange[3] + 1;
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] - 1;
        // Specify halo points
        ops_par_loop(KerSetNodeType, "KerSetNodeType", g_Block[blockIndex],
                     SPACEDIM, haloIterRng,
                     ops_arg_gbl(&immersedSolidType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                                 LOCALSTENCIL, "int", OPS_WRITE),
                     ops_arg_gbl(&compoId, 1, "int", OPS_READ));

        iterRange = BlockIterRng(blockIndex, IterRngKmax());
        // Specify general boundary type
        ops_par_loop(KerSetNodeType, "KerSetNodeType", g_Block[blockIndex],
                     SPACEDIM, iterRange,
                     ops_arg_gbl(&boundaryType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                                 LOCALSTENCIL, "int", OPS_WRITE),
                     ops_arg_gbl(&compoId, 1, "int", OPS_READ));
        haloIterRng[0] = iterRange[0] - 1;
        haloIterRng[1] = iterRange[1] + 1;
        haloIterRng[2] = iterRange[2] - 1;
        haloIterRng[3] = iterRange[3] + 1;
        haloIterRng[4] = iterRange[4] + 1;
        haloIterRng[5] = iterRange[5] + 1;
        // Specify halo points
        ops_par_loop(KerSetNodeType, "KerSetNodeType", g_Block[blockIndex],
                     SPACEDIM, haloIterRng,
                     ops_arg_gbl(&immersedSolidType, 1, "int", OPS_READ),
                     ops_arg_dat(g_NodeType[blockIndex], NUMCOMPONENTS,
                                 LOCALSTENCIL, "int", OPS_WRITE),
                     ops_arg_gbl(&compoId, 1, "int", OPS_READ));
    }
    FreeArrayMemory(haloIterRng);
}

void DefineBlocks(const SizeType blockNum,
                  const std::vector<SizeType>& blockSize, const Real meshSize,
                  const std::vector<Real>& startPos) {
    SetBlockNum(blockNum);
    SetBlockSize(blockSize);
    COORDINATES.resize(blockNum);
    SizeType numBlockStartPos;
    numBlockStartPos = startPos.size();
    if (numBlockStartPos == blockNum * SPACEDIM) {
        for (int blockIndex = 0; blockIndex < blockNum; blockIndex++) {
            COORDINATES.at(blockIndex).resize(SPACEDIM);
            for (int coordIndex = 0; coordIndex < SPACEDIM; coordIndex++) {
                int numOfGridPoints =
                    BlockSize(blockIndex)[SPACEDIM * blockIndex + coordIndex];
                COORDINATES.at(blockIndex)
                    .at(coordIndex)
                    .resize(numOfGridPoints);
                for (int nodeIndex = 0; nodeIndex < numOfGridPoints;
                     nodeIndex++) {
                    COORDINATES.at(blockIndex).at(coordIndex).at(nodeIndex) =
                        startPos.at(coordIndex) + nodeIndex * meshSize;
                }
            }
        }
    } else {
        ops_printf(
            "Error! Expected %i coordinates of three starting points %i, but "
            "received only =%i \n",
            SPACEDIM * blockNum, numBlockStartPos);
        assert(numBlockStartPos == blockNum * SPACEDIM);
    }

    g_Block = new ops_block[BLOCKNUM];
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
    for (int blockIndex = 0; blockIndex < BLOCKNUM; blockIndex++) {
        std::string label(std::to_string(blockIndex));
        std::string blockName("Block_" + label);
        // The name parameter is not properly typed in the definition of
        // ops_decl_block, so there is a minor warning here.
        int* size = new int[SPACEDIM];  // size of the dat
        for (int cordIdx = 0; cordIdx < SPACEDIM; cordIdx++) {
            size[cordIdx] = BlockSize(blockIndex)[cordIdx];
        }
        g_Block[blockIndex] =
            ops_decl_block(SPACEDIM, (char*)blockName.c_str());
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
        delete[] size;
    }
    //TODO temporialy here, later we shall have a define geometry
    g_NodeType.CreateFieldFromScratch();
    g_CoordinateXYZ.SetDataDim(SPACEDIM);
    g_CoordinateXYZ.CreateFieldFromScratch();
    g_GeometryProperty.CreateFieldFromScratch();
}

void AssignCoordinates(int blockIndex,
                       const std::vector<std::vector<Real>>& blockCoordinates) {
#ifdef OPS_2D
    if (SPACEDIM == 2) {
        int* range = BlockIterRng(blockIndex, IterRngWhole());
        ops_par_loop(KerSetCoordinates, "KerSetCoordinates",
                     g_Block[blockIndex], SPACEDIM, range,
                     ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM,
                                 LOCALSTENCIL, "double", OPS_WRITE),
                     ops_arg_idx(),
                     ops_arg_gbl(blockCoordinates.at(0).data(),
                                 BlockSize(blockIndex)[0], "double", OPS_READ),
                     ops_arg_gbl(blockCoordinates.at(1).data(),
                                 BlockSize(blockIndex)[1], "double", OPS_READ));
    }
#endif

#ifdef OPS_3D
    if (SPACEDIM == 3) {
        int* range = BlockIterRng(blockIndex, IterRngWhole());
        ops_par_loop(KerSetCoordinates3D, "KerSetCoordinates3D",
                     g_Block[blockIndex], SPACEDIM, range,
                     ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM,
                                 LOCALSTENCIL, "double", OPS_WRITE),
                     ops_arg_idx(),
                     ops_arg_gbl(blockCoordinates.at(0).data(),
                                 BlockSize(blockIndex)[0], "double", OPS_READ),
                     ops_arg_gbl(blockCoordinates.at(1).data(),
                                 BlockSize(blockIndex)[1], "double", OPS_READ),
                     ops_arg_gbl(blockCoordinates.at(2).data(),
                                 BlockSize(blockIndex)[2], "double", OPS_READ)

        );
    }
#endif
}

void PrepareFlowField() {
    ops_printf("The coordinates are assigned!\n");
    for (int blockId = 0; blockId < BlockNum(); blockId++) {
        SetBlockGeometryProperty(blockId);
        ops_printf("The geometry property for Block %i is set!\n", blockId);
        for (int compoId = 0; compoId < NUMCOMPONENTS; compoId++) {
            SetBulkandHaloNodesType(blockId, compoId);
            ops_printf(
                "The bulk and halo node property are set for Component %i at "
                "Block %i\n",
                compoId, blockId);
        }
        AssignCoordinates(blockId, COORDINATES.at(blockId));
    }
}