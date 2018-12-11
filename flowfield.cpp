// Copyright 2017 the MPLB team. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

/*! @brief   Implementing functions related to the flow field
  * @author  Jianping Meng
  * @details Implementing functions related to create the flow
  * field (allocate memory), set up the geometry and the boundary
  * property, and deallocate the memory.
  */

#include "flowfield.h"
std::string CASENAME;
int g_BlockNum = 1;
int SPACEDIM = 2;
int g_HaloDepth = 0;
ops_block* g_Block;
ops_dat* g_f;
ops_dat* g_fStage;
ops_dat* g_feq;
ops_dat* g_MacroVars;
ops_dat* g_MacroVarsCopy;
Real* g_ResidualError;
ops_reduction* g_ResidualErrorHandle;
ops_dat* g_Bodyforce;
Real g_dt = 1;
Real* KN;
ops_dat* g_Tau;
ops_dat* g_DiscreteConvectionTerm;
ops_dat* g_CoordinateXYZ;
/*!
 *metrics for 2D: 0 xi_x 1 xi_y  2 eta_x  3 eta_y
 *metrics for 3D:
 */
ops_dat* g_Metrics;
ops_dat* g_NodeType;
ops_dat* g_GeometryProperty;
int g_HaloNum = 0;
ops_halo* g_Halos;
ops_halo_group g_HaloGroups;
int* g_BlockIterRngWhole;
int* g_BlockIterRngJmin;
int* g_BlockIterRngJmax;
int* g_BlockIterRngImin;
int* g_BlockIterRngImax;
int* g_BlockIterRngKmax;
int* g_BlockIterRngKmin;
int* g_BlockIterRngBulk;
int* g_BlockSize;

const int HaloPtNum() { return std::max(SchemeHaloNum(), BoundaryHaloNum()); }
/*!
 * Manually setting up all the variables necessary for the simulation
 * This function has not been updated for 3D problems. 
*/
void DefineVariables() {
    void* temp = NULL;
    g_Block = new ops_block[g_BlockNum];
    g_f = new ops_dat[g_BlockNum];
    g_Bodyforce = new ops_dat[g_BlockNum];
    g_fStage = new ops_dat[g_BlockNum];
    g_feq = new ops_dat[g_BlockNum];
    g_MacroVars = new ops_dat[g_BlockNum];
    g_Tau = new ops_dat[g_BlockNum];
    g_CoordinateXYZ = new ops_dat[g_BlockNum];
    g_BlockIterRngWhole = new int[g_BlockNum * 2 * SPACEDIM];
    g_BlockIterRngJmin = new int[g_BlockNum * 2 * SPACEDIM];
    g_BlockIterRngJmax = new int[g_BlockNum * 2 * SPACEDIM];
    g_BlockIterRngImax = new int[g_BlockNum * 2 * SPACEDIM];
    g_BlockIterRngImin = new int[g_BlockNum * 2 * SPACEDIM];
    g_BlockIterRngBulk = new int[g_BlockNum * 2 * SPACEDIM];
    // if steady flow
    g_MacroVarsCopy = new ops_dat[g_BlockNum];
    g_ResidualErrorHandle = new ops_reduction[MacroVarsNum()];
    g_ResidualError = new Real[2 * MacroVarsNum()];
    // end if steady flow
    int haloDepth = HaloDepth();  // zero for the current case
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
    //g_Metrics= new ops_dat[g_BlockNum];
    // calculate the g_Metrics;
    // if cutting cell method
    g_NodeType = new ops_dat[g_BlockNum];
    g_GeometryProperty = new ops_dat[g_BlockNum];
    for (int blockIndex = 0; blockIndex < g_BlockNum; blockIndex++) {
        std::string label(std::to_string(blockIndex));
        std::string blockName("Block_" + label);
        // The name parameter is not properly typed in the definition of
        // ops_decl_block, so there is a minor warning here.
        g_Block[blockIndex] =
            ops_decl_block(SPACEDIM, (char*)blockName.c_str());
        int size[2] = {BlockSize(blockIndex)[0],
                       BlockSize(blockIndex)[1]};  // size of the dat

        g_BlockIterRngWhole[blockIndex * 2 * SPACEDIM] = 0;
        g_BlockIterRngWhole[blockIndex * 2 * SPACEDIM + 1] = size[0];
        g_BlockIterRngWhole[blockIndex * 2 * SPACEDIM + 2] = 0;
        g_BlockIterRngWhole[blockIndex * 2 * SPACEDIM + 3] = size[1];

        g_BlockIterRngBulk[blockIndex * 2 * SPACEDIM] = 1;
        g_BlockIterRngBulk[blockIndex * 2 * SPACEDIM + 1] = size[0] - 1;
        g_BlockIterRngBulk[blockIndex * 2 * SPACEDIM + 2] = 1;
        g_BlockIterRngBulk[blockIndex * 2 * SPACEDIM + 3] = size[1] - 1;

        g_BlockIterRngJmax[blockIndex * 2 * SPACEDIM] = 0;
        g_BlockIterRngJmax[blockIndex * 2 * SPACEDIM + 1] = size[0];
        g_BlockIterRngJmax[blockIndex * 2 * SPACEDIM + 2] = size[1] - 1;
        g_BlockIterRngJmax[blockIndex * 2 * SPACEDIM + 3] = size[1];

        g_BlockIterRngJmin[blockIndex * 2 * SPACEDIM] = 0;
        g_BlockIterRngJmin[blockIndex * 2 * SPACEDIM + 1] = size[0];
        g_BlockIterRngJmin[blockIndex * 2 * SPACEDIM + 2] = 0;
        g_BlockIterRngJmin[blockIndex * 2 * SPACEDIM + 3] = 1;

        g_BlockIterRngImax[blockIndex * 2 * SPACEDIM] = size[0] - 1;
        g_BlockIterRngImax[blockIndex * 2 * SPACEDIM + 1] = size[0];
        g_BlockIterRngImax[blockIndex * 2 * SPACEDIM + 2] = 0;
        g_BlockIterRngImax[blockIndex * 2 * SPACEDIM + 3] = size[1];

        g_BlockIterRngImin[blockIndex * 2 * SPACEDIM] = 0;
        g_BlockIterRngImin[blockIndex * 2 * SPACEDIM + 1] = 1;
        g_BlockIterRngImin[blockIndex * 2 * SPACEDIM + 2] = 0;
        g_BlockIterRngImin[blockIndex * 2 * SPACEDIM + 3] = size[1];

        std::string dataName("f_");
        dataName += label;
        g_f[blockIndex] =
            ops_decl_dat(g_Block[blockIndex], NUMXI, size, base, d_m, d_p,
                         (Real*)temp, RealC, dataName.c_str());
        dataName = "feq_" + label;
        g_feq[blockIndex] =
            ops_decl_dat(g_Block[blockIndex], NUMXI, size, base, d_m, d_p,
                         (Real*)temp, RealC, dataName.c_str());
        dataName = "fStage_" + label;
        g_fStage[blockIndex] =
            ops_decl_dat(g_Block[blockIndex], NUMXI, size, base, d_m, d_p,
                         (Real*)temp, RealC, dataName.c_str());
        dataName = "Bodyforce_" + label;
        g_Bodyforce[blockIndex] =
            ops_decl_dat(g_Block[blockIndex], NUMXI, size, base, d_m, d_p,
                         (Real*)temp, RealC, dataName.c_str());
        dataName = "MacroVars_" + label;
        g_MacroVars[blockIndex] =
            ops_decl_dat(g_Block[blockIndex], NUMMACROVAR, size, base, d_m, d_p,
                         (Real*)temp, RealC, dataName.c_str());
        dataName = "Tau_" + label;
        g_Tau[blockIndex] =
            ops_decl_dat(g_Block[blockIndex], NUMCOMPONENTS, size, base, d_m,
                         d_p, (Real*)temp, RealC, dataName.c_str());
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
                sizeof(double), "double", MACROVARNAME[localIdx].c_str());
        }
        // end if steady flow
    }
}
/*!
 * setting up all the variables necessary for the simulation from a HDF5 file
 * This function can be used for both 2D and 3D cases
*/
void DefineVariablesFromHDF5() {
    void* temp = NULL;
    g_Block = new ops_block[g_BlockNum];
    g_f = new ops_dat[g_BlockNum];
    g_Bodyforce = new ops_dat[g_BlockNum];
    g_fStage = new ops_dat[g_BlockNum];
    g_feq = new ops_dat[g_BlockNum];
    g_MacroVars = new ops_dat[g_BlockNum];
    g_Tau = new ops_dat[g_BlockNum];
    g_CoordinateXYZ = new ops_dat[g_BlockNum];
    g_BlockIterRngWhole = new int[g_BlockNum * 2 * SPACEDIM];
    g_BlockIterRngJmin = new int[g_BlockNum * 2 * SPACEDIM];
    g_BlockIterRngJmax = new int[g_BlockNum * 2 * SPACEDIM];
    g_BlockIterRngImax = new int[g_BlockNum * 2 * SPACEDIM];
    g_BlockIterRngImin = new int[g_BlockNum * 2 * SPACEDIM];
    if (3 == SPACEDIM) {
        g_BlockIterRngKmax = new int[g_BlockNum * 2 * SPACEDIM];
        g_BlockIterRngKmin = new int[g_BlockNum * 2 * SPACEDIM];
    }
    g_BlockIterRngBulk = new int[g_BlockNum * 2 * SPACEDIM];
    // if steady flow
    g_MacroVarsCopy = new ops_dat[g_BlockNum];
    g_ResidualErrorHandle = new ops_reduction[MacroVarsNum()];
    g_ResidualError = new Real[2 * MacroVarsNum()];
    // end if steady flow
    int haloDepth = HaloDepth();
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
    // g_Metrics= new ops_dat[g_BlockNum];
    // calculate the g_Metrics;
    // if cutting cell method
    g_NodeType = new ops_dat[g_BlockNum];
    g_GeometryProperty = new ops_dat[g_BlockNum];
    for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
        std::string label(std::to_string(blockIndex));
        std::string blockName("Block_" + label);
        std::string fileName = CASENAME + "_" + "Block" + label + ".h5";
        // The name parameter is not properly typed in the definition of
        // ops_decl_block, so there is a minor warning here.
        g_Block[blockIndex] =
            ops_decl_block(SPACEDIM, (char*)blockName.c_str());
        int* size = new int[SPACEDIM];  // size of the dat
        for (int cordIdx = 0; cordIdx < SPACEDIM; cordIdx++) {
            size[cordIdx] = BlockSize(blockIndex)[cordIdx];
        }
        // we assume a problem is at least 2D
        g_BlockIterRngWhole[blockIndex * 2 * SPACEDIM] = 0;
        g_BlockIterRngWhole[blockIndex * 2 * SPACEDIM + 1] = size[0];
        g_BlockIterRngWhole[blockIndex * 2 * SPACEDIM + 2] = 0;
        g_BlockIterRngWhole[blockIndex * 2 * SPACEDIM + 3] = size[1];
        g_BlockIterRngBulk[blockIndex * 2 * SPACEDIM] = 1;
        g_BlockIterRngBulk[blockIndex * 2 * SPACEDIM + 1] = size[0] - 1;
        g_BlockIterRngBulk[blockIndex * 2 * SPACEDIM + 2] = 1;
        g_BlockIterRngBulk[blockIndex * 2 * SPACEDIM + 3] = size[1] - 1;
        g_BlockIterRngJmax[blockIndex * 2 * SPACEDIM] = 0;
        g_BlockIterRngJmax[blockIndex * 2 * SPACEDIM + 1] = size[0];
        g_BlockIterRngJmax[blockIndex * 2 * SPACEDIM + 2] = size[1] - 1;
        g_BlockIterRngJmax[blockIndex * 2 * SPACEDIM + 3] = size[1];
        g_BlockIterRngJmin[blockIndex * 2 * SPACEDIM] = 0;
        g_BlockIterRngJmin[blockIndex * 2 * SPACEDIM + 1] = size[0];
        g_BlockIterRngJmin[blockIndex * 2 * SPACEDIM + 2] = 0;
        g_BlockIterRngJmin[blockIndex * 2 * SPACEDIM + 3] = 1;
        g_BlockIterRngImax[blockIndex * 2 * SPACEDIM] = size[0] - 1;
        g_BlockIterRngImax[blockIndex * 2 * SPACEDIM + 1] = size[0];
        g_BlockIterRngImax[blockIndex * 2 * SPACEDIM + 2] = 0;
        g_BlockIterRngImax[blockIndex * 2 * SPACEDIM + 3] = size[1];
        g_BlockIterRngImin[blockIndex * 2 * SPACEDIM] = 0;
        g_BlockIterRngImin[blockIndex * 2 * SPACEDIM + 1] = 1;
        g_BlockIterRngImin[blockIndex * 2 * SPACEDIM + 2] = 0;
        g_BlockIterRngImin[blockIndex * 2 * SPACEDIM + 3] = size[1];
        if (3 == SPACEDIM) {
            g_BlockIterRngWhole[blockIndex * 2 * SPACEDIM + 4] = 0;
            g_BlockIterRngWhole[blockIndex * 2 * SPACEDIM + 5] = size[2];
            g_BlockIterRngBulk[blockIndex * 2 * SPACEDIM + 4] = 1;
            g_BlockIterRngBulk[blockIndex * 2 * SPACEDIM + 5] = size[2] - 1;
            g_BlockIterRngJmax[blockIndex * 2 * SPACEDIM + 4] = 0;
            g_BlockIterRngJmax[blockIndex * 2 * SPACEDIM + 5] = size[2];
            g_BlockIterRngJmin[blockIndex * 2 * SPACEDIM + 4] = 0;
            g_BlockIterRngJmin[blockIndex * 2 * SPACEDIM + 5] = size[2];
            g_BlockIterRngImax[blockIndex * 2 * SPACEDIM + 4] = 0;
            g_BlockIterRngImax[blockIndex * 2 * SPACEDIM + 5] = size[2];
            g_BlockIterRngImin[blockIndex * 2 * SPACEDIM + 4] = 0;
            g_BlockIterRngImin[blockIndex * 2 * SPACEDIM + 5] = size[2];
            g_BlockIterRngKmax[blockIndex * 2 * SPACEDIM] = 0;
            g_BlockIterRngKmax[blockIndex * 2 * SPACEDIM + 1] = size[0];
            g_BlockIterRngKmax[blockIndex * 2 * SPACEDIM + 2] = 0;
            g_BlockIterRngKmax[blockIndex * 2 * SPACEDIM + 3] = size[1];
            g_BlockIterRngKmax[blockIndex * 2 * SPACEDIM + 4] = size[2] - 1;
            g_BlockIterRngKmax[blockIndex * 2 * SPACEDIM + 5] = size[2];
            g_BlockIterRngKmin[blockIndex * 2 * SPACEDIM] = 0;
            g_BlockIterRngKmin[blockIndex * 2 * SPACEDIM + 1] = size[0];
            g_BlockIterRngKmin[blockIndex * 2 * SPACEDIM + 2] = 0;
            g_BlockIterRngKmin[blockIndex * 2 * SPACEDIM + 3] = size[1];
            g_BlockIterRngKmin[blockIndex * 2 * SPACEDIM + 4] = 0;
            g_BlockIterRngKmin[blockIndex * 2 * SPACEDIM + 5] = 1;
        }
        std::string dataName("f_");
        dataName += label;
        g_f[blockIndex] =
            ops_decl_dat(g_Block[blockIndex], NUMXI, size, base, d_m, d_p,
                         (Real*)temp, RealC, dataName.c_str());
        dataName = "feq_" + label;
        g_feq[blockIndex] =
            ops_decl_dat(g_Block[blockIndex], NUMXI, size, base, d_m, d_p,
                         (Real*)temp, RealC, dataName.c_str());
        dataName = "fStage_" + label;
        g_fStage[blockIndex] =
            ops_decl_dat(g_Block[blockIndex], NUMXI, size, base, d_m, d_p,
                         (Real*)temp, RealC, dataName.c_str());
        dataName = "Bodyforce_" + label;
        g_Bodyforce[blockIndex] =
            ops_decl_dat(g_Block[blockIndex], NUMXI, size, base, d_m, d_p,
                         (Real*)temp, RealC, dataName.c_str());
        dataName = "MacroVars_" + label;
        g_MacroVars[blockIndex] =
            ops_decl_dat_hdf5(g_Block[blockIndex], NUMMACROVAR, "double",
                              dataName.c_str(), fileName.c_str());
        dataName = "Tau_" + label;
        g_Tau[blockIndex] =
            ops_decl_dat(g_Block[blockIndex], NUMCOMPONENTS, size, base, d_m,
                         d_p, (Real*)temp, RealC, dataName.c_str());
        dataName = "Nodetype_" + label;
        // problem specific -- cut cell method
        g_NodeType[blockIndex] = ops_decl_dat_hdf5(
            g_Block[blockIndex], 1, "int", dataName.c_str(), fileName.c_str());
        dataName = "GeometryProperty_" + label;
        g_GeometryProperty[blockIndex] = ops_decl_dat_hdf5(
            g_Block[blockIndex], 1, "int", dataName.c_str(), fileName.c_str());
        dataName = "CoordinateXYZ_" + label;
        g_CoordinateXYZ[blockIndex] =
            ops_decl_dat_hdf5(g_Block[blockIndex], SPACEDIM, RealC,
                              dataName.c_str(), fileName.c_str());
        // if steady flow
        // in the future, we may consider to add an option for the "if"
        dataName = "MacroVars_Copy" + label;
        g_MacroVarsCopy[blockIndex] =
            ops_decl_dat(g_Block[blockIndex], NUMMACROVAR, size, base, d_m, d_p,
                         (Real*)temp, RealC, dataName.c_str());
        for (int localIdx = 0; localIdx < MacroVarsNum(); localIdx++) {
            g_ResidualErrorHandle[localIdx] = ops_decl_reduction_handle(
                // this is double
                sizeof(double), "double", MACROVARNAME[localIdx].c_str());
        }
        // end if steady flow
        delete[] size;
    }
    delete[] d_p;
    delete[] d_m;
    delete[] base;
}
/*!
 * Manually define the halo relation between blocks.
 * When using
 */
void DefineHaloTransfer(){
    /*! @brief Defining the halo relationship
     *  @details Currently we need to manually define them,
     *  will be modified to read CGNF format in the future
     **/
    //     g_HaloNum = 8;
    //     g_Halos = new ops_halo[g_HaloNum];
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
    //         g_Halos[0] = ops_decl_halo ( g_f[0], g_f[0], halo_iter,
    //         base_from,
    //                                      base_to, dir, dir );
    //         base_from[0] = nx - 1; // need to be changed
    //         base_to[0] = d_m[0];
    //         g_Halos[1] = ops_decl_halo ( g_f[0], g_f[0], halo_iter,
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
    //         g_Halos[2] = ops_decl_halo ( g_f[0], g_f[0], halo_iter,
    //         base_from,
    //                                      base_to, dir, dir );
    //         base_from[1] = ny - 1; //need to be changed
    //         base_to[1] = d_m[1];
    //         g_Halos[3] = ops_decl_halo ( g_f[0], g_f[0], halo_iter,
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
    //         g_Halos[4] = ops_decl_halo ( g_f[0], g_f[0], halo_iter,
    //         base_from,
    //                                      base_to, dir, dir );
    //         base_from[0] = nx - 1; // need to be changed
    //         base_from[1] = ny - 1;
    //         base_to[0] = d_m[0];
    //         base_to[1] = d_m[1];
    //         g_Halos[5] = ops_decl_halo ( g_f[0], g_f[0], halo_iter,
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
    //         g_Halos[6] = ops_decl_halo ( g_f[0], g_f[0], halo_iter,
    //         base_from,
    //                                      base_to, dir, dir );
    //         base_from[0] = nx - 1; // need to be changed
    //         base_from[1] = 0;
    //         base_to[0] = d_m[0];
    //         base_to[1] = ny;
    //         g_Halos[7] = ops_decl_halo ( g_f[0], g_f[0], halo_iter,
    //         base_from,
    //                                      base_to, dir, dir );
    //     }
    //     g_HaloGroups = ops_decl_halo_group ( g_HaloNum, g_Halos );

    g_HaloNum = 2;
    g_Halos = new ops_halo[g_HaloNum];
    int haloDepth = HaloDepth();
    // max halo depths for the dat in the positive direction
    int d_p[2] = {haloDepth, haloDepth};
    // max halo depths for the dat in the negative direction
    int d_m[2] = {-haloDepth, -haloDepth};
    // The domain size in the Block 0
    int nx = BlockSize(0)[0];
    int ny = BlockSize(0)[1];
    int dir[] = {1, 2};
    {
        int halo_iter[] = {nx+d_p[0]-d_m[0], 1};
        int base_from[] = {d_m[0], 0};
        int base_to[] = {d_m[0], ny};
        g_Halos[0] = ops_decl_halo(g_f[0], g_f[0], halo_iter, base_from,
                                   base_to, dir, dir);
        base_from[1] = ny - 1;  // need to be changed
        base_to[1] = d_m[1];
        g_Halos[1] = ops_decl_halo(g_f[0], g_f[0], halo_iter, base_from,
                                   base_to, dir, dir);
    }

    g_HaloGroups = ops_decl_halo_group(g_HaloNum, g_Halos);
}
/*
 * We need a name to specify which file to input
 * To be decided: a single filename or an array of filenames
 */

void WriteFlowfieldToHdf5(const long timeStep) {
    for (int blockIndex = 0; blockIndex < g_BlockNum; blockIndex++) {
        std::string blockName("Block_");
        std::string label(std::to_string(blockIndex));
        std::string time(std::to_string(timeStep));
        blockName += (label + "_" + time);
        std::string fileName =CASENAME+"_"+blockName + ".h5";
        ops_fetch_block_hdf5_file(g_Block[blockIndex], fileName.c_str());
        ops_fetch_dat_hdf5_file(g_MacroVars[blockIndex], fileName.c_str());
        ops_fetch_dat_hdf5_file(g_Tau[blockIndex], fileName.c_str());
        ops_fetch_dat_hdf5_file(g_CoordinateXYZ[blockIndex],fileName.c_str());
    }
}

void WriteDistributionsToHdf5(const long timeStep) {
    for (int blockIndex = 0; blockIndex < g_BlockNum; blockIndex++) {
        std::string blockName("Block_");
        std::string label(std::to_string(blockIndex));
        std::string time(std::to_string(timeStep));
        blockName += (label + "_" + time);
        std::string fileName =CASENAME+"_"+blockName + ".h5";
        ops_fetch_block_hdf5_file(g_Block[blockIndex], fileName.c_str());
        ops_fetch_dat_hdf5_file(g_f[blockIndex], fileName.c_str());
        ops_fetch_dat_hdf5_file(g_feq[blockIndex], fileName.c_str());
        ops_fetch_dat_hdf5_file(g_fStage[blockIndex], fileName.c_str());
        ops_fetch_dat_hdf5_file(g_Bodyforce[blockIndex], fileName.c_str());
    }
}


void WriteNodePropertyToHdf5(const long timeStep) {
    for (int blockIndex = 0; blockIndex < g_BlockNum; blockIndex++) {
        std::string blockName("Block_");
        std::string label(std::to_string(blockIndex));
        std::string time(std::to_string(timeStep));
        blockName += (label + "_" + time);
        std::string fileName =CASENAME+"_"+blockName + ".h5";
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
void SetupFlowfieldfromHdf5() {
    CASENAME = "SWECircularDamBreak";  // Input parameter
    SPACEDIM = 2;
    g_BlockNum = 1;  // Input parameter
    g_BlockSize = new int[g_BlockNum * SPACEDIM];
    g_BlockSize[0] = 401;  // Input parameters
    g_BlockSize[1] = 401;  // Input parameters
    // g_BlockSize[2] = 11;  // Input parameters

    KN = new Real[ComponentNum()];
    KN[0] = 0.001;  // Input parameters
    // All above parameters should be written down by the python script
    Real minDx{2. / 400};  // Input parameters at this moment
    Real minDy{2. / 400};  // Input parameters at this moment
    // g_dt = 0.01 * fmin(minDx, minDy) / MaximumSpeed();  // finite difference
    // scheme
    g_dt = 0.0001414;
    //g_dt = fmin(fmin(minDx, minDy) / MaximumSpeed(),
    //              0.5 * KN[0]);  // finite difference scheme
    // g_dt = minDx / SoundSpeed();  // stream-collision
    g_HaloDepth = HaloPtNum();
    DefineVariablesFromHDF5();
    // DefineHaloTransfer();
    // above calls must be before the ops_partition call
    ops_partition((char*)"LBM");
}

void DestroyFlowfield() {
    delete[] g_f;
    delete[] g_fStage;
    delete[] g_feq;
    delete[] g_Bodyforce;
    delete[] g_Block;
    delete[] g_MacroVars;
    delete[] g_Tau;
    delete[] KN;
    delete[] g_CoordinateXYZ;
    if (g_HaloNum > 0) delete[] g_Halos;
    delete[] g_NodeType;
    delete[] g_GeometryProperty;
    delete[] g_BlockIterRngWhole;
    delete[] g_BlockIterRngBulk;
    delete[] g_BlockIterRngImax;
    delete[] g_BlockIterRngImin;
    delete[] g_BlockIterRngJmax;
    delete[] g_BlockIterRngJmin;
    delete[] g_BlockSize;
    // if steady flow
    delete[] g_MacroVarsCopy;
    delete[] g_ResidualErrorHandle;
    delete[] g_ResidualError;
	if (3 == SPACEDIM) {
		delete[] g_BlockIterRngKmax;
		delete[] g_BlockIterRngKmin;
	}
    // end if steady flow
    // delete[] halos;
}

Real TotalMeshSize(){
    Real size=1;
    for (int blockIdx=0;blockIdx<BlockNum();blockIdx++)
    {
        for (int cordIdx=0;cordIdx<SPACEDIM;cordIdx++ ){
            size*=BlockSize(blockIdx)[cordIdx];
        }
    }
    return size;
}
