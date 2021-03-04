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
#include <type_traits>
#include "ops_seq_v2.h"
#include "block.h"
#include "field.h"
#include "model.h"

std::string CASENAME;
bool TRANSIENT{false};
/*!
 * SPACEDIM=2 for 2D 3 for three 3D
 */
#ifdef OPS_3D
int SPACEDIM{3};
#endif  //  OPS_3D
#ifdef OPS_2D
int SPACEDIM{2};
#endif // ops_2D
BlockGroup BLOCKS;
RealField f{"f"};
RealField fStage{"fStage"};
RealField MacroVars{"MacroVars"};
RealField MacroVarsCopy{"MacroVars_Copy"};
Real* g_ResidualError{nullptr};
ops_reduction* g_ResidualErrorHandle{nullptr};
RealField MacroBodyforce{"MacroBodyForce"};

const BlockGroup& g_Block() { return BLOCKS; };
RealField& g_f() { return f; };
RealField& g_fStage() { return fStage; };
RealField& g_MacroVars() { return MacroVars; };
RealField& g_MacroVarsCopy() { return MacroVarsCopy; };
RealField& g_MacroBodyforce() { return MacroBodyforce; };
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
RealField CoordinateXYZ{"CoordinateXYZ"};
const RealField& g_CoordinateXYZ() { return CoordinateXYZ; };
std::vector<std::vector<std::vector<Real>>> COORDINATES;

IntField NodeType{"NodeType"};
IntField GeometryProperty{"GeometryProperty"};
IntField& g_NodeType() { return NodeType; };
IntField& g_GeometryProperty() { return GeometryProperty; };

/*!
 * Formal collection of halo relations required by the OPS library
 */
std::map<std::string,ops_halo_group> HALOGROUPS;

void DefineCase(const std::string& caseName, const int spaceDim,
                const bool transient) {
    if (SPACEDIM != spaceDim) {
        ops_printf("Error! The SPACEDIM here is inconsistent with the\n");
        assert(SPACEDIM == spaceDim);
    }
    SetCaseName(caseName);
    ops_decl_const("SPACEDIM", 1, "int", &SPACEDIM);
    TRANSIENT = transient;
}

bool IsTransient() { return TRANSIENT; }

//This version is for the distribution function
//Other version to follow.
//TODO this function only suitable for single halo
#ifdef OPS_3D
void DefinePeriodicHaloPair3D(const std::map<int, std::string>& haloPair) {
    if (BLOCKS.size() > 1) {
        ops_printf(
            "Periodic boundary conditions are only valid for single-block "
            "applications!");
        assert(BLOCKS.size() == 1);
    }
    const SizeType blockId{BLOCKS.begin()->first};
    DefinePeriodicHaloPair3D(haloPair, f[blockId], f.HaloDepth());
}

void DefinePeriodicHaloPair3D(const std::map<int, std::string>& haloPair,
                              RealField& data) {
    if (BLOCKS.size() > 1) {
        ops_printf(
            "Periodic boundary conditions are only valid for single-block "
            "applications!");
        assert(BLOCKS.size() == 1);
    }
    const SizeType blockId{BLOCKS.begin()->first};
#ifdef OPS_3D
    DefinePeriodicHaloPair3D(haloPair, data[blockId], data.HaloDepth());
#endif  // OPS_3D
}

void DefinePeriodicHaloPair3D(const std::map<int, std::string>& haloPair,
                              IntField& data) {
    if (BLOCKS.size() > 1) {
        ops_printf(
            "Periodic boundary conditions are only valid for single-block "
            "applications!");
        assert(BLOCKS.size() == 1);
    }
    const SizeType blockId{BLOCKS.begin()->first};
#ifdef OPS_3D
    DefinePeriodicHaloPair3D(haloPair, data[blockId], data.HaloDepth());
#endif  // OPS_3D
}

void DefinePeriodicHaloPair3D(const std::map<int,std::string>& haloPair, ops_dat dat,
                              const int haloDepth) {

    if (dat == nullptr) {
        ops_printf("Data must be allocated before defining halo relations!");
        assert(dat == nullptr);
    }
    // max halo depths for the dat in the positive direction
    int d_p[3] = {haloDepth, haloDepth, haloDepth};
    // max halo depths for the dat in the negative direction
    int d_m[3] = {-haloDepth, -haloDepth, -haloDepth};
    // The domain size in the Block 0
    const Block& block{BLOCKS.begin()->second};
    int nx{(int)block.Size().at(0)};
    int ny{(int)block.Size().at(1)};
    int nz{(int)block.Size().at(2)};
    int dir[] = {1, 2, 3};
    for (auto& pair : haloPair) {
        if (0 == pair.first) {
            // left and right pair
            int halo_iter[] = {haloDepth, ny + d_p[1] - d_m[1],
                               nz + d_p[2] - d_m[2]};
            int base_from[] = {0, d_m[1], d_m[2]};
            int base_to[] = {nx, d_m[1], d_m[2]};
            ops_halo leftToRight = ops_decl_halo(dat, dat, halo_iter, base_from,
                                                 base_to, dir, dir);
            base_from[0] = nx + d_m[0];
            base_to[0] = d_m[0];
            ops_halo rightToLeft = ops_decl_halo(dat, dat, halo_iter, base_from,
                                                 base_to, dir, dir);
            ops_halo group[]{leftToRight, rightToLeft};
            HALOGROUPS.emplace(pair.second,ops_decl_halo_group(2, group));
        }

        if (1 == pair.first) {
            // top and bottom pair
            int halo_iter[] = {nx + d_p[0] - d_m[0], haloDepth,
                               nz + d_p[2] - d_m[2]};
            int base_from[] = {d_m[0], 0, d_m[2]};
            int base_to[] = {d_m[0], ny, d_m[2]};
            ops_halo botToTop = ops_decl_halo(dat, dat, halo_iter, base_from,
                                              base_to, dir, dir);
            base_from[1] = ny + d_m[1];
            base_to[1] = d_m[1];
            ops_halo topToBot = ops_decl_halo(dat, dat, halo_iter, base_from,
                                              base_to, dir, dir);
            ops_halo group[]{botToTop, topToBot};
            HALOGROUPS.emplace(pair.second,ops_decl_halo_group(2, group));
        }

        if (2 == pair.first) {
            // front and back pair
            int halo_iter[] = {nx + d_p[0] - d_m[0], ny + d_p[1] - d_m[1],
                               haloDepth};
            int base_from[] = {d_m[0], d_m[1], 0};
            int base_to[] = {d_m[0], d_m[1], nz};
            ops_halo backToFront = ops_decl_halo(dat, dat, halo_iter, base_from,
                                                 base_to, dir, dir);
            base_from[2] = nz + d_m[2];
            base_to[2] = d_m[2];
            ops_halo frontToBack = ops_decl_halo(dat, dat, halo_iter, base_from,
                                                 base_to, dir, dir);
            ops_halo group[]{backToFront, frontToBack};
            HALOGROUPS.emplace(pair.second, ops_decl_halo_group(2, group));
        }
    }
}

#endif //OPS_3D

void Partition() {
    ops_partition((char*)"LBM Solver");
    PrepareFlowField();
}

/*
 * We need a name to specify which file to input
 * To be decided: a single filename or an array of filenames
 */

void WriteFlowfieldToHdf5(const SizeType timeStep) {
    MacroVars.WriteToHDF5(CASENAME, timeStep);
    CoordinateXYZ.WriteToHDF5(CASENAME, timeStep);
    MacroBodyforce.WriteToHDF5(CASENAME, timeStep);
}

void WriteDistributionsToHdf5(const SizeType timeStep) {
    f.WriteToHDF5(CASENAME, timeStep);

}

void WriteNodePropertyToHdf5(const SizeType timeStep) {
    GeometryProperty.WriteToHDF5(CASENAME, timeStep);
    NodeType.WriteToHDF5(CASENAME, timeStep);
}

const std::string& CaseName() { return CASENAME; }
void SetCaseName(const std::string& caseName) { CASENAME = caseName; }
void setCaseName(const char* caseName) {
    std::string tmp(caseName);
    CASENAME = tmp;
}
int SpaceDim() { return SPACEDIM; }

void DestroyFlowfield() {

    FreeArrayMemory(TAUREF);
    // if steady flow
    // FreeArrayMemory(g_MacroVarsCopy);
    FreeArrayMemory(g_ResidualErrorHandle);
    FreeArrayMemory(g_ResidualError);
    // end if steady flow
    // delete[] halos;
}





Real TotalMeshSize() {

}

Real TimeStep() { return DT; }
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

const std::map<std::string,ops_halo_group>& HaloGroups() { return HALOGROUPS; }

void TransferHalos() {
    if (HALOGROUPS.size() > 0) {
        for (auto haloGroup : HALOGROUPS){
             ops_halo_transfer(haloGroup.second);
        }
    }
}

void TransferHalos(const std::string key) {
    ops_halo_transfer(HALOGROUPS[key]);
}

void TransferHalos(const std::vector<std::string> keys) {
    if (keys.size() > 0) {
        for (auto key : keys) {
            ops_halo_transfer(HALOGROUPS[key]);
        }
    }
}

void  SetBlockGeometryProperty(Block& block) {
    int geometryProperty = (int)VG_Fluid;
    // int* iterRange = BlockIterRng(blockIndex, IterRngBulk());
    std::vector<int> iterRange;
    iterRange.assign(block.BulkRange().begin(), block.BulkRange().end());
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty", block.Get(),
                 SPACEDIM, iterRange.data(),
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(GeometryProperty[block.ID()], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    // specify halo points
    geometryProperty = VG_ImmersedSolid;
    iterRange.assign(block.JminRange().begin(), block.JminRange().end());
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
                 block.Get(), SPACEDIM, haloIterRng,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(GeometryProperty[block.ID()], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    iterRange.assign(block.JmaxRange().begin(), block.JmaxRange().end());
    haloIterRng[0] = iterRange[0] - 1;
    haloIterRng[1] = iterRange[1] + 1;
    haloIterRng[2] = iterRange[2] + 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                 block.Get(), SPACEDIM, haloIterRng,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(GeometryProperty[block.ID()], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    iterRange.assign(block.IminRange().begin(), block.IminRange().end());
    haloIterRng[0] = iterRange[0] - 1;
    haloIterRng[1] = iterRange[1] - 1;
    haloIterRng[2] = iterRange[2] - 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                block.Get(), SPACEDIM, haloIterRng,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(GeometryProperty[block.ID()], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    iterRange.assign(block.ImaxRange().begin(), block.ImaxRange().end());
    haloIterRng[0] = iterRange[0] + 1;
    haloIterRng[1] = iterRange[1] + 1;
    haloIterRng[2] = iterRange[2] - 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                block.Get(), SPACEDIM, haloIterRng,
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(GeometryProperty[block.ID()], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    if (3 == SPACEDIM) {
        iterRange.assign(block.KminRange().begin(), block.KminRange().end());
        haloIterRng[0] = iterRange[0] - 1;
        haloIterRng[1] = iterRange[1] + 1;
        haloIterRng[2] = iterRange[2] - 1;
        haloIterRng[3] = iterRange[3] + 1;
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] - 1;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, haloIterRng,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        iterRange.assign(block.KmaxRange().begin(), block.KmaxRange().end());
        haloIterRng[0] = iterRange[0] - 1;
        haloIterRng[1] = iterRange[1] + 1;
        haloIterRng[2] = iterRange[2] - 1;
        haloIterRng[3] = iterRange[3] + 1;
        haloIterRng[4] = iterRange[4] + 1;
        haloIterRng[5] = iterRange[5] + 1;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, haloIterRng,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
    }

    // specify domain
    geometryProperty = VG_JP;
    iterRange.assign(block.JminRange().begin(), block.JminRange().end());
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                block.Get(), SPACEDIM, iterRange.data(),
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(GeometryProperty[block.ID()], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    geometryProperty = VG_JM;
    iterRange.assign(block.JmaxRange().begin(), block.JmaxRange().end());
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                block.Get(), SPACEDIM, iterRange.data(),
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(GeometryProperty[block.ID()], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    geometryProperty = VG_IP;
    iterRange.assign(block.IminRange().begin(), block.IminRange().end());
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                block.Get(), SPACEDIM, iterRange.data(),
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(GeometryProperty[block.ID()], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    geometryProperty = VG_IM;
    iterRange.assign(block.ImaxRange().begin(), block.ImaxRange().end());
    ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                block.Get(), SPACEDIM, iterRange.data(),
                 ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                 ops_arg_dat(GeometryProperty[block.ID()], 1, LOCALSTENCIL,
                             "int", OPS_WRITE));
    if (3 == SPACEDIM) {
        geometryProperty = VG_KP;
        iterRange.assign(block.KminRange().begin(), block.KminRange().end());
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, iterRange.data(),
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        geometryProperty = VG_KM;
        iterRange.assign(block.KmaxRange().begin(), block.KmaxRange().end());
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, iterRange.data(),
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
    }

    const int nx{(int)block.Size().at(0)};
    const int ny{(int)block.Size().at(1)};
    // 2D Domain corner points four types
    if (2 == SPACEDIM) {
        int iminjmin[]{0, 1, 0, 1};
        geometryProperty = VG_IPJP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, iminjmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminjmax[] = {0, 1, ny - 1, ny};
        geometryProperty = VG_IPJM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, iminjmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmax[] = {nx - 1, nx, ny - 1, ny};
        geometryProperty = VG_IMJM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, imaxjmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmin[] = {nx - 1, nx, 0, 1};
        geometryProperty = VG_IMJP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, imaxjmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
    }

    if (3 == SPACEDIM) {
        const int nz{(int)block.Size().at(2)};
        // 3D Domain edges 12 types
        int iminjmin[]{0, 1, 0, 1, 0, nz};
        geometryProperty = VG_IPJP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, iminjmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminjmax[]{0, 1, ny - 1, ny, 0, nz};
        geometryProperty = VG_IPJM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, iminjmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmax[]{nx - 1, nx, ny - 1, ny, 0, nz};
        geometryProperty = VG_IMJM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, imaxjmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmin[]{nx - 1, nx, 0, 1, 0, nz};
        geometryProperty = VG_IMJP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, imaxjmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));

        int iminkmin[]{0, 1, 0, ny, 0, 1};
        geometryProperty = VG_IPKP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, iminkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminkmax[]{0, 1, 0, ny, nz - 1, nz};
        geometryProperty = VG_IPKM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, iminkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxkmax[]{nx - 1, nx, 0, ny, nz - 1, nz};
        geometryProperty = VG_IMKM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, imaxkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxkmin[]{nx - 1, nx, 0, ny, 0, 1};
        geometryProperty = VG_IMKP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, imaxkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));

        int jminkmin[]{0, nx, 0, 1, 0, 1};
        geometryProperty = VG_JPKP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, jminkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int jminkmax[]{0, nx, 0, 1, nz - 1, nz};
        geometryProperty = VG_JPKM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, jminkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int jmaxkmax[]{0, nx, ny - 1, ny, nz - 1, nz};
        geometryProperty = VG_JMKM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, jmaxkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int jmaxkmin[]{0, nx, ny - 1, ny, 0, 1};
        geometryProperty = VG_JMKP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, jmaxkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));

        // 3D domain corners 8 types
        int iminjminkmin[]{0, 1, 0, 1, 0, 1};
        geometryProperty = VG_IPJPKP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, iminjminkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminjminkmax[]{0, 1, 0, 1, nz - 1, nz};
        geometryProperty = VG_IPJPKM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, iminjminkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminjmaxkmin[]{0, 1, ny - 1, ny, 0, 1};
        geometryProperty = VG_IPJMKP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, iminjmaxkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int iminjmaxkmax[]{0, 1, ny - 1, ny, nz - 1, nz};
        geometryProperty = VG_IPJMKM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, iminjmaxkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjminkmin[]{nx - 1, nx, 0, 1, 0, 1};
        geometryProperty = VG_IMJPKP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, imaxjminkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjminkmax[]{nx - 1, nx, 0, 1, nz - 1, nz};
        geometryProperty = VG_IMJPKM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, imaxjminkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmaxkmin[]{nx - 1, nx, ny - 1, ny, 0, 1};
        geometryProperty = VG_IMJMKP_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, imaxjmaxkmin,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
        int imaxjmaxkmax[]{nx - 1, nx, ny - 1, ny, nz - 1, nz};
        geometryProperty = VG_IMJMKM_I;
        ops_par_loop(KerSetGeometryProperty, "KerSetGeometryProperty",
                    block.Get(), SPACEDIM, imaxjmaxkmax,
                     ops_arg_gbl(&geometryProperty, 1, "int", OPS_READ),
                     ops_arg_dat(GeometryProperty[block.ID()], 1,
                                 LOCALSTENCIL, "int", OPS_WRITE));
    }
}

void SetBoundaryNodeType(){
    for (auto boundary : BlockBoundaries()){
    const int boundaryType{(int)boundary.boundaryType};
    const Block& block{BLOCKS.at(boundary.blockIndex)};
    std::vector<int> iterRange{
        BoundarySurfaceRange(block, boundary.boundarySurface)};
    const SizeType compoId{boundary.componentID};
    // Specify general boundary type
    ops_par_loop(KerSetNodeType, "KerSetNodeType", block.Get(),
                 SPACEDIM, iterRange.data(),
                 ops_arg_gbl(&boundaryType, 1, "int", OPS_READ),
                 ops_arg_dat(NodeType[block.ID()], NUMCOMPONENTS,
                             LOCALSTENCIL, "int", OPS_WRITE),
                 ops_arg_gbl(&compoId, 1, "int", OPS_READ));
    }
}
void SetBulkandHaloNodesType(Block& block , int compoId) {
    const int fluidType{(int)VertexType::Fluid};
    const int immersedSolidType{(int)VertexType::ImmersedSolid};

    std::vector<int> iterRange;
    iterRange.assign(block.BulkRange().begin(), block.BulkRange().end());
    ops_par_loop(KerSetNodeType, "KerSetNodeType",block.Get(),
                 SPACEDIM, iterRange.data(),
                 ops_arg_gbl(&fluidType, 1, "int", OPS_READ),
                 ops_arg_dat(NodeType[block.ID()], NUMCOMPONENTS,
                             LOCALSTENCIL, "int", OPS_WRITE),
                 ops_arg_gbl(&compoId, 1, "int", OPS_READ));

    iterRange.assign(block.JminRange().begin(), block.JminRange().end());
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
    ops_par_loop(KerSetNodeType, "KerSetNodeType",block.Get(),
                 SPACEDIM, haloIterRng,
                 ops_arg_gbl(&immersedSolidType, 1, "int", OPS_READ),
                 ops_arg_dat(NodeType[block.ID()], NUMCOMPONENTS,
                             LOCALSTENCIL, "int", OPS_WRITE),
                 ops_arg_gbl(&compoId, 1, "int", OPS_READ));

    iterRange.assign(block.JmaxRange().begin(), block.JmaxRange().end());
    haloIterRng[0] = iterRange[0] - 1;
    haloIterRng[1] = iterRange[1] + 1;
    haloIterRng[2] = iterRange[2] + 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    // Specify halo points
    ops_par_loop(KerSetNodeType, "KerSetNodeType",block.Get(),
                 SPACEDIM, haloIterRng,
                 ops_arg_gbl(&immersedSolidType, 1, "int", OPS_READ),
                 ops_arg_dat(NodeType[block.ID()], NUMCOMPONENTS,
                             LOCALSTENCIL, "int", OPS_WRITE),
                 ops_arg_gbl(&compoId, 1, "int", OPS_READ));

    iterRange.assign(block.IminRange().begin(), block.IminRange().end());
    haloIterRng[0] = iterRange[0] - 1;
    haloIterRng[1] = iterRange[1] - 1;
    haloIterRng[2] = iterRange[2] - 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    // Specify halo points
    ops_par_loop(KerSetNodeType, "KerSetNodeType",block.Get(),
                 SPACEDIM, haloIterRng,
                 ops_arg_gbl(&immersedSolidType, 1, "int", OPS_READ),
                 ops_arg_dat(NodeType[block.ID()], NUMCOMPONENTS,
                             LOCALSTENCIL, "int", OPS_WRITE),
                 ops_arg_gbl(&compoId, 1, "int", OPS_READ));

    iterRange.assign(block.ImaxRange().begin(), block.ImaxRange().end());
    haloIterRng[0] = iterRange[0] + 1;
    haloIterRng[1] = iterRange[1] + 1;
    haloIterRng[2] = iterRange[2] - 1;
    haloIterRng[3] = iterRange[3] + 1;
    if (3 == SPACEDIM) {
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] + 1;
    }
    // Specify halo points
    ops_par_loop(KerSetNodeType, "KerSetNodeType",block.Get(),
                 SPACEDIM, haloIterRng,
                 ops_arg_gbl(&immersedSolidType, 1, "int", OPS_READ),
                 ops_arg_dat(NodeType[block.ID()], NUMCOMPONENTS,
                             LOCALSTENCIL, "int", OPS_WRITE),
                 ops_arg_gbl(&compoId, 1, "int", OPS_READ));

    if (3 == SPACEDIM) {
        iterRange.assign(block.KminRange().begin(), block.KminRange().end());
        haloIterRng[0] = iterRange[0] - 1;
        haloIterRng[1] = iterRange[1] + 1;
        haloIterRng[2] = iterRange[2] - 1;
        haloIterRng[3] = iterRange[3] + 1;
        haloIterRng[4] = iterRange[4] - 1;
        haloIterRng[5] = iterRange[5] - 1;
        // Specify halo points
        ops_par_loop(KerSetNodeType, "KerSetNodeType",block.Get(),
                     SPACEDIM, haloIterRng,
                     ops_arg_gbl(&immersedSolidType, 1, "int", OPS_READ),
                     ops_arg_dat(NodeType[block.ID()], NUMCOMPONENTS,
                                 LOCALSTENCIL, "int", OPS_WRITE),
                     ops_arg_gbl(&compoId, 1, "int", OPS_READ));

        iterRange.assign(block.KmaxRange().begin(), block.KmaxRange().end());
        haloIterRng[0] = iterRange[0] - 1;
        haloIterRng[1] = iterRange[1] + 1;
        haloIterRng[2] = iterRange[2] - 1;
        haloIterRng[3] = iterRange[3] + 1;
        haloIterRng[4] = iterRange[4] + 1;
        haloIterRng[5] = iterRange[5] + 1;
        // Specify halo points
        ops_par_loop(KerSetNodeType, "KerSetNodeType",block.Get(),
                     SPACEDIM, haloIterRng,
                     ops_arg_gbl(&immersedSolidType, 1, "int", OPS_READ),
                     ops_arg_dat(NodeType[block.ID()], NUMCOMPONENTS,
                                 LOCALSTENCIL, "int", OPS_WRITE),
                     ops_arg_gbl(&compoId, 1, "int", OPS_READ));
    }
    FreeArrayMemory(haloIterRng);
}

// void DefineBlocks(const SizeType blockNum,
//                   const std::vector<SizeType>& blockSize, const Real meshSize,
//                   const std::vector<Real>& startPos) {
//     SetBlockNum(blockNum);
//     SetBlockSize(blockSize);
//     COORDINATES.resize(blockNum);
//     SizeType numBlockStartPos;
//     numBlockStartPos = startPos.size();
//     if (numBlockStartPos == blockNum * SPACEDIM) {
//         for (int blockIndex = 0; blockIndex < blockNum; blockIndex++) {
//             COORDINATES.at(blockIndex).resize(SPACEDIM);
//             for (int coordIndex = 0; coordIndex < SPACEDIM; coordIndex++) {
//                 int numOfGridPoints =
//                     BlockSize(blockIndex)[SPACEDIM * blockIndex + coordIndex];
//                 COORDINATES.at(blockIndex)
//                     .at(coordIndex)
//                     .resize(numOfGridPoints);
//                 for (int nodeIndex = 0; nodeIndex < numOfGridPoints;
//                      nodeIndex++) {
//                     COORDINATES.at(blockIndex).at(coordIndex).at(nodeIndex) =
//                         startPos.at(coordIndex) + nodeIndex * meshSize;
//                 }
//             }
//         }
//     } else {
//         ops_printf(
//             "Error! Expected %i coordinates of three starting points %i, but "
//             "received only =%i \n",
//             SPACEDIM * blockNum, numBlockStartPos);
//         assert(numBlockStartPos == blockNum * SPACEDIM);
//     }

//     g_Block = new ops_block[BLOCKNUM];
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
//     for (int blockIndex = 0; blockIndex < BLOCKNUM; blockIndex++) {
//         std::string label(std::to_string(blockIndex));
//         std::string blockName("Block_" + label);
//         // The name parameter is not properly typed in the definition of
//         // ops_decl_block, so there is a minor warning here.
//         int* size = new int[SPACEDIM];  // size of the dat
//         for (int cordIdx = 0; cordIdx < SPACEDIM; cordIdx++) {
//             size[cordIdx] = BlockSize(blockIndex)[cordIdx];
//         }
//        block.Get() =
//             ops_decl_block(SPACEDIM, (char*)blockName.c_str());
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
//         delete[] size;
//     }
//     //TODO temporialy here, later we shall have a define geometry
//     g_NodeType.CreateFieldFromScratch();
//     g_CoordinateXYZ.SetDataDim(SPACEDIM);
//     g_CoordinateXYZ.CreateFieldFromScratch();
//     g_GeometryProperty.CreateFieldFromScratch();
// }

void AssignCoordinates(Block& block,
                       const std::vector<std::vector<Real>>& blockCoordinates) {
#ifdef OPS_2D
    if (SPACEDIM == 2) {
        std::vector<int> range(2 * SPACEDIM);
        range.assign(block.WholeRange().begin(), block.WholeRange().end());
        const Real* coordinateX{blockCoordinates.at(0).data()};
        const SizeType sizeX{block.Size().at(0)};
        const Real* coordinateY{blockCoordinates.at(1).data()};
        const SizeType sizeY{block.Size().at(1)};
        ops_par_loop(KerSetCoordinates, "KerSetCoordinates",
                    block.Get(), SPACEDIM, range.data(),
                     ops_arg_dat(g_CoordinateXYZ[blockIndex], SPACEDIM,
                                 LOCALSTENCIL, "double", OPS_WRITE),
                     ops_arg_idx(),
                     ops_arg_gbl(coordinateX, sizeX, "double", OPS_READ),
                     ops_arg_gbl(coordinateY, sizeY, "double", OPS_READ));
    }
#endif

#ifdef OPS_3D
    if (SPACEDIM == 3) {
        std::vector<int> range(2 * SPACEDIM);
        range.assign(block.WholeRange().begin(), block.WholeRange().end());
        const Real* coordinateX{blockCoordinates.at(0).data()};
        const SizeType sizeX{block.Size().at(0)};
        const Real* coordinateY{blockCoordinates.at(1).data()};
        const SizeType sizeY{block.Size().at(1)};
        const Real* coordinateZ{blockCoordinates.at(2).data()};
        const SizeType sizeZ{block.Size().at(2)};
        ops_par_loop(KerSetCoordinates3D, "KerSetCoordinates3D", block.Get(),
                     SPACEDIM, range.data(),
                     ops_arg_dat(CoordinateXYZ[block.ID()], SPACEDIM,
                                 LOCALSTENCIL, "double", OPS_WRITE),
                     ops_arg_idx(),
                     ops_arg_gbl(coordinateX, sizeX, "double", OPS_READ),
                     ops_arg_gbl(coordinateY, sizeY, "double", OPS_READ),
                     ops_arg_gbl(coordinateZ, sizeZ, "double", OPS_READ)

        );
    }
#endif
}

void DefineBlocks(const std::vector<SizeType>& blockIds,
                  const std::vector<std::string>& blockNames,
                  const std::vector<SizeType>& blockSizes) {
    const SizeType blockNum{blockIds.size()};
    if (blockNum != blockNames.size()) {
        ops_printf(
            "Error! The size of blockIds %i is inconsistent with the size of "
            "blockNames %i!\n",
            blockNum, blockNames.size());
        assert(blockNum == blockNames.size());
    }
    if ((SPACEDIM * blockNum) != blockSizes.size()) {
        ops_printf(
            "Error! The size of blockIds %i is inconsistent with the size of "
            "blockSize %i!\n",
            blockNum, blockSizes.size());
        assert(blockNum == blockSizes.size());
    }
    for (SizeType i = 0; i < blockNum; i++) {
        const SizeType blockId{blockIds.at(i)};
        std::vector<SizeType> blockSize(SPACEDIM);
        for (int j = 0; j < SPACEDIM; j++) {
            blockSize.at(j) = blockSizes.at(i * SPACEDIM + j);
        }
        Block block(blockId, blockNames[i], blockSize);
        BLOCKS.emplace(blockId, block);
    }
}

void PrepareFlowField() {
    ops_printf("The coordinates are assigned!\n");
    for (auto& idBlock: BLOCKS) {
        Block& block{idBlock.second};
        const SizeType blockId{idBlock.first};
        SetBlockGeometryProperty(block);
        ops_printf("The geometry property for Block %i is set!\n", blockId);
        for (int compoId = 0; compoId < NUMCOMPONENTS; compoId++) {
            SetBulkandHaloNodesType(block, compoId);
            ops_printf(
                "The bulk and halo node property are set for Component %i at "
                "Block %i\n",
                compoId, blockId);
        }
        AssignCoordinates(block, COORDINATES.at(blockId));
    }
    SetBoundaryNodeType();
}