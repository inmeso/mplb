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

/*! @brief   Head files for boundary conditions
 * @author  Jianping Meng
 * @details Declaring functions related to boundary conditions.
 */
//#include "boundary.h"
#include "flowfield.h"
#include "type.h"
#include "boundary.h"
#include <cassert>
#include "model.h"
/*!
 * boundaryHaloPt: the halo point needed by the boundary condition
 * In general, the periodic boundary conditions will need one halo point
 * even such as the scheme-collision scheme need no halo points.
 * Other physical boundary schemes usually don't need halo point.
 *
 */

// Vector to assemble all boundary conditions so as to use
// in TreatDomainBoundary().

std::vector<BlockBoundary> blockBoundaries;

int boundaryHaloPt{1};

const std::vector<BlockBoundary>& BlockBoundaries() { return blockBoundaries; }


void DefineBlockBoundary(int blockIndex, int componentID,
                         BoundarySurface boundarySurface,
                         const VertexType boundaryType) {
    if (boundaryType != VertexType::VirtualBoundary) {
        ops_printf(
            "Error: This routine is specially for defining virtual "
            "boundary!\n");
    }

    BlockBoundary blockBoundary;
    blockBoundary.blockIndex = blockIndex;
    blockBoundary.componentID = componentID;
    blockBoundary.givenVars = std::vector<Real>();
    blockBoundary.boundarySurface = boundarySurface;
    blockBoundary.boundaryScheme = BoundaryScheme::None;
    blockBoundary.boundaryType = boundaryType;
    blockBoundaries.push_back(blockBoundary);
    ops_printf(
        "The scheme %i is adopted for Component %i at Surface %i, boundary "
        "type %i and Block %i\n",
        blockBoundary.boundaryScheme, blockBoundary.componentID,
        blockBoundary.boundarySurface, blockBoundary.boundaryType,
        blockBoundary.blockIndex);
}

void DefineBlockBoundary(int blockIndex, int componentID,
                         BoundarySurface boundarySurface,
                         BoundaryScheme boundaryScheme,
                         const std::vector<VariableTypes>& macroVarTypes,
                         const std::vector<Real>& macroVarValues,
                         const VertexType boundaryType) {
    const SizeType numMacroVarTypes{macroVarTypes.size()};
    const SizeType numMacroVarValues{macroVarValues.size()};
    if (numMacroVarTypes != numMacroVarValues){
        ops_printf(
            "Error: Please specify same number of variable types and values\n");
    }

    if ((boundaryScheme == BoundaryScheme::EQMDiffuseRefl ||
         boundaryScheme == BoundaryScheme::EQNNoSlip ||
         boundaryScheme == BoundaryScheme::KineticDiffuseWall ||
         boundaryScheme == BoundaryScheme::BounceBack) &&
        (boundaryType == VertexType::Wall)) {
        if (macroVarTypes.at(0) != VariableTypes::Variable_U ||
            macroVarTypes.at(1) != VariableTypes::Variable_V
#ifdef OPS_3D
            || macroVarTypes.at(2) != VariableTypes::Variable_W
#endif
        ) {
            ops_printf(
                "Error: Please specify the velocity in the order of "
                "(u,v,w)!\n");
        }
    }

    BlockBoundary blockBoundary;
    blockBoundary.blockIndex = blockIndex;
    blockBoundary.componentID = componentID;
    blockBoundary.givenVars = macroVarValues;
    blockBoundary.boundarySurface = boundarySurface;
    blockBoundary.boundaryScheme = boundaryScheme;
    blockBoundary.boundaryType = boundaryType;
    blockBoundaries.push_back(blockBoundary);
    ops_printf(
        "The scheme %i is adopted for Component %i at Surface %i, boundary "
        "type %i and Block %i\n",
        blockBoundary.boundaryScheme, blockBoundary.componentID,
        blockBoundary.boundarySurface, blockBoundary.boundaryType,
        blockBoundary.blockIndex);
}

int BoundaryHaloNum() { return boundaryHaloPt; }

void SetBoundaryHaloNum(const int boundaryHaloNum) {
    boundaryHaloPt = boundaryHaloNum;
}

void BoundaryNormal3D(const VertexGeometryType vg, int* unitNormal) {
    if (nullptr != unitNormal) {
        switch (vg) {
            case VG_IP: {
                unitNormal[0] = 1;
                unitNormal[1] = 0;
                unitNormal[2] = 0;
            } break;
            case VG_IM: {
                unitNormal[0] = -1;
                unitNormal[1] = 0;
                unitNormal[2] = 0;
            } break;
            case VG_JP: {
                unitNormal[0] = 0;
                unitNormal[1] = 1;
                unitNormal[2] = 0;
            } break;
            case VG_JM: {
                unitNormal[0] = 0;
                unitNormal[1] = -1;
                unitNormal[2] = 0;
            } break;
            case VG_KP: {
                unitNormal[0] = 0;
                unitNormal[1] = 0;
                unitNormal[2] = 1;
            } break;
            case VG_KM: {
                unitNormal[0] = 0;
                unitNormal[1] = 0;
                unitNormal[2] = -1;
            } break;
            case VG_IPJP_I: {
                unitNormal[0] = 1 / sqrt(2);
                unitNormal[1] = 1 / sqrt(2);
                unitNormal[2] = 0;
            } break;
            case VG_IPJM_I: {
                unitNormal[0] = 1 / sqrt(2);
                unitNormal[1] = -1 / sqrt(2);
                unitNormal[2] = 0;
            } break;
            case VG_IMJP_I: {
                unitNormal[0] = -1 / sqrt(2);
                unitNormal[1] = 1 / sqrt(2);
                unitNormal[2] = 0;
            } break;
            case VG_IMJM_I: {
                unitNormal[0] = -1 / sqrt(2);
                unitNormal[1] = -1 / sqrt(2);
                unitNormal[2] = 0;
            } break;
            case VG_IPKP_I: {
                unitNormal[0] = 1 / sqrt(2);
                unitNormal[1] = 0;
                unitNormal[2] = 1 / sqrt(2);
            } break;
            case VG_IPKM_I: {
                unitNormal[0] = 1 / sqrt(2);
                unitNormal[1] = 0;
                unitNormal[2] = -1 / sqrt(2);
            } break;
            case VG_IMKP_I: {
                unitNormal[0] = -1 / sqrt(2);
                unitNormal[1] = 0;
                unitNormal[2] = 1 / sqrt(2);
            } break;
            case VG_IMKM_I: {
                unitNormal[0] = -1 / sqrt(2);
                unitNormal[1] = 0;
                unitNormal[2] = -1 / sqrt(2);
            } break;
            case VG_JPKP_I: {
                unitNormal[0] = 0;
                unitNormal[1] = 1 / sqrt(2);
                unitNormal[2] = 1 / sqrt(2);
            } break;
            case VG_JPKM_I: {
                unitNormal[0] = 0;
                unitNormal[1] = 1 / sqrt(2);
                unitNormal[2] = -1 / sqrt(2);
            } break;
            case VG_JMKP_I: {
                unitNormal[0] = 0;
                unitNormal[1] = -1 / sqrt(2);
                unitNormal[2] = 1 / sqrt(2);
            } break;
            case VG_JMKM_I: {
                unitNormal[0] = 0;
                unitNormal[1] = -1 / sqrt(2);
                unitNormal[2] = -1 / sqrt(2);
            } break;

            case VG_IPJP_O: {
                unitNormal[0] = 1 / sqrt(2);
                unitNormal[1] = 1 / sqrt(2);
                unitNormal[2] = 0;
            } break;
            case VG_IPJM_O: {
                unitNormal[0] = 1 / sqrt(2);
                unitNormal[1] = -1 / sqrt(2);
                unitNormal[2] = 0;
            } break;
            case VG_IMJP_O: {
                unitNormal[0] = -1 / sqrt(2);
                unitNormal[1] = 1 / sqrt(2);
                unitNormal[2] = 0;
            } break;
            case VG_IMJM_O: {
                unitNormal[0] = -1 / sqrt(2);
                unitNormal[1] = -1 / sqrt(2);
                unitNormal[2] = 0;
            } break;
            case VG_IPKP_O: {
                unitNormal[0] = 1 / sqrt(2);
                unitNormal[1] = 0;
                unitNormal[2] = 1 / sqrt(2);
            } break;
            case VG_IPKM_O: {
                unitNormal[0] = 1 / sqrt(2);
                unitNormal[1] = 0;
                unitNormal[2] = -1 / sqrt(2);
            } break;
            case VG_IMKP_O: {
                unitNormal[0] = -1 / sqrt(2);
                unitNormal[1] = 0;
                unitNormal[2] = 1 / sqrt(2);
            } break;
            case VG_IMKM_O: {
                unitNormal[0] = -1 / sqrt(2);
                unitNormal[1] = 0;
                unitNormal[2] = -1 / sqrt(2);
            } break;
            case VG_JPKP_O: {
                unitNormal[0] = 0;
                unitNormal[1] = 1 / sqrt(2);
                unitNormal[2] = 1 / sqrt(2);
            } break;
            case VG_JPKM_O: {
                unitNormal[0] = 0;
                unitNormal[1] = 1 / sqrt(2);
                unitNormal[2] = -1 / sqrt(2);
            } break;
            case VG_JMKP_O: {
                unitNormal[0] = 0;
                unitNormal[1] = -1 / sqrt(2);
                unitNormal[2] = 1 / sqrt(2);
            } break;
            case VG_JMKM_O: {
                unitNormal[0] = 0;
                unitNormal[1] = -1 / sqrt(2);
                unitNormal[2] = -1 / sqrt(2);
            } break;

            case VG_IPJPKP_I: {
                unitNormal[0] = 1 / sqrt(3);
                unitNormal[1] = 1 / sqrt(3);
                unitNormal[2] = 1 / sqrt(3);
            } break;
            case VG_IPJPKM_I: {
                unitNormal[0] = 1 / sqrt(3);
                unitNormal[1] = 1 / sqrt(3);
                unitNormal[2] = -1 / sqrt(3);
            } break;
            case VG_IPJMKP_I: {
                unitNormal[0] = 1 / sqrt(3);
                unitNormal[1] = -1 / sqrt(3);
                unitNormal[2] = 1 / sqrt(3);
            } break;
            case VG_IPJMKM_I: {
                unitNormal[0] = 1 / sqrt(3);
                unitNormal[1] = -1 / sqrt(3);
                unitNormal[2] = -1 / sqrt(3);
            } break;
            case VG_IMJPKP_I: {
                unitNormal[0] = -1 / sqrt(3);
                unitNormal[1] = 1 / sqrt(3);
                unitNormal[2] = 1 / sqrt(3);
            } break;
            case VG_IMJPKM_I: {
                unitNormal[0] = -1 / sqrt(3);
                unitNormal[1] = 1 / sqrt(3);
                unitNormal[2] = -1 / sqrt(3);
            } break;
            case VG_IMJMKP_I: {
                unitNormal[0] = -1 / sqrt(3);
                unitNormal[1] = -1 / sqrt(3);
                unitNormal[2] = 1 / sqrt(3);
            } break;
            case VG_IMJMKM_I: {
                unitNormal[0] = -1 / sqrt(3);
                unitNormal[1] = -1 / sqrt(3);
                unitNormal[2] = -1 / sqrt(3);
            } break;
            case VG_IPJPKP_O: {
                unitNormal[0] = 1 / sqrt(3);
                unitNormal[1] = 1 / sqrt(3);
                unitNormal[2] = 1 / sqrt(3);
            } break;
            case VG_IPJPKM_O: {
                unitNormal[0] = 1 / sqrt(3);
                unitNormal[1] = 1 / sqrt(3);
                unitNormal[2] = -1 / sqrt(3);
            } break;
            case VG_IPJMKP_O: {
                unitNormal[0] = 1 / sqrt(3);
                unitNormal[1] = -1 / sqrt(3);
                unitNormal[2] = 1 / sqrt(3);
            } break;
            case VG_IPJMKM_O: {
                unitNormal[0] = 1 / sqrt(3);
                unitNormal[1] = -1 / sqrt(3);
                unitNormal[2] = -1 / sqrt(3);
            } break;
            case VG_IMJPKP_O: {
                unitNormal[0] = -1 / sqrt(3);
                unitNormal[1] = 1 / sqrt(3);
                unitNormal[2] = 1 / sqrt(3);
            } break;
            case VG_IMJPKM_O: {
                unitNormal[0] = -1 / sqrt(3);
                unitNormal[1] = 1 / sqrt(3);
                unitNormal[2] = -1 / sqrt(3);
            } break;
            case VG_IMJMKP_O: {
                unitNormal[0] = -1 / sqrt(3);
                unitNormal[1] = -1 / sqrt(3);
                unitNormal[2] = 1 / sqrt(3);
            } break;
            case VG_IMJMKM_O: {
                unitNormal[0] = -1 / sqrt(3);
                unitNormal[1] = -1 / sqrt(3);
                unitNormal[2] = -1 / sqrt(3);
            } break;
            default:
                break;
        }
    }
}
#ifdef OPS_3D
void ImplementBoundary3D() {
    for (const auto& boundary : BlockBoundaries()) {
        const Block& block{g_Block().at(boundary.blockIndex)};
        TreatBlockBoundary3D(block, boundary.componentID,
                             boundary.givenVars.data(), boundary.boundaryScheme,
                             boundary.boundarySurface);
    }
}
#endif

#ifdef OPS_2D
void ImplementBoundary() {
    for (const auto& boundary : BlockBoundaries()) {
        const Block& block{g_Block().at(boundary.blockIndex)};
        TreatBlockBoundary(block, boundary.componentID,
                             boundary.givenVars.data(), boundary.boundaryScheme,
                             boundary.boundarySurface);
    }
}
#endif
