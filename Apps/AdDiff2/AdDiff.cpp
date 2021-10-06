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

/** @brief An example main source code of stimulating 3D lid-driven cavity flow
 *  @author Jianping Meng
 **/
#include <cmath>
#include <iostream>
#include <ostream>
#include <string>
#include "mplb.h"
#include "ops_seq_v2.h"
#include "AdDiff_kernel.inc"
//Provide macroscopic initial conditions

//////////////////////////////////////////////////////////////////////
////////////////                  2D                  ////////////////
//////////////////////////////////////////////////////////////////////

#ifdef OPS_2D
void UpdateConcentration() {
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        const Real* pdt{pTimeStep()};
        const Real* ttt{g_Time()};

            //const Component& compo{idCompo.second};
            auto it = g_Components().begin();
            const auto& Iter1{*it};
            const Component& compoVel{Iter1.second};
            std::advance(it, 1);
            const auto& Iter2{*it};
            const Component& compoRho{Iter2.second};
            const int compoRhoId{compoRho.id};
            //std::cout << compoRhoId << "\n";
            //g_MacroVars().at(compoRho.macroVars.at(Variable_Rho).id).CreateHalos();
            ops_par_loop(
                KerCalcConcentration, "KerCalcConcentration", block.Get(),
                SpaceDim(), iterRng.data(),
                ops_arg_dat(g_MacroVars()
                                    .at(compoRho.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex), 1,
                                     LOCALSTENCIL, "Real", OPS_RW),
                ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                     LOCALSTENCIL, "Real", OPS_READ),
                ops_arg_gbl(ttt, 1, "Real", OPS_READ),
                ops_arg_dat(g_MacroVars().at(compoVel.uId).at(blockIndex),
                                     1, LOCALSTENCIL, "Real", OPS_READ),
                ops_arg_dat(g_MacroVars().at(compoVel.vId).at(blockIndex),
                                     1, LOCALSTENCIL, "Real", OPS_READ),
                ops_arg_gbl(pdt, 1, "double", OPS_READ));
            //break;
        
        
    }
}



void Calcphi2Gradients() {
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        const int order{2};
            auto it = g_Components().begin();
            std::advance(it, 1);
            const auto& Iter2{*it};
            const Component& compoRho{Iter2.second};
            const int compoRhoId{compoRho.id};
            ops_par_loop(
                KerCalcGradients, "KerCalcGradients", block.Get(), SpaceDim(),
                iterRng.data(),
                ops_arg_dat(g_phiGrad()[blockIndex], SpaceDim(), LOCALSTENCIL, "double",
                            OPS_RW),
                ops_arg_dat(g_MacroVars()
                                    .at(compoRho.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex), 1,
                                     ONEPTREGULARSTENCIL, "Real", OPS_READ),
                ops_arg_dat(g_NodeType().at(Iter2.first).at(blockIndex), 1,
                            LOCALSTENCIL, "int", OPS_READ),
                ops_arg_dat(g_GeometryProperty()[blockIndex], 1,
                            LOCALSTENCIL, "int", OPS_READ),
                ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                     LOCALSTENCIL, "Real", OPS_READ),
                ops_arg_gbl(Iter2.second.index, 2, "int", OPS_READ),
                ops_arg_gbl(&order, 1, "int", OPS_READ),
                ops_arg_idx());
        
    }
}


void CalcmuGradients() {
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        const int order{1};
            auto it = g_Components().begin();
            std::advance(it, 1);
            const auto& Iter2{*it};
            ops_par_loop(
                KerCalcGradients, "KerCalcGradients", block.Get(), SpaceDim(),
                iterRng.data(),
                ops_arg_dat(g_muGrad()[blockIndex], SpaceDim(), LOCALSTENCIL, "double",
                            OPS_RW),
                ops_arg_dat(g_mu()[blockIndex], 1, ONEPTREGULARSTENCIL, "double",
                            OPS_READ),
                ops_arg_dat(g_NodeType().at(Iter2.first).at(blockIndex), 1,
                            LOCALSTENCIL, "int", OPS_READ),
                ops_arg_dat(g_GeometryProperty()[blockIndex], 1,
                            LOCALSTENCIL, "int", OPS_READ),
                ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                     LOCALSTENCIL, "Real", OPS_READ),
                ops_arg_gbl(Iter2.second.index, 2, "int", OPS_READ),
                ops_arg_gbl(&order, 1, "int", OPS_READ),
                ops_arg_idx());

        
    }
} 

void CalcMu() {
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
            auto it = g_Components().begin();
            std::advance(it, 1);
            const auto& Iter2{*it};
            const Component& compoRho{Iter2.second};
            const int compoRhoId{compoRho.id};
            ops_par_loop(
                KerCalcMu, "KerCalcMu", block.Get(), SpaceDim(),
                iterRng.data(),
                ops_arg_dat(g_mu()[blockIndex], 1, LOCALSTENCIL, "double",
                            OPS_RW),
                ops_arg_dat(g_MacroVars()
                                    .at(compoRho.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex), 1,
                                     LOCALSTENCIL, "Real", OPS_READ),
                ops_arg_dat(g_phiGrad()[blockIndex], SpaceDim(),
                            LOCALSTENCIL, "double", OPS_READ),
                ops_arg_dat(g_NodeType().at(compoRho.id).at(blockIndex), 1,
                                     LOCALSTENCIL, "int", OPS_READ),
                ops_arg_dat(g_GeometryProperty().at(blockIndex), 1,
                                     LOCALSTENCIL, "int", OPS_READ),
                ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                     LOCALSTENCIL, "Real", OPS_READ));
            

        
    }
} 


void SetInitialMacrosVars() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIdx{block.ID()};
        auto it = g_Components().begin();
        const auto& Iter{*it};
            //const Component& compo{idCompo.second};
            //const int rhoId{compo.macroVars.at(Variable_Rho).id};
            ops_par_loop(KerSetInitialMacroVars, "KerSetInitialMacroVars",
                         block.Get(), SpaceDim(), iterRng.data(),
                         ops_arg_dat(g_MacroVars().at(Iter.second.macroVars.at(Variable_Rho).id).at(blockIdx), 1,
                                     LOCALSTENCIL, "Real", OPS_RW),
                         ops_arg_dat(g_MacroVars().at(Iter.second.uId).at(blockIdx),
                                     1, LOCALSTENCIL, "Real", OPS_RW),
                         ops_arg_dat(g_MacroVars().at(Iter.second.vId).at(blockIdx),
                                     1, LOCALSTENCIL, "Real", OPS_RW),
                         ops_arg_dat(g_CoordinateXYZ()[blockIdx], SpaceDim(),
                                     LOCALSTENCIL, "Real", OPS_READ),
                         ops_arg_idx());
        
    }
}

void CalcPhiWetting() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        auto compo = g_Components().at(1);
            //const Component& compo{idCompo.second};
            //const int rhoId{compo.macroVars.at(Variable_Rho).id};
            ops_par_loop(KerUpdateRhoWetting, "KerUpdateRhoWetting",
                         block.Get(), SpaceDim(), iterRng.data(),
                         ops_arg_dat(g_MacroVars().at(compo.macroVars.at(Variable_Rho).id).at(blockIndex), 1,
                                     ONEPTLATTICESTENCIL, "Real", OPS_RW),
                         ops_arg_dat(g_GeometryProperty().at(blockIndex), 1,
                                     ONEPTREGULARSTENCIL, "int", OPS_READ),
                         ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                                     LOCALSTENCIL, "int", OPS_READ),
                         ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                     LOCALSTENCIL, "Real", OPS_READ));
        
    }
}

void PrintPhi() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
            auto it = g_Components().begin();
            const auto& Iter1{*it};
            std::advance(it, 1);
            const auto& Iter2{*it};
            const Component& compoRho{Iter1.second};
            const int compoRhoId{compoRho.id};
            //std::cout << compoRhoId << "\n";

            ops_par_loop(KerPrintPhi, "KerPrintPhi",
                         block.Get(), SpaceDim(), iterRng.data(),
                         ops_arg_dat(g_f()[blockIndex], 1, LOCALSTENCIL,
                                        "double", OPS_READ));
    }
}

//Provide macroscopic body-force term
void UpdateMacroscopicBodyForce(const Real time) {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIdx{block.ID()};
            auto it = g_Components().begin();
            const auto& Iter1{*it};
            std::advance(it, 1);
            const auto& Iter2{*it};
            const Component& compoRho{Iter2.second};

            ops_par_loop(KerUpdateMacroBodyForce, "KerUpdateMacroBodyForce",
                         block.Get(), SpaceDim(), iterRng.data(),
                         ops_arg_dat(g_MacroBodyforce().at(Iter1.second.id).at(blockIdx), SpaceDim(),
                                     LOCALSTENCIL, "Real", OPS_RW),
                         ops_arg_dat(g_MacroVars()
                                    .at(compoRho.macroVars.at(Variable_Rho).id)
                                    .at(blockIdx), 1,
                                     LOCALSTENCIL, "Real", OPS_READ),
                         ops_arg_dat(g_muGrad()[blockIdx], SpaceDim(),
                            LOCALSTENCIL, "double", OPS_READ),
                         ops_arg_dat(g_GeometryProperty().at(blockIdx), 1,
                                     LOCALSTENCIL, "int", OPS_READ),
                         ops_arg_dat(g_NodeType().at(compoRho.id).at(blockIdx), 1,
                                     LOCALSTENCIL, "int", OPS_READ),
                         ops_arg_dat(g_CoordinateXYZ()[blockIdx], SpaceDim(),
                                     LOCALSTENCIL, "Real", OPS_READ),
                         ops_arg_idx());

        
    }
}

void SetSolid() {
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        const int sizeX{block.Size().at(0)};
        const int sizeY{block.Size().at(1)};
        for (const auto& idCompo : g_Components()) {
        const Component& compo{idCompo.second};
        ops_par_loop(KerSetSolid, "KerSetSolid",
                     block.Get(), SpaceDim(), iterRng.data(),
                     ops_arg_dat(g_GeometryProperty().at(blockIndex), 1,
                                 LOCALSTENCIL, "int", OPS_WRITE),
                     ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                                 LOCALSTENCIL, "int", OPS_WRITE),
                     ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                 LOCALSTENCIL, "Real", OPS_READ),
                     ops_arg_gbl(&sizeX, 1, "double", OPS_READ),
                     ops_arg_gbl(&sizeY, 1, "double", OPS_READ));
        }
    }
}

void SetEmbeddedBodyGeometry() {
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        for (const auto& idCompo : g_Components()) {
        const Component& compo{idCompo.second};
        ops_par_loop(KerSetEmbeddedBodyGeometry, "KerSetEmbeddedBodyGeometry",
                     block.Get(), SpaceDim(), iterRng.data(),
                     ops_arg_dat(g_GeometryProperty().at(blockIndex), 1,
                                 LOCALSTENCIL, "int", OPS_WRITE),
                     ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                                 ONEPTLATTICESTENCIL, "int", OPS_READ));
        }
    }
}

#endif
//////////////////////////////////////////////////////////////////////
////////////////                  3D                  ////////////////
//////////////////////////////////////////////////////////////////////
#ifdef OPS_3D
void UpdateConcentration3D() {
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        const Real* pdt{pTimeStep()};
        const Real* ttt{g_Time()};

            //const Component& compo{idCompo.second};
            auto it = g_Components().begin();
            const auto& Iter1{*it};
            const Component& compoVel{Iter1.second};
            std::advance(it, 1);
            const auto& Iter2{*it};
            const Component& compoRho{Iter2.second};
            const int compoRhoId{compoRho.id};
            //std::cout << compoRhoId << "\n";
            //g_MacroVars().at(compoRho.macroVars.at(Variable_Rho).id).CreateHalos();
            ops_par_loop(
                KerCalcConcentration3D, "KerCalcConcentration", block.Get(),
                SpaceDim(), iterRng.data(),
                ops_arg_dat(g_MacroVars()
                                    .at(compoRho.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex), 1,
                                     LOCALSTENCIL, "Real", OPS_RW),
                ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                     LOCALSTENCIL, "Real", OPS_READ),
                ops_arg_gbl(ttt, 1, "Real", OPS_READ),
                ops_arg_dat(g_MacroVars().at(compoVel.uId).at(blockIndex),
                                     1, LOCALSTENCIL, "Real", OPS_READ),
                ops_arg_dat(g_MacroVars().at(compoVel.vId).at(blockIndex),
                                     1, LOCALSTENCIL, "Real", OPS_READ),
                ops_arg_gbl(pdt, 1, "double", OPS_READ));
            //break;
        
        
    }
}

void Calcphi2Gradients3D() {
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        const int order{2};
            auto it = g_Components().begin();
            std::advance(it, 1);
            const auto& Iter2{*it};
            const Component& compoRho{Iter2.second};
            const int compoRhoId{compoRho.id};
            ops_par_loop(
                KerCalcGradients3D, "KerCalcGradients", block.Get(), SpaceDim(),
                iterRng.data(),
                ops_arg_dat(g_phiGrad()[blockIndex], SpaceDim(), LOCALSTENCIL, "double",
                            OPS_RW),
                ops_arg_dat(g_MacroVars()
                                    .at(compoRho.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex), 1,
                                     ONEPTREGULARSTENCIL, "Real", OPS_READ),
                ops_arg_dat(g_NodeType().at(Iter2.first).at(blockIndex), 1,
                            LOCALSTENCIL, "int", OPS_READ),
                ops_arg_dat(g_GeometryProperty()[blockIndex], 1,
                            LOCALSTENCIL, "int", OPS_READ),
                ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                     LOCALSTENCIL, "Real", OPS_READ),
                ops_arg_gbl(Iter2.second.index, 2, "int", OPS_READ),
                ops_arg_gbl(&order, 1, "int", OPS_READ),
                ops_arg_idx());
        
    }
}

void CalcmuGradients3D() {
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        const int order{1};
            auto it = g_Components().begin();
            std::advance(it, 1);
            const auto& Iter2{*it};
            ops_par_loop(
                KerCalcGradients3D, "KerCalcGradients", block.Get(), SpaceDim(),
                iterRng.data(),
                ops_arg_dat(g_muGrad()[blockIndex], SpaceDim(), LOCALSTENCIL, "double",
                            OPS_RW),
                ops_arg_dat(g_mu()[blockIndex], 1, ONEPTREGULARSTENCIL, "double",
                            OPS_READ),
                ops_arg_dat(g_NodeType().at(Iter2.first).at(blockIndex), 1,
                            LOCALSTENCIL, "int", OPS_READ),
                ops_arg_dat(g_GeometryProperty()[blockIndex], 1,
                            LOCALSTENCIL, "int", OPS_READ),
                ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                     LOCALSTENCIL, "Real", OPS_READ),
                ops_arg_gbl(Iter2.second.index, 2, "int", OPS_READ),
                ops_arg_gbl(&order, 1, "int", OPS_READ),
                ops_arg_idx());

        
    }
} 

void CalcMu3D() {
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
            auto it = g_Components().begin();
            std::advance(it, 1);
            const auto& Iter2{*it};
            const Component& compoRho{Iter2.second};
            const int compoRhoId{compoRho.id};
            ops_par_loop(
                KerCalcMu3D, "KerCalcMu", block.Get(), SpaceDim(),
                iterRng.data(),
                ops_arg_dat(g_mu()[blockIndex], 1, LOCALSTENCIL, "double",
                            OPS_RW),
                ops_arg_dat(g_MacroVars()
                                    .at(compoRho.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex), 1,
                                     LOCALSTENCIL, "Real", OPS_READ),
                ops_arg_dat(g_phiGrad()[blockIndex], SpaceDim(),
                            LOCALSTENCIL, "double", OPS_READ),
                         ops_arg_dat(g_NodeType().at(compoRho.id).at(blockIndex), 1,
                                     LOCALSTENCIL, "int", OPS_READ));
            

        
    }
} 


void SetInitialMacrosVars3D() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIdx{block.ID()};
        auto it = g_Components().begin();
        const auto& Iter{*it};
            //const Component& compo{idCompo.second};
            //const int rhoId{compo.macroVars.at(Variable_Rho).id};
            ops_par_loop(KerSetInitialMacroVars3D, "KerSetInitialMacroVars",
                         block.Get(), SpaceDim(), iterRng.data(),
                         ops_arg_dat(g_MacroVars().at(Iter.second.macroVars.at(Variable_Rho).id).at(blockIdx), 1,
                                     LOCALSTENCIL, "Real", OPS_RW),
                         ops_arg_dat(g_MacroVars().at(Iter.second.uId).at(blockIdx),
                                     1, LOCALSTENCIL, "Real", OPS_RW),
                         ops_arg_dat(g_MacroVars().at(Iter.second.vId).at(blockIdx),
                                     1, LOCALSTENCIL, "Real", OPS_RW),
                         ops_arg_dat(g_MacroVars().at(Iter.second.wId).at(blockIdx),
                                     1, LOCALSTENCIL, "Real", OPS_RW),
                         ops_arg_dat(g_CoordinateXYZ()[blockIdx], SpaceDim(),
                                     LOCALSTENCIL, "Real", OPS_READ),
                         ops_arg_idx());
        
    }
}

void CalcPhiWetting3D() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        auto compo = g_Components().at(1);
            //const Component& compo{idCompo.second};
            //const int rhoId{compo.macroVars.at(Variable_Rho).id};
            ops_par_loop(KerUpdateRhoWetting3D, "KerUpdateRhoWetting3D",
                         block.Get(), SpaceDim(), iterRng.data(),
                         ops_arg_dat(g_MacroVars().at(compo.macroVars.at(Variable_Rho).id).at(blockIndex), 1,
                                     ONEPTLATTICESTENCIL, "Real", OPS_RW),
                         ops_arg_dat(g_GeometryProperty().at(blockIndex), 1,
                                     ONEPTREGULARSTENCIL, "int", OPS_READ),
                         ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                                     LOCALSTENCIL, "int", OPS_READ),
                         ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                     LOCALSTENCIL, "Real", OPS_READ));
        
    }
}


void PrintPhi3D() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
            auto it = g_Components().begin();
            const auto& Iter1{*it};
            std::advance(it, 1);
            const auto& Iter2{*it};
            const Component& compoRho{Iter2.second};
            const int compoRhoId{compoRho.id};
            //::cout << compoRhoId << "\n";

            ops_par_loop(KerPrintPhi3D, "KerPrintPhi",
                         block.Get(), SpaceDim(), iterRng.data(),
                         ops_arg_dat(g_MacroVars()
                                    .at(compoRho.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex), 1,
                                     ONEPTLATTICESTENCIL, "Real", OPS_READ),
                         ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                     LOCALSTENCIL, "Real", OPS_READ));
    }
}

//Provide macroscopic body-force term
void UpdateMacroscopicBodyForce3D(const Real time) {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIdx{block.ID()};
            auto it = g_Components().begin();
            const auto& Iter1{*it};
            std::advance(it, 1);
            const auto& Iter2{*it};
            const Component& compoRho{Iter2.second};

            ops_par_loop(KerUpdateMacroBodyForce3D, "KerUpdateMacroBodyForce",
                         block.Get(), SpaceDim(), iterRng.data(),
                         ops_arg_dat(g_MacroBodyforce().at(Iter1.second.id).at(blockIdx), SpaceDim(),
                                     LOCALSTENCIL, "Real", OPS_RW),
                         ops_arg_dat(g_MacroVars()
                                    .at(compoRho.macroVars.at(Variable_Rho).id)
                                    .at(blockIdx), 1,
                                     LOCALSTENCIL, "Real", OPS_READ),
                         ops_arg_dat(g_muGrad()[blockIdx], SpaceDim(),
                                     LOCALSTENCIL, "Real", OPS_READ),
                         ops_arg_dat(g_GeometryProperty().at(blockIdx), 1,
                                     LOCALSTENCIL, "int", OPS_READ),
                         ops_arg_dat(g_NodeType().at(compoRho.id).at(blockIdx), 1,
                                     LOCALSTENCIL, "int", OPS_READ),
                         ops_arg_dat(g_CoordinateXYZ()[blockIdx], SpaceDim(),
                                     LOCALSTENCIL, "Real", OPS_READ),
                         ops_arg_idx());

        
    }
}

void SetSolid3D() {
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        const int sizeX{block.Size().at(0)};
        const int sizeY{block.Size().at(1)};
        const int sizeZ{block.Size().at(2)};
        for (const auto& idCompo : g_Components()) {
        const Component& compo{idCompo.second};
        ops_par_loop(KerSetSolid3D, "KerSetSolid",
                     block.Get(), SpaceDim(), iterRng.data(),
                     ops_arg_dat(g_GeometryProperty().at(blockIndex), 1,
                                 LOCALSTENCIL, "int", OPS_WRITE),
                     ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                                 LOCALSTENCIL, "int", OPS_WRITE),
                     ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                 LOCALSTENCIL, "Real", OPS_READ),
                     ops_arg_gbl(&sizeX, 1, "double", OPS_READ),
                     ops_arg_gbl(&sizeY, 1, "double", OPS_READ),
                     ops_arg_gbl(&sizeZ, 1, "double", OPS_READ));
        }
    }
}

void SetEmbeddedBodyGeometry3D() {
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        for (const auto& idCompo : g_Components()) {
        const Component& compo{idCompo.second};
        ops_par_loop(KerSetEmbeddedBodyGeometry3D, "KerSetEmbeddedBodyGeometry",
                     block.Get(), SpaceDim(), iterRng.data(),
                     ops_arg_dat(g_GeometryProperty().at(blockIndex), 1,
                                 LOCALSTENCIL, "int", OPS_WRITE),
                     ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                                 ONEPTLATTICESTENCIL, "int", OPS_READ));
        }
    }
}

#endif
//////////////////////////////////////////////////////////////////////
////////////////               Simulate               ////////////////
//////////////////////////////////////////////////////////////////////
#ifdef OPS_2D
void simulate() {

    std::string caseName{"Advection_Diffusion_DS"}; //DS means diffusive scaling
    SizeType spaceDim{2};
    DefineCase(caseName, spaceDim);
    //Define one block this application
    std::vector<int> blockIds{0};
    std::vector<std::string> blockNames{"Cavity"};
    std::vector<int> blockSize{100, 100};
    Real meshSize{1.};
    std::map<int, std::vector<Real>> startPos{{0, {0.0, 0.0}}};
    DefineBlocks(blockIds, blockNames, blockSize, meshSize, startPos);

    std::vector<int> fromBlockIds{0, 0, 0, 0};
    std::vector<int> toBlockIds{0, 0, 0, 0};

    std::vector<BoundarySurface> fromBoundarySurface{BoundarySurface::Top,
        BoundarySurface::Bottom,BoundarySurface::Left,
        BoundarySurface::Right};
    std::vector<BoundarySurface> toBoundarySurface{BoundarySurface::Bottom,
        BoundarySurface::Top,BoundarySurface::Right,
        BoundarySurface::Left};

    std::vector<VertexType> blockConnectionType{
        VertexType::MDPeriodic, VertexType::MDPeriodic,
        VertexType::MDPeriodic, VertexType::MDPeriodic};

    DefineBlockConnection(fromBlockIds, fromBoundarySurface, toBlockIds,
                          toBoundarySurface, blockConnectionType);


    std::vector<std::string> compoNames{"Fluid","Diffusion"};
    std::vector<int> compoid{0,1};
    std::vector<std::string> lattNames{"d2q9_diffusive","d2q9_diffusive"};
    std::vector<Real> tauRef{1,1};
    DefineComponents(compoNames, compoid, lattNames, tauRef);

    std::vector<VariableTypes> marcoVarTypes{Variable_Rho, Variable_U,
                                             Variable_V, Variable_Rho};
    std::vector<std::string> macroVarNames{"rho", "u", "v", "C"};
    std::vector<int> macroVarId{0, 1, 2, 3};
    std::vector<int> macroCompoId{0, 0, 0, 1};
    DefineMacroVars(marcoVarTypes, macroVarNames, macroVarId, macroCompoId);

    std::vector<CollisionType> collisionTypes{Collision_BGKADF,Collision_BGKADF};
    std::vector<int> collisionCompoId{0,1};
    DefineCollision(collisionTypes, collisionCompoId);

    std::vector<BodyForceType> bodyForceTypes{BodyForce_None,BodyForce_None};
    std::vector<SizeType> bodyForceCompoId{0,1};
    DefineBodyForce(bodyForceTypes, bodyForceCompoId);

    SchemeType scheme{Scheme_StreamCollision};
    DefineScheme(scheme);

    // Setting boundary conditions
    SizeType componentId{0};
    std::vector<VariableTypes> macroVarTypesatBoundary{
        Variable_U, Variable_V};
    std::vector<Real> noSlipStationaryWall{0, 0};

    // Periodic Boundary Conditions
    
    DefineBlockBoundary(0, componentId, BoundarySurface::Top,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);
    DefineBlockBoundary(0, componentId, BoundarySurface::Bottom,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);
    DefineBlockBoundary(0, componentId, BoundarySurface::Left,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);
    DefineBlockBoundary(0, componentId, BoundarySurface::Right,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);
    //DefineBlockBoundary(0, componentId, BoundarySurface::None, BoundaryScheme::EQMDiffuseReflF, macroVarTypesatBoundary,
    //                    noSlipStationaryWall);
    /*DefineBlockBoundary(0, componentId, BoundarySurface::Left,
                        BoundaryScheme::EQMDiffuseReflF, macroVarTypesatBoundary,
                        noSlipStationaryWall);
    DefineBlockBoundary(0, componentId, BoundarySurface::Right,
                        BoundaryScheme::EQMDiffuseReflF, macroVarTypesatBoundary,
                        noSlipStationaryWall);*/

    /*DefineBlockBoundary(0, 1, BoundarySurface::Top,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);
    DefineBlockBoundary(0, 1, BoundarySurface::Bottom,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);*/
    //DefineBlockBoundary(0, 1, BoundarySurface::None, BoundaryScheme::ZouHeVelocity, macroVarTypesatBoundary,
    //                    noSlipStationaryWall);
    /*DefineBlockBoundary(0, 1, BoundarySurface::Left,
                        BoundaryScheme::EQMDiffuseReflG, macroVarTypesatBoundary,
                        noSlipStationaryWall);
    DefineBlockBoundary(0, 1, BoundarySurface::Right,
                        BoundaryScheme::EQMDiffuseReflG, macroVarTypesatBoundary,
                        noSlipStationaryWall);*/
    DefineBlockBoundary(0, 1, BoundarySurface::Left,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);
    DefineBlockBoundary(0, 1, BoundarySurface::Right,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);
    DefineBlockBoundary(0, 1, BoundarySurface::Top,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);
    DefineBlockBoundary(0, 1, BoundarySurface::Bottom,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);


    std::vector<InitialType> initType{Initial_BGKFeq2ndAD,Initial_BGKFeq2ndAD};
    std::vector<SizeType> initalCompoId{0,1};

    //g_mu().CreateHalos();
    //const auto& compoC = g_Components().at(1);
    //g_MacroVars().at(compoC.macroVars.at(Variable_Rho).id).CreateHalos();
    DefineInitialCondition(initType, initalCompoId);

    Partition();
    //SetSolid();
    //SetEmbeddedBodyGeometry();
    ops_diagnostic_output();

    SetInitialMacrosVars();
    //std::cout << "test\n";
    UpdateConcentration();
    //CalcPhiWetting();
    //std::cout << "test1\n";
    //PreDefinedInitialCondition();
    
    //PreDefinedInitialConditionADF();
    //PrintPhi();
    PreDefinedInitialConditionAD();
    UpdateMacroVars();
    
    //PrintPhi();
    std::cout<<"test2\n";
    SetTimeStep(1);

    const Real convergenceCriteria{-1E-7};
    const SizeType checkPeriod{125};
    Real iter{0};
    WriteFlowfieldToHdf5(iter);
    WriteDistributionsToHdf5(iter);
    WriteNodePropertyToHdf5(iter);
    Iterate(StreamCollision,convergenceCriteria, checkPeriod);
}
#endif
#ifdef OPS_3D
void simulate() {

    std::string caseName{"Advection_Diffusion"};
    SizeType spaceDim{3};
    DefineCase(caseName, spaceDim);
    std::vector<int> blockIds{0};
    std::vector<std::string> blockNames{"Cavity"};
    std::vector<int> blockSize{40, 40, 40};
    Real meshSize{1.};
    std::map<int, std::vector<Real>> startPos{{0, {0.0, 0.0, 0.0}}};
    DefineBlocks(blockIds, blockNames, blockSize, meshSize, startPos);

    std::vector<int> fromBlockIds{0, 0, 0, 0};
    std::vector<int> toBlockIds{0, 0, 0, 0};

    std::vector<BoundarySurface> fromBoundarySurface{BoundarySurface::Front,
        BoundarySurface::Back,BoundarySurface::Left,
        BoundarySurface::Right};
    std::vector<BoundarySurface> toBoundarySurface{BoundarySurface::Back,
        BoundarySurface::Front,BoundarySurface::Right,
        BoundarySurface::Left};

    std::vector<VertexType> blockConnectionType{
        VertexType::MDPeriodic, VertexType::MDPeriodic,VertexType::MDPeriodic, VertexType::MDPeriodic};

    DefineBlockConnection(fromBlockIds, fromBoundarySurface, toBlockIds,
                          toBoundarySurface, blockConnectionType);


    std::vector<std::string> compoNames{"Fluid","Diffusion"};
    std::vector<int> compoid{0,1};
    std::vector<std::string> lattNames{"d3q19_diffusive","d3q19_diffusive"};
    std::vector<Real> tauRef{1,1};
    DefineComponents(compoNames, compoid, lattNames, tauRef);

    std::vector<VariableTypes> marcoVarTypes{Variable_Rho, Variable_U_Force,
                                             Variable_V_Force, Variable_W_Force, Variable_Rho};
    std::vector<std::string> macroVarNames{"rho", "u", "v", "w", "C"};
    std::vector<int> macroVarId{0, 1, 2, 3, 4};
    std::vector<int> macroCompoId{0, 0, 0, 0, 1};
    DefineMacroVars(marcoVarTypes, macroVarNames, macroVarId, macroCompoId);

    std::vector<CollisionType> collisionTypes{Collision_BGKADF,Collision_BGKADG};
    std::vector<int> collisionCompoId{0,1};
    DefineCollision(collisionTypes, collisionCompoId);

    std::vector<BodyForceType> bodyForceTypes{GuoForce,BodyForce_None};
    std::vector<SizeType> bodyForceCompoId{0,1};
    DefineBodyForce(bodyForceTypes, bodyForceCompoId);

    SchemeType scheme{Scheme_StreamCollision};
    DefineScheme(scheme);

    // Setting boundary conditions
    SizeType componentId{0};
    std::vector<VariableTypes> macroVarTypesatBoundary{
        Variable_U_Force, Variable_V_Force, Variable_W_Force};
    std::vector<Real> noSlipStationaryWall{0, 0, 0};

    // Periodic Boundary Conditions for component Fluid 0
    /*
    DefineBlockBoundary(0, componentId, BoundarySurface::Top,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);
    DefineBlockBoundary(0, componentId, BoundarySurface::Bottom,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);
    DefineBlockBoundary(0, componentId, BoundarySurface::Left,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);
    DefineBlockBoundary(0, componentId, BoundarySurface::Right,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);
    DefineBlockBoundary(0, componentId, BoundarySurface::Front,
                        BoundaryScheme::EQMDiffuseReflF, macroVarTypesatBoundary,
                        noSlipStationaryWall);
    DefineBlockBoundary(0, componentId, BoundarySurface::Back,
                        BoundaryScheme::EQMDiffuseReflF, macroVarTypesatBoundary,
                        noSlipStationaryWall);

    DefineBlockBoundary(0, 1, BoundarySurface::Top,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);
    DefineBlockBoundary(0, 1, BoundarySurface::Bottom,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);
    DefineBlockBoundary(0, 1, BoundarySurface::Left,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);
    DefineBlockBoundary(0, 1, BoundarySurface::Right,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);
    DefineBlockBoundary(0, 1, BoundarySurface::Front,
                        BoundaryScheme::EQMDiffuseReflG, macroVarTypesatBoundary,
                        noSlipStationaryWall);
    DefineBlockBoundary(0, 1, BoundarySurface::Back,
                        BoundaryScheme::EQMDiffuseReflG, macroVarTypesatBoundary,
                        noSlipStationaryWall);
    */
    DefineBlockBoundary(0, componentId, BoundarySurface::Front,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);
    DefineBlockBoundary(0, componentId, BoundarySurface::Back,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);
    DefineBlockBoundary(0, 1, BoundarySurface::Front,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);
    DefineBlockBoundary(0, 1, BoundarySurface::Back,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);
    /*
    DefineBlockBoundary(0, componentId, BoundarySurface::Left,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);
    DefineBlockBoundary(0, componentId, BoundarySurface::Right,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);
    DefineBlockBoundary(0, 1, BoundarySurface::Left,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);
    DefineBlockBoundary(0, 1, BoundarySurface::Right,
                        BoundaryScheme::MDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::MDPeriodic);
    */
    DefineBlockBoundary(0, componentId, BoundarySurface::None, BoundaryScheme::EQMDiffuseReflF, macroVarTypesatBoundary,
                        noSlipStationaryWall);
    DefineBlockBoundary(0, 1, BoundarySurface::None, BoundaryScheme::ZouHeVelocity, macroVarTypesatBoundary,
                        noSlipStationaryWall);

    std::vector<InitialType> initType{Initial_BGKFeq2ndAD,Initial_BGKGeq2ndAD};
    std::vector<SizeType> initalCompoId{0,1};
    g_mu().CreateHalos();
    const auto& compoC = g_Components().at(1);
    g_MacroVars().at(compoC.macroVars.at(Variable_Rho).id).CreateHalos();
    DefineInitialCondition(initType,initalCompoId);
    Partition();
    SetSolid3D();
    SetEmbeddedBodyGeometry3D();
    ops_diagnostic_output();
    SetInitialMacrosVars3D();
    
    std::cout << "test\n";
    UpdateConcentration3D();
    CalcPhiWetting3D();
    std::cout << "test1\n";
    PrintPhi3D();
    PreDefinedInitialConditionAD3D();
    UpdateMacroVars3D();
    
    std::cout<<"test2\n";
    SetTimeStep(1);
    PrintPhi3D();
    const Real convergenceCriteria{-1E-7};
    const SizeType checkPeriod{50};
    Real iter{0};
    WriteFlowfieldToHdf5(iter);
    WriteDistributionsToHdf5(iter);
    WriteNodePropertyToHdf5(iter);
    Iterate(StreamCollision,convergenceCriteria, checkPeriod);
}
#endif

void simulate(const Configuration & config, const SizeType timeStep=0) {
    // DefineCase(config.caseName, config.spaceDim);
    // DefineBlocks(config.blockNum, config.blockSize, config.meshSize,
    //                     config.startPos);
    // if (timeStep == 0) {
    //     DefineComponents(config.compoNames, config.compoIds, config.lattNames);
    //     DefineMacroVars(config.macroVarTypes, config.macroVarNames,
    //                     config.macroVarIds, config.macroCompoIds);
    // } else {
    //     // restart from a time step
    //     DefineComponents(config.compoNames, config.compoIds, config.lattNames,
    //                      timeStep);
    //     DefineMacroVars(config.macroVarTypes, config.macroVarNames,
    //                     config.macroVarIds, config.macroCompoIds,timeStep);
    // }

    // DefineCollision(config.CollisionTypes, config.CollisionCompoIds);
    // DefineBodyForce(config.bodyForceTypes, config.bodyForceCompoIds);
    // DefineScheme(config.schemeType);
    // DefineInitialCondition(config.initialTypes,config.initialConditionCompoId);
    // for (auto& bcConfig : config.blockBoundaryConfig) {
    //     DefineBlockBoundary(bcConfig.blockIndex, bcConfig.componentID,
    //                         bcConfig.boundarySurface, bcConfig.boundaryScheme,
    //                         bcConfig.macroVarTypesatBoundary,
    //                         bcConfig.givenVars, bcConfig.boundaryType);
    // }
    // Partition();
    // if (timeStep == 0) {
    //     SetInitialMacrosVars();
    //     PreDefinedInitialCondition3D();
    // } else{
    //     //Help function for restart a steady simulation
    //     //Mainly make the residual calculation correct at first iteration.
    //     RestartMacroVars4SteadySim();
    // }

    // SetTauRef(config.tauRef);
    // SetTimeStep(config.meshSize / SoundSpeed());
    // if (config.transient){
    //     Iterate(config.timeSteps, config.checkPeriod,timeStep);
    // } else{
    //     Iterate(config.convergenceCriteria, config.checkPeriod,timeStep);
    // }
}

int main(int argc, const char** argv) {
    // OPS initialisation
    ops_init(argc, argv, 4);
    double ct0, ct1, et0, et1;
    ops_timers(&ct0, &et0);
    // start a simulation by hard-coding
    if (argc <= 1) {
        simulate();
    }
    // start a new simulaton from a configuration file
    if (argc>1 && argc <=2){
        std::string configFileName(argv[1]);
        ReadConfiguration(configFileName);
        simulate(Config());
    }
    // restart from the time step specified by argv[2]
    if (argc>2 && argc <=3){
        std::string configFileName(argv[1]);
        ReadConfiguration(configFileName);
        const SizeType timeStep{static_cast<SizeType>(std::stoi(argv[2]))};
        simulate(Config(),timeStep);
    }

    ops_timers(&ct1, &et1);
    ops_printf("\nTotal Wall time %lf\n", et1 - et0);
    //Print OPS performance details to output stream
    ops_timing_output(std::cout);
    ops_exit();
}