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
#include "FreeEnergy_kernel.inc"
//Provide macroscopic initial conditions

//////////////////////////////////////////////////////////////////////
////////////////                  2D                  ////////////////
//////////////////////////////////////////////////////////////////////

#ifdef OPS_2D
//Subroutine to update the order parameter at each lattice point
void UpdateConcentration() {
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        const Real* pdt{pTimeStep()};
        const Real* ttt{g_Time()};
        /*
        auto it = g_Components().begin();
        const auto& Iter1{*it};
        const Component& compoVel{Iter1.second};
        std::advance(it, 1);
        const auto& Iter2{*it};
        */
        auto compoVel = g_Components().at(0);
        auto compoRho = g_Components().at(1);
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
        
        
    }
}

//Subroutine to print a given macroscopic variable at every lattice point
void PrintPhi() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        auto compo = g_Components().at(1);

        ops_par_loop(KerPrintPhi, "KerPrintPhi",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_f()[blockIndex], 1, LOCALSTENCIL,
                                    "double", OPS_READ));
    }
}
//Subroutine to set initial macroscopic variables in the system
void SetInitialMacrosVars() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIdx{block.ID()};
        /*
        auto it = g_Components().begin();
        const auto& Iter{*it};
        */
        auto compo = g_Components().at(0);
        ops_par_loop(KerSetInitialMacroVars, "KerSetInitialMacroVars",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_MacroVars().at(compo.macroVars.at(Variable_Rho).id).at(blockIdx), 1,
                                    LOCALSTENCIL, "Real", OPS_RW),
                        ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIdx),
                                    1, LOCALSTENCIL, "Real", OPS_RW),
                        ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIdx),
                                    1, LOCALSTENCIL, "Real", OPS_RW),
                        ops_arg_dat(g_CoordinateXYZ()[blockIdx], SpaceDim(),
                                    LOCALSTENCIL, "Real", OPS_READ),
                        ops_arg_idx());
        
    }
}


//Provide macroscopic body-force term
void UpdateMacroscopicBodyForce(const Real time) {}
//Subroutine to implement complex geometry
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

//Subroutine to cover solid boundary in surface nodes with normal direction information
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
//Subroutine to update the order parameter at each lattice point
void UpdateConcentration3D() {
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        const Real* pdt{pTimeStep()};
        const Real* ttt{g_Time()};

        /*
        auto it = g_Components().begin();
        const auto& Iter1{*it};
        const Component& compoVel{Iter1.second};
        std::advance(it, 1);
        const auto& Iter2{*it};
        const Component& compoRho{Iter2.second};
        const int compoRhoId{compoRho.id};
        */
        auto compoVel = g_Components().at(0);
        auto compoRho = g_Components().at(1);

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

//Subroutine to print a given macroscopic variable at every lattice point
void PrintPhi3D() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};

        auto compoRho = g_Components().at(1);

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
void UpdateMacroscopicBodyForce(const Real time) {}
//Subroutine to set initial macroscopic variables in the system
void SetInitialMacrosVars3D() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIdx{block.ID()};
        auto it = g_Components().begin();
        const auto& Iter{*it};
        auto compo = g_Components().at(0);
        ops_par_loop(KerSetInitialMacroVars3D, "KerSetInitialMacroVars",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_MacroVars().at(compo.macroVars.at(Variable_Rho).id).at(blockIdx), 1,
                                    LOCALSTENCIL, "Real", OPS_RW),
                        ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIdx),
                                    1, LOCALSTENCIL, "Real", OPS_RW),
                        ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIdx),
                                    1, LOCALSTENCIL, "Real", OPS_RW),
                        ops_arg_dat(g_MacroVars().at(compo.wId).at(blockIdx),
                                    1, LOCALSTENCIL, "Real", OPS_RW),
                        ops_arg_dat(g_CoordinateXYZ()[blockIdx], SpaceDim(),
                                    LOCALSTENCIL, "Real", OPS_READ),
                        ops_arg_idx());
        
    }
}

//Subroutine to implement complex geometry
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

//Subroutine to cover solid boundary in surface nodes with normal direction information
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
    ModelType modeltypes{Free_Energy};
    DefineModelType(modeltypes);
    //Define one block this application
    std::vector<int> blockIds{0};
    std::vector<std::string> blockNames{"Cavity"};
    std::vector<int> blockSize{360, 80};
    Real meshSize{1.};
    std::map<int, std::vector<Real>> startPos{{0, {0.0, 0.0}}};
    DefineBlocks(blockIds, blockNames, blockSize, meshSize, startPos);

    std::vector<int> fromBlockIds{0, 0, 0, 0,0,0,0,0};
    std::vector<int> toBlockIds{0, 0, 0, 0,0,0,0,0};

    std::vector<BoundarySurface> fromBoundarySurface{BoundarySurface::Top,
        BoundarySurface::Bottom,BoundarySurface::Left,
        BoundarySurface::Right,BoundarySurface::LeftTop,
        BoundarySurface::LeftBottom,BoundarySurface::RightTop,
        BoundarySurface::RightBottom};
    std::vector<BoundarySurface> toBoundarySurface{BoundarySurface::Bottom,
        BoundarySurface::Top,BoundarySurface::Right,
        BoundarySurface::Left,BoundarySurface::RightBottom,
        BoundarySurface::RightTop,BoundarySurface::LeftBottom,
        BoundarySurface::LeftTop};

    std::vector<VertexType> blockConnectionType{
        VertexType::MDPeriodic, VertexType::MDPeriodic, VertexType::MDPeriodic, VertexType::MDPeriodic,VertexType::FDPeriodic, VertexType::FDPeriodic, VertexType::FDPeriodic, VertexType::FDPeriodic};

    DefineBlockConnection(fromBlockIds, fromBoundarySurface, toBlockIds,
                          toBoundarySurface, blockConnectionType);


    std::vector<std::string> compoNames{"Fluid","Diffusion"};
    std::vector<int> compoid{0,1};
    std::vector<std::string> lattNames{"d2q9_diffusive","d2q9_diffusive"};
    std::vector<Real> tauRef{1.0,1.0};

    SetFEParams(M_PI/3,0.02,0.02);

    DefineComponents(compoNames, compoid, lattNames, tauRef);

    std::vector<VariableTypes> marcoVarTypes{Variable_Rho, Variable_U_Force,
                                             Variable_V_Force, Variable_Rho};
    std::vector<std::string> macroVarNames{"rho", "u", "v", "C"};
    std::vector<int> macroVarId{0, 1, 2, 3};
    std::vector<int> macroCompoId{0, 0, 0, 1};
    DefineMacroVars(marcoVarTypes, macroVarNames, macroVarId, macroCompoId);

    std::vector<CollisionType> collisionTypes{Collision_BGKFEF,Collision_BGKFEG};
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
        Variable_U_Force, Variable_V_Force};
    std::vector<Real> noSlipStationaryWall{0, 0};

    // Periodic Boundary Conditions

    
    DefineBlockBoundary(0,componentId , BoundarySurface::Top,
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
    /*
    DefineBlockBoundary(0,componentId , BoundarySurface::LeftTop,
                        BoundaryScheme::FDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::FDPeriodic);
    DefineBlockBoundary(0, componentId, BoundarySurface::RightBottom,
                        BoundaryScheme::FDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::FDPeriodic);
    DefineBlockBoundary(0, componentId, BoundarySurface::LeftBottom,
                        BoundaryScheme::FDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::FDPeriodic);
    DefineBlockBoundary(0, componentId, BoundarySurface::RightTop,
                        BoundaryScheme::FDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::FDPeriodic);
    */

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
    /*
    DefineBlockBoundary(0, 1, BoundarySurface::LeftTop,
                        BoundaryScheme::FDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::FDPeriodic);
    DefineBlockBoundary(0, 1, BoundarySurface::RightBottom,
                        BoundaryScheme::FDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::FDPeriodic);
    DefineBlockBoundary(0, 1, BoundarySurface::RightTop,
                        BoundaryScheme::FDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::FDPeriodic);
    DefineBlockBoundary(0, 1, BoundarySurface::LeftBottom,
                        BoundaryScheme::FDPeriodic, macroVarTypesatBoundary,
                        noSlipStationaryWall, VertexType::FDPeriodic);
    */
    
    DefineBlockBoundary(0, componentId, BoundarySurface::None, BoundaryScheme::ZouHeVelocity, macroVarTypesatBoundary,
                        noSlipStationaryWall);
    DefineBlockBoundary(0, 1, BoundarySurface::None, BoundaryScheme::ZouHeVelocityG, macroVarTypesatBoundary,
                        noSlipStationaryWall);

    std::vector<InitialType> initType{Initial_BGKFeq2ndFE,Initial_BGKGeq2ndFE};
    std::vector<SizeType> initalCompoId{0,1};

    g_mu().CreateHalos();
    const auto& compoC = g_Components().at(1);
    g_MacroVars().at(compoC.macroVars.at(Variable_Rho).id).CreateHalos();
    for (auto& pair : g_NodeType()) {
            pair.second.TransferHalos();
    }
    DefineInitialCondition(initType, initalCompoId);

    Partition();
    SetSolid();
    SetEmbeddedBodyGeometry();
    ops_diagnostic_output();

    SetInitialMacrosVars();
    UpdateConcentration();
    CalcPhiWetting();

    PreDefinedInitialConditionAD();
    TransferHalos();
    UpdateMacroVars();
    
    SetTimeStep(1);
    
    const Real convergenceCriteria{-1E-7};
    const SizeType checkPeriod{5000};
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
    ModelType modeltypes{Free_Energy};
    DefineModelType(modeltypes);
    std::vector<int> blockIds{0};
    std::vector<std::string> blockNames{"Cavity"};
    std::vector<int> blockSize{80, 80, 80};
    Real meshSize{1.};
    SetFEParams(2*M_PI/3,0.02,0.02);
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

    std::vector<CollisionType> collisionTypes{Collision_BGKFEF,Collision_BGKFEG};
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
    /*
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
    */
    DefineBlockBoundary(0, componentId, BoundarySurface::None, BoundaryScheme::ZouHeVelocity, macroVarTypesatBoundary,
                        noSlipStationaryWall);
    DefineBlockBoundary(0, 1, BoundarySurface::None, BoundaryScheme::ZouHeVelocityG, macroVarTypesatBoundary,
                        noSlipStationaryWall);

    std::vector<InitialType> initType{Initial_BGKFeq2ndFE,Initial_BGKGeq2ndFE};
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
    

    UpdateConcentration3D();
    CalcPhiWetting3D();

    PreDefinedInitialConditionAD3D();
    UpdateMacroVars3D();

    SetTimeStep(1);

    const Real convergenceCriteria{-1E-7};
    const SizeType checkPeriod{1};
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