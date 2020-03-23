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
#ifndef FLOWFIELD_H
#define FLOWFIELD_H
#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <typeinfo>
#include "boundary.h"
#include "model.h"
#include "scheme.h"
#include "type.h"
/*!
 * This module is set for defining blocks and variables defined on a block
 * including distribution functions, macroscopic variables, node properties,
 * and relevant parameters.
 * The responsibilities including:
 * 1. Create all variables from files or annually written subroutines
 * 2. Initialise the required macroscopic variables and thereby the
 *    distribution functions.
 * 3. Provide some tools for accessing variables.
 */
extern int SPACEDIM;
extern ops_block* g_Block;
// Cutting cell
int* IterRngWhole();
int* IterRngJmin();
int* IterRngJmax();
int* IterRngImin();
int* IterRngImax();
int* IterRngBulk();
int* IterRngKmax();
int* IterRngKmin();
/*!
 *Get the pointer pointing to the starting position of IterRng of this block
 *No NULL check for efficiency
 *Note: it looks that ops_par_loop call does not support const point.
 */
inline int* BlockIterRng(const int blockId, int* iterRng) {
    return &iterRng[blockId * 2 * SPACEDIM];
}
/*!
 * Return the starting position of memory in which we store the size of each
 * block
 */
const int* BlockSize(const int blockId);
const int BlockNum();
const int SpaceDim();
const int HaloDepth();
const Real TimeStep();
const Real* pTimeStep();
const Real* TauRef();
// const int* GetBlockNum();
const std::string CaseName();
const int HaloPtNum();
Real TotalMeshSize();
const std::map<std::string,ops_halo_group>& HaloGroups();
void SetTimeStep(Real dt);
void SetCaseName(const std::string caseName);
void setCaseName(const char* caseName);
void SetTauRef(const std::vector<Real> tauRef);
void SetBlockSize(const std::vector<SizeType> blockSize);
void SetBlockNum(const SizeType blockNum);
template <typename T>
class Field {
   private:
    std::map<SizeType, ops_dat> data;
    std::string name;
    SizeType dim{1};
    SizeType haloDepth{1};
    std::string type;

   public:
    Field(std::string& varName, const SizeType dataDim = 1,
          const SizeType halo = 1);

    Field(const char* varName, const SizeType dataDim = 1,
          const SizeType halo = 1);

    void CreateFieldFromScratch();

    void CreateFieldFromCheckPoint(const SizeType timeStep);

    void CreateFieldFromFile(const std::vector<std::string>& fileNames);

    void SetDataDim(const SizeType dataDim) { dim = dataDim; };

    void SetDataHalo(const SizeType halo) { haloDepth = halo; };

    const SizeType HaloDepth() const { return haloDepth; };

    const SizeType DataDim() const { return dim; };

    ~Field(){};

    ops_dat at(SizeType blockIdx) { return data.at(blockIdx); };
    ops_dat operator[](SizeType blockIdx) { return this->at(blockIdx); };
};

template <typename T>
Field<T>::Field(std::string& varName, const SizeType dataDim,
                const SizeType halo) {
    name = varName;
    dim = dataDim;
    haloDepth = halo;
    std::string name{typeid(T).name()};
    if (name == "i") {
        type = "int";
    }
    if (name == "f") {
        type = "float";
    }
    if (name == "d") {
        type = "double";
    }
}
template <typename T>
Field<T>::Field(const char* varName, const SizeType dataDim,
                const SizeType halo) {
    name = std::string{varName};
    dim = dataDim;
    haloDepth = halo;
    std::string name{typeid(T).name()};
    if (name == "i") {
        type = "int";
    }
    if (name == "f") {
        type = "float";
    }
    if (name == "d") {
        type = "double";
    }
}

template <typename T>
void Field<T>::CreateFieldFromScratch() {
    T* temp = NULL;
    int* size = new int[SPACEDIM];
    int* d_p = new int[SPACEDIM];
    int* d_m = new int[SPACEDIM];
    int* base = new int[SPACEDIM];
    for (int cordIdx = 0; cordIdx < SPACEDIM; cordIdx++) {
        d_p[cordIdx] = haloDepth;
        d_m[cordIdx] = -haloDepth;
        base[cordIdx] = 0;
    }
    for (SizeType blockIdx = 0; blockIdx < BlockNum(); blockIdx++) {
        std::string dataName{name + "_" + std::to_string(blockIdx)};
        for (int cordIdx = 0; cordIdx < SPACEDIM; cordIdx++) {
            size[cordIdx] = BlockSize(blockIdx)[cordIdx];
        }
        ops_dat localDat =
            ops_decl_dat(g_Block[blockIdx], dim, size, base, d_m, d_p, temp,
                         type.c_str(), dataName.c_str());
        data.emplace(blockIdx, localDat);
    }
    delete[] size;
    delete[] d_p;
    delete[] d_m;
    delete[] base;
}

template <typename T>
void Field<T>::CreateFieldFromCheckPoint(const SizeType timeStep) {
    for (SizeType blockIdx = 0; blockIdx < BlockNum(); blockIdx++) {
        std::string dataName{name + "_" + std::to_string(blockIdx)};
        std::string label{std::to_string(blockIdx)};
        std::string stepStr{std::to_string(timeStep)};
        std::string blockName{"Block_" + label + "_" + stepStr};
        std::string fileName{CaseName() + "_" + blockName + ".h5"};
        ops_dat localDat =
            ops_decl_dat_hdf5(g_Block[blockIdx], dim, type.c_str(),
                              dataName.c_str(), fileName.c_str());
        data.emplace(blockIdx, localDat);
    }
}
template <typename T>
void Field<T>::CreateFieldFromFile(const std::vector<std::string>& fileNames) {
    for (SizeType blockIdx = 0; blockIdx < BlockNum(); blockIdx++) {
        std::string dataName{name + "_" + std::to_string(blockIdx)};
        ops_dat localDat =
            ops_decl_dat_hdf5(g_Block[blockIdx], dim, type.c_str(),
                              dataName.c_str(), fileNames.at(blockIdx).c_str());
        data.emplace(blockIdx, localDat);
    }
}

/*!
 * The size of g_f in each node will be determined by the employed quadrature
 * and the model. For example, if we are simulating a two-phase flow, then the
 * size will be the product of NUMXI and NUMCOMPONENTS.
 */
extern Field<Real> g_f;
/*! might be changed to a local temporary variable
 * if we use some control routine in the main.cpp
 */
extern Field<Real> g_fStage;
/*!
 * Macroscopic bodyforce
 */
extern Field<Real> g_MacroBodyforce;
/*!
 * g_MacroVars: for storing the macroscopic variables, to reduce
 * the complexity of calculating equilibrium, it will has a specific order
 */
extern Field<Real> g_MacroVars;
/*!
 * Save the macroscopic variables at the previous step
 * Typically used for steady flow.
 */
extern Field<Real> g_MacroVarsCopy;
/*!
 * the residual error for steady flows
 * for each macroscopic variable, there are two values: the absolute
 * and relative
 * for each component of a vector, two values are allocated
 */
extern Real* g_ResidualError;
extern ops_reduction* g_ResidualErrorHandle;
extern Field<int> g_NodeType;
/*!
 * immersed solid? or the end point of the body.
 */
extern Field<int> g_GeometryProperty;
/*!
 * Coordinate
 */
extern Field<Real> g_CoordinateXYZ;

/*!
 * Manually setup the flow field.
 */
void SetupFlowfield();
void SetupFlowfieldfromHdf5();
// TODO This function is temporary, will be removed in the near future
//void AllocateMemory();
void WriteFlowfieldToHdf5(const  SizeType timeStep);
void WriteDistributionsToHdf5(const  SizeType timeStep);
void WriteNodePropertyToHdf5(const  SizeType timeStep);
template <typename T>
//Write a Field variable to HDF5 file
void WriteVariableToHdf5(const SizeType timeStep,Field<T>& var ) {
    for (int blockIndex = 0; blockIndex < BlockNum(); blockIndex++) {
        std::string blockName("Block_");
        std::string label(std::to_string(blockIndex));
        std::string time(std::to_string(timeStep));
        blockName += (label + "_" + time);
        std::string fileName = CaseName() + "_" + blockName + ".h5";
        ops_fetch_block_hdf5_file(g_Block[blockIndex], fileName.c_str());
        ops_fetch_dat_hdf5_file(var[blockIndex], fileName.c_str());
    }
}
void DestroyFlowfield();
//Define the halo pair of implementing periodic boundary condition for
//distribution function itself
void DefinePeriodicHaloPair3D(const std::map<int, std::string>& haloPair);
//Define the halo pair of implementing periodic boundary condition for a
//specfied ops_dat variable
void DefinePeriodicHaloPair3D(const std::map<int, std::string>& haloPair,
                              ops_dat dat, const int haloDepth);
// Transfer all defined halo pair at once
void TransferHalos();
// Transfer the halo pair specified by the key
void TransferHalos(const std::string key);
// Transfer a group of halo pairs specified by the keys
void TransferHalos(const std::vector<std::string> keys);
void SetHaloDepth(const int haloDepth);
void Partition();
void PrepareFlowField();
// caseName: case name
// spaceDim: 2D or 3D application
void DefineCase(const std::string& caseName, const int spaceDim,
                const bool transient=false);
Real GetMaximumResidual(const SizeType checkPeriod);
// blockNum: total number if blocks.
// blockSize: array of integers specifying the block blocksize.
// meshSize: The size of mesh i.e. dx (At present dx = dy = dz).
// startPos: Starting position of each block.
void DefineBlocks(const SizeType blockNum,
                  const std::vector<SizeType>& blockSize, const Real meshSize,
                  const std::vector<Real>& startPos);

bool IsTransient();
template <typename T>
//Define the halo pair of implementing periodic boundary condition for a
//specfied Field variable
void DefinePeriodicHaloPair3D(const std::map<int,std::string>& haloPair,
                              Field<T>& data) {
#ifdef OPS_3D
    DefinePeriodicHaloPair3D(haloPair, data[0], data.HaloDepth());
#endif  // OPS_3D
}
#endif
