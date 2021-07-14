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
#ifndef FIELD_H
#define FIELD_H
#include <list>
#include <map>
#include <string>
#include <vector>
#include "block.h"
#include "ops_lib_core.h"
#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif
#include "type.h"
template <typename T>
class Field {
   private:
    std::map<int, ops_dat> data;
    std::map<int, const Block&> dataBlock;
    std::string name;
    int dim{1};
    int haloDepth{1};
    ops_halo_group haloGroup{nullptr};
#ifdef OPS_3D
    int spaceDim{3};
#endif
#ifdef OPS_2D
    int spaceDim{2};
#endif
    std::string type;

   public:
    Field(const std::string& varName, const int dataDim = 1,
          const int halo = 1);
    Field(const char* varName, const int dataDim = 1,
          const int halo = 1);
    void CreateFieldFromScratch(const BlockGroup& blocks);
    void CreateFieldFromScratch(const Block& block);
    void CreateFieldFromFile(const std::string& fileName, const Block& block);
    void CreateFieldFromFile(const std::string& caseName, const Block& block,
                             const SizeType timeStep);
    void CreateFieldFromFile(const std::string& caseName,
                             const BlockGroup& blocks,
                             const SizeType timeStep);
    void SetDataDim(const int dataDim) { dim = dataDim; };
    void SetDataHalo(const int halo) { haloDepth = halo; };
    void WriteToHDF5(const std::string& caseName, const SizeType timeStep) const;
    int HaloDepth() const { return haloDepth; };
    int DataDim() const { return dim; };
    ~Field(){};
    ops_dat& at(int blockIdx) { return data.at(blockIdx); };
    const ops_dat& at(int blockIdx) const { return data.at(blockIdx); };
    ops_dat& operator[](int blockIdx) { return this->at(blockIdx); };
    const ops_dat& operator[](int blockIdx) const { return this->at(blockIdx); };
    void CreateHalos();
    void TransferHalos();
};
template <typename T>
void Field<T>::TransferHalos() {
    if (haloGroup != nullptr) {
        ops_halo_transfer(haloGroup);
    }
};
template <typename T>
Field<T>::Field(const std::string& varName, const int dataDim, const int halo) {
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
Field<T>::Field(const char* varName, const int dataDim,
                const int halo) {
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
void Field<T>::CreateFieldFromScratch(const Block& block) {
    T* temp{nullptr};
    int* d_p = new int[spaceDim];
    int* d_m = new int[spaceDim];
    int* base = new int[spaceDim];
    for (int cordIdx = 0; cordIdx < spaceDim; cordIdx++) {
        d_p[cordIdx] = haloDepth;
        d_m[cordIdx] = -haloDepth;
        base[cordIdx] = 0;
    }
    const int blockId{block.ID()};
    std::string dataName{name + "_" + block.Name()};
    std::vector<int> size{block.Size()};
    ops_dat localDat =
        ops_decl_dat(block.Get(), dim, size.data(), base,
                     d_m, d_p, temp, type.c_str(), dataName.c_str());
    data.emplace(blockId, localDat);
    dataBlock.emplace(blockId, block);
    delete[] d_p;
    delete[] d_m;
    delete[] base;
}

template <typename T>
void Field<T>::CreateFieldFromScratch(const BlockGroup& blocks) {
    for (const auto& idBlock : blocks) {
        const Block& block{idBlock.second};
        CreateFieldFromScratch(block);
    }
}

template <typename T>
void Field<T>::CreateFieldFromFile(const std::string& fileName,
                                   const Block& block) {
    std::string dataName{name + "_" + block.Name()};
    ops_dat localDat = ops_decl_dat_hdf5(block.Get(), dim, type.c_str(),
                                         dataName.c_str(), fileName.c_str());
    data.emplace(block.ID(), localDat);
    dataBlock.emplace(block.ID(), block);
}

template <typename T>
void Field<T>::CreateFieldFromFile(const std::string& caseName,
                                   const Block& block,
                                   const SizeType timeStep) {
    std::string fileName{caseName + "_" + block.Name() + "_T" +
                         std::to_string(timeStep) + ".h5"};
    CreateFieldFromFile(fileName, block);
}

template <typename T>
void Field<T>::CreateFieldFromFile(const std::string& caseName,
                                   const BlockGroup& blocks,
                                   const SizeType timeStep) {
    for (const auto& idBlock : blocks) {
        const Block& block{idBlock.second};
        CreateFieldFromFile(caseName, block, timeStep);
    }
}
template <typename T>
void Field<T>::WriteToHDF5(const std::string& caseName,
                           const SizeType timeStep) const {
    for (const auto& idData : data) {
        const int blockId{idData.first};
        const Block& block{dataBlock.at(blockId)};
        std::string fileName = caseName + "_" + block.Name() + "_T" +
                               std::to_string(timeStep) + ".h5";
        ops_fetch_block_hdf5_file(block.Get(), fileName.c_str());
        ops_fetch_dat_hdf5_file(idData.second, fileName.c_str());
    }
}
/**
 * @brief This method creates all halos for communicating between blocks.
 *
 * The method works under the assumption that blocks are connecteed in an exact
 * point-to-point fashine without rotating coordinates.
 *
 * The periodic boundary can be treated by setting the neighbor to the
 * block itself.
 */
template <typename T>
void Field<T>::CreateHalos() {
    std::vector<ops_halo> halos;
    for (const auto& idBlock : dataBlock) {
        const int id{idBlock.first};
        const Block& block{idBlock.second};
        const int nx{(int)block.Size().at(0)};
        const int ny{(int)block.Size().at(1)};
#ifdef OPS_3D
        int nz{(int)block.Size().at(2)};
        const int d_p[]{haloDepth, haloDepth, haloDepth};
        const int d_m[]{-haloDepth, -haloDepth, -haloDepth};
        int dir[]{1, 2, 3};
#endif
#ifdef OPS_2D
        const int d_p[]{haloDepth, haloDepth};
        const int d_m[]{-haloDepth, -haloDepth};
        int dir[]{1, 2};
#endif
        for (const auto& surfaceNeighbor : block.Neighbors()) {
            const Neighbor& neighbor{surfaceNeighbor.second};
            const BoundarySurface surface{surfaceNeighbor.first};
            const VertexType type{neighbor.type};
            switch (surface) {
                case BoundarySurface::Right: {
                    int disp{0};
                    if (type == VertexType::FDPeriodic ||
                        type == VertexType::MDPeriodic) {
                        disp = d_m[0];
                    }
                    if (type == VertexType::VirtualBoundary) {
                        disp = d_m[0] - 1;
                    }
#ifdef OPS_3D
                    int halo_iter[] = {haloDepth, ny + d_p[1] - d_m[1],
                                       nz + d_p[2] - d_m[2]};
                    int base_from[] = {nx + disp, d_m[1], d_m[2]};
                    int base_to[] = {d_m[0], d_m[1], d_m[2]};
#endif
#ifdef OPS_2D
                    int halo_iter[] = {haloDepth, ny + d_p[1] - d_m[1]};
                    int base_from[] = {nx + disp, d_m[1]};
                    int base_to[] = {d_m[0], d_m[1]};
#endif
                    ops_halo rightToLeft = ops_decl_halo(
                        data.at(id), data.at(neighbor.blockId),
                        halo_iter, base_from, base_to, dir, dir);
                    halos.push_back(rightToLeft);
                } break;
                case BoundarySurface::Left: {
                    const int neighborBase{
                        dataBlock.at(neighbor.blockId).Size().at(0)};
                    int disp{0};
                    if (type == VertexType::FDPeriodic ||
                        type == VertexType::MDPeriodic) {
                        disp = 0;
                    }
                    if (type == VertexType::VirtualBoundary) {
                        disp = 1;
                    }
#ifdef OPS_3D
                    int halo_iter[] = {haloDepth, ny + d_p[1] - d_m[1],
                                       nz + d_p[2] - d_m[2]};
                    int base_from[] = {disp, d_m[1], d_m[2]};
                    int base_to[] = {neighborBase, d_m[1], d_m[2]};
#endif
#ifdef OPS_2D
                    int halo_iter[] = {haloDepth, ny + d_p[1] - d_m[1]};
                    int base_from[] = {disp, d_m[1]};
                    int base_to[] = {neighborBase, d_m[1]};
#endif
                    ops_halo leftToRight = ops_decl_halo(
                        data.at(id), data.at(neighbor.blockId),
                        halo_iter, base_from, base_to, dir, dir);
                    halos.push_back(leftToRight);
                } break;
                case BoundarySurface::Bottom: {
                    const int neighborBase{
                        dataBlock.at(neighbor.blockId).Size().at(1)};
                    int disp{0};
                    if (type == VertexType::FDPeriodic ||
                        type == VertexType::MDPeriodic) {
                        disp = 0;
                    }
                    if (type == VertexType::VirtualBoundary) {
                        disp = 1;
                    }
#ifdef OPS_3D
                    int halo_iter[] = {nx + d_p[0] - d_m[0], haloDepth,
                                       nz + d_p[2] - d_m[2]};
                    int base_from[] = {d_m[0], disp, d_m[2]};
                    int base_to[] = {d_m[0], neighborBase, d_m[2]};
#endif
#ifdef OPS_2D
                    int halo_iter[] = {nx + d_p[0] - d_m[0], haloDepth};
                    int base_from[] = {d_m[0], disp};
                    int base_to[] = {d_m[0], neighborBase};
#endif
                    ops_halo botToTop = ops_decl_halo(
                        data.at(id), data.at(neighbor.blockId),
                        halo_iter, base_from, base_to, dir, dir);
                    halos.push_back(botToTop);
                } break;

                case BoundarySurface::Top: {
                    int disp{0};
                    if (type == VertexType::FDPeriodic ||
                        type == VertexType::MDPeriodic) {
                        disp = d_m[1];
                    }
                    if (type == VertexType::VirtualBoundary) {
                        disp = d_m[1] - 1;
                    }
#ifdef OPS_3D
                    int halo_iter[] = {nx + d_p[0] - d_m[0], haloDepth,
                                       nz + d_p[2] - d_m[2]};
                    int base_from[] = {d_m[0], ny + disp, d_m[2]};
                    int base_to[] = {d_m[0], d_m[1], d_m[2]};
#endif
#ifdef OPS_2D
                    int halo_iter[] = {nx + d_p[0] - d_m[0], haloDepth};
                    int base_from[] = {d_m[0], ny + disp};
                    int base_to[] = {d_m[0], d_m[1]};
#endif
                    ops_halo topToBot = ops_decl_halo(
                        data.at(id), data.at(neighbor.blockId),
                        halo_iter, base_from, base_to, dir, dir);
                    halos.push_back(topToBot);
                } break;
#ifdef OPS_3D
                case BoundarySurface::Back: {
                    const int neighborBase{
                        dataBlock.at(neighbor.blockId).Size().at(2)};
                    int disp{0};
                    if (type == VertexType::FDPeriodic ||
                        type == VertexType::MDPeriodic) {
                        disp = 0;
                    }
                    if (type == VertexType::VirtualBoundary) {
                        disp = 1;
                    }
                    int halo_iter[] = {nx + d_p[0] - d_m[0],
                                       ny + d_p[1] - d_m[1], haloDepth};
                    int base_from[] = {d_m[0], d_m[1], disp};
                    int base_to[] = {d_m[0], d_m[1], neighborBase};

                    ops_halo backToFront = ops_decl_halo(
                        data.at(id), data.at(neighbor.blockId),
                        halo_iter, base_from, base_to, dir, dir);
                    halos.push_back(backToFront);
                } break;

                case BoundarySurface::Front: {
                    int disp{0};
                    if (type == VertexType::FDPeriodic ||
                        type == VertexType::MDPeriodic) {
                        disp = d_m[2];
                    }
                    if (type == VertexType::VirtualBoundary) {
                        disp = d_m[2] - 1;
                    }
                    int halo_iter[] = {nx + d_p[0] - d_m[0],
                                       ny + d_p[1] - d_m[1], haloDepth};
                    int base_from[] = {d_m[0], d_m[1], nz + disp};
                    int base_to[] = {d_m[0], d_m[1], d_m[2]};

                    ops_halo frontToBack = ops_decl_halo(
                        data.at(id), data.at(neighbor.blockId),
                        halo_iter, base_from, base_to, dir, dir);
                    halos.push_back(frontToBack);
                } break;
#endif
                default:
                    break;
            }
        }
    }
    if (halos.size()>=1){
        haloGroup = ops_decl_halo_group(halos.size(), halos.data());
    }
}

using RealField = Field<Real>;
using IntField = Field<int>;
using IntFieldGroup = std::map<int, IntField>;
using RealFieldGroup = std::map<int, RealField>;
#endif