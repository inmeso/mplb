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

#ifndef BLOCK_H
#define BLOCK_H
#include <map>
#include <string>
#include <vector>
#include <cassert>

#include "ops_seq_v2.h"
#include "type.h"
class Block {
   private:
    ops_block block;
#ifdef OPS_3D
    const SizeType spaceDim{3};
#endif
#ifdef OPS_2D
    const SizeType spaceDim{2};
#endif
    std::string name;
    SizeType id;
    std::vector<SizeType> size;
    std::vector<SizeType> wholeRange;
    std::vector<SizeType> bulkRange;
    std::vector<SizeType> iminRange;
    std::vector<SizeType> imaxRange;
    std::vector<SizeType> jminRange;
    std::vector<SizeType> jmaxRange;
#ifdef OPS_3D
    std::vector<SizeType> kminRange;
    std::vector<SizeType> kmaxRange;
#endif
   public:
    Block(const SizeType blockId, const std::string& blockName,
          const std::vector<SizeType> & blockSize);
    SizeType ID() const { return id; };
    const std::string& Name() const { return name; };
    const ops_block& Get() const { return block; };
    ops_block& Get() { return block; };
    const std::vector<SizeType>& Size() const { return size; };
    const SizeType* pSize() const { return size.data(); };
    const std::vector<SizeType>& WholeRange() const { return wholeRange; };
    const std::vector<SizeType>& BulkRange() const { return bulkRange; };
    const std::vector<SizeType>& ImaxRange() const { return imaxRange; };
    const std::vector<SizeType>& IminRange() const { return iminRange; };
    const std::vector<SizeType>& JmaxRange() const { return jmaxRange; };
    const std::vector<SizeType>& JminRange() const { return jminRange; };
#ifdef OPS_3D
    const std::vector<SizeType>& KmaxRange() const { return kmaxRange; };
    const std::vector<SizeType>& KminRange() const { return kminRange; };
#endif
};

typedef std::map<SizeType, Block> BlockGroup;

Block::Block(const SizeType blockId, const std::string& blockName,
             const std::vector<SizeType>& blockSize) {
    id = blockId;
    name = blockName;
    if (spaceDim != blockSize.size()) {
        ops_printf(
            "Error! The SPACEDIM here is inconsistent with the block size\n");
        assert(spaceDim == blockSize.size());
    }
    size = blockSize;
    block = ops_decl_block(spaceDim, blockName.c_str());
    wholeRange.reserve(2 * spaceDim);
    wholeRange.at(0) = 0;
    wholeRange.at(1) = size.at(0);
    wholeRange.at(2) = 0;
    wholeRange.at(3) = size.at(1);
#ifdef OPS_3D
    wholeRange.at(4) = 0;
    wholeRange.at(5) = size.at(2);
#endif
    bulkRange.reserve(2 * spaceDim);
    bulkRange.at(0) = 1;
    bulkRange.at(1) = size.at(0)-1;
    bulkRange.at(2) = 1;
    bulkRange.at(3) = size.at(1)-1;
#ifdef OPS_3D
    bulkRange.at(4) = 1;
    bulkRange.at(5) = size.at(2)-1;
#endif
    iminRange.reserve(2 * spaceDim);
    iminRange.at(0) = 0;
    iminRange.at(1) = 1;
    iminRange.at(2) = 0;
    iminRange.at(3) = size.at(1);
#ifdef OPS_3D
    iminRange.at(4) = 0;
    iminRange.at(5) = size.at(2);
#endif
    imaxRange.reserve(2 * spaceDim);
    imaxRange.at(0) = size.at(0) - 1;
    imaxRange.at(1) = size.at(0);
    imaxRange.at(2) = 0;
    imaxRange.at(3) = size.at(1);
#ifdef OPS_3D
    imaxRange.at(4) = 0;
    imaxRange.at(5) = size.at(2);
#endif
    jminRange.reserve(2 * spaceDim);
    jminRange.at(0) = 0;
    jminRange.at(1) = size.at(0);
    jminRange.at(2) = 0;
    jminRange.at(3) = 1;
#ifdef OPS_3D
    jminRange.at(4) = 0;
    jminRange.at(5) = size.at(2);
#endif
    jmaxRange.reserve(2 * spaceDim);
    jmaxRange.at(0) = 0;
    jmaxRange.at(1) = size.at(0);
    jmaxRange.at(2) = size.at(1)-1;
    jmaxRange.at(3) = size.at(1);
#ifdef OPS_3D
    jmaxRange.at(4) = 0;
    jmaxRange.at(5) = size.at(2);
#endif
#ifdef OPS_3D
    kminRange.reserve(2 * spaceDim);
    kminRange.at(0) = 0;
    kminRange.at(1) = size.at(0);
    kminRange.at(2) = 0;
    kminRange.at(3) = size.at(1);
    kminRange.at(4) = 0;
    kminRange.at(5) = 1;
    kmaxRange.reserve(2 * spaceDim);
    kmaxRange.at(0) = 0;
    kmaxRange.at(1) = size.at(0);
    kmaxRange.at(2) = 0;
    kmaxRange.at(3) = size.at(1);
    kmaxRange.at(4) = size.at(2) - 1;
    kmaxRange.at(5) = size.at(2);
#endif
}

#endif  // BLOCK_H