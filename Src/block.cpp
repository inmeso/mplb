#include "block.h"
#include <string>
#include <vector>
#include <cassert>
#include "ops_lib_core.h"
#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif

Block::Block(const SizeType blockId, const std::string& blockName,
             const std::vector<int>& blockSize) {
    id = blockId;
    name = blockName;
    if (spaceDim != blockSize.size()) {
        ops_printf(
            "Error! The SPACEDIM here is inconsistent with the block size\n");
        assert(spaceDim == blockSize.size());
    }
    size = blockSize;
    block = ops_decl_block(spaceDim, blockName.c_str());
    wholeRange.resize(2 * spaceDim);
    wholeRange.at(0) = 0;
    wholeRange.at(1) = size.at(0);
    wholeRange.at(2) = 0;
    wholeRange.at(3) = size.at(1);
#ifdef OPS_3D
    wholeRange.at(4) = 0;
    wholeRange.at(5) = size.at(2);
#endif
    bulkRange.resize(2 * spaceDim);
    bulkRange.at(0) = 1;
    bulkRange.at(1) = size.at(0)-1;
    bulkRange.at(2) = 1;
    bulkRange.at(3) = size.at(1)-1;
#ifdef OPS_3D
    bulkRange.at(4) = 1;
    bulkRange.at(5) = size.at(2)-1;
#endif
    iminRange.resize(2 * spaceDim);
    iminRange.at(0) = 0;
    iminRange.at(1) = 1;
    iminRange.at(2) = 0;
    iminRange.at(3) = size.at(1);
#ifdef OPS_3D
    iminRange.at(4) = 0;
    iminRange.at(5) = size.at(2);
#endif
    imaxRange.resize(2 * spaceDim);
    imaxRange.at(0) = size.at(0) - 1;
    imaxRange.at(1) = size.at(0);
    imaxRange.at(2) = 0;
    imaxRange.at(3) = size.at(1);
#ifdef OPS_3D
    imaxRange.at(4) = 0;
    imaxRange.at(5) = size.at(2);
#endif
    jminRange.resize(2 * spaceDim);
    jminRange.at(0) = 0;
    jminRange.at(1) = size.at(0);
    jminRange.at(2) = 0;
    jminRange.at(3) = 1;
#ifdef OPS_3D
    jminRange.at(4) = 0;
    jminRange.at(5) = size.at(2);
#endif
    jmaxRange.resize(2 * spaceDim);
    jmaxRange.at(0) = 0;
    jmaxRange.at(1) = size.at(0);
    jmaxRange.at(2) = size.at(1)-1;
    jmaxRange.at(3) = size.at(1);
#ifdef OPS_3D
    jmaxRange.at(4) = 0;
    jmaxRange.at(5) = size.at(2);
#endif
#ifdef OPS_3D
    kminRange.resize(2 * spaceDim);
    kminRange.at(0) = 0;
    kminRange.at(1) = size.at(0);
    kminRange.at(2) = 0;
    kminRange.at(3) = size.at(1);
    kminRange.at(4) = 0;
    kminRange.at(5) = 1;
    kmaxRange.resize(2 * spaceDim);
    kmaxRange.at(0) = 0;
    kmaxRange.at(1) = size.at(0);
    kmaxRange.at(2) = 0;
    kmaxRange.at(3) = size.at(1);
    kmaxRange.at(4) = size.at(2) - 1;
    kmaxRange.at(5) = size.at(2);
#endif
}