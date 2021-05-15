#include "block.h"
#include <string>
#include <vector>
#include <map>
#include <cassert>
#include "ops_lib_core.h"
#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif
#include "type.h"
#include "boundary_host_device.h"
#include "boundary.h"

int Block::RangeStart(const int axis) { return 0; }
int Block::RangeEnd(const int axis) { return size.at(axis); }

int Block::RangeStart(const int axis, const BoundarySurface surface) {
    int start{RangeStart(axis)};
    if (axis == xaxis) {
        if (surface == BoundarySurface::Right ||
            surface == BoundarySurface::RightTop ||
            surface == BoundarySurface::RightBottom
#ifdef OPS_3D
            || surface == BoundarySurface::RightFront ||
            surface == BoundarySurface::RightBack
#endif
        ) {
            start = size.at(axis) - 1;
        }
    }
    if (axis == yaxis) {
        if (surface == BoundarySurface::Top ||
            surface == BoundarySurface::LeftTop ||
            surface == BoundarySurface::RightTop
#ifdef OP_3D
            || surface == BoundarySurface::TopBack ||
            surface == BoundarySurface::TopFront
#endif

        ) {
            start = size.at(axis) - 1;
        }
    }
#ifdef OPS_3D
    if (axis == zaxis) {
        if (surface == BoundarySurface::Front ||
            surface == BoundarySurface::LeftFront ||
            surface == BoundarySurface::RightFront ||
            surface == BoundarySurface::TopFront ||
            surface == BoundarySurface::BottomFront) {
            start = size.at(axis) - 1;
        }
    }
#endif
    return start;
}

int Block::RangeEnd(const int axis, const BoundarySurface surface) {
    int end{RangeEnd(axis)};
    if (axis == xaxis) {
        if (surface == BoundarySurface::Left ||
            surface == BoundarySurface::LeftBottom ||
            surface == BoundarySurface::LeftTop
#ifdef OPS_3D
            || surface == BoundarySurface::LeftFront ||
            surface == BoundarySurface::LeftBack
#endif
        ) {
            end = 1;
        }
    }
    if (axis == yaxis) {
        if (surface == BoundarySurface::Bottom ||
            surface == BoundarySurface::LeftBottom ||
            surface == BoundarySurface::RightBottom
#ifdef OPS_3D

            || surface == BoundarySurface::BottomBack ||
            surface == BoundarySurface::BottomFront
#endif
        ) {
            end = 1;
        }
    }
#ifdef OPS_3D
    if (axis == zaxis) {
        if (surface == BoundarySurface::Back ||
            surface == BoundarySurface::LeftBack ||
            surface == BoundarySurface::RightBack ||
            surface == BoundarySurface::TopBack ||
            surface == BoundarySurface::BottomBack) {
            end = 1;
        }
    }
#endif
    return end;
}

Block::Block(const int blockId, const std::string& blockName,
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
    bulkRange.at(1) = size.at(0) - 1;
    bulkRange.at(2) = 1;
    bulkRange.at(3) = size.at(1) - 1;
#ifdef OPS_3D
    bulkRange.at(4) = 1;
    bulkRange.at(5) = size.at(2) - 1;
#endif
    for (const auto surface : AllBoundarySurface) {
        boundarySurfaceRange[surface] = {RangeStart(xaxis, surface),
                                         RangeEnd(xaxis, surface),
                                         RangeStart(yaxis, surface),
                                         RangeEnd(yaxis, surface)
#ifdef OPS_3D
                                             ,
                                         RangeStart(zaxis, surface),
                                         RangeEnd(zaxis, surface)
#endif
        };
    }
}

void Block::AddNeighbor(BoundarySurface surface, const Neighbor& neighbor) {
    if (neighbors.find(surface) != neighbors.end()) {
        ops_printf("Error! There is already a neighbor defined for Block $s!\n",
                   name.c_str());
        assert(neighbors.find(surface) == neighbors.end());
    }
    neighbors.emplace(surface, neighbor);
}