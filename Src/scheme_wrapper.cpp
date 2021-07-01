#include <map>
#include <vector>
#include "flowfield.h"
#include "scheme.h"
#include "ops_seq_v2.h"
#include "scheme_kernel.inc"
#ifdef OPS_3D
void Stream3D() {
#ifdef OPS_3D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        for (const auto& compo : g_Components()) {
            ops_par_loop(
                KerStream3D, "KerStream3D", block.Get(), SpaceDim(),
                iterRng.data(),
                ops_arg_dat(g_f().at(blockIndex), NUMXI, LOCALSTENCIL, "double",
                            OPS_RW),
                ops_arg_dat(g_fStage().at(blockIndex), NUMXI,
                            ONEPTLATTICESTENCIL, "double", OPS_READ),
                ops_arg_dat(g_NodeType().at(compo.first).at(blockIndex), 1,
                            LOCALSTENCIL, "int", OPS_READ),
                ops_arg_dat(g_GeometryProperty().at(blockIndex), 1,
                            LOCALSTENCIL, "int", OPS_READ),
                ops_arg_gbl(compo.second.index, 2, "int", OPS_READ));
        }
    }
#endif // OPS_3Ds
}
#endif  // OPS_3D

#ifdef OPS_2D
void Stream() {
#ifdef OPS_2D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        for (const auto& compo : g_Components()) {
            ops_par_loop(
                KerStream, "KerStream", block.Get(), SpaceDim(),
                iterRng.data(),
                ops_arg_dat(g_f().at(blockIndex), NUMXI, LOCALSTENCIL, "double",
                            OPS_RW),
                ops_arg_dat(g_fStage().at(blockIndex), NUMXI,
                            ONEPTLATTICESTENCIL, "double", OPS_READ),
                ops_arg_dat(g_NodeType().at(compo.first).at(blockIndex), 1,
                            LOCALSTENCIL, "int", OPS_READ),
                ops_arg_dat(g_GeometryProperty().at(blockIndex), 1,
                            LOCALSTENCIL, "int", OPS_READ),
                ops_arg_gbl(compo.second.index, 2, "int", OPS_READ));
        }
    }
#endif // OPS_2D
}
#endif  // OPS_2D