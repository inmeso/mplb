#include <map>
#include <vector>
#include "flowfield.h"
#include "scheme.h"
#include "ops_seq_v2.h"
#include "scheme_kernel.inc"
#ifdef OPS_3D
void Stream3D() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const SizeType blockIndex{block.ID()};
        ops_par_loop(KerStream3D, "KerStream3D", block.Get(), SpaceDim(),
                     iterRng.data(),
                     ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                 "double", OPS_RW),
                     ops_arg_dat(g_fStage()[blockIndex], NUMXI,
                                 ONEPTLATTICESTENCIL, "double", OPS_READ),
                     ops_arg_dat(g_NodeType()[blockIndex], NUMCOMPONENTS,
                                 LOCALSTENCIL, "int", OPS_READ),
                     ops_arg_dat(g_GeometryProperty()[blockIndex], 1,
                                 LOCALSTENCIL, "int", OPS_READ));
    }
}
#endif // OPS_3D