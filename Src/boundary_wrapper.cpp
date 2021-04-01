#include <vector>
#include "type.h"
#include "flowfield.h"
#include "flowfield_host_device.h"
#include "boundary.h"
#include "model.h"
#include "boundary_host_device.h"
#include "scheme.h"
#include "ops_seq_v2.h"
#include "boundary_kernel.inc"
#ifdef OPS_3D
void TreatBlockBoundary3D(Block& block, const int componentID,
                          const Real* givenVars, int* range,
                          const BoundaryScheme boundaryScheme,
                          const BoundarySurface boundarySurface) {
    const int surface{(int)boundarySurface};
    const SizeType blockIndex{block.ID()};
    switch (boundaryScheme) {
        case BoundaryScheme::ExtrapolPressure1ST: {
            ops_par_loop(
                KerCutCellExtrapolPressure1ST3D,
                "KerCutCellExtrapolPressure1ST3D", block.Get(), SpaceDim(),
                range,
                ops_arg_dat(g_f()[blockIndex], NUMXI, ONEPTREGULARSTENCIL,
                            "double", OPS_RW),
                ops_arg_dat(g_NodeType()[blockIndex], NUMCOMPONENTS,
                            ONEPTREGULARSTENCIL, "int", OPS_READ),
                ops_arg_dat(g_GeometryProperty()[blockIndex], 1, LOCALSTENCIL,
                            "int", OPS_READ),
                ops_arg_gbl(givenVars, NUMMACROVAR, "double", OPS_READ),
                ops_arg_gbl(&componentID, 1, "int", OPS_READ),
                ops_arg_gbl(&surface, 1, "int", OPS_READ));
        } break;
        case BoundaryScheme::EQMDiffuseRefl: {
            ops_par_loop(
                KerCutCellEQMDiffuseRefl3D, "KerCutCellEQMDiffuseRefl3D",
                block.Get(), SpaceDim(), range,
                ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL, "double",
                            OPS_RW),
                ops_arg_dat(g_NodeType()[blockIndex], NUMCOMPONENTS,
                            LOCALSTENCIL, "int", OPS_READ),
                ops_arg_dat(g_GeometryProperty()[blockIndex], 1, LOCALSTENCIL,
                            "int", OPS_READ),
                ops_arg_gbl(givenVars, NUMMACROVAR, "double", OPS_READ),
                ops_arg_gbl(&componentID, 1, "int", OPS_READ));
        } break;
        case BoundaryScheme::EQN: {
            ops_par_loop(
                KerCutCellNoslipEQN3D, "KerCutCellNoslipEQN3D", block.Get(),
                SpaceDim(), range,
                ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL, "double",
                            OPS_RW),
                ops_arg_dat(g_NodeType()[blockIndex], NUMCOMPONENTS,
                            LOCALSTENCIL, "int", OPS_READ),
                ops_arg_gbl(givenVars, NUMMACROVAR, "double", OPS_READ),
                ops_arg_gbl(&componentID, 1, "int", OPS_READ));
        } break;
        case BoundaryScheme::FDPeriodic: {
            ops_par_loop(KerCutCellPeriodic3D, "KerCutCellPeriodic3D",
                         block.Get(), SpaceDim(), range,
                         ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                     "double", OPS_RW),
                         ops_arg_dat(g_NodeType()[blockIndex], NUMCOMPONENTS,
                                     LOCALSTENCIL, "int", OPS_READ),
                         ops_arg_dat(g_GeometryProperty()[blockIndex], 1,
                                     LOCALSTENCIL, "int", OPS_READ),
                         ops_arg_gbl(&componentID, 1, "int", OPS_READ),
                         ops_arg_gbl(&surface, 1, "int", OPS_READ));
        } break;
        default:
            break;
    }
}
#endif //OPS_3D
