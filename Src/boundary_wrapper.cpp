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
void TreatBlockBoundary3D(const Block& block, const int componentID,
                          const Real* givenVars,
                          const BoundaryScheme boundaryScheme,
                          const BoundarySurface boundarySurface) {
    const int surface{(int)boundarySurface};
    const int blockIndex{block.ID()};
    std::vector<int> range(2 * SpaceDim());
    range.assign(block.BoundarySurfaceRange().at(boundarySurface).begin(),
                     block.BoundarySurfaceRange().at(boundarySurface).end());
    switch (boundaryScheme) {
        case BoundaryScheme::ExtrapolPressure1ST: {
            ops_par_loop(
                KerCutCellExtrapolPressure1ST3D,
                "KerCutCellExtrapolPressure1ST3D", block.Get(), SpaceDim(),
                range.data(),
                ops_arg_dat(g_f()[blockIndex], NUMXI, ONEPTREGULARSTENCIL,
                            "double", OPS_RW),
                ops_arg_dat(g_NodeType().at(componentID).at(blockIndex), 1,
                            ONEPTREGULARSTENCIL, "int", OPS_READ),
                ops_arg_dat(g_GeometryProperty()[blockIndex], 1, LOCALSTENCIL,
                            "int", OPS_READ),
                ops_arg_gbl(givenVars, 1, "double", OPS_READ),
                ops_arg_gbl(&surface, 1, "int", OPS_READ),
                ops_arg_gbl(g_Components().at(componentID).index, 2, "int",
                            OPS_READ));
        } break;
        case BoundaryScheme::EQMDiffuseRefl: {
            ops_par_loop(
                KerCutCellEQMDiffuseRefl3D, "KerCutCellEQMDiffuseRefl3D",
                block.Get(), SpaceDim(), range.data(),
                ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL, "double",
                            OPS_RW),
                ops_arg_dat(g_NodeType().at(componentID).at(blockIndex), 1,
                            LOCALSTENCIL, "int", OPS_READ),
                ops_arg_dat(g_GeometryProperty()[blockIndex], 1, LOCALSTENCIL,
                            "int", OPS_READ),
                ops_arg_gbl(givenVars, 3, "double", OPS_READ),
                ops_arg_gbl(g_Components().at(componentID).index, 2, "int",
                            OPS_READ));
        } break;
        case BoundaryScheme::FDPeriodic: {
            ops_par_loop(
                KerCutCellPeriodic3D, "KerCutCellPeriodic3D", block.Get(),
                SpaceDim(), range.data(),
                ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL, "double",
                            OPS_RW),
                ops_arg_dat(g_NodeType().at(componentID).at(blockIndex), 1,
                            LOCALSTENCIL, "int", OPS_READ),
                ops_arg_dat(g_GeometryProperty()[blockIndex], 1, LOCALSTENCIL,
                            "int", OPS_READ),
                ops_arg_gbl(g_Components().at(componentID).index, 2, "int",
                            OPS_READ),
                ops_arg_gbl(&surface, 1, "int", OPS_READ));
        } break;
        default:
            break;
    }
}
#endif //OPS_3D

#ifdef OPS_2D
void TreatBlockBoundary(const Block& block, const int componentID,
                          const Real* givenVars,
                          const BoundaryScheme boundaryScheme,
                          const BoundarySurface boundarySurface) {
    const int surface{(int)boundarySurface};
    const int blockIndex{block.ID()};
    std::vector<int> range(2 * SpaceDim());
    range.assign(block.BoundarySurfaceRange().at(boundarySurface).begin(),
                     block.BoundarySurfaceRange().at(boundarySurface).end());
    switch (boundaryScheme) {
        case BoundaryScheme::ExtrapolPressure1ST: {
            ops_par_loop(
                KerCutCellExtrapolPressure1ST,
                "KerCutCellExtrapolPressure1ST", block.Get(), SpaceDim(),
                range.data(),
                ops_arg_dat(g_f()[blockIndex], NUMXI, ONEPTREGULARSTENCIL,
                            "double", OPS_RW),
                ops_arg_dat(g_NodeType().at(componentID).at(blockIndex), 1,
                            ONEPTREGULARSTENCIL, "int", OPS_READ),
                ops_arg_dat(g_GeometryProperty()[blockIndex], 1, LOCALSTENCIL,
                            "int", OPS_READ),
                ops_arg_gbl(givenVars, 1, "double", OPS_READ),
                ops_arg_gbl(&surface, 1, "int", OPS_READ),
                ops_arg_gbl(g_Components().at(componentID).index, 2, "int",
                            OPS_READ));
        } break;
        case BoundaryScheme::EQMDiffuseRefl: {
            ops_par_loop(
                KerCutCellEQMDiffuseRefl, "KerCutCellEQMDiffuseRefl",
                block.Get(), SpaceDim(), range.data(),
                ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL, "double",
                            OPS_RW),
                ops_arg_dat(g_NodeType().at(componentID).at(blockIndex), 1,
                            LOCALSTENCIL, "int", OPS_READ),
                ops_arg_dat(g_GeometryProperty()[blockIndex], 1, LOCALSTENCIL,
                            "int", OPS_READ),
                ops_arg_gbl(givenVars, 2, "double", OPS_READ),
                ops_arg_gbl(g_Components().at(componentID).index, 2, "int",
                            OPS_READ));
        } break;
        case BoundaryScheme::FDPeriodic: {
            ops_par_loop(
                KerCutCellPeriodic, "KerCutCellPeriodic", block.Get(),
                SpaceDim(), range.data(),
                ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL, "double",
                            OPS_RW),
                ops_arg_dat(g_NodeType().at(componentID).at(blockIndex), 1,
                            LOCALSTENCIL, "int", OPS_READ),
                ops_arg_dat(g_GeometryProperty()[blockIndex], 1, LOCALSTENCIL,
                            "int", OPS_READ),
                ops_arg_gbl(g_Components().at(componentID).index, 2, "int",
                            OPS_READ),
                ops_arg_gbl(&surface, 1, "int", OPS_READ));
        } break;
        default:
            break;
    }
}
#endif //OPS_2D
