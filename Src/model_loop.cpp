#include <vector>
#include <map>
#include "flowfield.h"
#include "model.h"
#include "scheme.h"
#include "ops_seq_v2.h"
#include "model_kernel.inc"
#ifdef OPS_3D
void PreDefinedCollision3D() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const SizeType blockIndex{block.ID()};
        for (auto& pair : CollisionTerms()) {
            const SizeType compoId{pair.first};
            const CollisionType collisionType{pair.second};
            const Real tau{TauRef()[compoId]};
            const Real* pdt{pTimeStep()};

            switch (collisionType) {
                case Collision_BGKIsothermal2nd:
                    ops_par_loop(
                        KerCollideBGKIsothermal3D, "KerCollideBGKIsothermal3D",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_READ),
                        ops_arg_dat(g_MacroVars()[blockIndex], NUMMACROVAR,
                                    LOCALSTENCIL, "double", OPS_RW),
                        ops_arg_dat(g_NodeType()[blockIndex], NUMCOMPONENTS,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_gbl(&tau, 1, "double", OPS_READ),
                        ops_arg_gbl(pdt, 1, "double", OPS_READ),
                        ops_arg_gbl(&compoId, 1, "int", OPS_READ));
                    break;
                case Collision_BGKThermal4th:
                    ops_par_loop(
                        KerCollideBGKThermal3D, "KerCollideBGKThermal3D",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_READ),
                        ops_arg_dat(g_MacroVars()[blockIndex], NUMMACROVAR,
                                    LOCALSTENCIL, "double", OPS_RW),
                        ops_arg_dat(g_NodeType()[blockIndex], NUMCOMPONENTS,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_gbl(&tau, 1, "double", OPS_READ),
                        ops_arg_gbl(pdt, 1, "double", OPS_READ),
                        ops_arg_gbl(&compoId, 1, "int", OPS_READ));
                    break;
                default:
                    ops_printf(
                        "The specified collision type is not implemented!\n");
                    break;
            }
        }
    }
}

void UpdateMacroVars3D() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const SizeType blockIndex{block.ID()};
        const Real* pdt{pTimeStep()};
        ops_par_loop(KerCalcMacroVars3D, "KerCalcMacroVars3D", block.Get(),
                     SpaceDim(), iterRng.data(),
                     ops_arg_dat(g_MacroVars()[blockIndex], NUMMACROVAR,
                                 LOCALSTENCIL, "double", OPS_RW),
                     ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                 "double", OPS_READ),
                     ops_arg_dat(g_NodeType()[blockIndex], NUMCOMPONENTS,
                                 LOCALSTENCIL, "int", OPS_READ),
                     ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                 LOCALSTENCIL, "double", OPS_READ),
                     ops_arg_dat(g_MacroBodyforce()[blockIndex],
                                 SpaceDim() * NUMCOMPONENTS, LOCALSTENCIL,
                                 "double", OPS_READ),
                     ops_arg_gbl(pdt, 1, "double", OPS_READ));
    }
}

void PreDefinedBodyForce3D() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const SizeType blockIndex{block.ID()};
        for (auto& pair : BodyForceTerms()) {
            const SizeType compoId{pair.first};
            const BodyForceType forceType{pair.second};
            switch (forceType) {
                case BodyForce_1st:
                    ops_par_loop(
                        KerCalcBodyForce1ST3D, "KerCalcBodyForce1ST3D",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(g_MacroBodyforce()[blockIndex],
                                    SpaceDim() * NUMCOMPONENTS, LOCALSTENCIL,
                                    "double", OPS_READ),
                        ops_arg_dat(g_MacroVars()[blockIndex], NUMMACROVAR,
                                    LOCALSTENCIL, "double", OPS_RW),
                        ops_arg_dat(g_NodeType()[blockIndex], NUMCOMPONENTS,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_gbl(&compoId, 1, "int", OPS_READ));
                    break;
                case BodyForce_None:
                    ops_par_loop(
                        KerCalcBodyForceNone3D, "KerCalcBodyForceNone",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(g_MacroBodyforce()[blockIndex],
                                    SpaceDim() * NUMCOMPONENTS, LOCALSTENCIL,
                                    "double", OPS_READ),
                        ops_arg_dat(g_MacroVars()[blockIndex], NUMMACROVAR,
                                    LOCALSTENCIL, "double", OPS_RW),
                        ops_arg_dat(g_NodeType()[blockIndex], NUMCOMPONENTS,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_gbl(&compoId, 1, "int", OPS_READ));
                    break;
                default:
                    ops_printf(
                        "The specified force type is not implemented!\n");
                    break;
            }
        }
    }
}


// TODO to be updated according to the new idea
void PreDefinedInitialCondition3D() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const SizeType blockIndex{block.ID()};
        for (auto& pair : InitialTerms()) {
            const SizeType compoId{pair.first};
            const InitialType initialType{pair.second};
            switch (initialType) {
                case Initial_BGKFeq2nd:
                    ops_par_loop(
                        KerInitialiseBGK2nd3D, "KerInitialiseBGK2nd3D",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(g_MacroVars()[blockIndex], NUMMACROVAR,
                                    LOCALSTENCIL, "double", OPS_RW),
                        ops_arg_dat(g_NodeType()[blockIndex], NUMCOMPONENTS,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_gbl(&compoId, 1, "int", OPS_READ));
                    break;
                default:
                    ops_printf(
                        "The specified initial type is not implemented!\n");
                    break;
            }
        }
    }
    //TODO this may be better arranged.
    if (!IsTransient()) {
        CopyCurrentMacroVar();
    }
}
#endif //OPS_3D