// Copyright 2017 the MPLB team. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

/*! @brief   Define kernel functions related to discrete velocity model
  * @author  Jianping Meng
  * @details Define kernel functions for calculating the equilibrium functions,
  * the relaxation time, body force term, and macroscopic variables.
  */

#ifndef MODEL_KERNEL_H
#define MODEL_KERNEL_H
#include "model.h"

/*!
 * We assume that the layout of MacroVars is rho, u, v, w, T, ...
 * In the macroVars, all variables are conserved, i.e., with density.
 * @todo how to deal with overflow in a kernel function? in particular, GPU
 */

#ifdef OPS_2D

//two dimensional code
void KerCalcFeq(const int* nodeType, const Real* macroVars, Real* feq) {
    VertexTypes vt{(VertexTypes)nodeType[OPS_ACC0(0, 0)]};
    if (vt != Vertex_ImmersedSolid) {
        for (int compoIndex = 0; compoIndex < NUMCOMPONENTS; compoIndex++) {
            EquilibriumType equilibriumType{
                (EquilibriumType)EQUILIBRIUMTYPE[compoIndex]};
            const int startPos{VARIABLECOMPPOS[2*compoIndex]};
            if (Equilibrium_BGKIsothermal2nd == equilibriumType) {
                Real rho{macroVars[OPS_ACC_MD1(startPos, 0, 0)]};
                Real u{macroVars[OPS_ACC_MD1(startPos+1, 0, 0)]};
                Real v{macroVars[OPS_ACC_MD1(startPos+2, 0, 0)]};
                const Real T{1};
                const int polyOrder{2};
                for (int xiIndex = COMPOINDEX[2 * compoIndex];
                     xiIndex <= COMPOINDEX[2 * compoIndex + 1]; xiIndex++) {
                    feq[OPS_ACC_MD2(xiIndex, 0, 0)] =
                        CalcBGKFeq(xiIndex, rho, u, v, T, polyOrder);
                }
            }
            if (Equilibrium_BGKThermal4th == equilibriumType) {
                Real rho{macroVars[OPS_ACC_MD1(startPos, 0, 0)]};
                Real u{macroVars[OPS_ACC_MD1(startPos + 1, 0, 0)]};
                Real v{macroVars[OPS_ACC_MD1(startPos + 2, 0, 0)]};
                Real T{macroVars[OPS_ACC_MD1(startPos + 3, 0, 0)]};
                const int polyOrder{4};
                for (int xiIndex = COMPOINDEX[2 * compoIndex];
                     xiIndex <= COMPOINDEX[2 * compoIndex + 1]; xiIndex++) {
                    feq[OPS_ACC_MD2(xiIndex, 0, 0)] =
                        CalcBGKFeq(xiIndex, rho, u, v, T, polyOrder);
                }
            }
            if (Equilibrium_BGKSWE4th == equilibriumType) {
                Real h{macroVars[OPS_ACC_MD1(startPos, 0, 0)]};
                Real u{macroVars[OPS_ACC_MD1(startPos + 1, 0, 0)]};
                Real v{macroVars[OPS_ACC_MD1(startPos + 2, 0, 0)]};
                const int polyOrder{4};
                for (int xiIndex = COMPOINDEX[2 * compoIndex];
                     xiIndex <= COMPOINDEX[2 * compoIndex + 1]; xiIndex++) {
                    feq[OPS_ACC_MD2(xiIndex, 0, 0)] =
                        CalcSWEFeq(xiIndex, h, u, v, polyOrder);
                }
            }
        }
    }
}


/*!
 * For most of cases, in particular, tau may not change actually
 * But for compressible flows, it will change
 * Due to the property of LBM, even modelling a liquid,there will be a
 * 'Knudsen' number
 */
void KerCalcTau(const int* nodeType, const Real* tauRef, const Real* macroVars,
                Real* tau) {  
    VertexTypes vt = (VertexTypes)nodeType[OPS_ACC0(0, 0)];    
    if (vt != Vertex_ImmersedSolid) {
        for (int compoIndex = 0; compoIndex < NUMCOMPONENTS; compoIndex++) {
            EquilibriumType equilibriumType{
                (EquilibriumType)EQUILIBRIUMTYPE[compoIndex]};
            const int startPos{VARIABLECOMPPOS[2*compoIndex]};
            if (Equilibrium_BGKIsothermal2nd == equilibriumType) {
                Real rho{macroVars[OPS_ACC_MD2(startPos, 0, 0)]};                
                const Real T{1};               
                tau[OPS_ACC_MD3(compoIndex, 0, 0)] =
                    tauRef[compoIndex] / (rho * sqrt(T));
            }
            if (Equilibrium_BGKThermal4th == equilibriumType) {
                Real rho{macroVars[OPS_ACC_MD2(startPos, 0, 0)]};
                Real T{macroVars[OPS_ACC_MD2(startPos+3, 0, 0)]};
                tau[OPS_ACC_MD3(compoIndex, 0, 0)] =
                    tauRef[compoIndex] / (rho * sqrt(T));
            }
            if (Equilibrium_BGKSWE4th == equilibriumType) {
                Real h{macroVars[OPS_ACC_MD2(startPos, 0, 0)]};
                tau[OPS_ACC_MD3(compoIndex, 0, 0)] = tauRef[compoIndex] / h;
            }           
        }
    }
}
/*!
 * If a Newton-Cotes quadrature is used, it can be converted to the way
 * similar to the Gauss-Hermite quadrature *
 */
void KerCalcMacroVars(const int* nodeType, const Real* f, Real* macroVars) {
    VertexTypes vt = (VertexTypes)nodeType[OPS_ACC0(0, 0)];
    if (vt != Vertex_ImmersedSolid) {
        bool rhoCalculated{false};
        bool veloCalculated[LATTDIM];
        for (int lattIdx = 0; lattIdx < LATTDIM; lattIdx++) {
            veloCalculated[lattIdx] = false;
        }
        for (int m = 0; m < NUMMACROVAR; m++) {
            macroVars[OPS_ACC_MD2(m, 0, 0)] = 0;
            VariableTypes varType = (VariableTypes)VARIABLETYPE[m];
            switch (varType) {
                case Variable_Rho: {
                    rhoCalculated = true;
                    for (int xiIdx = 0; xiIdx < NUMXI; xiIdx++) {
                        macroVars[OPS_ACC_MD2(m, 0, 0)] +=
                            f[OPS_ACC_MD1(xiIdx, 0, 0)];
                    }
#ifdef debug
                    Real rho{macroVars[OPS_ACC_MD2(m, 0, 0)]};
                    if (isnan(rho) || rho <= 0 || isinf(rho)) {
                        ops_printf(
                            "%sDensity=%f\n",
                            "Density becomes invalid! Maybe something wrong...",
                            rho);
                    }
#endif
                } break;
                case Variable_U: {
                    if (rhoCalculated) {
                        veloCalculated[0] = true;
                        for (int xiIdx = 0; xiIdx < NUMXI; xiIdx++) {
                            macroVars[OPS_ACC_MD2(m, 0, 0)] +=
                                CS * XI[xiIdx * LATTDIM] *
                                f[OPS_ACC_MD1(xiIdx, 0, 0)];
                        }
                        macroVars[OPS_ACC_MD2(m, 0, 0)] /=
                            macroVars[OPS_ACC_MD2(0, 0, 0)];
#ifdef debug
                        Real u{macroVars[OPS_ACC_MD2(m, 0, 0)]};
                        if (isnan(u) || isinf(u)) {
                            ops_printf("%sU=%f\n",
                                       "Velocity becomes invalid! Maybe "
                                       "something wrong...",
                                       u);
                        }

#endif
                    } else {
                        ops_printf("%s\n",
                                   "Density has not been calculated before "
                                   "calculating U!");
                    }
                } break;
                case Variable_V: {
                    if (rhoCalculated) {
                        veloCalculated[1] = true;
                        for (int xiIdx = 0; xiIdx < NUMXI; xiIdx++) {
                            macroVars[OPS_ACC_MD2(m, 0, 0)] +=
                                CS * XI[xiIdx * LATTDIM + 1] *
                                f[OPS_ACC_MD1(xiIdx, 0, 0)];
                        }
                        macroVars[OPS_ACC_MD2(m, 0, 0)] /=
                            macroVars[OPS_ACC_MD2(0, 0, 0)];
#ifdef debug
                        Real v{macroVars[OPS_ACC_MD2(m, 0, 0)]};
                        if (isnan(v) || isinf(v)) {
                            ops_printf("%sV=%f\n",
                                       "Velocity becomes invalid! Maybe "
                                       "something wrong...",
                                       v);
                        }
#endif
                    } else {
                        ops_printf("%s\n",
                                   "Density has not been calculated before "
                                   "calculating V!");
                    }
                } break;
                case Variable_W: {
                    if (rhoCalculated) {
                        veloCalculated[2] = true;
                        for (int xiIdx = 0; xiIdx < NUMXI; xiIdx++) {
                            macroVars[OPS_ACC_MD2(m, 0, 0)] +=
                                CS * XI[xiIdx * LATTDIM + 2] *
                                f[OPS_ACC_MD1(xiIdx, 0, 0)];
                        }
                        macroVars[OPS_ACC_MD2(m, 0, 0)] /=
                            macroVars[OPS_ACC_MD2(0, 0, 0)];
#ifdef debug
                        Real w{macroVars[OPS_ACC_MD2(m, 0, 0)]};
                        if (isnan(w) || isinf(w)) {
                            ops_printf("%sW=%f\n",
                                       "Velocity becomes invalid! Maybe "
                                       "something wrong...",
                                       w);
                        }
#endif
                    } else {
                        ops_printf("%s\n",
                                   "Density has not been calculated before "
                                   "calculating W!");
                    }
                } break;
                case Variable_Qx: {
                    bool ifCalc = true;
                    Real velo[LATTDIM];
                    for (int d = 0; d < LATTDIM; d++) {
                        ifCalc = ifCalc && veloCalculated[d];
                        velo[d] = macroVars[OPS_ACC_MD2(d + 1, 0, 0)];
                    }
                    if (ifCalc) {
                        for (int xiIdx = 0; xiIdx < NUMXI; xiIdx++) {
                            Real T = 0;
                            for (int d = 0; d < LATTDIM; d++) {
                                T += (CS * XI[xiIdx * LATTDIM + d] - velo[d]) *
                                     (CS * XI[xiIdx * LATTDIM + d] - velo[d]) *
                                     f[OPS_ACC_MD1(xiIdx, 0, 0)];
                            }
                            macroVars[OPS_ACC_MD2(m, 0, 0)] +=
                                (0.5 * (CS * XI[xiIdx * LATTDIM] - velo[0]) *
                                 T);
                        }
#ifdef debug
                        Real qx{macroVars[OPS_ACC_MD2(m, 0, 0)]};
                        if (isnan(qx) || isinf(qx)) {
                            ops_printf("%sQx=%f\n",
                                       "Heatflux becomes invalid! Maybe "
                                       "something wrong...",
                                       qx);
                        }
#endif
                    } else {
                        ops_printf("%s\n",
                                   "The macroscopic velocity have "
                                   "not been calculated before calculating the "
                                   "temperature!");
                    }
                } break;
                case Variable_Qy: {
                    bool ifCalc = true;
                    Real velo[LATTDIM];
                    for (int d = 0; d < LATTDIM; d++) {
                        ifCalc = ifCalc && veloCalculated[d];
                        velo[d] = macroVars[OPS_ACC_MD2(d + 1, 0, 0)];
                    }
                    if (ifCalc) {
                        for (int xiIdx = 0; xiIdx < NUMXI; xiIdx++) {
                            Real T = 0;
                            for (int d = 0; d < LATTDIM; d++) {
                                T += (CS * XI[xiIdx * LATTDIM + d] - velo[d]) *
                                     (CS * XI[xiIdx * LATTDIM + d] - velo[d]) *
                                     f[OPS_ACC_MD1(xiIdx, 0, 0)];
                            }
                            macroVars[OPS_ACC_MD2(m, 0, 0)] +=
                                (0.5 *
                                 (CS * XI[xiIdx * LATTDIM + 1] - velo[1]) * T);
                        }
#ifdef debug
                        Real qy{macroVars[OPS_ACC_MD2(m, 0, 0)]};
                        if (isnan(qy) || isinf(qy)) {
                            ops_printf("%sQx=%f\n",
                                       "Heatflux becomes invalid! Maybe "
                                       "something wrong...",
                                       qy);
                        }
#endif
                    } else {
                        ops_printf("%s\n",
                                   "The density and macroscopic velocity have "
                                   "not been calculated before calculating the "
                                   "temperature!");
                    }
                } break;
                case Variable_Qz: {
                    bool ifCalc = true;
                    Real velo[LATTDIM];
                    for (int d = 0; d < LATTDIM; d++) {
                        ifCalc = ifCalc && veloCalculated[d];
                        velo[d] = macroVars[OPS_ACC_MD2(d + 1, 0, 0)];
                    }
                    if (ifCalc) {
                        for (int xiIdx = 0; xiIdx < NUMXI; xiIdx++) {
                            Real T = 0;
                            for (int d = 0; d < LATTDIM; d++) {
                                T += (CS * XI[xiIdx * LATTDIM + d] - velo[d]) *
                                     (CS * XI[xiIdx * LATTDIM + d] - velo[d]) *
                                     f[OPS_ACC_MD1(xiIdx, 0, 0)];
                            }
                            macroVars[OPS_ACC_MD2(m, 0, 0)] +=
                                (0.5 *
                                 (CS * XI[xiIdx * LATTDIM + 2] - velo[2]) * T);
                        }
#ifdef debug
                        Real qz{macroVars[OPS_ACC_MD2(m, 0, 0)]};
                        if (isnan(qz) || isinf(qz)) {
                            ops_printf("%sQx=%f\n",
                                       "Heatflux becomes invalid! Maybe "
                                       "something wrong...",
                                       qz);
                        }
#endif
                    } else {
                        ops_printf("%s\n",
                                   "The density and macroscopic velocity have "
                                   "not been calculated before calculating the "
                                   "temperature!");
                    }
                } break;
                case Variable_T: {
                    bool ifCalc = rhoCalculated;
                    Real rho = macroVars[OPS_ACC_MD2(0, 0, 0)];
                    Real velo[LATTDIM];
                    for (int d = 0; d < LATTDIM; d++) {
                        ifCalc = ifCalc && veloCalculated[d];
                        velo[d] = macroVars[OPS_ACC_MD2(d + 1, 0, 0)];
                    }
                    if (ifCalc) {
                        for (int xiIdx = 0; xiIdx < NUMXI; xiIdx++) {
                            for (int d = 0; d < LATTDIM; d++) {
                                macroVars[OPS_ACC_MD2(m, 0, 0)] +=
                                    (CS * XI[xiIdx * LATTDIM + d] - velo[d]) *
                                    (CS * XI[xiIdx * LATTDIM + d] - velo[d]) *
                                    f[OPS_ACC_MD1(xiIdx, 0, 0)];
                            }
                        }
                        macroVars[OPS_ACC_MD2(m, 0, 0)]/=(rho*LATTDIM);
#ifdef debug
                        Real T{macroVars[OPS_ACC_MD2(m, 0, 0)]};
                        if (isnan(T) || isinf(T)) {
                            ops_printf("%sW=%f\n",
                                       "Temperature becomes invalid! Maybe "
                                       "something wrong...",
                                       T);
                        }
#endif
                    } else {
                        ops_printf("%s\n",
                                   "The density and macroscopic velocity have "
                                   "not been calculated before calculating the "
                                   "temperature!");
                    }
                } break;
                default:
                    break;
            }//Switch
        }  // m
    } // isVertex
}
#endif
#ifdef OPS_3D
void KerCalcFeq3D(const int* nodeType, const Real* macroVars, Real* feq) {
    VertexTypes vt = (VertexTypes)nodeType[OPS_ACC0(0, 0, 0)];
    if (vt != Vertex_ImmersedSolid) {
        for (int compoIndex = 0; compoIndex < NUMCOMPONENTS; compoIndex++) {
            EquilibriumType equilibriumType{
                (EquilibriumType)EQUILIBRIUMTYPE[compoIndex]};
            const int startPos{VARIABLECOMPPOS[2*compoIndex]};
            if (Equilibrium_BGKIsothermal2nd == equilibriumType) {
                Real rho{macroVars[OPS_ACC_MD1(startPos, 0, 0, 0)]};
                Real u{macroVars[OPS_ACC_MD1(startPos + 1, 0, 0, 0)]};
                Real v{macroVars[OPS_ACC_MD1(startPos + 2, 0, 0, 0)]};
                Real w{macroVars[OPS_ACC_MD1(startPos + 3, 0, 0, 0)]};
                const Real T{1};
                const int polyOrder{2};
                for (int xiIndex = COMPOINDEX[2 * compoIndex];
                     xiIndex <= COMPOINDEX[2 * compoIndex + 1]; xiIndex++) {
                    feq[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                        CalcBGKFeq(xiIndex, rho, u, v, w, T, polyOrder);
                }
            }
            if (Equilibrium_BGKThermal4th == equilibriumType) {
                Real rho{macroVars[OPS_ACC_MD1(startPos, 0, 0, 0)]};
                Real u{macroVars[OPS_ACC_MD1(startPos + 1, 0, 0, 0)]};
                Real v{macroVars[OPS_ACC_MD1(startPos + 2, 0, 0, 0)]};
                Real w{macroVars[OPS_ACC_MD1(startPos + 3, 0, 0, 0)]};
                Real T{macroVars[OPS_ACC_MD1(startPos + 4, 0, 0, 0)]};
                const int polyOrder{4};
                for (int xiIndex = COMPOINDEX[2 * compoIndex];
                     xiIndex <= COMPOINDEX[2 * compoIndex + 1]; xiIndex++) {
                    feq[OPS_ACC_MD2(xiIndex, 0, 0, 0)] =
                        CalcBGKFeq(xiIndex, rho, u, v, w, T, polyOrder);
                }
            }
        }
    }
}

void KerCalcBodyForce3D(const Real* time, const int* nodeType,
                        const Real* coordinates, const Real* macroVars,
                        Real* bodyForce) {
    // here we assume the force is constant
    // user may introduce a function of g(r,t) for the body force
    const Real g[]{0.0001, 0, 0};
    VertexTypes vt = (VertexTypes)nodeType[OPS_ACC1(0, 0, 0)];
    if (vt != Vertex_ImmersedSolid) {
        for (int compoIndex = 0; compoIndex < NUMCOMPONENTS; compoIndex++) {
            BodyForceType forceType{(BodyForceType)FORCETYPE[compoIndex]};
            const int startPos{VARIABLECOMPPOS[2 * compoIndex]};
            if (BodyForce_1st == forceType) {
                Real rho{macroVars[OPS_ACC_MD3(startPos, 0, 0, 0)]};
                //ops_printf("rho=%f g=%f\n",rho,g[0]);
                for (int xiIndex = COMPOINDEX[2 * compoIndex];
                     xiIndex <= COMPOINDEX[2 * compoIndex + 1]; xiIndex++) {
                    bodyForce[OPS_ACC_MD4(xiIndex, 0, 0, 0)] =
                        CalcBodyForce(xiIndex, rho, g);
                }
            }
            if (BodyForce_None == forceType) {
                for (int xiIndex = COMPOINDEX[2 * compoIndex];
                     xiIndex <= COMPOINDEX[2 * compoIndex + 1]; xiIndex++) {
                    bodyForce[OPS_ACC_MD4(xiIndex, 0, 0, 0)] = 0;
                }
            }
        }
    }
}

void KerCalcTau3D(const int* nodeType, const Real* tauRef,
                  const Real* macroVars, Real* tau) {
    /*
     *@note: multicomponent ready for incompressible flows.
     */
    VertexTypes vt = (VertexTypes)nodeType[OPS_ACC0(0, 0, 0)];
    if (vt != Vertex_ImmersedSolid) {
        for (int compoIndex = 0; compoIndex < NUMCOMPONENTS; compoIndex++) {
            EquilibriumType equilibriumType{
                (EquilibriumType)EQUILIBRIUMTYPE[compoIndex]};
            const int startPos{VARIABLECOMPPOS[2*compoIndex]};
            if (Equilibrium_BGKIsothermal2nd == equilibriumType) {
                Real rho{macroVars[OPS_ACC_MD2(startPos, 0, 0, 0)]};
                const Real T{1};
                tau[OPS_ACC_MD3(compoIndex, 0, 0, 0)] =
                    tauRef[compoIndex] / (rho * sqrt(T));
            }
            if (Equilibrium_BGKThermal4th == equilibriumType) {
                Real rho{macroVars[OPS_ACC_MD2(startPos, 0, 0, 0)]};
                Real T{macroVars[OPS_ACC_MD2(startPos + 3, 0, 0, 0)]};
                tau[OPS_ACC_MD3(compoIndex, 0, 0, 0)] =
                    tauRef[compoIndex] / (rho * sqrt(T));
            }
        }
    }
}
/*!
 * If a Newton-Cotes quadrature is used, it can be converted to the way
 * similar to the Gauss-Hermite quadrature
 *
 */
void KerCalcMacroVars3D(const Real* dt, const int* nodeType, const Real* coordinates, const Real* f, Real* macroVars) {
    VertexTypes vt = (VertexTypes)nodeType[OPS_ACC1(0, 0, 0)];
    Real* acceleration = new Real[LATTDIM * NUMCOMPONENTS];
    for (int compoIndex = 0; compoIndex < NUMCOMPONENTS; compoIndex++) {
        for (int i = 0; i < LATTDIM; i++) {
            acceleration[compoIndex+i] = 0;
        }
    }
    acceleration[0] = 0.0001;
    if (vt != Vertex_ImmersedSolid) {
        for (int compoIndex = 0; compoIndex < NUMCOMPONENTS; compoIndex++) {
            bool rhoCalculated{false};
            Real rho{0};
            bool* veloCalculated = new bool[LATTDIM];
            Real* velo = new Real[LATTDIM];
            for (int lattIdx = 0; lattIdx < LATTDIM; lattIdx++) {
                veloCalculated[lattIdx] = false;
                velo[lattIdx] = 0;
            }
            for (int m = VARIABLECOMPPOS[2 * compoIndex];
                 m <= VARIABLECOMPPOS[2 * compoIndex + 1]; m++) {
                macroVars[OPS_ACC_MD4(m, 0, 0, 0)] = 0;
                VariableTypes varType = (VariableTypes)VARIABLETYPE[m];
                switch (varType) {
                    case Variable_Rho: {
                        rhoCalculated = true;
                        for (int xiIdx = COMPOINDEX[2 * compoIndex];
                             xiIdx <= COMPOINDEX[2 * compoIndex + 1]; xiIdx++) {
                            macroVars[OPS_ACC_MD4(m, 0, 0, 0)] +=
                                f[OPS_ACC_MD3(xiIdx, 0, 0, 0)];
                        }
                        rho = macroVars[OPS_ACC_MD4(m, 0, 0, 0)];
#ifdef debug
                        if (isnan(rho) || rho <= 0 || isinf(rho)) {
                            ops_printf("%sDensity=%f\n",
                                       "Density becomes invalid! Maybe "
                                       "something wrong...",
                                       rho);
                        }
#endif
                    } break;
                    case Variable_U: {
                        if (rhoCalculated) {
                            veloCalculated[0] = true;
                            for (int xiIdx = COMPOINDEX[2 * compoIndex];
                                 xiIdx <= COMPOINDEX[2 * compoIndex + 1];
                                 xiIdx++) {
                                macroVars[OPS_ACC_MD4(m, 0, 0, 0)] +=
                                    CS * XI[xiIdx * LATTDIM] *
                                    f[OPS_ACC_MD3(xiIdx, 0, 0, 0)];
                            }
                            macroVars[OPS_ACC_MD4(m, 0, 0, 0)] /=
                                macroVars[OPS_ACC_MD4(0, 0, 0, 0)];
                            velo[0] = macroVars[OPS_ACC_MD4(m, 0, 0, 0)];
#ifdef debug
                            if (isnan(velo[0]) || isinf(velo[0])) {
                                ops_printf("%sU=%f\n",
                                           "Velocity becomes invalid! Maybe "
                                           "something wrong...",
                                           u);
                            }
#endif
                        } else {
                            ops_printf("%s\n",
                                       "Density has not been calculated before "
                                       "calculating U!");
                        }
                    } break;
                    case Variable_V: {
                        if (rhoCalculated) {
                            veloCalculated[1] = true;
                            for (int xiIdx = COMPOINDEX[2 * compoIndex];
                                 xiIdx <= COMPOINDEX[2 * compoIndex + 1];
                                 xiIdx++) {
                                macroVars[OPS_ACC_MD4(m, 0, 0, 0)] +=
                                    CS * XI[xiIdx * LATTDIM + 1] *
                                    f[OPS_ACC_MD3(xiIdx, 0, 0, 0)];
                            }
                            macroVars[OPS_ACC_MD4(m, 0, 0, 0)] /=
                                macroVars[OPS_ACC_MD4(0, 0, 0, 0)];
                            velo[1] = macroVars[OPS_ACC_MD4(m, 0, 0, 0)];
#ifdef debug
                            if (isnan(velo[1]) || isinf(velo[1])) {
                                ops_printf("%sV=%f\n",
                                           "Velocity becomes invalid! Maybe "
                                           "something wrong...",
                                           v);
                            }
#endif
                        } else {
                            ops_printf("%s\n",
                                       "Density has not been calculated before "
                                       "calculating V!");
                        }
                    } break;
                    case Variable_W: {
                        if (rhoCalculated) {
                            veloCalculated[2] = true;
                            for (int xiIdx = COMPOINDEX[2 * compoIndex];
                                 xiIdx <= COMPOINDEX[2 * compoIndex + 1];
                                 xiIdx++) {
                                macroVars[OPS_ACC_MD4(m, 0, 0, 0)] +=
                                    CS * XI[xiIdx * LATTDIM + 2] *
                                    f[OPS_ACC_MD3(xiIdx, 0, 0, 0)];
                            }
                            macroVars[OPS_ACC_MD4(m, 0, 0, 0)] /=
                                macroVars[OPS_ACC_MD4(0, 0, 0, 0)];
                            velo[2] = macroVars[OPS_ACC_MD4(m, 0, 0, 0)];
#ifdef debug
                            if (isnan(velo[2]) || isinf(velo[2])) {
                                ops_printf("%sW=%f\n",
                                           "Velocity becomes invalid! Maybe "
                                           "something wrong...",
                                           w);
                            }
#endif
                        } else {
                            ops_printf("%s\n",
                                       "Density has not been calculated before "
                                       "calculating W!");
                        }
                    } break;
                    case Variable_Qx: {
                        bool ifCalc = true;
                        for (int d = 0; d < LATTDIM; d++) {
                            ifCalc = ifCalc && veloCalculated[d];
                        }
                        if (ifCalc) {
                            for (int xiIdx = COMPOINDEX[2 * compoIndex];
                                 xiIdx <= COMPOINDEX[2 * compoIndex + 1];
                                 xiIdx++) {
                                Real T = 0;
                                for (int d = 0; d < LATTDIM; d++) {
                                    T += (CS * XI[xiIdx * LATTDIM + d] -
                                          velo[d]) *
                                         (CS * XI[xiIdx * LATTDIM + d] -
                                          velo[d]) *
                                         f[OPS_ACC_MD3(xiIdx, 0, 0, 0)];
                                }
                                macroVars[OPS_ACC_MD4(m, 0, 0, 0)] +=
                                    (0.5 *
                                     (CS * XI[xiIdx * LATTDIM] - velo[0]) * T);
                            }
#ifdef debug
                            Real qx{macroVars[OPS_ACC_MD4(m, 0, 0, 0)]};
                            if (isnan(qx) || isinf(qx)) {
                                ops_printf("%sQx=%f\n",
                                           "Heatflux becomes invalid! Maybe "
                                           "something wrong...",
                                           qx);
                            }
#endif
                        } else {
                            ops_printf(
                                "%s\n",
                                "The macroscopic velocity have "
                                "not been calculated before calculating the "
                                "temperature!");
                        }
                    } break;
                    case Variable_Qy: {
                        bool ifCalc = true;
                        for (int d = 0; d < LATTDIM; d++) {
                            ifCalc = ifCalc && veloCalculated[d];
                        }
                        if (ifCalc) {
                            for (int xiIdx = COMPOINDEX[2 * compoIndex];
                                 xiIdx <= COMPOINDEX[2 * compoIndex + 1];
                                 xiIdx++) {
                                Real T = 0;
                                for (int d = 0; d < LATTDIM; d++) {
                                    T += (CS * XI[xiIdx * LATTDIM + d] -
                                          velo[d]) *
                                         (CS * XI[xiIdx * LATTDIM + d] -
                                          velo[d]) *
                                         f[OPS_ACC_MD3(xiIdx, 0, 0, 0)];
                                }
                                macroVars[OPS_ACC_MD4(m, 0, 0, 0)] +=
                                    (0.5 *
                                     (CS * XI[xiIdx * LATTDIM + 1] - velo[1]) *
                                     T);
                            }
#ifdef debug
                            Real qy{macroVars[OPS_ACC_MD4(m, 0, 0, 0)]};
                            if (isnan(qy) || isinf(qy)) {
                                ops_printf("%sQx=%f\n",
                                           "Heatflux becomes invalid! Maybe "
                                           "something wrong...",
                                           qy);
                            }
#endif
                        } else {
                            ops_printf(
                                "%s\n",
                                "The density and macroscopic velocity have "
                                "not been calculated before calculating the "
                                "temperature!");
                        }
                    } break;
                    case Variable_Qz: {
                        bool ifCalc = true;
                        for (int d = 0; d < LATTDIM; d++) {
                            ifCalc = ifCalc && veloCalculated[d];
                        }
                        if (ifCalc) {
                            for (int xiIdx = COMPOINDEX[2 * compoIndex];
                                 xiIdx <= COMPOINDEX[2 * compoIndex + 1];
                                 xiIdx++) {
                                Real T = 0;
                                for (int d = 0; d < LATTDIM; d++) {
                                    T += (CS * XI[xiIdx * LATTDIM + d] -
                                          velo[d]) *
                                         (CS * XI[xiIdx * LATTDIM + d] -
                                          velo[d]) *
                                         f[OPS_ACC_MD3(xiIdx, 0, 0, 0)];
                                }
                                macroVars[OPS_ACC_MD4(m, 0, 0, 0)] +=
                                    (0.5 *
                                     (CS * XI[xiIdx * LATTDIM + 2] - velo[2]) *
                                     T);
                            }
#ifdef debug
                            Real qz{macroVars[OPS_ACC_MD4(m, 0, 0, 0)]};
                            if (isnan(qz) || isinf(qz)) {
                                ops_printf("%sQx=%f\n",
                                           "Heatflux becomes invalid! Maybe "
                                           "something wrong...",
                                           qz);
                            }
#endif
                        } else {
                            ops_printf(
                                "%s\n",
                                "The density and macroscopic velocity have "
                                "not been calculated before calculating the "
                                "temperature!");
                        }
                    } break;
                    case Variable_T: {
                        bool ifCalc = rhoCalculated;
                        for (int d = 0; d < LATTDIM; d++) {
                            ifCalc = ifCalc && veloCalculated[d];
                        }
                        if (ifCalc) {
                            for (int xiIdx = COMPOINDEX[2 * compoIndex];
                                 xiIdx <= COMPOINDEX[2 * compoIndex + 1];
                                 xiIdx++) {
                                for (int d = 0; d < LATTDIM; d++) {
                                    macroVars[OPS_ACC_MD4(m, 0, 0, 0)] +=
                                        (CS * XI[xiIdx * LATTDIM + d] -
                                         velo[d]) *
                                        (CS * XI[xiIdx * LATTDIM + d] -
                                         velo[d]) *
                                        f[OPS_ACC_MD3(xiIdx, 0, 0, 0)];
                                }
                            }
                            macroVars[OPS_ACC_MD4(m, 0, 0, 0)] /=
                                (rho * LATTDIM);
#ifdef debug
                            Real T{macroVars[OPS_ACC_MD4(m, 0, 0, 0)]};
                            if (isnan(T) || isinf(T)) {
                                ops_printf("%sW=%f\n",
                                           "Temperature becomes invalid! Maybe "
                                           "something wrong...",
                                           T);
                            }
#endif
                        } else {
                            ops_printf(
                                "%s\n",
                                "The density and macroscopic velocity have "
                                "not been calculated before calculating the "
                                "temperature!");
                        }
                    } break;
                    case Variable_U_Force: {
                        if (rhoCalculated) {
                            veloCalculated[0] = true;
                            for (int xiIdx = COMPOINDEX[2 * compoIndex];
                                 xiIdx <= COMPOINDEX[2 * compoIndex + 1];
                                 xiIdx++) {
                                macroVars[OPS_ACC_MD4(m, 0, 0, 0)] +=
                                    CS * XI[xiIdx * LATTDIM] *
                                    f[OPS_ACC_MD3(xiIdx, 0, 0, 0)];
                            }
                            macroVars[OPS_ACC_MD4(m, 0, 0, 0)] /=
                                macroVars[OPS_ACC_MD4(0, 0, 0, 0)];
                            if (Vertex_Fluid == vt) {
                                macroVars[OPS_ACC_MD4(m, 0, 0, 0)] +=
                                    ((*dt) *
                                     acceleration[compoIndex * LATTDIM] / 2);
                            }
                            velo[0] = macroVars[OPS_ACC_MD4(m, 0, 0, 0)];
#ifdef debug
                            if (isnan(velo[0]) || isinf(velo[0])) {
                                ops_printf("%sU=%f\n",
                                           "Velocity becomes invalid! Maybe "
                                           "something wrong...",
                                           u);
                            }
#endif
                        } else {
                            ops_printf("%s\n",
                                       "Density has not been calculated before "
                                       "calculating U!");
                        }
                    } break;
                    case Variable_V_Force: {
                        if (rhoCalculated) {
                            veloCalculated[1] = true;
                            for (int xiIdx = COMPOINDEX[2 * compoIndex];
                                 xiIdx <= COMPOINDEX[2 * compoIndex + 1];
                                 xiIdx++) {
                                macroVars[OPS_ACC_MD4(m, 0, 0, 0)] +=
                                    CS * XI[xiIdx * LATTDIM + 1] *
                                    f[OPS_ACC_MD3(xiIdx, 0, 0, 0)];
                            }
                            macroVars[OPS_ACC_MD4(m, 0, 0, 0)] /=
                                macroVars[OPS_ACC_MD4(0, 0, 0, 0)];
                            if (Vertex_Fluid == vt) {
                                macroVars[OPS_ACC_MD4(m, 0, 0, 0)] +=
                                    ((*dt) *
                                     acceleration[compoIndex * LATTDIM + 1] /
                                     2);
                            }
                            velo[1] = macroVars[OPS_ACC_MD4(m, 0, 0, 0)];
#ifdef debug
                            if (isnan(velo[1]) || isinf(velo[1])) {
                                ops_printf("%sV=%f\n",
                                           "Velocity becomes invalid! Maybe "
                                           "something wrong...",
                                           v);
                            }
#endif
                        } else {
                            ops_printf("%s\n",
                                       "Density has not been calculated before "
                                       "calculating V!");
                        }
                    } break;
                    case Variable_W_Force: {
                        if (rhoCalculated) {
                            veloCalculated[2] = true;
                            for (int xiIdx = COMPOINDEX[2 * compoIndex];
                                 xiIdx <= COMPOINDEX[2 * compoIndex + 1];
                                 xiIdx++) {
                                macroVars[OPS_ACC_MD4(m, 0, 0, 0)] +=
                                    CS * XI[xiIdx * LATTDIM + 2] *
                                    f[OPS_ACC_MD3(xiIdx, 0, 0, 0)];
                            }
                            macroVars[OPS_ACC_MD4(m, 0, 0, 0)] /=
                                macroVars[OPS_ACC_MD4(0, 0, 0, 0)];
                            if (Vertex_Fluid == vt) {
                                macroVars[OPS_ACC_MD4(m, 0, 0, 0)] +=
                                    ((*dt) *
                                     acceleration[compoIndex * LATTDIM + 2] /
                                     2);
                            }
                            velo[2] = macroVars[OPS_ACC_MD4(m, 0, 0, 0)];
#ifdef debug
                            if (isnan(velo[2]) || isinf(velo[2])) {
                                ops_printf("%sW=%f\n",
                                           "Velocity becomes invalid! Maybe "
                                           "something wrong...",
                                           w);
                            }
#endif
                        } else {
                            ops_printf("%s\n",
                                       "Density has not been calculated before "
                                       "calculating W!");
                        }
                    } break;
                    default:
                        break;
                }  // Switch
            }      // m
            delete[] veloCalculated;
            delete[] velo;
        }  // compoIdx
    }  // isVertex
    delete[] acceleration;
}
#endif
#endif  // MODEL_KERNEL_H
