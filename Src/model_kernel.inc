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

/*! @brief   Define kernel functions related to discrete velocity model
 * @author  Jianping Meng
 * @details Define kernel functions for calculating the equilibrium functions,
 * the relaxation time, body force term, and macroscopic variables.
 */

#ifndef MODEL_KERNEL_INC
#define MODEL_KERNEL_INC
#include "ops_lib_core.h"

#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif
#include "type.h"
#include "model.h"
#include "flowfield_host_device.h"
#include "model_host_device.h"

/*!
 * We assume that the layout of MacroVars is rho, u, v, w, T, ...
 * In the macroVars, all variables are conserved, i.e., with density.
 * @todo how to deal with overflow in a kernel function? in particular, GPU
 */
#ifdef OPS_2D

//Function to compute the product of two nxn matrices
std::vector<std::vector<Real> > KerMultMatrix(std::vector<std::vector<double> > matA, std::vector<std::vector<double> > matB){
    Real S{0};
    std::vector<std::vector<Real> > mA{matA};
    std::vector<std::vector<Real> > mB{matB};
    int rowsA =  mA.size();
    int colsA = mA[0].size();
    int rowsB =  mB.size();
    int colsB = mB[0].size();
    std::vector<std::vector<Real> > mult(rowsA, std::vector<Real>(colsB));
    int i,j,k;
    Real prod;
    //std::cout<<colsA<<" "<<rowsA<<" "<<colsB<<" "<<rowsB<<"   ";
    //std::cout<<mA.size()<<" "<<mB[0].size()<<" "<<mA[0].size()<<"  ";
    /*
    for(i = 0; i < 9; ++i)
        for(j = 0; j < 9; ++j)
        {
            std::cout<<mA[i][j]<<"   ";
        }
    */
    for(i = 0; i < rowsA; ++i)
            for(j = 0; j < colsB; ++j)
            {
                mult[i][j]=0.0;
            }
    if (colsA==rowsB){


        for(i = 0; i < rowsA; ++i)
            for(j = 0; j < colsB; ++j)
                for(k = 0; k < colsA; ++k)
                {
                    prod = mA[i][k] * mB[k][j];
                    mult[i][j] += prod;
                }

        }
    else{
        std::cout<<"Number of columns in the second matrix must equal the number of rows of the first matrix!\n";

    }
    return mult;
}

//Calculates first-order or second-order finite difference gradients
void KerCalcGradients(ACC<Real>&grad, const ACC<Real>&f, const ACC<int>&nodeType,
                      const ACC<int>&geometry, const ACC<Real>& coordinates,
                      const int *lattIdx, const int *Order, const int *idx){
    VertexType vt = (VertexType)nodeType(0, 0);
    int order={*Order};
        //First order
        if (order == 1){
            if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic ||
                vt == VertexType::VirtualBoundary || vt == VertexType::Wall) {
                grad(0,0,0) = 0.5*(f(1, 0)-f(-1, 0));
                grad(1,0,0) = 0.5*(f(0, 1)-f(0, -1));
            }
        }
        //Second order
        else if (order == 2){
            if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic ||
                vt == VertexType::VirtualBoundary || vt == VertexType::Wall) {
                grad(0,0,0) = f(-1, 0)-2*f(0, 0)+f(1, 0);
                grad(1,0,0) = f(0, -1)-2*f(0, 0)+f(0, 1);

            }
        }
        else{
            std::cout << "Only First and Second order finite difference gradients"
                      << " are supported\n";
        }

}

//Wetting boundary implementation
//TODO: this kernel name is repeated in the AdDiff app
void KerUpdateRhoWetting(ACC<Real>& rho,const ACC<int>& geometry,
                         const ACC<int>&nodeType,const ACC<Real>& coordinates,
                         const Real *theta, const Real *Aconst, const Real *Kappa) {
    const Real teta={*theta};
    const Real alpha={acos(sin(teta)*sin(teta))};
    const Real Omga={2.0*((teta <= M_PI/2.0) - (M_PI/2.0 < teta))*sqrt(cos(alpha/3.0)*
                     (1.0-cos(alpha/3.0)))};
    const Real K={*Kappa};
    const Real A={*Aconst};
    const Real h={Omga*sqrt(K*A/2.0)};
    VertexType vt = (VertexType)nodeType(0, 0);
    VertexGeometryType vg = (VertexGeometryType)geometry(0, 0);
    VertexGeometryType vg1 = (VertexGeometryType)geometry(1, 0);
    VertexGeometryType vg2 = (VertexGeometryType)geometry(-1, 0);
    VertexGeometryType vg3 = (VertexGeometryType)geometry(0, 1);
    VertexGeometryType vg4 = (VertexGeometryType)geometry(0, -1);
    if (vt == VertexType::Wall) {
        switch (vg) {
            //Planar walls
            //+ve X normal
            case VG_IP: {
                if(vg3!=VG_IPJP_O&&vg4!=VG_IPJM_O){
                    rho(-1,0)=rho(1,0)+2*h/K;
                }
            } break;
            //-ve X normal
            case VG_IM: {
                if(vg3!=VG_IMJP_O&&vg4!=VG_IMJM_O){
                    rho(1,0)=rho(-1,0)+2*h/K;
                }
            } break;
            //+ve Y normal
            case VG_JP: {
                if(vg1!=VG_IPJP_O&&vg2!=VG_IMJP_O){
                    rho(0,-1)=rho(0,1)+2*h/K;
                }
            } break;
            //-ve Y normal
            case VG_JM: {
                if(vg1!=VG_IPJM_O&&vg2!=VG_IMJM_O){
                    rho(0,1)=rho(0,-1)+2*h/K;
                }
            } break;

            //Outer corners
            //+ve X, +ve Y normal
            case VG_IPJP_O: {
                rho(-1,-1)=(rho(1,0)+rho(0,1)+rho(1,-1)+rho(-1,1))/4+2*h/K;

            } break;
            //+ve X, -ve Y normal
            case VG_IPJM_O: {
                rho(-1,1)=(rho(1,0)+rho(0,-1)+rho(1,1)+rho(-1,-1))/4+2*h/K;

            } break;
            //-ve X, +ve Y normal
            case VG_IMJP_O: {
                rho(1,-1)=(rho(-1,0)+rho(0,1)+rho(1,1)+rho(-1,-1))/4+2*h/K;

            } break;
            //-ve X, -ve Y normal
            case VG_IMJM_O: {
                rho(1,1)=(rho(-1,0)+rho(0,-1)+rho(1,-1)+rho(-1,1))/4+2*h/K;
            } break;

            //Inner corners
            //+ve X, +ve Y normal
            case VG_IPJP_I: {
                if(vg3!=VG_IPJP_O){
                    rho(-1,0)=rho(1,0)+2*h/K;
                }
                if(vg1!=VG_IPJP_O){
                    rho(0,-1)=rho(0,1)+2*h/K;
                }
            } break;
            //+ve X, -ve Y normal
            case VG_IPJM_I: {
                if(vg4!=VG_IPJM_O){
                    rho(-1,0)=rho(1,0)+2*h/K;
                }
                if(vg1!=VG_IPJM_O){
                    rho(0,1)=rho(0,-1)+2*h/K;
                }
            } break;
            //-ve X, +ve Y normal
            case VG_IMJP_I: {
                if(vg3!=VG_IMJP_O){
                    rho(1,0)=rho(-1,0)+2*h/K;
                }
                if(vg2!=VG_IMJP_O){
                    rho(0,-1)=rho(0,1)+2*h/K;
                }
            } break;
            //-ve X, -ve Y normal
            case VG_IMJM_I: {
                if(vg4!=VG_IMJM_O){
                    rho(1,0)=rho(-1,0)+2*h/K;
                }
                if(vg2!=VG_IMJM_O){
                    rho(0,1)=rho(0,-1)+2*h/K;
                }
            } break;
        }
    }
}

//Calculates chemical potential in the Free Energy model
//TODO: This kernel name is repeated in the AdDiff App
void KerCalcMu(ACC<Real>&mu,const ACC<Real>&Phi,const ACC<Real>&d2phi,
               const ACC<int>& geometry, const ACC<int>&nodeType,
               const ACC<Real>& coordinates, const Real *Aconst, const Real *Kappa){
    Real phi{Phi(0,0)};
    Real phi3 = {phi*phi*phi};
    Real d2phix = {d2phi(0,0,0)};
    Real d2phiy = {d2phi(1,0,0)};
    const Real K={*Kappa};
    const Real A={*Aconst};
    //std::cout<<A<<" "<<K<<" "<<phi<<" "<<d2phix<<" "<<d2phiy<<"   ";
    mu(0,0)=-A*phi+A*phi3-K*(d2phix+d2phiy);
}



//Update surface tension force between fluid components
//TODO: this kernel name is reserved for an app. Also, this function should be named more close to multiphase
void KerUpdateMacroBodyForce(ACC<Real>& force, const ACC<Real>& Phi,
                             const ACC<Real>& dmu, const ACC<int>& geometry,
                             const ACC<int>&nodeType,
                             const ACC<Real>& coordinates, const int* idx) {
    Real phi = {Phi(0,0)};
    Real dmux = {dmu(0,0,0)};
    Real dmuy = {dmu(1,0,0)};
    VertexType vt = (VertexType)nodeType(0, 0);
    VertexGeometryType vg = (VertexGeometryType)geometry(0, 0);
    //Set gradient of chemical potential perpendicular to wall = 0
    if (vt == VertexType::Wall) {

        switch (vg) {
            case VG_IP: {
                dmux=0;
            } break;
            case VG_IM: {
                dmux=0;
            } break;
            case VG_JP: {
                dmuy=0;
            } break;
            case VG_JM: {
                dmuy=0;
            } break;
            case VG_IPJP_O: {
                dmux=0;
                dmuy=0;
            } break;
            case VG_IPJM_O: {
                dmux=0;
                dmuy=0;
            } break;
            case VG_IMJP_O: {
                dmux=0;
                dmuy=0;
            } break;
            case VG_IMJM_O: {
                dmux=0;
                dmuy=0;
            } break;
            case VG_IPJP_I: {
                dmux=0;
                dmuy=0;
            } break;
            case VG_IPJM_I: {
                dmux=0;
                dmuy=0;
            } break;
            case VG_IMJP_I: {
                dmux=0;
                dmuy=0;
            } break;
            case VG_IMJM_I: {
                dmux=0;
                dmuy=0;
            } break;
        }
    }

    force(0,0,0) = -phi*dmux;
    force(1,0,0) = -phi*dmuy;

}


void KerCalcDensity(ACC<Real>& Rho, const ACC<Real>& f,
                    const ACC<int>& nodeType, const int* lattIdx) {
#ifdef OPS_2D
    VertexType vt = (VertexType)nodeType(0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real rho{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {

            rho += f(xiIdx, 0, 0);
        }

#ifdef CPU
        if (isnan(rho) || rho <= -10 || isinf(rho)) {
            ops_printf(
                "Error! Density %f becomes invalid！Something "
                "wrong...",
                rho);
            assert(!(isnan(rho) || rho <= -10 || isinf(rho)));
        }
#endif
        Rho(0, 0) = rho;
    }
#endif  // OPS_2D
}

void KerCalcU(ACC<Real>& U, const ACC<Real>& f, const ACC<int>& nodeType,
              const ACC<Real>& Rho, const int* lattIdx) {
#ifdef OPS_2D
    VertexType vt = (VertexType)nodeType(0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real u{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            u += CS * XI[xiIdx * LATTDIM] * f(xiIdx, 0, 0);
        }
        u /= Rho(0, 0);
#ifdef CPU
        if (isnan(u) || isinf(u)) {
            ops_printf(
                "Error! Velocity U=%f becomes invalid! "
                "Maybe something wrong...\n",
                u);
            assert(!(isnan(u) || isinf(u)));
        }
#endif
        U(0, 0) = u;
    }
#endif  // OPS_2D
}

void KerCalcV(ACC<Real>& V, const ACC<Real>& f, const ACC<int>& nodeType,
              const ACC<Real>& Rho, const int* lattIdx) {
#ifdef OPS_2D
    VertexType vt = (VertexType)nodeType(0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real v{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            v += CS * XI[xiIdx * LATTDIM + 1] * f(xiIdx, 0, 0);
        }
        v /= Rho(0, 0);
#ifdef CPU
        if (isnan(v) || isinf(v)) {
            ops_printf(
                "Error! Velocity V=%f becomes invalid! "
                "Maybe something wrong...\n",
                v);
            assert(!(isnan(v) || isinf(v)));
        }
#endif
        V(0, 0) = v;
    }
#endif  // OPS_2D
}

void KerCalcUForce(ACC<Real>& U, const ACC<Real>& f, const ACC<int>& nodeType,
                   const ACC<Real>& coordinates, const ACC<Real>& acceleration,
                   const ACC<Real>& Rho, const Real* dt, const int* lattIdx) {
#ifdef OPS_2D
    const Real x{coordinates(0, 0, 0)};
    const Real y{coordinates(1, 0, 0)};
    VertexType vt = (VertexType)nodeType(0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real u{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            u += CS * XI[xiIdx * LATTDIM] * f(xiIdx, 0, 0);
        }
        u /= Rho(0, 0);
        if (VertexType::Fluid == vt || VertexType::MDPeriodic == vt) {
            u += ((*dt) * acceleration(0, 0, 0) / (2*Rho(0,0)));
        }
#ifdef CPU
        if (isnan(u) || isinf(u)) {
            ops_printf(
                "Error! Velocity U=%f becomes invalid! Maybe something wrong "
                "at x=%f y=%f, a=%f, rho=%f\n",
                u, x, y, acceleration(0,0,0), Rho(0,0));
            assert(!(isnan(u) || isinf(u)));
        }
#endif
        U(0, 0) = u;
    }
#endif  // OPS_2D
}

void KerCalcVForce(ACC<Real>& V, const ACC<Real>& f, const ACC<int>& nodeType,
                   const ACC<Real>& coordinates, const ACC<Real>& acceleration,
                   const ACC<Real>& Rho, const Real* dt, const int* lattIdx) {
#ifdef OPS_2D
    const Real x{coordinates(0, 0, 0)};
    const Real y{coordinates(1, 0, 0)};
    VertexType vt = (VertexType)nodeType(0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real v{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            v += CS * XI[xiIdx * LATTDIM + 1] * f(xiIdx, 0, 0);
        }
        v /= Rho(0, 0);
        if (VertexType::Fluid == vt || VertexType::MDPeriodic == vt) {
            v += ((*dt) * acceleration(1, 0, 0) / (2*Rho(0,0)));
        }
#ifdef CPU
        if (isnan(v) || isinf(v)) {
            ops_printf(
                "Error! Velocity V=%f becomes invalid! Maybe something wrong "
                "at x=%f y=%f\n",
                v, x, y);
            assert(!(isnan(v) || isinf(v)));
        }
#endif
        V(0, 0) = v;
    }
#endif  // OPS_2D
}

/*!
 * If a Newton-Cotes quadrature is used, it can be converted to the way
 * similar to the Gauss-Hermite quadrature *
 */

void KerInitialiseBGK2nd(ACC<Real>& f, const ACC<int>& nodeType,
                         const ACC<Real>& Rho, const ACC<Real>& U,
                         const ACC<Real>& V, const int* lattIdx) {
#ifdef OPS_2D
    VertexType vt = (VertexType)nodeType(0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real rho{Rho(0, 0)};
        Real u{U(0, 0)};
        Real v{V(0, 0)};
        const Real T{1};
        const int polyOrder{2};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            f(xiIdx, 0, 0) = CalcBGKFeq(xiIdx, rho, u, v, T, polyOrder);
#ifdef CPU
            const Real res{f(xiIdx, 0, 0)};
            if (isnan(res) || res <= 0 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function %f becomes invalid at the "
                    "lattice %i\n",
                    res, xiIdx);
                assert(!(isnan(res) || res <= 0 || isinf(res)));
            }
#endif  // CPU
        }
    }
#endif  // OPS_2D
}

void KerCollideBGKIsothermal(ACC<Real>& fStage, const ACC<Real>& f,
                             const ACC<Real>& coordinates,
                             const ACC<int>& nodeType, const ACC<Real>& Rho,
                             const ACC<Real>& U, const ACC<Real>& V,
                              const Real* tauRef,
                             const Real* dt, const int* lattIdx) {
#ifdef OPS_2D
    VertexType vt = (VertexType)nodeType(0, 0);
    // collisionRequired: means if collision is required at boundary
    // e.g., the ZouHe boundary condition explicitly requires collision
    bool collisionRequired = (vt != VertexType::ImmersedSolid);
    if (collisionRequired) {
        Real rho{Rho(0, 0)};
        Real u{U(0, 0)};
        Real v{V(0, 0)};
        const Real T{1};
        const int polyOrder{2};
        Real tau = (*tauRef);
        Real dtOvertauPlusdt = (*dt) / (tau + 0.5 * (*dt));
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            const Real feq{CalcBGKFeq(xiIndex, rho, u, v, T, polyOrder)};
            if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic) {
                //fStage(xiIndex, 0, 0) =
                //    feq + (1 - dtOvertauPlusdt) * (f(xiIndex, 0, 0) - feq) +
                //    tau * dtOvertauPlusdt * fStage(xiIndex, 0, 0);
                fStage(xiIndex, 0, 0) =
                    feq + (1 - dtOvertauPlusdt) * (f(xiIndex, 0, 0) - feq) +
                    tau * dtOvertauPlusdt * fStage(xiIndex, 0, 0);
            } else {
                fStage(xiIndex, 0, 0) =
                    feq + (1 - dtOvertauPlusdt) * (f(xiIndex, 0, 0) - feq);
            }
#ifdef CPU
            const Real res{fStage(xiIndex, 0, 0)};
            if (isnan(res) || res <= 0 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function = %e becomes invalid at  "
                    "the lattice %i where feq=%e and rho=%e u=%e v=%e at "
                    "x=%e y=%e\n",
                    res, xiIndex, feq, rho, u, v, coordinates(0, 0, 0),
                    coordinates(1, 0, 0));
                assert(!(isnan(res) || res <= 0 || isinf(res)));
            }
#endif  // CPU
        }
    }
#endif  // OPS_2D
}

void KerInitialiseBGKADF(ACC<Real>& g, const ACC<int>& nodeType,
                         ACC<Real>& Conc, const ACC<Real>& U,
                         const ACC<Real>& V,const ACC<Real>& acceleration,const Real* dt, const int* lattIdx) {
#ifdef OPS_2D
    VertexType vt = (VertexType)nodeType(0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real C{Conc(0,0)};
        Real u{U(0, 0)};
        Real v{V(0, 0)};
        const Real T{1};
        const int polyOrder{2};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            g(xiIdx, 0, 0) = CalcBGKFeqAD(xiIdx, C, u, v, T, polyOrder);

#ifdef CPU
            const Real res{g(xiIdx, 0, 0)};
            if (isnan(res) || res <= -100.01 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function %f becomes invalid at the "
                    "lattice %i\n",
                    res, xiIdx);
                assert(!(isnan(res) || res <= -10.01 || isinf(res)));
            }
#endif  // CPU
        }
    }
#endif  // OPS_2D
}

void KerInitialiseBGKFEF(ACC<Real>& g, const ACC<int>& nodeType,
                         ACC<Real>& Conc, const ACC<Real>& U,
                         const ACC<Real>& V,const ACC<Real>& acceleration,const Real* dt, const int* lattIdx) {
#ifdef OPS_2D
    VertexType vt = (VertexType)nodeType(0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real C{Conc(0,0)};
        Real u{U(0, 0)};
        Real v{V(0, 0)};
        Real ff[]{acceleration(0, 0, 0), acceleration(1, 0, 0)};
        const Real T{1};
        const int polyOrder{2};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            g(xiIdx, 0, 0) = CalcBGKFeqFE(xiIdx, C, u, v, ff, dt, T, polyOrder);

#ifdef CPU
            const Real res{g(xiIdx, 0, 0)};
            if (isnan(res) || res <= -100.01 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function %f becomes invalid at the "
                    "lattice %i\n",
                    res, xiIdx);
                assert(!(isnan(res) || res <= -10.01 || isinf(res)));
            }
#endif  // CPU
        }
    }
#endif  // OPS_2D
}

void KerInitialiseBGKFEG(ACC<Real>& g, const ACC<int>& nodeType,const ACC<int>& geometry,
                         ACC<Real>& Conc, const ACC<Real>& U,
                         const ACC<Real>& V, const ACC<Real>& Mu, const Real* dt, const int* lattIdx) {
#ifdef OPS_2D
    VertexType vt = (VertexType)nodeType(0, 0);
    if (vt != VertexType::ImmersedSolid) {

        Real C{Conc(0,0)};
        Real u{U(0, 0)};
        Real v{V(0, 0)};
        Real mu{Mu(0,0)};
        const int polyOrder{2};
        Real S{0.0};

        for (int xiIdx = lattIdx[0]+1; xiIdx <= lattIdx[1]; xiIdx++) {
            g(xiIdx, 0, 0) = CalcBGKGeqFE(xiIdx, C, u, v, mu, dt, polyOrder);
            S+=g(xiIdx,0,0);
        g(lattIdx[0],0,0)=C-S;



#ifdef CPU
            const Real res{g(xiIdx, 0, 0)};
            if (isnan(res) || res <= -1.01 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function %f becomes invalid at the "
                    "lattice %i\n",
                    res, xiIdx);
                assert(!(isnan(res) || res <= -1.01 || isinf(res)));
            }
#endif  // CPU
        }
    }
#endif  // OPS_2D
}


void KerCollideBGKADF(ACC<Real>& gStage, const ACC<Real>& g,
                             const ACC<Real>& coordinates,
                             const ACC<int>& nodeType, const ACC<Real>& Conc,
                             const ACC<Real>& U, const ACC<Real>& V,
                             const ACC<Real>& acceleration,
                              const Real* tauRef,
                             const Real* dt, const int* lattIdx) {
#ifdef OPS_2D
    VertexType vt = (VertexType)nodeType(0, 0);
    // collisionRequired: means if collision is required at boundary
    // e.g., the ZouHe boundary condition explicitly requires collision
    bool collisionRequired = (vt != VertexType::ImmersedSolid);
    if (collisionRequired) {
        Real u{U(0, 0)};
        Real v{V(0, 0)};
        Real C{Conc(0,0)};

        const Real T{1};
        const int polyOrder{2};
        Real tau = (*tauRef);
        Real dtOvertauPlusdt = (*dt) / (tau);
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            const Real geq{CalcBGKFeqAD(xiIndex, C, u, v, T, polyOrder)};
            if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic) {

                gStage(xiIndex, 0, 0) =
                    geq - (dtOvertauPlusdt-1) * (g(xiIndex, 0, 0) - geq) +
                    (1.0-1.0/(2.0*tau)) * gStage(xiIndex, 0, 0);
            } else {
                gStage(xiIndex, 0, 0) =
                    geq - (dtOvertauPlusdt-1) * (g(xiIndex, 0, 0) - geq);
            }
#ifdef CPU
            const Real res{gStage(xiIndex, 0, 0)};
            if (isnan(res) || res <= -100.01 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function = %e becomes invalid at  "
                    "the lattice %i where feq=%e and rho=%e u=%e v=%e at "
                    "x=%e y=%e\n",
                    res, xiIndex, geq, C, u, v, coordinates(0, 0, 0),
                    coordinates(1, 0, 0));
                std::cout << C << "Conc\n";
                assert(!(isnan(res) || res <= -10.01 || isinf(res)));
            }
#endif  // CPU
        }

    }
#endif  // OPS_2D
}



void KerCollideBGKFEF(ACC<Real>& gStage, const ACC<Real>& g,
                             const ACC<Real>& coordinates,
                             const ACC<int>& nodeType, const ACC<Real>& Conc,
                             const ACC<Real>& Phi,
                             const ACC<Real>& U, const ACC<Real>& V,
                             const ACC<Real>& acceleration,
                              const Real* tauRef,
                             const Real* dt, const int* lattIdx) {
#ifdef OPS_2D
    VertexType vt = (VertexType)nodeType(0, 0);
    // collisionRequired: means if collision is required at boundary
    // e.g., the ZouHe boundary condition explicitly requires collision
    bool collisionRequired = (vt != VertexType::ImmersedSolid);
    if (collisionRequired) {
        Real u{U(0, 0)};
        Real v{V(0, 0)};
        Real C{Conc(0,0)};
        Real phi{Phi(0,0)};
        Real ff[]{acceleration(0, 0, 0), acceleration(1, 0, 0)};
        const Real T{1};
        const int polyOrder{2};
        //Real tau = (*tauRef);
        //Real dtOvertauPlusdt = (*dt) / (tau);
        Real css=1/sqrt(3);
        Real taul=1.0;
        Real taug=1.0;
        Real Dt{*dt};
        Real vl=css*css*(taul-Dt/2.0);
        Real vg=css*css*(taug-Dt/2.0);
        //Real vl=1.17;
        //Real vg=0.05;
        Real V=vg+(phi+1.0)/2.0*(vl-vg);
        Real W=2.0/(6.0*V+1.0);
        //Real dtOvertauPlusdt = (*dt) / (W);
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            const Real geq{CalcBGKFeqFE(xiIndex, C, u, v, ff, dt, T, polyOrder)};
            if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic) {

                gStage(xiIndex, 0, 0) =
                    geq - (W-1) * (g(xiIndex, 0, 0) - g                gStage(xiIndex, 0, 0) =
                    geq - (W-1) * (g(xiIndex, 0, 0) - geq);eq) +
                    (*dt)*(1.0-(*dt*W)/(2.0)) * gStage(xiIndex, 0, 0);
            } else {

            }
#ifdef CPU
            const Real res{gStage(xiIndex, 0, 0)};
            if (isnan(res) || res <= -100.01 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function = %e becomes invalid at  "
                    "the lattice %i where feq=%e and rho=%e u=%e v=%e at "
                    "x=%e y=%e force=%e,%e\n",
                    res, xiIndex, geq, C, u, v, coordinates(0, 0, 0),
                    coordinates(1, 0, 0),ff[0],ff[1]);
                std::cout << C << "Conc\n";
                assert(!(isnan(res) || res <= -10.01 || isinf(res)));
            }
#endif  // CPU
        }

    }
#endif  // OPS_2D
}


void KerCollideBGKFEG(ACC<Real>& gStage, const ACC<Real>& g,
                             const ACC<Real>& coordinates,
                             const ACC<int>& nodeType,const ACC<int>& geometry, const ACC<Real>& Conc,
                             const ACC<Real>& U, const ACC<Real>& V, const ACC<Real>& Mu,
                              const Real* tauRef,
                             const Real* dt, const int* lattIdx) {
#ifdef OPS_2D
    VertexType vt = (VertexType)nodeType(0, 0);
    // collisionRequired: means if collision is required at boundary
    // e.g., the ZouHe boundary condition explicitly requires collision
    bool collisionRequired = (vt != VertexType::ImmersedSolid);

    if (collisionRequired) {
        Real u{U(0, 0)};
        Real v{V(0, 0)};
        Real C{Conc(0,0)};
        Real mu{Mu(0,0)};
        const int polyOrder{2};
        Real tau = (*tauRef);
        Real dtOvertauPlusdt = (*dt) / (tau);
        Real S{0.0};

        for (int xiIndex = lattIdx[0]+1; xiIndex <= lattIdx[1]; xiIndex++) {
            const Real geq{CalcBGKGeqFE(xiIndex, C, u, v, mu, dt, polyOrder)};
            if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic) {
                gStage(xiIndex, 0, 0) =
                    geq - (dtOvertauPlusdt-1) * (g(xiIndex, 0, 0) - geq);
            } else {
                gStage(xiIndex, 0, 0) =
                    geq - (dtOvertauPlusdt-1) * (g(xiIndex, 0, 0) - geq);
            }
            S+=gStage(xiIndex,0,0);

#ifdef CPU
            const Real res{gStage(xiIndex, 0, 0)};
            if (isnan(res) || res <= -100.01 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function = %e becomes invalid at  "
                    "the lattice %i where feq=%e and rho=%e u=%e v=%e at "
                    "x=%e y=%e\n",
                    res, xiIndex, geq, Conc, u, v, coordinates(0, 0, 0),
                    coordinates(1, 0, 0));
                std::cout << C << "Conc\n";
                assert(!(isnan(res) || res <= -1.01 || isinf(res)));
            }
#endif  // CPU
        }

        gStage(lattIdx[0],0,0)=C-S;
        //std::cout<<gStage(lattIdx[0],0,0)<<"   ";
    }
#endif  // OPS_2D
}



void KerCollideMRTFEF(ACC<Real>& gStage, const ACC<Real>& g,
                             const ACC<Real>& coordinates,
                             const ACC<int>& nodeType, const ACC<Real>& Conc,
                             const ACC<Real>& U, const ACC<Real>& V,
                             const ACC<Real>& acceleration,
                              const Real* tauRef,
                             const Real* dt, const int* lattIdx) {
#ifdef OPS_2D
    VertexType vt = (VertexType)nodeType(0, 0);
    // collisionRequired: means if collision is required at boundary
    // e.g., the ZouHe boundary condition explicitly requires collision
    bool collisionRequired = (vt != VertexType::ImmersedSolid);
    if (collisionRequired) {
        Real u{U(0, 0)};
        Real v{V(0, 0)};
        Real C{Conc(0,0)};
        Real ff[]{acceleration(0, 0, 0), acceleration(1, 0, 0)};
        //STORE SEPERATELY


        ///////
        Real taul=1.0;
        Real taug=1.0;
        Real Dt{*dt};
        Real vl=CS*CS*(taul-Dt/2);
        Real vg=CS*CS*(taug-Dt/2);
        Real V=vg+(C+1.0)/2.0*(vl-vg);
        Real W=2.0/(6.0*V+1.0);

        std::vector<std::vector<Real> > S = {
            {0, 0, 0, 0, 0, 0, 0, 0, 0} ,
            {0, Dt, 0, 0, 0, 0, 0, 0, 0} ,
            {0, 0, Dt, 0, 0, 0, 0, 0, 0} ,
            {0, 0, 0, 0, 0, 0, 0, 0, 0} ,
            {0, 0, 0, 0, Dt, 0, 0, 0, 0} ,
            {0, 0, 0, 0, 0, 0, 0, 0, 0} ,
            {0, 0, 0, 0, 0, 0, Dt, 0, 0} ,
            {0, 0, 0, 0, 0, 0, 0, Dt*W, 0} ,
            {0, 0, 0, 0, 0, 0, 0, 0, Dt*W} ,
        };
        const Real T{1};
        /*
        Real G1{g(lattIdx[0], 0, 0)};
        Real G2{};
        Real G3{};
        Real G4{};
        Real G2{};
        Real G2{};
        Real G2{};
        Real G2{};
        Real G2{};
        */
        const int polyOrder{2};
        Real tau = (*tauRef);
        Real dtOvertauPlusdt = (Dt) / (tau);
        std::vector<std::vector<Real> > gm = {
            {(Real)g(lattIdx[0], 0, 0)-(Real)CalcBGKFeqFE(lattIdx[0], C, u, v, ff, dt, T, polyOrder)},
            {(Real)g(lattIdx[0]+1, 0, 0)-(Real)CalcBGKFeqFE(lattIdx[0]+1, C, u, v, ff, dt, T, polyOrder)},
            {(Real)g(lattIdx[0]+2, 0, 0)-(Real)CalcBGKFeqFE(lattIdx[0]+2, C, u, v, ff, dt, T, polyOrder)},
            {(Real)g(lattIdx[0]+3, 0, 0)-(Real)CalcBGKFeqFE(lattIdx[0]+3, C, u, v, ff, dt, T, polyOrder)},
            {(Real)g(lattIdx[0]+4, 0, 0)-(Real)CalcBGKFeqFE(lattIdx[0]+4, C, u, v, ff, dt, T, polyOrder)},
            {(Real)g(lattIdx[0]+5, 0, 0)-(Real)CalcBGKFeqFE(lattIdx[0]+5, C, u, v, ff, dt, T, polyOrder)},
            {(Real)g(lattIdx[0]+6, 0, 0)-(Real)CalcBGKFeqFE(lattIdx[0]+6, C, u, v, ff, dt, T, polyOrder)},
            {(Real)g(lattIdx[0]+7, 0, 0)-(Real)CalcBGKFeqFE(lattIdx[0]+7, C, u, v, ff, dt, T, polyOrder)},
            {(Real)g(lattIdx[0]+8, 0, 0)-(Real)CalcBGKFeqFE(lattIdx[0]+8, C, u, v, ff, dt, T, polyOrder)}};
        std::vector<std::vector<Real> > Force = {
            {0.5*gStage(lattIdx[0],0,0)},
            {0.5*gStage(lattIdx[0]+1,0,0)},
            {0.5*gStage(lattIdx[0]+2,0,0)},
            {0.5*gStage(lattIdx[0]+3,0,0)},
            {0.5*gStage(lattIdx[0]+4,0,0)},
            {0.5*gStage(lattIdx[0]+5,0,0)},
            {0.5*gStage(lattIdx[0]+6,0,0)},
            {0.5*gStage(lattIdx[0]+7,0,0)},
            {0.5*gStage(lattIdx[0]+8,0,0)},
        };
        std::vector<std::vector<Real> > Mgm=KerMultMatrix(M,gm);
        std::vector<std::vector<Real> > SMgm=KerMultMatrix(S,Mgm);
        std::vector<std::vector<Real> > MF=KerMultMatrix(M,Force);
        std::vector<std::vector<Real> > SMF=KerMultMatrix(S,MF);
        std::vector<std::vector<Real> > MinvSMgm=KerMultMatrix(Minv,SMgm);
        std::vector<std::vector<Real> > MinvSMF=KerMultMatrix(Minv,SMF);


        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic) {
                std::cout<<Mgm[0][xiIndex]<<" "<<SMgm[0][xiIndex]<<" "<<MinvSMgm[0][xiIndex]<<" ";
                gStage(xiIndex, 0, 0) =
                    g(xiIndex, 0, 0) + gStage(xiIndex, 0, 0) - MinvSMgm[0][xiIndex] - MinvSMF[0][xiIndex];

            } else {
                gStage(xiIndex, 0, 0) =
                    g(xiIndex, 0, 0) - MinvSMgm[0][xiIndex];
            }
#ifdef CPU
            const Real res{gStage(xiIndex, 0, 0)};
            if (isnan(res) || res <= -100.01 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function = %e becomes invalid at  "
                    "the lattice %i where rho=%e u=%e v=%e at "
                    "x=%e y=%e force=%e,%e\n",
                    res, xiIndex, C, u, v, coordinates(0, 0, 0),
                    coordinates(1, 0, 0),ff[0],ff[1]);
                std::cout << C << "Conc\n";
                //assert(!(isnan(res) || res <= -10.01 || isinf(res)));
            }
#endif  // CPU
        }

    }
#endif  // OPS_2D
}



void KerCollideBGKThermal(ACC<Real>& fStage, const ACC<Real>& f,
                          const ACC<int>& nodeType, const ACC<Real>& Rho,
                          const ACC<Real>& U, const ACC<Real>& V,
                          const ACC<Real>& Temperature, const Real* tauRef,
                          const Real* dt, const int* lattIdx) {
#ifdef OPS_2D
    VertexType vt = (VertexType)nodeType(0, 0);
    // collisionRequired: means if collision is required at boundary
    // e.g., the ZouHe boundary condition explicitly requires collision
    bool collisionRequired = (vt != VertexType::ImmersedSolid);
    if (collisionRequired) {
        Real rho{Rho(0, 0)};
        Real u{U(0, 0)};
        Real v{V(0, 0)};
        Real T{Temperature(0, 0)};
        const int polyOrder{4};
        Real tau = (*tauRef) / (rho * sqrt(T));
        Real dtOvertauPlusdt = (*dt) / (tau + 0.5 * (*dt));
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            const Real feq{CalcBGKFeq(xiIndex, rho, u, v, T, polyOrder)};
            if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic) {
                fStage(xiIndex, 0, 0) =
                    f(xiIndex, 0, 0) -
                    dtOvertauPlusdt * (f(xiIndex, 0, 0) - feq) +
                    tau * dtOvertauPlusdt * fStage(xiIndex, 0, 0);
            } else {
                fStage(xiIndex, 0, 0) =
                    f(xiIndex, 0, 0) -
                    dtOvertauPlusdt * (f(xiIndex, 0, 0) - feq);
            }
#ifdef CPU
            const Real res{fStage(xiIndex, 0, 0)};
            if (isnan(res) || res <= 0 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function %f becomes invalid at the "
                    "lattice %i\n",
                    res, xiIndex);
                assert(!(isnan(res) || res <= 0 || isinf(res)));
            }
#endif  // CPU
        }
    }
#endif  // OPS_2D
}

void KerCalcBodyForce1ST(ACC<Real>& fStage, const ACC<Real>& acceration,
                         const ACC<Real>& Rho, const ACC<int>& nodeType,
                         const int* lattIdx) {
#ifdef OPS_2D

    VertexType vt = (VertexType)nodeType(0, 0);
    if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic) {
        Real rho{Rho(0, 0)};
        Real g[]{acceration(0, 0, 0), acceration(1, 0, 0)};
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            const Real bodyForce{CalcBodyForce(xiIndex, rho, g)};
#ifdef CPU
            if (isnan(bodyForce) || isinf(bodyForce)) {
                ops_printf(
                    "Error! Body force  %f becomes invalid  at  the lattice "
                    "%i\n",
                    bodyForce, xiIndex);
                assert(!(isnan(bodyForce) || isinf(bodyForce)));
            }
#endif
            fStage(xiIndex, 0, 0) = bodyForce;
        }
    }
#endif  // OPS_2D
}

void KerCalcGuoForce(ACC<Real>& fStage, const ACC<Real>& acceration,
                         const ACC<Real>& Rho, const ACC<Real>& U, const ACC<Real>& V,
                         const Real* dt, const ACC<int>& nodeType,
                         const int* lattIdx) {
#ifdef OPS_2D

    VertexType vt = (VertexType)nodeType(0, 0);
    if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic) {
        Real u{U(0,0)};
        Real v{V(0,0)};
        Real rho{Rho(0, 0)};
        Real g[]{acceration(0, 0, 0), acceration(1, 0, 0)};
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            const Real bodyForce{CalcGuoForce(xiIndex, rho, g, u, v, dt)};
#ifdef CPU
            if (isnan(bodyForce) || isinf(bodyForce)) {
                ops_printf(
                    "Error! Body force  %f becomes invalid  at  the lattice "
                    "%i\n",
                    bodyForce, xiIndex);
                assert(!(isnan(bodyForce) || isinf(bodyForce)));
            }
#endif
            if (bodyForce>=0){
                fStage(xiIndex, 0, 0) = bodyForce;
            }
            else {
                fStage(xiIndex, 0, 0) = bodyForce;
            }

        }
    }
#endif  // OPS_2D
}

void KerCalcBodyForceNone(ACC<Real>& fStage, const ACC<Real>& acceration,
                          const ACC<int>& nodeType, const int* lattIdx) {
#ifdef OPS_2D
    VertexType vt = (VertexType)nodeType(0, 0);
    if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic) {
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            fStage(xiIndex, 0, 0) = 0;
        }
    }
#endif  // OPS_2D
}
#endif // OPS_2D outter

#ifdef OPS_3D
void KerInitialiseBGK2nd3D(ACC<Real>& f, const ACC<int>& nodeType,
                           const ACC<Real>& Rho, const ACC<Real>& U,
                           const ACC<Real>& V, const ACC<Real>& W,
                           const int* lattIdx) {
#ifdef OPS_3D
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real rho{Rho(0, 0, 0)};
        Real u{U(0, 0, 0)};
        Real v{V(0, 0, 0)};
        Real w{W(0, 0, 0)};
        const Real T{1};
        const int polyOrder{2};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            f(xiIdx, 0, 0, 0) = CalcBGKFeq(xiIdx, rho, u, v, w, T, polyOrder);
#ifdef CPU
            const Real res{f(xiIdx, 0, 0, 0)};
            if (isnan(res) || res <= 0 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function %f becomes invalid at the "
                    "lattice %i\n",
                    res, xiIdx);
                assert(!(isnan(res) || res <= 0 || isinf(res)));
            }
#endif  // CPU
        }
    }
#endif  // OPS_3D
}

void KerCollideBGKIsothermal3D(ACC<Real>& fStage, const ACC<Real>& f,
                               const ACC<Real>& coordinates,
                               const ACC<int>& nodeType, const ACC<Real>& Rho,
                               const ACC<Real>& U, const ACC<Real>& V,
                               const ACC<Real>& W, const Real* tauRef,
                               const Real* dt, const int* lattIdx) {
#ifdef OPS_3D
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    // collisionRequired: means if collision is required at boundary
    // e.g., the ZouHe boundary condition explicitly requires collision
    bool collisionRequired = (vt != VertexType::ImmersedSolid);
    if (collisionRequired) {
        Real rho{Rho(0, 0, 0)};
        Real u{U(0, 0, 0)};
        Real v{V(0, 0, 0)};
        Real w{W(0, 0, 0)};
        const Real T{1};
        const int polyOrder{2};
        Real tau = (*tauRef);
        Real dtOvertauPlusdt = (*dt) / (tau + 0.5 * (*dt));
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            const Real feq{CalcBGKFeq(xiIndex, rho, u, v, w, T, polyOrder)};
            if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic) {
                fStage(xiIndex, 0, 0, 0) =
                    feq + (1 - dtOvertauPlusdt) * (f(xiIndex, 0, 0, 0) - feq) +
                    tau * dtOvertauPlusdt * fStage(xiIndex, 0, 0, 0);
            } else {
                fStage(xiIndex, 0, 0, 0) =
                    feq + (1 - dtOvertauPlusdt) * (f(xiIndex, 0, 0, 0) - feq);
            }
#ifdef CPU
            const Real res{fStage(xiIndex, 0, 0, 0)};
            if (isnan(res) || res <= 0 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function = %e becomes invalid at  "
                    "the lattice %i where feq=%e and rho=%e u=%e v=%e w=%e at "
                    "x=%e y=%e z=%e\n",
                    res, xiIndex, feq, rho, u, v, w, coordinates(0, 0, 0, 0),
                    coordinates(1, 0, 0, 0), coordinates(2, 0, 0, 0));
                assert(!(isnan(res) || res <= 0 || isinf(res)));
            }
#endif  // CPU
        }
    }
#endif  // OPS_3D
}

void KerCollideBGKThermal3D(ACC<Real>& fStage, const ACC<Real>& f,
                            const ACC<int>& nodeType, const ACC<Real>& Rho,
                            const ACC<Real>& U, const ACC<Real>& V,
                            const ACC<Real>& W, const ACC<Real>& Temperature,
                            const Real* tauRef, const Real* dt,
                            const int* lattIdx) {
#ifdef OPS_3D
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    // collisionRequired: means if collision is required at boundary
    // e.g., the ZouHe boundary condition explicitly requires collision
    bool collisionRequired = (vt != VertexType::ImmersedSolid);
    if (collisionRequired) {
        Real rho{Rho(0, 0, 0)};
        Real u{U(0, 0, 0)};
        Real v{V(0, 0, 0)};
        Real w{W(0, 0, 0)};
        Real T{Temperature(0, 0, 0)};
        const int polyOrder{4};
        Real tau = (*tauRef) / (rho * sqrt(T));
        Real dtOvertauPlusdt = (*dt) / (tau + 0.5 * (*dt));
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            const Real feq{CalcBGKFeq(xiIndex, rho, u, v, w, T, polyOrder)};
            if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic) {
                fStage(xiIndex, 0, 0, 0) =
                    f(xiIndex, 0, 0, 0) -
                    dtOvertauPlusdt * (f(xiIndex, 0, 0, 0) - feq) +
                    tau * dtOvertauPlusdt * fStage(xiIndex, 0, 0, 0);
            } else {
                fStage(xiIndex, 0, 0, 0) =
                    f(xiIndex, 0, 0, 0) -
                    dtOvertauPlusdt * (f(xiIndex, 0, 0, 0) - feq);
            }
#ifdef CPU
            const Real res{fStage(xiIndex, 0, 0, 0)};
            if (isnan(res) || res <= 0 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function %f becomes invalid at the "
                    "lattice %i\n",
                    res, xiIndex);
                assert(!(isnan(res) || res <= 0 || isinf(res)));
            }
#endif  // CPU
        }
    }
#endif  // OPS_3D
}


void KerInitialiseBGKADF3D(ACC<Real>& g, const ACC<int>& nodeType,
                         ACC<Real>& Conc, const ACC<Real>& U,
                         const ACC<Real>& V,const ACC<Real>& W,const ACC<Real>& acceleration,const Real* dt, const int* lattIdx) {
#ifdef OPS_3D
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real C{Conc(0,0,0)};
        Real u{U(0, 0, 0)};
        Real v{V(0, 0, 0)};
        Real w{W(0, 0, 0)};

        const Real T{1};
        const int polyOrder{2};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            g(xiIdx, 0, 0, 0) = CalcBGKFeqAD3D(xiIdx, C, u, v, w, T, polyOrder);
            //std::cout<<xiIdx<<" "<<g(xiIdx,0,0,0)<<" "<<u<<"   ";

#ifdef CPU
            const Real res{g(xiIdx, 0, 0, 0)};
            if (isnan(res) || res <= -100.01 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function %f becomes invalid at the "
                    "lattice %i\n",
                    res, xiIdx);
                assert(!(isnan(res) || res <= -10.01 || isinf(res)));
            }
#endif  // CPU
        }
    }
#endif  // OPS_3D
}


void KerInitialiseBGKFEF3D(ACC<Real>& g, const ACC<int>& nodeType,
                         ACC<Real>& Conc, const ACC<Real>& U,
                         const ACC<Real>& V,const ACC<Real>& W,const ACC<Real>& acceleration,const Real* dt, const int* lattIdx) {
#ifdef OPS_3D
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real C{Conc(0,0,0)};
        Real u{U(0, 0, 0)};
        Real v{V(0, 0, 0)};
        Real w{W(0, 0, 0)};
        Real ff[]{acceleration(0, 0, 0, 0), acceleration(1, 0, 0, 0), acceleration(2, 0, 0, 0)};
        const Real T{1};
        const int polyOrder{2};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            g(xiIdx, 0, 0, 0) = CalcBGKFeqFE3D(xiIdx, C, u, v, w, ff, dt, T, polyOrder);

#ifdef CPU
            const Real res{g(xiIdx, 0, 0, 0)};
            if (isnan(res) || res <= -100.01 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function %f becomes invalid at the "
                    "lattice %i\n",
                    res, xiIdx);
                assert(!(isnan(res) || res <= -10.01 || isinf(res)));
            }
#endif  // CPU
        }
    }
#endif  // OPS_3D
}

void KerInitialiseBGKFEG3D(ACC<Real>& g, const ACC<int>& nodeType,
                         ACC<Real>& Conc, const ACC<Real>& U,
                         const ACC<Real>& V, const ACC<Real>& W, const ACC<Real>& Mu, const Real* dt, const int* lattIdx) {
#ifdef OPS_3D
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real C{Conc(0,0,0)};
        Real u{U(0, 0, 0)};
        Real v{V(0, 0, 0)};
        Real w{W(0, 0, 0)};
        Real mu{Mu(0,0,0)};
        const int polyOrder{2};
        Real S{0.0};

        for (int xiIdx = lattIdx[0]+1; xiIdx <= lattIdx[1]; xiIdx++) {
            g(xiIdx, 0, 0, 0) = CalcBGKGeqFE3D(xiIdx, C, u, v, w, mu, dt, polyOrder);
            S+=g(xiIdx,0,0,0);
        g(lattIdx[0],0,0,0)=C-S;

#ifdef CPU
            const Real res{g(xiIdx, 0, 0, 0)};
            if (isnan(res) || res <= -1.01 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function %f becomes invalid at the "
                    "lattice %i\n",
                    res, xiIdx);
                assert(!(isnan(res) || res <= -1.01 || isinf(res)));
            }
#endif  // CPU
        }
    }
#endif  // OPS_3D
}


void KerCollideBGKADF3D(ACC<Real>& gStage, const ACC<Real>& g,
                             const ACC<Real>& coordinates,
                             const ACC<int>& nodeType, const ACC<Real>& Conc,
                             const ACC<Real>& U, const ACC<Real>& V, const ACC<Real>& W,
                             const ACC<Real>& acceleration,
                              const Real* tauRef,
                             const Real* dt, const int* lattIdx) {
#ifdef OPS_3D
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    // collisionRequired: means if collision is required at boundary
    // e.g., the ZouHe boundary condition explicitly requires collision
    //std::cout<<"test ";
    bool collisionRequired = (vt != VertexType::ImmersedSolid);
    if (collisionRequired) {
        Real u{U(0, 0, 0)};
        Real v{V(0, 0, 0)};
        Real w{W(0, 0, 0)};
        Real C{Conc(0,0, 0)};

        const Real T{1};
        const int polyOrder{2};
        Real tau = (*tauRef);
        Real test=0;
        Real dtOvertauPlusdt = (*dt) / (tau);

        //Real dtOvertauPlusdt = (*dt) / (tau + 0.5 * (*dt));
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            const Real geq{CalcBGKFeqAD3D(xiIndex, C, u, v, w, T, polyOrder)};
            test=gStage(xiIndex, 0, 0, 0);
            if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic) {

                gStage(xiIndex, 0, 0, 0) =
                    geq + (1 - dtOvertauPlusdt) * (g(xiIndex, 0, 0, 0) - geq) +
                    //tau * dtOvertauPlusdt * gStage(xiIndex, 0, 0, 0);
                    (1-1/(2*tau)) * gStage(xiIndex, 0, 0, 0);
                //std::cout<<xiIndex<<" "<<g(xiIndex,0,0,0)<<" "<<u<<"   ";
            } else {
                gStage(xiIndex, 0, 0, 0) =
                    geq + (1-dtOvertauPlusdt) * (g(xiIndex, 0, 0, 0) - geq);
            }
#ifdef CPU
            const Real res{gStage(xiIndex, 0, 0, 0)};
            if (isnan(res) || res <= -200.01 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function = %e becomes invalid at  "
                    "the lattice %i where feq=%e and rho=%e u=%e v=%e at "
                    "x=%e y=%e z=%e force=%e\n",
                    res, xiIndex, geq, C, u, v, coordinates(0, 0, 0, 0),
                    coordinates(1, 0, 0, 0),coordinates(2, 0, 0, 0),test);
                std::cout << C << "Conc\n";
                assert(!(isnan(res) || res <= 0 || isinf(res)));
            }
#endif  // CPU
        }

    }
#endif  // OPS_3D
}


void KerCollideBGKFEF3D(ACC<Real>& gStage, const ACC<Real>& g,
                             const ACC<Real>& coordinates,
                             const ACC<int>& nodeType, const ACC<Real>& Conc,
                             const ACC<Real>& U, const ACC<Real>& V, const ACC<Real>& W,
                             const ACC<Real>& acceleration,
                              const Real* tauRef,
                             const Real* dt, const int* lattIdx) {
#ifdef OPS_3D
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    // collisionRequired: means if collision is required at boundary
    // e.g., the ZouHe boundary condition explicitly requires collision
    bool collisionRequired = (vt != VertexType::ImmersedSolid);
    if (collisionRequired) {
        Real u{U(0, 0, 0)};
        Real v{V(0, 0, 0)};
        Real w{W(0, 0, 0)};
        Real C{Conc(0,0, 0)};
        Real ff[]{acceleration(0, 0, 0, 0), acceleration(1, 0, 0, 0), acceleration(2, 0, 0, 0)};
        const Real T{1};
        const int polyOrder{2};
        Real tau = (*tauRef);
        //Real dtOvertauPlusdt = (*dt) / (tau + 0.5 * (*dt));
        Real dtOvertauPlusdt = (*dt) / (tau);
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            const Real geq{CalcBGKFeqFE3D(xiIndex, C, u, v, w, ff, dt, T, polyOrder)};
            if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic) {

                gStage(xiIndex, 0, 0, 0) =
                    geq + (1 - dtOvertauPlusdt) * (g(xiIndex, 0, 0, 0) - geq) +
                    (1-1/(2*tau)) * gStage(xiIndex, 0, 0, 0);
            } else {
                gStage(xiIndex, 0, 0, 0) =
                    geq + (1 - dtOvertauPlusdt) * (g(xiIndex, 0, 0, 0) - geq);
            }
#ifdef CPU
            const Real res{gStage(xiIndex, 0, 0, 0)};
            if (isnan(res) || res <= -200.01 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function = %e becomes invalid at  "
                    "the lattice %i where feq=%e and rho=%e u=%e v=%e at "
                    "x=%e y=%e z=%e force=%e\n",
                    res, xiIndex, geq, C, u, v, coordinates(0, 0, 0, 0),
                    coordinates(1, 0, 0, 0),coordinates(2, 0, 0, 0),gStage(xiIndex,0,0,0));
                std::cout << C << "Conc\n";
                assert(!(isnan(res) || res <= 0 || isinf(res)));
            }
#endif  // CPU
        }

    }
#endif  // OPS_3D
}

void KerCollideBGKFEG3D(ACC<Real>& gStage, const ACC<Real>& g,
                             const ACC<Real>& coordinates,
                             const ACC<int>& nodeType, const ACC<Real>& Conc,
                             const ACC<Real>& U, const ACC<Real>& V, const ACC<Real>& W,
                             const ACC<Real>& Mu, const Real* tauRef,
                             const Real* dt, const int* lattIdx) {
#ifdef OPS_3D
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    // collisionRequired: means if collision is required at boundary
    // e.g., the ZouHe boundary condition explicitly requires collision
    bool collisionRequired = (vt != VertexType::ImmersedSolid);
    if (collisionRequired) {
        Real u{U(0, 0, 0)};
        Real v{V(0, 0, 0)};
        Real w{W(0, 0, 0)};
        Real C{Conc(0,0,0)};
        Real mu{Mu(0,0,0)};
        const int polyOrder{2};
        Real tau = (*tauRef);
        //Real dtOvertauPlusdt = (*dt) / (tau + 0.5 * (*dt));
        Real dtOvertauPlusdt = (*dt) / (tau);
        Real S{0.0};
        for (int xiIndex = lattIdx[0]+1; xiIndex <= lattIdx[1]; xiIndex++) {
            const Real geq{CalcBGKGeqFE3D(xiIndex, C, u, v, w, mu, dt, polyOrder)};
            if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic) {
                gStage(xiIndex, 0, 0, 0) =
                    geq + (1 - dtOvertauPlusdt) * (g(xiIndex, 0, 0, 0) - geq);
            } else {
                gStage(xiIndex, 0, 0, 0) =
                    geq + (1 - dtOvertauPlusdt) * (g(xiIndex, 0, 0, 0) - geq);
            }
            S+=gStage(xiIndex,0,0,0);

#ifdef CPU
            const Real res{gStage(xiIndex, 0, 0,0)};
            if (isnan(res) || res <= -10.01 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function = %e becomes invalid at  "
                    "the lattice %i where feq=%e and rho=%e u=%e v=%e at "
                    "x=%e y=%e z=%e\n",
                    res, xiIndex, geq, Conc, u, v, coordinates(0, 0, 0, 0),
                    coordinates(1, 0, 0, 0),coordinates(2, 0, 0, 0));
                std::cout << C << "Conc\n";
                assert(!(isnan(res) || res <= -10.01 || isinf(res)));
            }
#endif  // CPU
        }
        gStage(lattIdx[0],0,0,0)=C-S;
    }
#endif  // OPS_3D
}

void KerCalcGuoForce3D(ACC<Real>& fStage, const ACC<Real>& acceleration,
                         const ACC<Real>& Rho, const ACC<Real>& U, const ACC<Real>& V,
                         const ACC<Real>& W, const Real* dt, const ACC<int>& nodeType,
                         const int* lattIdx) {
#ifdef OPS_3D

    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic) {
        Real u{U(0,0,0)};
        Real v{V(0,0,0)};
        Real w{W(0,0,0)};
        Real rho{Rho(0, 0, 0)};
        Real g[]{acceleration(0, 0, 0, 0), acceleration(1, 0, 0, 0), acceleration(2, 0, 0, 0)};
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            const Real bodyForce{CalcGuoForce3D(xiIndex, rho, g, u, v, w, dt)};
#ifdef CPU
            if (isnan(bodyForce) || isinf(bodyForce)) {
                ops_printf(
                    "Error! Body force  %f becomes invalid  at  the lattice "
                    "%i\n",
                    bodyForce, xiIndex);
                assert(!(isnan(bodyForce) || isinf(bodyForce)));
            }
#endif
            if (bodyForce>=0){
                fStage(xiIndex, 0, 0, 0) = bodyForce;
            }
            else {
                fStage(xiIndex, 0, 0, 0) = bodyForce;
            }

        }
    }
#endif  // OPS_3D
}

void KerCalcBodyForce1ST3D(ACC<Real>& fStage, const ACC<Real>& acceration,
                           const ACC<Real>& Rho, const ACC<int>& nodeType,
                           const int* lattIdx) {
#ifdef OPS_3D

    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic) {
        Real rho{Rho(0, 0, 0)};
        Real g[]{acceration(0, 0, 0, 0), acceration(1, 0, 0, 0),
                 acceration(2, 0, 0, 0)};
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            const Real bodyForce{CalcBodyForce(xiIndex, rho, g)};
#ifdef CPU
            if (isnan(bodyForce) || isinf(bodyForce)) {
                ops_printf(
                    "Error! Body force  %f becomes invalid  at  the lattice "
                    "%i\n",
                    bodyForce, xiIndex);
                assert(!(isnan(bodyForce) || isinf(bodyForce)));
            }
#endif
            fStage(xiIndex, 0, 0, 0) = bodyForce;
        }
    }
#endif  // OPS_3D
}

void KerCalcBodyForceNone3D(ACC<Real>& fStage, const ACC<Real>& acceration,
                            const ACC<int>& nodeType, const int* lattIdx) {
#ifdef OPS_3D
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic) {
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            fStage(xiIndex, 0, 0, 0) = 0;
        }
    }
#endif  // OPS_3D
}

void KerCalcDensity3D(ACC<Real>& Rho, const ACC<Real>& f,
                      const ACC<int>& nodeType, const int* lattIdx) {
#ifdef OPS_3D
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real rho{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            rho += f(xiIdx, 0, 0, 0);
            //if (rho<=-10){
            //std::cout<<f(xiIdx,0,0,0)<<" "<<xiIdx<<"  ";
            //}
        }
#ifdef CPU
        if (isnan(rho) || rho <= -100000 || isinf(rho)) {
            ops_printf(
                "Error! Density %f becomes invalid！Something "
                "wrong...",
                rho);
            //assert(!(isnan(rho) || rho <= -1000 || isinf(rho)));
        }
#endif
        Rho(0, 0, 0) = rho;
        //if (rho<=-2){
        //std::cout<<rho<<" ";
        //}
    }
#endif // OPS_3D
}

void KerCalcU3D(ACC<Real>& U, const ACC<Real>& f, const ACC<int>& nodeType,
                const ACC<Real>& Rho, const int* lattIdx) {
#ifdef OPS_3D
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real u{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            u += CS * XI[xiIdx * LATTDIM] * f(xiIdx, 0, 0, 0);
            //std::cout<<f(xiIdx,0,0,0)<<" "<<u<<" ";
        }
        u /= Rho(0, 0, 0);
#ifdef CPU
        if (isnan(u) || isinf(u)) {
            ops_printf(
                "Error! Velocity U=%f becomes invalid! "
                "Maybe something wrong...\n",
                u);
            assert(!(isnan(u) || isinf(u)));
        }
#endif
        U(0, 0, 0) = u;
    }
#endif // OPS_3D
}

void KerCalcV3D(ACC<Real>& V, const ACC<Real>& f, const ACC<int>& nodeType,
                const ACC<Real>& Rho, const int* lattIdx) {
#ifdef OPS_3D
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real v{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            v += CS * XI[xiIdx * LATTDIM + 1] * f(xiIdx, 0, 0, 0);
        }
        v /= Rho(0, 0, 0);
#ifdef CPU
        if (isnan(v) || isinf(v)) {
            ops_printf(
                "Error! Velocity V=%f becomes invalid! "
                "Maybe something wrong...\n",
                v);
            assert(!(isnan(v) || isinf(v)));
        }
#endif
        V(0, 0, 0) = v;
    }
#endif // OPS_3D
}

void KerCalcW3D(ACC<Real>& W, const ACC<Real>& f, const ACC<int>& nodeType,
                const ACC<Real>& Rho, const int* lattIdx) {
#ifdef OPS_3D
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real w{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            w += CS * XI[xiIdx * LATTDIM + 2] * f(xiIdx, 0, 0, 0);
        }
        w /= Rho(0, 0, 0);
#ifdef CPU
        if (isnan(w) || isinf(w)) {
            ops_printf(
                "Error! Velocity W=%f becomes invalid! "
                "Maybe something wrong...\n",
                w);
            assert(!(isnan(w) || isinf(w)));
        }
#endif
        W(0, 0, 0) = w;
    }
#endif // OPS_3D
}

void KerCalcUForce3D(ACC<Real>& U, const ACC<Real>& f, const ACC<int>& nodeType,
                     const ACC<Real>& coordinates,
                     const ACC<Real>& acceleration, const ACC<Real>& Rho,
                     const Real* dt, const int* lattIdx) {
#ifdef OPS_3D
    const Real x{coordinates(0, 0, 0, 0)};
    const Real y{coordinates(1, 0, 0, 0)};
    const Real z{coordinates(2, 0, 0, 0)};
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real u{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            u += CS * XI[xiIdx * LATTDIM] * f(xiIdx, 0, 0, 0);
        }
        u /= Rho(0, 0, 0);
        if (VertexType::Fluid == vt || VertexType::MDPeriodic == vt) {
            u += ((*dt) * acceleration(0, 0, 0, 0) / 2);
        }
#ifdef CPU
        if (isnan(u) || isinf(u)) {
            ops_printf(
                "Error! Velocity U=%f becomes invalid! Maybe something wrong "
                "at x=%f y=%f z=%f, a=%f, rho=%f, f=%f,\n",
                u, x, y, z,acceleration(0,0,0,0),Rho(0, 0, 0));
            for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
                std::cout<<f(xiIdx, 0, 0, 0)<<"  ";
            }
            assert(!(isnan(u) || isinf(u)));
        }
#endif
        U(0, 0, 0) = u;
    }
#endif // OPS_3D
}

void KerCalcVForce3D(ACC<Real>& V, const ACC<Real>& f, const ACC<int>& nodeType,
                     const ACC<Real>& coordinates,
                     const ACC<Real>& acceleration, const ACC<Real>& Rho,
                     const Real* dt, const int* lattIdx) {
#ifdef OPS_3D
    const Real x{coordinates(0, 0, 0, 0)};
    const Real y{coordinates(1, 0, 0, 0)};
    const Real z{coordinates(2, 0, 0, 0)};
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real v{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            v += CS * XI[xiIdx * LATTDIM + 1] * f(xiIdx, 0, 0, 0);
        }
        v /= Rho(0, 0, 0);
        if (VertexType::Fluid == vt || VertexType::MDPeriodic == vt) {
            v += ((*dt) * acceleration(1, 0, 0, 0) / 2);
        }
#ifdef CPU
        if (isnan(v) || isinf(v)) {
            ops_printf(
                "Error! Velocity V=%f becomes invalid! Maybe something wrong "
                "at x=%f y=%f z=%f\n",
                v, x, y, z);
            assert(!(isnan(v) || isinf(v)));
        }
#endif
        V(0, 0, 0) = v;
    }
#endif // OPS_3D
}

void KerCalcWForce3D(ACC<Real>& W, const ACC<Real>& f, const ACC<int>& nodeType,
                     const ACC<Real>& coordinates,
                     const ACC<Real>& acceleration, const ACC<Real>& Rho,
                     const Real* dt, const int* lattIdx) {
#ifdef OPS_3D
    const Real x{coordinates(0, 0, 0, 0)};
    const Real y{coordinates(1, 0, 0, 0)};
    const Real z{coordinates(2, 0, 0, 0)};
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt != VertexType::ImmersedSolid) {
        Real w{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            w += CS * XI[xiIdx * LATTDIM + 2] * f(xiIdx, 0, 0, 0);
        }
        w /= Rho(0, 0, 0);
        if (VertexType::Fluid == vt || VertexType::MDPeriodic == vt) {
            w += ((*dt) * acceleration(2, 0, 0, 0) / 2);
        }
#ifdef CPU
        if (isnan(w) || isinf(w)) {
            ops_printf(
                "Error! Velocity W=%f becomes invalid! Maybe something wrong "
                "at x=%f y=%f z=%f\n",
                w, x, y, z);
            assert(!(isnan(w) || isinf(w)));
        }
#endif
        W(0, 0, 0) = w;
    }
#endif //OPS_3D
}


//Calculates first-order or second-order finite difference gradients
void KerCalcGradients3D(ACC<Real>&grad, const ACC<Real>&f,
                        const ACC<int>&nodeType, const ACC<int>&geometry,
                        const ACC<Real>& coordinates, const int *lattIdx,
                        const int *Order, const int *idx){
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    int order={*Order};
        //First order
        if (order == 1){
            if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic ||
                vt == VertexType::VirtualBoundary || vt==VertexType::Wall) {
                grad(0,0,0,0) = 0.5*(f(1, 0, 0)-f(-1, 0, 0));
                grad(1,0,0,0) = 0.5*(f(0, 1, 0)-f(0, -1, 0));
                grad(2,0,0,0) = 0.5*(f(0, 0, 1)-f(0, 0, -1));
            }
        }
        //Second order
        else if (order == 2){
            if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic ||
                vt == VertexType::VirtualBoundary || vt==VertexType::Wall) {
                grad(0,0,0,0) = f(-1, 0, 0)-2*f(0, 0, 0)+f(1, 0, 0);
                grad(1,0,0,0) = f(0, -1, 0)-2*f(0, 0, 0)+f(0, 1, 0);
                grad(2,0,0,0) = f(0, 0, -1)-2*f(0, 0, 0)+f(0, 0, 1);
            }
        }
        else{
            std::cout << "Only First and Second order finite difference gradients"
                      << " are supported\n";
        }

}

//Wetting boundary implementation
//TODO: this kernel name is repeated in the AdDiff app
void KerUpdateRhoWetting3D(ACC<Real>& rho,const ACC<int>& geometry,
                           const ACC<int>&nodeType,
                           const ACC<Real>& coordinates,
                           const Real *theta, const Real *Aconst, const Real *Kappa) {
    const Real teta={*theta};
    const Real alpha={acos(sin(teta)*sin(teta))};
    const Real Omga={2*((teta <= M_PI/2.0) - (M_PI/2.0 < teta))*sqrt(cos(alpha/3.0)
                     *(1-cos(alpha/3.0)))};
    const Real K={*Kappa};
    const Real A={*Aconst};
    const Real h={Omga*sqrt(K*A/2.0)};
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    VertexGeometryType vg = (VertexGeometryType)geometry(0, 0, 0);
    VertexGeometryType vg0 = (VertexGeometryType)geometry(1, 0, 0);
    VertexGeometryType vg1 = (VertexGeometryType)geometry(-1, 0, 0);
    VertexGeometryType vg2 = (VertexGeometryType)geometry(0, 1, 0);
    VertexGeometryType vg3 = (VertexGeometryType)geometry(0, -1, 0);
    VertexGeometryType vg4 = (VertexGeometryType)geometry(0, 0, 1);
    VertexGeometryType vg5 = (VertexGeometryType)geometry(0, 0, -1);
    if (vt == VertexType::Wall) {

        switch (vg) {
            case VG_IP: {
                if(vg2!=VG_IPJP_O&&vg3!=VG_IPJM_O
                   &&vg4!=VG_IPKP_O&&vg5!=VG_IPKM_O){
                    rho(-1,0,0)=rho(1,0,0)+2*h/K;
                }
            } break;
            case VG_IM: {
                if(vg2!=VG_IMJP_O&&vg3!=VG_IMJM_O
                   &&vg4!=VG_IMKP_O&&vg5!=VG_IMKM_O){
                    rho(1,0,0)=rho(-1,0,0)+2*h/K;
                }
            } break;
            case VG_JP: {
                if(vg0!=VG_IPJP_O&&vg1!=VG_IMJP_O
                   &&vg4!=VG_JPKP_O&&vg5!=VG_JPKM_O){
                    rho(0,-1,0)=rho(0,1,0)+2*h/K;
                }
            } break;
            case VG_JM: {
                if(vg0!=VG_IPJM_O&&vg1!=VG_IMJM_O
                   &&vg4!=VG_JMKP_O&&vg5!=VG_JMKM_O){
                    rho(0,1,0)=rho(0,-1,0)+2*h/K;
                }
            } break;
            case VG_KP: {
                if(vg0!=VG_IPKP_O&&vg1!=VG_IMKP_O
                   &&vg2!=VG_JPKP_O&&vg3!=VG_JMKP_O){
                    rho(0,0,-1)=rho(0,0,1)+2*h/K;
                }
            } break;
            case VG_KM: {
                if(vg0!=VG_IPKM_O&&vg1!=VG_IMKM_O
                   &&vg2!=VG_JPKM_O&&vg3!=VG_JMKM_O){
                    rho(0,0,1)=rho(0,0,-1)+2*h/K;
                }
            } break;
            case VG_IPJP_O: {
                if(vg4!=VG_IPJPKP_O&&vg5!=VG_IPJPKM_O){
                    rho(-1,-1,0)=(rho(1,0,0)+rho(0,1,0)+rho(1,-1,0)
                    +rho(-1,1,0))/4+2*h/K;
                }
            } break;
            case VG_IPJM_O: {
                if(vg4!=VG_IPJMKP_O&&vg5!=VG_IPJMKM_O){
                    rho(-1,1,0)=(rho(1,0,0)+rho(0,-1,0)+rho(-1,-1,0)
                    +rho(1,1,0))/4+2*h/K;
                }
            } break;
            case VG_IMJP_O: {
                if(vg4!=VG_IMJPKP_O&&vg5!=VG_IMJPKM_O){
                    rho(1,-1,0)=(rho(-1,0,0)+rho(0,1,0)+rho(-1,-1,0)
                    +rho(1,1,0))/4+2*h/K;
                }
            } break;
            case VG_IMJM_O: {
                if(vg4!=VG_IMJMKP_O&&vg5!=VG_IMJMKM_O){
                    rho(1,1,0)=(rho(-1,0,0)+rho(0,-1,0)+rho(1,-1,0)
                    +rho(-1,1,0))/4+2*h/K;
                }
            } break;
            case VG_IPKP_O: {
                if(vg2!=VG_IPJPKP_O&&vg3!=VG_IPJMKP_O){
                    rho(-1,0,-1)=(rho(1,0,0)+rho(0,0,1)+rho(1,0,-1)
                    +rho(1,0,-1))/4+2*h/K;
                }
            } break;
            case VG_IPKM_O: {
                if(vg2!=VG_IPJPKM_O&&vg3!=VG_IPJMKM_O){
                    rho(-1,0,1)=(rho(1,0,0)+rho(0,0,-1)+rho(-1,0,-1)
                    +rho(1,0,1))/4+2*h/K;
                }
            } break;
            case VG_IMKP_O: {
                if(vg2!=VG_IMJPKP_O&&vg3!=VG_IMJMKP_O){
                    rho(1,0,-1)=(rho(-1,0,0)+rho(0,0,1)+rho(-1,0,-1)
                    +rho(1,0,1))/4+2*h/K;
                }
            } break;
            case VG_IMKM_O: {
                if(vg2!=VG_IMJPKM_O&&vg3!=VG_IMJMKM_O){
                    rho(1,0,1)=(rho(-1,0,0)+rho(0,0,-1)+rho(1,0,-1)
                    +rho(1,0,-1))/4+2*h/K;
                }
            } break;
            case VG_JMKP_O: {
                if(vg0!=VG_IPJMKP_O&&vg1!=VG_IMJMKP_O){
                    rho(0,1,-1)=(rho(0,-1,0)+rho(0,0,1)+rho(0,-1,-1)
                    +rho(0,1,1))/4+2*h/K;
                }
            } break;
            case VG_JMKM_O: {
                if(vg0!=VG_IPJMKM_O&&vg1!=VG_IMJMKM_O){
                    rho(0,1,1)=(rho(0,-1,0)+rho(0,0,-1)+rho(0,-1,1)
                    +rho(0,1,-1))/4+2*h/K;
                }
            } break;
            case VG_JPKP_O: {
                if(vg0!=VG_IPJPKP_O&&vg1!=VG_IMJPKP_O){
                    rho(0,-1,-1)=(rho(0,1,0)+rho(0,0,1)+rho(0,-1,1)
                    +rho(0,1,-1))/4+2*h/K;
                }
            } break;
            case VG_JPKM_O: {
                if(vg0!=VG_IPJPKM_O&&vg1!=VG_IMJPKM_O){
                    rho(0,-1,1)=(rho(0,1,0)+rho(0,0,-1)+rho(0,-1,1)
                    +rho(0,1,-1))/4+2*h/K;
                }
            } break;
            case VG_IPJPKP_O: {
                rho(-1,-1,-1)=(rho(1,0,0)+rho(0,1,0)+rho(0,0,1)+rho(1,-1,-1)
                +rho(-1,1,-1)+rho(-1,-1,1)+rho(1,0,-1)+rho(0,1,-1)+rho(0,-1,1)
                +rho(1,-1,0)+rho(-1,1,0)+rho(-1,0,1))/12+2*h/K;


            } break;
            case VG_IPJMKP_O: {
                rho(-1,1,-1)=(rho(1,0,0)+rho(0,-1,0)+rho(0,0,1)+rho(1,1,-1)
                +rho(-1,-1,-1)+rho(-1,1,1)+rho(1,0,-1)+rho(0,-1,-1)+rho(0,1,1)
                +rho(1,1,0)+rho(-1,-1,0)+rho(-1,0,1))/12+2*h/K;


            } break;
            case VG_IMJPKP_O: {
                rho(1,-1,-1)=(rho(-1,0,0)+rho(0,1,0)+rho(0,0,1)+rho(-1,-1,-1)
                +rho(1,1,-1)+rho(1,-1,1)+rho(-1,0,-1)+rho(0,1,-1)+rho(0,-1,1)
                +rho(-1,-1,0)+rho(1,1,0)+rho(1,0,1))/12+2*h/K;


            } break;
            case VG_IMJMKP_O: {
                rho(1,1,-1)=(rho(-1,0,0)+rho(0,-1,0)+rho(0,0,1)+rho(-1,1,-1)
                +rho(1,-1,-1)+rho(1,1,1)+rho(-1,0,-1)+rho(0,-1,-1)+rho(0,1,1)
                +rho(-1,1,0)+rho(1,-1,0)+rho(1,0,1))/12+2*h/K;

            } break;
            case VG_IPJPKM_O: {
                rho(-1,-1,1)=(rho(1,0,0)+rho(0,1,0)+rho(0,0,-1)+rho(1,-1,1)
                +rho(-1,1,1)+rho(-1,-1,-1)+rho(1,0,1)+rho(0,1,1)+rho(0,-1,-1)
                +rho(1,-1,0)+rho(-1,1,0)+rho(-1,0,-1))/12+2*h/K;


            } break;
            case VG_IPJMKM_O: {
                rho(-1,1,1)=(rho(1,0,0)+rho(0,-1,0)+rho(0,0,-1)+rho(1,1,1)
                +rho(-1,-1,1)+rho(-1,1,-1)+rho(1,0,1)+rho(0,-1,1)+rho(0,1,-1)
                +rho(1,1,0)+rho(-1,-1,0)+rho(-1,0,-1))/12+2*h/K;


            } break;
            case VG_IMJPKM_O: {
                rho(1,-1,1)=(rho(-1,0,0)+rho(0,1,0)+rho(0,0,-1))+rho(-1,-1,1)
                +rho(1,1,1)+rho(1,-1,-1)+rho(-1,0,1)+rho(0,1,1)+rho(0,-1,-1)
                +rho(-1,-1,0)+rho(1,1,0)+rho(1,0,-1)/12+2*h/K;


            } break;
            case VG_IMJMKM_O: {
                rho(1,1,1)=(rho(-1,0,0)+rho(0,-1,0)+rho(0,0,-1)+rho(-1,1,1)
                +rho(1,-1,1)+rho(1,1,-1)+rho(-1,0,1)+rho(0,-1,1)+rho(0,1,-1)
                +rho(-1,1,0)+rho(1,-1,0)+rho(1,0,-1))/12+2*h/K;

            } break;
            case VG_IPJP_I: {
                if(vg2!=VG_IPJP_O){
                    rho(-1,0,0)=rho(1,0,0)+2*h/K;
                }
                if(vg0!=VG_IPJP_O){
                    rho(0,-1,0)=rho(0,1,0)+2*h/K;
                }
            } break;
            case VG_IPJM_I: {
                if(vg3!=VG_IPJM_O){
                    rho(-1,0,0)=rho(1,0,0)+2*h/K;
                }
                if(vg0!=VG_IPJM_O){
                    rho(0,1,0)=rho(0,-1,0)+2*h/K;
                }
            } break;
            case VG_IMJP_I: {
                if(vg2!=VG_IMJP_O){
                    rho(1,0,0)=rho(-1,0,0)+2*h/K;
                }
                if(vg1!=VG_IMJP_O){
                    rho(0,-1,0)=rho(0,1,0)+2*h/K;
                }
            } break;
            case VG_IMJM_I: {
                if(vg3!=VG_IMJM_O){
                    rho(1,0,0)=rho(-1,0,0)+2*h/K;
                }
                if(vg1!=VG_IMJM_O){
                    rho(0,1,0)=rho(0,-1,0)+2*h/K;
                }
            } break;
            case VG_IPKP_I: {
                if(vg4!=VG_IPKP_O){
                    rho(-1,0,0)=rho(1,0,0)+2*h/K;
                }
                if(vg0!=VG_IPKP_O){
                    rho(0,0,-1)=rho(0,0,1)+2*h/K;
                }
            } break;
            case VG_IPKM_I: {
                if(vg5!=VG_IPKM_O){
                    rho(-1,0,0)=rho(1,0,0)+2*h/K;
                }
                if(vg0!=VG_IPKM_O){
                    rho(0,0,1)=rho(0,0,-1)+2*h/K;
                }
            } break;
            case VG_IMKP_I: {
                if(vg4!=VG_IMKP_O){
                    rho(1,0,0)=rho(-1,0,0)+2*h/K;
                }
                if(vg1!=VG_IMKP_O){
                    rho(0,0,-1)=rho(0,0,1)+2*h/K;
                }
            } break;
            case VG_IMKM_I: {
                if(vg5!=VG_IMKM_O){
                    rho(1,0,0)=rho(-1,0,0)+2*h/K;
                }
                if(vg1!=VG_IMKM_O){
                    rho(0,0,1)=rho(0,0,-1)+2*h/K;
                }
            } break;
            case VG_JPKP_I: {
                if(vg4!=VG_JPKP_O){
                    rho(0,-1,0)=rho(0,1,0)+2*h/K;
                }
                if(vg2!=VG_JPKP_O){
                    rho(0,0,-1)=rho(0,0,1)+2*h/K;
                }
            } break;
            case VG_JPKM_I: {
                if(vg5!=VG_JPKM_O){
                    rho(0,-1,0)=rho(0,1,0)+2*h/K;
                }
                if(vg2!=VG_JPKM_O){
                    rho(0,0,1)=rho(0,0,-1)+2*h/K;
                }
            } break;
            case VG_JMKP_I: {
                if(vg4!=VG_JMKP_O){
                    rho(0,1,0)=rho(0,-1,0)+2*h/K;
                }
                if(vg3!=VG_JMKP_O){
                    rho(0,0,-1)=rho(0,0,1)+2*h/K;
                }
            } break;
            case VG_JMKM_I: {
                if(vg5!=VG_JMKM_O){
                    rho(0,1,0)=rho(0,-1,0)+2*h/K;
                }
                if(vg3!=VG_JMKM_O){
                    rho(0,0,1)=rho(0,0,-1)+2*h/K;
                }
            } break;
            case VG_IPJPKP_I: {
                if(vg4!=VG_IPJPKP_O){
                    rho(0,-1,0)=rho(0,1,0)+2*h/K;
                }
                if(vg2!=VG_IPJPKP_O){
                    rho(0,0,-1)=rho(0,0,1)+2*h/K;
                }
                if(vg0!=VG_IPJPKP_O){
                    rho(-1,0,0)=rho(1,0,0)+2*h/K;
                }
            } break;
            case VG_IPJPKM_I: {
                if(vg5!=VG_IPJPKM_O){
                    rho(0,-1,0)=rho(0,1,0)+2*h/K;
                }
                if(vg2!=VG_IPJPKM_O){
                    rho(0,0,1)=rho(0,0,-1)+2*h/K;
                }
                if(vg0!=VG_IPJPKM_O){
                    rho(-1,0,0)=rho(1,0,0)+2*h/K;
                }
            } break;
            case VG_IPJMKP_I: {
                if(vg4!=VG_IPJMKP_O){
                    rho(0,1,0)=rho(0,-1,0)+2*h/K;
                }
                if(vg3!=VG_IPJMKP_O){
                    rho(0,0,-1)=rho(0,0,1)+2*h/K;
                }
                if(vg0!=VG_IPJMKP_O){
                    rho(-1,0,0)=rho(1,0,0)+2*h/K;
                }
            } break;
            case VG_IPJMKM_I: {
                if(vg5!=VG_IPJMKM_O){
                    rho(0,1,0)=rho(0,-1,0)+2*h/K;
                }
                if(vg3!=VG_IPJMKM_O){
                    rho(0,0,1)=rho(0,0,-1)+2*h/K;
                }
                if(vg0!=VG_IPJMKM_O){
                    rho(-1,0,0)=rho(1,0,0)+2*h/K;
                }
            } break;
            case VG_IMJPKP_I: {
                if(vg4!=VG_IMJPKP_O){
                    rho(0,-1,0)=rho(0,1,0)+2*h/K;
                }
                if(vg2!=VG_IMJPKP_O){
                    rho(0,0,-1)=rho(0,0,1)+2*h/K;
                }
                if(vg1!=VG_IMJPKP_O){
                    rho(1,0,0)=rho(-1,0,0)+2*h/K;
                }
            } break;
            case VG_IMJPKM_I: {
                if(vg5!=VG_IMJPKM_O){
                    rho(0,-1,0)=rho(0,1,0)+2*h/K;
                }
                if(vg2!=VG_IMJPKM_O){
                    rho(0,0,1)=rho(0,0,-1)+2*h/K;
                }
                if(vg1!=VG_IMJPKM_O){
                    rho(1,0,0)=rho(-1,0,0)+2*h/K;
                }
            } break;
            case VG_IMJMKP_I: {
                if(vg4!=VG_IMJMKP_O){
                    rho(0,1,0)=rho(0,-1,0)+2*h/K;
                }
                if(vg3!=VG_IMJMKP_O){
                    rho(0,0,-1)=rho(0,0,1)+2*h/K;
                }
                if(vg1!=VG_IMJMKP_O){
                    rho(1,0,0)=rho(-1,0,0)+2*h/K;
                }
            } break;
            case VG_IMJMKM_I: {
                if(vg5!=VG_IMJMKM_O){
                    rho(0,1,0)=rho(0,-1,0)+2*h/K;
                }
                if(vg3!=VG_IMJMKM_O){
                    rho(0,0,1)=rho(0,0,-1)+2*h/K;
                }
                if(vg1!=VG_IMJMKM_O){
                    rho(1,0,0)=rho(-1,0,0)+2*h/K;
                }
            } break;
            case VG_IPJP_O_KP: {
                if(vg0!=VG_IPKP_O&&vg2!=VG_JPKP_O){
                    rho(0,0,-1)=rho(0,0,1)+2*h/K;

                }
                if(vg4!=VG_IPJPKP_O){
                    rho(-1,-1,0)=(rho(1,0,0)+rho(0,1,0)+rho(1,-1,0)+rho(-1,1,0))/4
                    +2*h/K;

                }
            } break;
            case VG_IPJM_O_KP: {
                if(vg0!=VG_IPKP_O&&vg3!=VG_JMKP_O){
                    rho(0,0,-1)=rho(0,0,1)+2*h/K;

                }
                if(vg4!=VG_IPJMKP_O){
                    rho(-1,1,0)=(rho(1,0,0)+rho(0,-1,0)+rho(-1,-1,0)+rho(1,1,0))/4
                    +2*h/K;

                }
            } break;
            case VG_IMJP_O_KP: {
                if(vg1!=VG_IMKP_O&&vg2!=VG_JPKP_O){
                    rho(0,0,-1)=rho(0,0,1)+2*h/K;

                }
                if(vg4!=VG_IMJPKP_O){
                    rho(1,-1,0)=(rho(-1,0,0)+rho(0,1,0)+rho(-1,-1,0)+rho(1,1,0))/4
                    +2*h/K;

                }
            } break;
            case VG_IMJM_O_KP: {
                if(vg1!=VG_IMKP_O&&vg3!=VG_JMKP_O){
                    rho(0,0,-1)=rho(0,0,1)+2*h/K;

                }
                if(vg4!=VG_IMJMKP_O){
                    rho(1,1,0)=(rho(-1,0,0)+rho(0,-1,0)+rho(1,-1,0)+rho(-1,1,0))/4
                    +2*h/K;
                }
            } break;
            case VG_IPJP_O_KM: {
                if(vg0!=VG_IPKM_O&&vg2!=VG_JPKM_O){
                    rho(0,0,1)=rho(0,0,-1)+2*h/K;

                }
                if(vg5!=VG_IPJPKM_O){
                    rho(-1,-1,0)=(rho(1,0,0)+rho(0,1,0)+rho(1,-1,0)+rho(-1,1,0))/4
                    +2*h/K;

                }
            } break;
            case VG_IPJM_O_KM: {
                if(vg0!=VG_IPKM_O&&vg3!=VG_JMKM_O){
                    rho(0,0,1)=rho(0,0,-1)+2*h/K;

                }
                if(vg5!=VG_IPJMKM_O){
                    rho(-1,1,0)=(rho(1,0,0)+rho(0,-1,0)+rho(-1,-1,0)+rho(1,1,0))/4
                    +2*h/K;

                }
            } break;
            case VG_IMJP_O_KM: {
                if(vg1!=VG_IMKP_O&&vg2!=VG_JPKM_O){
                    rho(0,0,1)=rho(0,0,-1)+2*h/K;

                }
                if(vg5!=VG_IMJPKM_O){
                    rho(1,-1,0)=(rho(-1,0,0)+rho(0,1,0)+rho(-1,-1,0)+rho(1,1,0))/4
                    +2*h/K;

                }
            } break;
            case VG_IMJM_O_KM: {
                if(vg1!=VG_IMKM_O&&vg3!=VG_JMKM_O){
                    rho(0,0,1)=rho(0,0,-1)+2*h/K;

                }
                if(vg5!=VG_IMJMKM_O){
                    rho(1,1,0)=(rho(-1,0,0)+rho(0,-1,0)+rho(1,-1,0)+rho(-1,1,0))/4
                    +2*h/K;
                }
            } break;
            case VG_IPKP_O_JP: {
                if(vg0!=VG_IPJP_O&&vg4!=VG_JPKP_O){
                    rho(0,-1,0)=rho(0,1,0)+2*h/K;

                }
                if(vg2!=VG_IPJPKP_O){
                    rho(-1,0,-1)=(rho(1,0,0)+rho(0,0,1)+rho(1,0,-1)+rho(1,0,-1))/4
                    +2*h/K;
                }
            } break;
            case VG_IPKM_O_JP: {
                if(vg0!=VG_IPJP_O&&vg5!=VG_JPKM_O){
                    rho(0,-1,0)=rho(0,1,0)+2*h/K;

                }
                if(vg2!=VG_IPJPKM_O){
                    rho(-1,0,1)=(rho(1,0,0)+rho(0,0,-1)+rho(-1,0,-1)+rho(1,0,1))/4
                    +2*h/K;
                }
            } break;
            case VG_IMKP_O_JP: {
                if(vg1!=VG_IMJP_O&&vg4!=VG_JPKP_O){
                    rho(0,-1,0)=rho(0,1,0)+2*h/K;

                }
                if(vg2!=VG_IMJPKM_O){
                    rho(1,0,-1)=(rho(-1,0,0)+rho(0,0,1)+rho(-1,0,-1)+rho(1,0,1))/4
                    +2*h/K;
                }
            } break;
            case VG_IMKM_O_JP: {
                if(vg1!=VG_IMJP_O&&vg5!=VG_JPKM_O){
                    rho(0,-1,0)=rho(0,1,0)+2*h/K;

                }
                if(vg2!=VG_IMJPKM_O){
                    rho(1,0,1)=(rho(-1,0,0)+rho(0,0,-1)+rho(1,0,-1)+rho(1,0,-1))/4
                    +2*h/K;
                }
            } break;
            case VG_IPKP_O_JM: {
                if(vg0!=VG_IPJM_O&&vg4!=VG_JMKP_O){
                    rho(0,1,0)=rho(0,-1,0)+2*h/K;

                }
                if(vg3!=VG_IPJMKP_O){
                    rho(-1,0,-1)=(rho(1,0,0)+rho(0,0,1)+rho(1,0,-1)+rho(1,0,-1))/4
                    +2*h/K;
                }
            } break;
            case VG_IPKM_O_JM: {
                if(vg0!=VG_IPJM_O&&vg5!=VG_JMKM_O){
                    rho(0,1,0)=rho(0,-1,0)+2*h/K;

                }
                if(vg3!=VG_IPJMKM_O){
                    rho(-1,0,1)=(rho(1,0,0)+rho(0,0,-1)+rho(-1,0,-1)+rho(1,0,1))/4
                    +2*h/K;
                }
            } break;
            case VG_IMKP_O_JM: {
                if(vg1!=VG_IMJM_O&&vg4!=VG_JMKP_O){
                    rho(0,1,0)=rho(0,-1,0)+2*h/K;

                }
                if(vg3!=VG_IMJMKP_O){
                    rho(1,0,-1)=(rho(-1,0,0)+rho(0,0,1)+rho(-1,0,-1)+rho(1,0,1))/4
                    +2*h/K;
                }
            } break;
            case VG_IMKM_O_JM: {
                if(vg1!=VG_IMJM_O&&vg5!=VG_JMKM_O){
                    rho(0,1,0)=rho(0,-1,0)+2*h/K;

                }
                if(vg3!=VG_IMJMKM_O){
                    rho(1,0,1)=(rho(-1,0,0)+rho(0,0,-1)+rho(1,0,-1)+rho(1,0,-1))/4
                    +2*h/K;
                }
            } break;
            case VG_JPKP_O_IP: {
                if(vg2!=VG_IPJP_O&&vg4!=VG_IPKP_O){
                    rho(-1,0,0)=rho(1,0,0)+2*h/K;
                }
                if(vg0!=VG_IPJPKP_O){
                    rho(0,-1,-1)=(rho(0,1,0)+rho(0,0,1)+rho(0,-1,1)+rho(0,1,-1))/4
                    +2*h/K;
                }
            } break;
            case VG_JPKM_O_IP: {
                if(vg2!=VG_IPJP_O&&vg5!=VG_IPKM_O){
                    rho(-1,0,0)=rho(1,0,0)+2*h/K;
                }
                if(vg0!=VG_IPJPKM_O){
                    rho(0,-1,1)=(rho(0,1,0)+rho(0,0,-1)+rho(0,-1,1)+rho(0,1,-1))/4
                    +2*h/K;
                }
            } break;
            case VG_JMKP_O_IP: {
                if(vg3!=VG_IPJM_O&&vg4!=VG_IPKP_O){
                    rho(-1,0,0)=rho(1,0,0)+2*h/K;
                }
                if(vg0!=VG_IPJMKP_O){
                    rho(0,1,-1)=(rho(0,-1,0)+rho(0,0,1)+rho(0,-1,-1)+rho(0,1,1))/4
                    +2*h/K;
                }
            } break;
            case VG_JMKM_O_IP: {
                if(vg3!=VG_IPJM_O&&vg5!=VG_IPKM_O){
                    rho(-1,0,0)=rho(1,0,0)+2*h/K;
                }
                if(vg0!=VG_IPJMKM_O){
                    rho(0,1,1)=(rho(0,-1,0)+rho(0,0,-1)+rho(0,-1,1)+rho(0,1,-1))/4
                    +2*h/K;
                }
            } break;
            case VG_JPKP_O_IM: {
                if(vg2!=VG_IMJP_O&&vg4!=VG_IMKP_O){
                    rho(1,0,0)=rho(-1,0,0)+2*h/K;
                }
                if(vg1!=VG_IMJPKP_O){
                    rho(0,-1,-1)=(rho(0,1,0)+rho(0,0,1)+rho(0,-1,1)+rho(0,1,-1))/4
                    +2*h/K;
                }
            } break;
            case VG_JPKM_O_IM: {
                if(vg2!=VG_IMJP_O&&vg5!=VG_IMKM_O){
                    rho(1,0,0)=rho(-1,0,0)+2*h/K;
                }
                if(vg1!=VG_IMJPKM_O){
                    rho(0,-1,1)=(rho(0,1,0)+rho(0,0,-1)+rho(0,-1,1)+rho(0,1,-1))/4
                    +2*h/K;
                }
            } break;
            case VG_JMKP_O_IM: {
                if(vg3!=VG_IMJM_O&&vg4!=VG_IMKP_O){
                    rho(1,0,0)=rho(-1,0,0)+2*h/K;
                }
                if(vg1!=VG_IMJMKP_O){
                    rho(0,1,-1)=(rho(0,-1,0)+rho(0,0,1)+rho(0,-1,-1)+rho(0,1,1))/4
                    +2*h/K;
                }
            } break;
            case VG_JMKM_O_IM: {
                if(vg3!=VG_IMJM_O&&vg5!=VG_IMKM_O){
                    rho(1,0,0)=rho(-1,0,0)+2*h/K;
                }
                if(vg1!=VG_IMJMKM_O){
                    rho(0,1,1)=(rho(0,-1,0)+rho(0,0,-1)+rho(0,-1,1)+rho(0,1,-1))/4
                    +2*h/K;
                }
            } break;
        }
    }
}

//Calculate chemical potential within Free Energy Model
//TODO:This kernel name is repeated in the AdDiff App
void KerCalcMu3D(ACC<Real>&mu,const ACC<Real>&Phi,const ACC<Real>&d2phi,
                 const ACC<int>&nodeType, const Real *Aconst, const Real *Kappa){
    Real phi{Phi(0,0,0)};
    Real phi3 = {phi*phi*phi};
    const Real A={*Aconst};
    const Real kappa={*Kappa};
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt != VertexType::ImmersedSolid){
    mu(0,0,0)=-A*phi+A*phi3-kappa*(d2phi(0,0,0,0)+d2phi(1,0,0,0)+d2phi(2,0,0,0));
    //if (mu(0,0,0)>=10){
    //std::cout<<A<<" "<<kappa<<" "<<phi<<" "<<mu(0,0,0)<<"    ";
    //}
    }
}



//Update surface tension force between fluid components
void KerUpdateMacroBodyForce3D(ACC<Real>& force, const ACC<Real>& Phi,
                               const ACC<Real>& dmu, const ACC<int>& geometry,
                               const ACC<int>&nodeType,
                               const ACC<Real>& coordinates, const int* idx) {
    Real phi = {Phi(0,0,0)};
    Real dmux = {dmu(0,0,0,0)};
    Real dmuy = {dmu(1,0,0,0)};
    Real dmuz = {dmu(2,0,0,0)};

    VertexType vt = (VertexType)nodeType(0, 0, 0);
    VertexGeometryType vg = (VertexGeometryType)geometry(0, 0, 0);
    //Set gradient of chemical potential perpendicular to wall = 0
    if (vt == VertexType::Wall) {

        switch (vg) {
            case VG_IP: {
                dmux=0;
            } break;
            case VG_IM: {
                dmux=0;
            } break;
            case VG_JP: {
                dmuy=0;
            } break;
            case VG_JM: {
                dmuy=0;
            } break;
            case VG_KP: {
                dmuz=0;
            } break;
            case VG_KM: {
                dmuz=0;
            } break;
            case VG_IPJP_O: {
                dmux=0;
                dmuy=0;

            } break;
            case VG_IPJM_O: {
                dmux=0;
                dmuy=0;

            } break;
            case VG_IMJP_O: {
                dmux=0;
                dmuy=0;

            } break;
            case VG_IMJM_O: {
                dmux=0;
                dmuy=0;
            } break;
            case VG_IPKP_O: {
                dmux=0;
                dmuz=0;

            } break;
            case VG_IPKM_O: {
                dmux=0;
                dmuz=0;

            } break;
            case VG_IMKP_O: {
                dmux=0;
                dmuz=0;

            } break;
            case VG_IMKM_O: {
                dmux=0;
                dmuz=0;

            } break;
            case VG_JMKP_O: {
                dmuy=0;
                dmuz=0;

            } break;
            case VG_JMKM_O: {
                dmuy=0;
                dmuz=0;
            } break;
            case VG_JPKP_O: {
                dmuy=0;
                dmuz=0;

            } break;
            case VG_JPKM_O: {
                dmuy=0;
                dmuz=0;
            } break;


            case VG_IPJP_I: {
                dmux=0;
                dmuy=0;

            } break;
            case VG_IPJM_I: {
                dmux=0;
                dmuy=0;

            } break;
            case VG_IMJP_I: {
                dmux=0;
                dmuy=0;

            } break;
            case VG_IMJM_I: {
                dmux=0;
                dmuy=0;
            } break;
            case VG_IPKP_I: {
                dmux=0;
                dmuz=0;

            } break;
            case VG_IPKM_I: {
                dmux=0;
                dmuz=0;

            } break;
            case VG_IMKP_I: {
                dmux=0;
                dmuz=0;

            } break;
            case VG_IMKM_I: {
                dmux=0;
                dmuz=0;

            } break;
            case VG_JMKP_I: {
                dmuy=0;
                dmuz=0;

            } break;
            case VG_JMKM_I: {
                dmuy=0;
                dmuz=0;
            } break;
            case VG_JPKP_I: {
                dmuy=0;
                dmuz=0;

            } break;
            case VG_JPKM_I: {
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IPJPKP_O: {
                dmux=0;
                dmuy=0;
                dmuz=0;

            } break;
            case VG_IPJMKP_O: {
                dmux=0;
                dmuy=0;
                dmuz=0;

            } break;
            case VG_IMJPKP_O: {
                dmux=0;
                dmuy=0;
                dmuz=0;

            } break;
            case VG_IMJMKP_O: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IPJPKM_O: {
                dmux=0;
                dmuy=0;
                dmuz=0;

            } break;
            case VG_IPJMKM_O: {
                dmux=0;
                dmuy=0;
                dmuz=0;

            } break;
            case VG_IMJPKM_O: {
                dmux=0;
                dmuy=0;
                dmuz=0;

            } break;
            case VG_IMJMKM_O: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;

            case VG_IPJPKP_I: {
                dmux=0;
                dmuy=0;
                dmuz=0;

            } break;
            case VG_IPJMKP_I: {
                dmux=0;
                dmuy=0;
                dmuz=0;

            } break;
            case VG_IMJPKP_I: {
                dmux=0;
                dmuy=0;
                dmuz=0;

            } break;
            case VG_IMJMKP_I: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IPJPKM_I: {
                dmux=0;
                dmuy=0;
                dmuz=0;

            } break;
            case VG_IPJMKM_I: {
                dmux=0;
                dmuy=0;
                dmuz=0;

            } break;
            case VG_IMJPKM_I: {
                dmux=0;
                dmuy=0;
                dmuz=0;

            } break;
            case VG_IMJMKM_I: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IPJP_I_KP: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IPJM_I_KP: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IMJP_I_KP: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IMJM_I_KP: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IPJP_I_KM: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IPJM_I_KM: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IMJP_I_KM: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IMJM_I_KM: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IPKP_I_JP: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IPKM_I_JP: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IMKP_I_JP: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IMKM_I_JP: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IPKP_I_JM: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IPKM_I_JM: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IMKP_I_JM: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IMKM_I_JM: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_JPKP_I_IP: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_JPKM_I_IP: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_JMKP_I_IP: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_JMKM_I_IP: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_JPKP_I_IM: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_JPKM_I_IM: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_JMKP_I_IM: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_JMKM_I_IM: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IPJP_O_KP: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IPJM_O_KP: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IMJP_O_KP: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IMJM_O_KP: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IPJP_O_KM: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IPJM_O_KM: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IMJP_O_KM: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IMJM_O_KM: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IPKP_O_JP: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IPKM_O_JP: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IMKP_O_JP: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IMKM_O_JP: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IPKP_O_JM: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IPKM_O_JM: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IMKP_O_JM: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_IMKM_O_JM: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_JPKP_O_IP: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_JPKM_O_IP: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_JMKP_O_IP: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_JMKM_O_IP: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_JPKP_O_IM: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_JPKM_O_IM: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_JMKP_O_IM: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;
            case VG_JMKM_O_IM: {
                dmux=0;
                dmuy=0;
                dmuz=0;
            } break;

        }
    }
    if (vt != VertexType::ImmersedSolid){
        force(0,0,0,0) = -phi*dmux;
        force(1,0,0,0) = -phi*dmuy;
        force(2,0,0,0) = -phi*dmuz;
    }

}


#endif //OPS_3D outter

#endif //MODEL_KERNEL_INC