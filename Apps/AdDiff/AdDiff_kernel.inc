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

#ifndef CAVITY2D_KERNEL_INC
#define CAVITY2D_KERNEL_INC
#ifdef OPS_2D

void KerSetInitialMacroVars(ACC<Real>& rho, ACC<Real>& u, ACC<Real>& v,
                            const ACC<Real>& coordinates, const int* idx) {
    rho(0, 0) = 1;
    u(0, 0) = 0;
    v(0, 0) = 0;
}

//ADCase

//void KerCalcConcentration(ACC<Real>& concentration, const ACC<Real>& xiIndex, const Real* TIME, const ACC<Real>& U, const ACC<Real>& V, const Real* Dt) {
//    Real Sig02{16};
//    Real DiffCFact{1.5};
//    Real dt = {*Dt};
//    Real time = {*TIME};
//    Real DiffC{DiffCFact/dt};
//    Real x0[2] = {8.0,8.0};
//    Real Vel[2] = {U(0,0),V(0,0)};
//    Real Pref{Sig02/(Sig02)};
//    Real xiIndexVar1 = {xiIndex(0,0,0)};
//    Real xiIndexVar2 = {xiIndex(1,0,0)};
//    Real numeru{(xiIndexVar1-x0[0])};
//    Real numerv{(xiIndexVar2-x0[1])};
//    Real numer2{numeru*numeru+numerv*numerv};
//    Real Gauss{exp(-numer2/(2*(Sig02)))};
//
//    concentration(0,0) = Pref * Gauss;
//    //std::cout << concentration(0,0) << " " << xiIndex(0,0,0) << " " << U(0,0) << " " << V(0,0) << " " << numer2 << "\n";
//}

void KerCalcConcentration(ACC<Real>& concentration, const ACC<Real>& xiIndex, const Real* TIME, const ACC<Real>& U, const ACC<Real>& V, const Real* Dt) {
    Real x0[2] = {0.0,60.0};
    Real xiIndexVar1 = {xiIndex(0,0,0)};
    Real xiIndexVar2 = {xiIndex(1,0,0)};
    Real x=(xiIndexVar1-x0[0]);
    Real y=(xiIndexVar2-x0[1]);
    Real R=sqrt(x*x+y*y);
    Real r=30.0;
    
    if (R<=r){
        concentration(0,0) = 1;
    }
    else{
        concentration(0,0) = -1;
    }
    
   //if (x>30 && x<90){
   //    concentration(0,0) = 1;
   //}
   //else{
   //    concentration(0,0) = -1;
   //}

}



void KerCalcGradients(ACC<Real>&grad, const ACC<Real>&f, const ACC<int>&nodeType, const ACC<int>&geometry, const ACC<Real>& coordinates, const int *lattIdx, const int *Order, const int *idx){
    VertexType vt = (VertexType)nodeType(0, 0);
    int order={*Order};
        
        if (order == 1){
            if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic ||
                vt == VertexType::VirtualBoundary || vt == VertexType::Wall) {
                grad(0,0,0) = 0.5*(f(1, 0)-f(-1, 0));
                grad(1,0,0) = 0.5*(f(0, 1)-f(0, -1));
                //std::cout<<grad(0,0,0)<<" "<<grad(1,0,0)<<" ";
            }    
        }  
        else if (order == 2){     
            if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic ||
                vt == VertexType::VirtualBoundary || vt == VertexType::Wall) {
                grad(0,0,0) = f(-1, 0)-2*f(0, 0)+f(1, 0);
                grad(1,0,0) = f(0, -1)-2*f(0, 0)+f(0, 1);
                
            }
        }
        else{
            std::cout << order << "\n";
        }
    
}

void KerUpdateRhoWetting(ACC<Real>& rho,const ACC<int>& geometry, const ACC<int>&nodeType,const ACC<Real>& coordinates) {
    const Real teta={M_PI/2};
    const Real alpha={acos(sin(teta)*sin(teta))};
    const Real Omga={2*((teta <= M_PI/2) - (M_PI/2 < teta))*sqrt(cos(alpha/3)*(1-cos(alpha/3)))};
    const Real K={0.02};
    const Real A=0.02;
    const Real h={Omga*sqrt(K*A/2)};
    VertexType vt = (VertexType)nodeType(0, 0);
    VertexGeometryType vg = (VertexGeometryType)geometry(0, 0);
    if (vt == VertexType::Wall) {
        switch (vg) {
            case VG_IP: {
                rho(-1,0)=rho(1,0)+2*h/K;
            } break;
            case VG_IM: {
                rho(1,0)=rho(-1,0)+2*h/K;
            } break;
            case VG_JP: {
                rho(0,-1)=rho(0,1)+2*h/K;
            } break;
            case VG_JM: {
                rho(0,1)=rho(0,-1)+2*h/K;
            } break;
        }
    }
}

void KerCalcMu(ACC<Real>&mu,const ACC<Real>&Phi,const ACC<Real>&d2phi, const ACC<int>& geometry, const ACC<int>&nodeType,const ACC<Real>& coordinates){
    Real phi{Phi(0,0)};
    Real phi3 = {phi*phi*phi};
    Real d2phix = {d2phi(0,0,0)};
    Real d2phiy = {d2phi(1,0,0)};
    Real A=0.02;
    Real kappa=0.02;
    
    mu(0,0)=-A*phi+A*phi3-kappa*(d2phix+d2phiy);
}

void KerPrintPhi(const ACC<Real>& Phi){
    std::cout << Phi(0,0) << " ";
}

void KerUpdateMacroBodyForce(ACC<Real>& force, const ACC<Real>& Phi, const ACC<Real>& dmu, const ACC<int>& geometry, const ACC<int>&nodeType, const ACC<Real>& coordinates,
                            const int* idx) {
    Real phi = {Phi(0,0)};
    Real dmux = {dmu(0,0,0)};
    Real dmuy = {dmu(1,0,0)};
    VertexType vt = (VertexType)nodeType(0, 0);
    VertexGeometryType vg = (VertexGeometryType)geometry(0, 0);
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
        }
    }
    force(0,0,0) = -phi*dmux;
    force(1,0,0) = -phi*dmuy;
}

void KerSetSolid(ACC<int>& nodeType, ACC<int>& geometry, const ACC<Real>& coordinates, const int* xMax, const int* yMax) {
    const Real X{coordinates(0,0,0)};
    const Real Y{coordinates(1,0,0)};
    int XMax{*xMax};
    int YMax{*yMax};
    std::cout<<XMax<<" ";
    if (X>15||X<XMax-15){
        nodeType(0,0)=(int)VertexType::ImmersedSolid;
        geometry(0,0)=(int)VG_ImmersedSolid;
    }
}

void KerSetEmbeddedBodyGeometry(ACC<int>& geometryProperty,
                                ACC<int>& nodeType) {
#ifdef OPS_2D
    VertexType vt = (VertexType)nodeType(0, 0);
    if (VertexType::ImmersedSolid == vt) {
        VertexType neiborVertexType[8];
        /*
                        6*****2*****4
                        *     *     *
                        *     *     *
                        1***********0
                        *     *     *
                        *     *     *
                        5*****3*****7
        */
        neiborVertexType[0] = (VertexType)nodeType(1, 0);
        neiborVertexType[1] = (VertexType)nodeType(-1, 0);
        neiborVertexType[2] = (VertexType)nodeType(0, 1);
        neiborVertexType[3] = (VertexType)nodeType(0, -1);
        neiborVertexType[4] = (VertexType)nodeType(1, 1);
        neiborVertexType[5] = (VertexType)nodeType(-1, -1);
        neiborVertexType[6] = (VertexType)nodeType(-1, 1);
        neiborVertexType[7] = (VertexType)nodeType(1, -1);
        int fluidNeiborNum{0};
        for (int i = 0; i < 8; i++) {
            if (VertexType::ImmersedSolid != neiborVertexType[i]) {
                fluidNeiborNum++;
            }
        }
        int solidNeiborNumatCoord{0};
        for (int i = 0; i < 4; i++) {
            if (VertexType::ImmersedSolid == neiborVertexType[i]) {
                solidNeiborNumatCoord++;
            }
        }
        if (fluidNeiborNum > 0) {
            nodeType(0,0)=(int)VertexType::Wall;
            // outer corner
            if (2 == solidNeiborNumatCoord) {
                if ((VertexType::ImmersedSolid == neiborVertexType[0] &&
                     VertexType::ImmersedSolid == neiborVertexType[1]) ||
                    (VertexType::ImmersedSolid == neiborVertexType[0] &&
                     VertexType::ImmersedSolid == neiborVertexType[1])) {
                    ops_printf("%s\n",
                               "There appears to be hanged solid points,i.e., "
                               "the solid body may be too thin");
                }
                if (VertexType::ImmersedSolid == neiborVertexType[2] &&
                    VertexType::ImmersedSolid == neiborVertexType[1]) {
                    if (VertexType::ImmersedSolid == neiborVertexType[6]) {
                        geometryProperty(0, 0) = (int)VG_IPJM_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points");
                    }
                }

                if (VertexType::ImmersedSolid == neiborVertexType[3] &&
                    VertexType::ImmersedSolid == neiborVertexType[1]) {
                    if (VertexType::ImmersedSolid == neiborVertexType[5]) {
                        geometryProperty(0, 0) = (int)VG_IPJP_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points");
                    }
                }

                if (VertexType::ImmersedSolid == neiborVertexType[3] &&
                    VertexType::ImmersedSolid == neiborVertexType[0]) {
                    if (VertexType::ImmersedSolid == neiborVertexType[7]) {
                        geometryProperty(0, 0) = (int)VG_IMJP_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points");
                    }
                }
                if (VertexType::ImmersedSolid == neiborVertexType[2] &&
                    VertexType::ImmersedSolid == neiborVertexType[0]) {
                    if (VertexType::ImmersedSolid == neiborVertexType[4]) {
                        geometryProperty(0, 0) = (int)VG_IMJM_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points");
                    }
                }
            }
            // Planlar corner
            if (3 == solidNeiborNumatCoord) {
                if (VertexType::ImmersedSolid != neiborVertexType[0]) {
                    geometryProperty(0, 0) = (int)VG_IP;
                }
                if (VertexType::ImmersedSolid != neiborVertexType[1]) {
                    geometryProperty(0, 0) = (int)VG_IM;
                }
                if (VertexType::ImmersedSolid != neiborVertexType[2]) {
                    geometryProperty(0, 0) = (int)VG_JP;
                }
                if (VertexType::ImmersedSolid != neiborVertexType[3]) {
                    geometryProperty(0, 0) = (int)VG_JM;
                }
            }
            // Inner corner
            if (4 == solidNeiborNumatCoord) {
                if (1 == fluidNeiborNum) {
                    if (VertexType::ImmersedSolid != neiborVertexType[4]) {
                        geometryProperty(0, 0) = (int)VG_IPJP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[5]) {
                        geometryProperty(0, 0) = (int)VG_IMJM_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[6]) {
                        geometryProperty(0, 0) = (int)VG_IMJP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[7]) {
                        geometryProperty(0, 0) = (int)VG_IPJM_I;
                    }
                } else {
                    ops_printf("%s\n",
                               "There appears to be hanged solid points");
                }
            }
        }
    }
#endif
}
#endif
#endif // CAVITY2D_KERNEL_INC

#ifndef CAVITY3D_KERNEL_INC
#define CAVITY3D_KERNEL_INC
#ifdef OPS_3D
void KerSetInitialMacroVars3D(ACC<Real>& rho, ACC<Real>& u, ACC<Real>& v, ACC<Real>& w,
                            const ACC<Real>& coordinates, const int* idx) {
    rho(0, 0, 0) = 1;
    u(0, 0, 0) = 0;
    v(0, 0, 0) = 0;
    w(0, 0, 0) = 0;
}

//ADCase

//void KerCalcConcentration(ACC<Real>& concentration, const ACC<Real>& xiIndex, const Real* TIME, const ACC<Real>& U, const ACC<Real>& V, const Real* Dt) {
//    Real Sig02{16};
//    Real DiffCFact{1.5};
//    Real dt = {*Dt};
//    Real time = {*TIME};
//    Real DiffC{DiffCFact/dt};
//    Real x0[2] = {8.0,8.0};
//    Real Vel[2] = {U(0,0),V(0,0)};
//    Real Pref{Sig02/(Sig02)};
//    Real xiIndexVar1 = {xiIndex(0,0,0)};
//    Real xiIndexVar2 = {xiIndex(1,0,0)};
//    Real numeru{(xiIndexVar1-x0[0])};
//    Real numerv{(xiIndexVar2-x0[1])};
//    Real numer2{numeru*numeru+numerv*numerv};
//    Real Gauss{exp(-numer2/(2*(Sig02)))};
//
//    concentration(0,0) = Pref * Gauss;
//    //std::cout << concentration(0,0) << " " << xiIndex(0,0,0) << " " << U(0,0) << " " << V(0,0) << " " << numer2 << "\n";
//}

void KerCalcConcentration3D(ACC<Real>& concentration, const ACC<Real>& xiIndex, const Real* TIME, const ACC<Real>& U, const ACC<Real>& V, const Real* Dt) {
    Real x0[3] = {20.0,20.0,0.0};
    Real xiIndexVar1 = {xiIndex(0,0,0,0)};
    Real xiIndexVar2 = {xiIndex(1,0,0,0)};
    Real xiIndexVar3 = {xiIndex(2,0,0,0)};
    Real x=(xiIndexVar1-x0[0]);
    Real y=(xiIndexVar2-x0[1]);
    Real z=(xiIndexVar3-x0[2]);
    Real R=sqrt(x*x+y*y+z*z);
    Real r=15.0;
    if (R<=r){
        concentration(0,0,0) = 1;
    }
    else{
        concentration(0,0,0) = -1;
    }
    //std::cout << concentration(0,0) << "\n";
    
    //std::cout << concentration(0,0) << " " << xiIndex(0,0,0) << " " << U(0,0) << " " << V(0,0) << " " << numer2 << "\n";
}

void KerCalcGradients3D(ACC<Real>&grad, const ACC<Real>&f, const ACC<int>&nodeType, const ACC<int>&geometry, const ACC<Real>& coordinates, const int *lattIdx, const int *Order, const int *idx){
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    int order={*Order};
        
        if (order == 1){
            if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic ||
                vt == VertexType::VirtualBoundary || vt==VertexType::Wall) {
                grad(0,0,0,0) = 0.5*(f(1, 0, 0)-f(-1, 0, 0));
                grad(1,0,0,0) = 0.5*(f(0, 1, 0)-f(0, -1, 0));
                grad(2,0,0,0) = 0.5*(f(0, 0, 1)-f(0, 0, -1));
            }    
        }  
        else if (order == 2){     
            if (vt == VertexType::Fluid || vt == VertexType::MDPeriodic ||
                vt == VertexType::VirtualBoundary || vt==VertexType::Wall) {
                grad(0,0,0,0) = f(-1, 0, 0)-2*f(0, 0, 0)+f(1, 0, 0);
                grad(1,0,0,0) = f(0, -1, 0)-2*f(0, 0, 0)+f(0, 1, 0);
                grad(2,0,0,0) = f(0, 0, -1)-2*f(0, 0, 0)+f(0, 0, 1);
            }
        }
        else{
            std::cout << order << "\n";
        }
    
}

void KerUpdateRhoWetting3D(ACC<Real>& rho,const ACC<int>& geometry, const ACC<int>&nodeType, const ACC<Real>& coordinates) {
    const Real teta={M_PI/4};
    const Real alpha={acos(sin(teta)*sin(teta))};
    const Real Omga={2*((teta <= M_PI/2) - (M_PI/2 < teta))*sqrt(cos(alpha/3)*(1-cos(alpha/3)))};
    const Real K={0.01};
    const Real A=0.01;
    const Real h={Omga*sqrt(K*A/2)};
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    VertexGeometryType vg = (VertexGeometryType)geometry(0, 0, 0);
    if (vt == VertexType::Wall) {

        switch (vg) {
            case VG_IP: {
                rho(-1,0,0)=rho(1,0,0)+2*h/K;
            } break;
            case VG_IM: {
                rho(1,0,0)=rho(-1,0,0)+2*h/K;
            } break;
            case VG_JP: {
                rho(0,-1,0)=rho(0,1,0)+2*h/K;
            } break;
            case VG_JM: {
                rho(0,1,0)=rho(0,-1,0)+2*h/K;
            } break;
            case VG_KP: {
                rho(0,0,-1)=rho(0,0,1)+2*h/K;
            } break;
            case VG_KM: {
                rho(0,0,1)=rho(0,0,-1)+2*h/K;
            } break;
        }
    }
}

void KerCalcMu3D(ACC<Real>&mu,const ACC<Real>&Phi,const ACC<Real>&d2phi){
    Real phi{Phi(0,0,0)};
    Real phi3 = {phi*phi*phi};
    Real A=0.01;
    Real kappa=0.01;
    mu(0,0,0)=-A*phi+A*phi3-kappa*(d2phi(0,0,0,0)+d2phi(1,0,0,0)+d2phi(2,0,0,0));
}

void KerPrintPhi3D(const ACC<Real>& Phi){
    std::cout << Phi(0,0,0) << " ";
}

void KerUpdateMacroBodyForce3D(ACC<Real>& force, const ACC<Real>& Phi, const ACC<Real>& dmu, const ACC<int>& geometry, const ACC<int>&nodeType, const ACC<Real>& coordinates,
                            const int* idx) {
    Real phi = {Phi(0,0,0)};
    Real dmux = {dmu(0,0,0,0)};
    Real dmuy = {dmu(1,0,0,0)};
    Real dmuz = {dmu(2,0,0,0)};
    //std::cout << "Stuff: " << dphix << " " << dphiy << " " << dmux << " " << dmuy << " " << phi << " " << coordinates(0,0,0) << " " << coordinates(1,0,0) << "\n";
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    VertexGeometryType vg = (VertexGeometryType)geometry(0, 0, 0);
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
        }
    }
    
    force(0,0,0,0) = -phi*dmux;
    force(1,0,0,0) = -phi*dmuy;
    force(2,0,0,0) = -phi*dmuz;
    //std::cout<<"Force: "<<force(0,0,0)<<" "<<force(1,0,0)<<"\n";
}
#endif
#endif // CAVITY2D_KERNEL_INC