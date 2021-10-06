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

/**
 * The Following contains kernel functions for the Free Energy model for
 * multiphase/multicomponent flow in 2D and 3D.
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

//Update order parameter (or concentration in advection/diffusion implementation)
void KerCalcConcentration(ACC<Real>& concentration, const ACC<Real>& xiIndex, const Real* TIME, const ACC<Real>& U, const ACC<Real>& V, const Real* Dt) {
    Real x0[2] = {100.0,100.0};
    Real xiIndexVar1 = {xiIndex(0,0,0)};
    Real xiIndexVar2 = {xiIndex(1,0,0)};
    Real x=(xiIndexVar1-x0[0]);
    Real y=(xiIndexVar2-x0[1]);
    Real R=sqrt(x*x+y*y);
    Real r=70.0;
    
    if (R<=r){
        concentration(0,0) = 1;
    }
    else{
        concentration(0,0) = -1;
    }

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
    const Real teta={2*M_PI/3};
    const Real alpha={acos(sin(teta)*sin(teta))};
    const Real Omga={2*((teta <= M_PI/2) - (M_PI/2 < teta))*sqrt(cos(alpha/3)*(1-cos(alpha/3)))};
    const Real K={0.02};
    const Real A=0.02;
    const Real h={Omga*sqrt(K*A/2)};
    VertexType vt = (VertexType)nodeType(0, 0);
    VertexGeometryType vg = (VertexGeometryType)geometry(0, 0);
    VertexGeometryType vg1 = (VertexGeometryType)geometry(1, 0);
    VertexGeometryType vg2 = (VertexGeometryType)geometry(-1, 0);
    VertexGeometryType vg3 = (VertexGeometryType)geometry(0, 1);
    VertexGeometryType vg4 = (VertexGeometryType)geometry(0, -1);
    if (vt == VertexType::Wall) {
        switch (vg) {
            case VG_IP: {
                if(vg3!=VG_IPJP_O&&vg4!=VG_IPJM_O){
                    rho(-1,0)=rho(1,0)+2*h/K;
                }
            } break;
            case VG_IM: {
                if(vg3!=VG_IMJP_O&&vg4!=VG_IMJM_O){
                    rho(1,0)=rho(-1,0)+2*h/K;
                }
            } break;
            case VG_JP: {
                if(vg1!=VG_IPJP_O&&vg2!=VG_IMJP_O){
                    rho(0,-1)=rho(0,1)+2*h/K;
                }
            } break;
            case VG_JM: {
                if(vg1!=VG_IPJM_O&&vg2!=VG_IMJM_O){
                    rho(0,1)=rho(0,-1)+2*h/K;
                }
            } break;
            case VG_IPJP_O: {
                rho(-1,-1)=(rho(1,0)+rho(0,1)+rho(1,-1)+rho(-1,1))/4+2*h/K;

            } break;
            case VG_IPJM_O: {
                rho(-1,1)=(rho(1,0)+rho(0,-1)+rho(1,1)+rho(-1,-1))/4+2*h/K;

            } break;
            case VG_IMJP_O: {
                rho(1,-1)=(rho(-1,0)+rho(0,1)+rho(1,1)+rho(-1,-1))/4+2*h/K;

            } break;
            case VG_IMJM_O: {
                rho(1,1)=(rho(-1,0)+rho(0,-1)+rho(1,-1)+rho(-1,1))/4+2*h/K;
            } break;
            case VG_IPJP_I: {
                if(vg3!=VG_IPJP_O){
                    rho(-1,0)=rho(1,0)+2*h/K;
                }
                if(vg1!=VG_IPJP_O){
                    rho(0,-1)=rho(0,1)+2*h/K;
                }
            } break;
            case VG_IPJM_I: {
                if(vg4!=VG_IPJM_O){
                    rho(-1,0)=rho(1,0)+2*h/K;
                }
                if(vg1!=VG_IPJM_O){
                    rho(0,1)=rho(0,-1)+2*h/K;
                }
            } break;
            case VG_IMJP_I: {
                if(vg3!=VG_IMJP_O){
                    rho(1,0)=rho(-1,0)+2*h/K;
                }
                if(vg2!=VG_IMJP_O){
                    rho(0,-1)=rho(0,1)+2*h/K;
                }
            } break;
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

void KerCalcMu(ACC<Real>&mu,const ACC<Real>&Phi,const ACC<Real>&d2phi, const ACC<int>& geometry, const ACC<int>&nodeType,const ACC<Real>& coordinates){
    Real phi{Phi(0,0)};
    Real phi3 = {phi*phi*phi};
    Real d2phix = {d2phi(0,0,0)};
    Real d2phiy = {d2phi(1,0,0)};
    Real A=0.02;
    Real kappa=0.02;

    mu(0,0)=-A*phi+A*phi3-kappa*(d2phix+d2phiy);
    //if (coordinates(0,0,0)==14||coordinates(0,0,0)==117){
    //std::cout<<phi<<" "<<mu(0,0)<<" "<<d2phix<<" "<<d2phiy<<"   ";
    //}
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

void KerSetSolid(ACC<int>& nodeType, ACC<int>& geometry, const ACC<Real>& coordinates, const int* xMax, const int* yMax) {
    const Real X{coordinates(0,0,0)};
    const Real Y{coordinates(1,0,0)};
    int XMax{*xMax};
    int YMax{*yMax};
    //std::cout<<XMax<<" ";
    if (X==XMax-1||X==0||Y==YMax-1||Y==0){
        nodeType(0,0)=(int)VertexType::ImmersedSolid;
        //std::cout<<nodeType(0,0)<<" ";
        geometry(0,0)=(int)VG_ImmersedSolid;
    }
    //std::cout<<nodeType(0,0)<<" ";
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
            if (VertexType::ImmersedSolid != neiborVertexType[i]&&VertexType::Wall != neiborVertexType[i]) {
                fluidNeiborNum++;
            }
        }
        int solidNeiborNumatCoord{0};
        for (int i = 0; i < 4; i++) {
            if (VertexType::ImmersedSolid == neiborVertexType[i]||VertexType::Wall == neiborVertexType[i]) {
                solidNeiborNumatCoord++;
            }
        }
        if (fluidNeiborNum > 0) {
            
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
                if ((VertexType::ImmersedSolid == neiborVertexType[2]||VertexType::Wall == neiborVertexType[2]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]||VertexType::Wall == neiborVertexType[1])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[6]||VertexType::Wall == neiborVertexType[6]) {
                        nodeType(0,0)=(int)VertexType::Wall;
                        geometryProperty(0, 0) = (int)VG_IPJM_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points");
                    }
                }

                if ((VertexType::ImmersedSolid == neiborVertexType[3]||VertexType::Wall == neiborVertexType[3]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]||VertexType::Wall == neiborVertexType[1])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[5]||VertexType::Wall == neiborVertexType[5]) {
                        nodeType(0,0)=(int)VertexType::Wall;
                        geometryProperty(0, 0) = (int)VG_IPJP_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points");
                    }
                }

                if ((VertexType::ImmersedSolid == neiborVertexType[3]||VertexType::Wall == neiborVertexType[3]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[0]||VertexType::Wall == neiborVertexType[0])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[7]||VertexType::Wall == neiborVertexType[7]) {
                        nodeType(0,0)=(int)VertexType::Wall;
                        geometryProperty(0, 0) = (int)VG_IMJP_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points");
                    }
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[2]||VertexType::Wall == neiborVertexType[2]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[0]||VertexType::Wall == neiborVertexType[0])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[4]||VertexType::Wall == neiborVertexType[4]) {
                        nodeType(0,0)=(int)VertexType::Wall;
                        geometryProperty(0, 0) = (int)VG_IMJM_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points");
                    }
                }
            }
            // Planlar corner
            if (3 == solidNeiborNumatCoord) {
                nodeType(0,0)=(int)VertexType::Wall;
                if (VertexType::ImmersedSolid != neiborVertexType[0]&&VertexType::Wall != neiborVertexType[0]) {
                    geometryProperty(0, 0) = (int)VG_IP;
                }
                if (VertexType::ImmersedSolid != neiborVertexType[1]&&VertexType::Wall != neiborVertexType[1]) {
                    geometryProperty(0, 0) = (int)VG_IM;
                }
                if (VertexType::ImmersedSolid != neiborVertexType[2]&&VertexType::Wall != neiborVertexType[2]) {
                    geometryProperty(0, 0) = (int)VG_JP;
                }
                if (VertexType::ImmersedSolid != neiborVertexType[3]&&VertexType::Wall != neiborVertexType[3]) {
                    geometryProperty(0, 0) = (int)VG_JM;
                }
            }
            // Inner corner
            if (4 == solidNeiborNumatCoord) {
                nodeType(0,0)=(int)VertexType::Wall;
                if (1 == fluidNeiborNum) {
                    if (VertexType::ImmersedSolid != neiborVertexType[4]&&VertexType::Wall != neiborVertexType[4]) {
                        geometryProperty(0, 0) = (int)VG_IPJP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[5]&&VertexType::Wall != neiborVertexType[5]) {
                        geometryProperty(0, 0) = (int)VG_IMJM_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[6]&&VertexType::Wall != neiborVertexType[6]) {
                        geometryProperty(0, 0) = (int)VG_IMJP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[7]&&VertexType::Wall != neiborVertexType[7]) {
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

void KerCalcConcentration3D(ACC<Real>& concentration, const ACC<Real>& xiIndex, const Real* TIME, const ACC<Real>& U, const ACC<Real>& V, const Real* Dt) {
    Real x0[3] = {0.0,20.0,20.0};
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
    const Real teta={M_PI/3};
    const Real alpha={acos(sin(teta)*sin(teta))};
    const Real Omga={2*((teta <= M_PI/2) - (M_PI/2 < teta))*sqrt(cos(alpha/3)*(1-cos(alpha/3)))};
    const Real K={0.05};
    const Real A=0.05;
    const Real h={Omga*sqrt(K*A/2)};
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
                if(vg2!=VG_IPJP_O&&vg3!=VG_IPJM_O&&vg4!=VG_IPKP_O&&vg5!=VG_IPKM_O){
                    rho(-1,0,0)=rho(1,0,0)+2*h/K;
                }
            } break;
            case VG_IM: {
                if(vg2!=VG_IMJP_O&&vg3!=VG_IMJM_O&&vg4!=VG_IMKP_O&&vg5!=VG_IMKM_O){
                    rho(1,0,0)=rho(-1,0,0)+2*h/K;
                }
            } break;
            case VG_JP: {
                if(vg0!=VG_IPJP_O&&vg1!=VG_IMJP_O&&vg4!=VG_JPKP_O&&vg5!=VG_JPKM_O){
                    rho(0,-1,0)=rho(0,1,0)+2*h/K;
                }
            } break;
            case VG_JM: {
                if(vg0!=VG_IPJM_O&&vg1!=VG_IMJM_O&&vg4!=VG_JMKP_O&&vg5!=VG_JMKM_O){
                    rho(0,1,0)=rho(0,-1,0)+2*h/K;
                }
            } break;
            case VG_KP: {
                if(vg0!=VG_IPKP_O&&vg1!=VG_IMKP_O&&vg2!=VG_JPKP_O&&vg3!=VG_JMKP_O){
                    rho(0,0,-1)=rho(0,0,1)+2*h/K;
                }
            } break;
            case VG_KM: {
                if(vg0!=VG_IPKM_O&&vg1!=VG_IMKM_O&&vg2!=VG_JPKM_O&&vg3!=VG_JMKM_O){
                    rho(0,0,1)=rho(0,0,-1)+2*h/K;
                }
            } break;
            case VG_IPJP_O: {
                if(vg4!=VG_IPJPKP_O&&vg5!=VG_IPJPKM_O){
                    rho(-1,-1,0)=(rho(1,0,0)+rho(0,1,0)+rho(1,-1,0)+rho(-1,1,0))/4+2*h/K;
                }
            } break;
            case VG_IPJM_O: {
                if(vg4!=VG_IPJMKP_O&&vg5!=VG_IPJMKM_O){
                    rho(-1,1,0)=(rho(1,0,0)+rho(0,-1,0)+rho(-1,-1,0)+rho(1,1,0))/4+2*h/K;
                }
            } break;
            case VG_IMJP_O: {
                if(vg4!=VG_IMJPKP_O&&vg5!=VG_IMJPKM_O){
                    rho(1,-1,0)=(rho(-1,0,0)+rho(0,1,0)+rho(-1,-1,0)+rho(1,1,0))/4+2*h/K;
                }
            } break;
            case VG_IMJM_O: {
                if(vg4!=VG_IMJMKP_O&&vg5!=VG_IMJMKM_O){
                    rho(1,1,0)=(rho(-1,0,0)+rho(0,-1,0)+rho(1,-1,0)+rho(-1,1,0))/4+2*h/K;
                }
            } break;
            case VG_IPKP_O: {
                if(vg2!=VG_IPJPKP_O&&vg3!=VG_IPJMKP_O){
                    rho(-1,0,-1)=(rho(1,0,0)+rho(0,0,1)+rho(1,0,-1)+rho(1,0,-1))/4+2*h/K;
                }
            } break;
            case VG_IPKM_O: {
                if(vg2!=VG_IPJPKM_O&&vg3!=VG_IPJMKM_O){
                    rho(-1,0,1)=(rho(1,0,0)+rho(0,0,-1)+rho(-1,0,-1)+rho(1,0,1))/4+2*h/K;
                }
            } break;
            case VG_IMKP_O: {
                if(vg2!=VG_IMJPKP_O&&vg3!=VG_IMJMKP_O){
                    rho(1,0,-1)=(rho(-1,0,0)+rho(0,0,1)+rho(-1,0,-1)+rho(1,0,1))/4+2*h/K;
                }
            } break;
            case VG_IMKM_O: {
                if(vg2!=VG_IMJPKM_O&&vg3!=VG_IMJMKM_O){
                    rho(1,0,1)=(rho(-1,0,0)+rho(0,0,-1)+rho(1,0,-1)+rho(1,0,-1))/4+2*h/K;
                }
            } break;
            case VG_JMKP_O: {
                if(vg0!=VG_IPJMKP_O&&vg1!=VG_IMJMKP_O){
                    rho(0,1,-1)=(rho(0,-1,0)+rho(0,0,1)+rho(0,-1,-1)+rho(0,1,1))/4+2*h/K;
                }
            } break;
            case VG_JMKM_O: {
                if(vg0!=VG_IPJMKM_O&&vg1!=VG_IMJMKM_O){
                    rho(0,1,1)=(rho(0,-1,0)+rho(0,0,-1)+rho(0,-1,1)+rho(0,1,-1))/4+2*h/K;
                }
            } break;
            case VG_JPKP_O: {
                if(vg0!=VG_IPJPKP_O&&vg1!=VG_IMJPKP_O){
                    rho(0,-1,-1)=(rho(0,1,0)+rho(0,0,1)+rho(0,-1,1)+rho(0,1,-1))/4+2*h/K;
                }
            } break;
            case VG_JPKM_O: {
                if(vg0!=VG_IPJPKM_O&&vg1!=VG_IMJPKM_O){
                    rho(0,-1,1)=(rho(0,1,0)+rho(0,0,-1)+rho(0,-1,1)+rho(0,1,-1))/4+2*h/K;
                }
            } break;
            case VG_IPJPKP_O: {
                rho(-1,-1,-1)=(rho(1,0,0)+rho(0,1,0)+rho(0,0,1)+rho(1,-1,-1)+rho(-1,1,-1)+rho(-1,-1,1)+
                rho(1,0,-1)+rho(0,1,-1)+rho(0,-1,1)+rho(1,-1,0)+rho(-1,1,0)+rho(-1,0,1))/12+2*h/K;


            } break;
            case VG_IPJMKP_O: {
                rho(-1,1,-1)=(rho(1,0,0)+rho(0,-1,0)+rho(0,0,1)+rho(1,1,-1)+rho(-1,-1,-1)+rho(-1,1,1)+
                rho(1,0,-1)+rho(0,-1,-1)+rho(0,1,1)+rho(1,1,0)+rho(-1,-1,0)+rho(-1,0,1))/12+2*h/K;


            } break;
            case VG_IMJPKP_O: {
                rho(1,-1,-1)=(rho(-1,0,0)+rho(0,1,0)+rho(0,0,1)+rho(-1,-1,-1)+rho(1,1,-1)+rho(1,-1,1)+
                rho(-1,0,-1)+rho(0,1,-1)+rho(0,-1,1)+rho(-1,-1,0)+rho(1,1,0)+rho(1,0,1))/12+2*h/K;


            } break;
            case VG_IMJMKP_O: {
                rho(1,1,-1)=(rho(-1,0,0)+rho(0,-1,0)+rho(0,0,1)+rho(-1,1,-1)+rho(1,-1,-1)+rho(1,1,1)+
                rho(-1,0,-1)+rho(0,-1,-1)+rho(0,1,1)+rho(-1,1,0)+rho(1,-1,0)+rho(1,0,1))/12+2*h/K;

            } break;
            case VG_IPJPKM_O: {
                rho(-1,-1,1)=(rho(1,0,0)+rho(0,1,0)+rho(0,0,-1)+rho(1,-1,1)+rho(-1,1,1)+rho(-1,-1,-1)+
                rho(1,0,1)+rho(0,1,1)+rho(0,-1,-1)+rho(1,-1,0)+rho(-1,1,0)+rho(-1,0,-1))/12+2*h/K;


            } break;
            case VG_IPJMKM_O: {
                rho(-1,1,1)=(rho(1,0,0)+rho(0,-1,0)+rho(0,0,-1)+rho(1,1,1)+rho(-1,-1,1)+rho(-1,1,-1)+
                rho(1,0,1)+rho(0,-1,1)+rho(0,1,-1)+rho(1,1,0)+rho(-1,-1,0)+rho(-1,0,-1))/12+2*h/K;


            } break;
            case VG_IMJPKM_O: {
                rho(1,-1,1)=(rho(-1,0,0)+rho(0,1,0)+rho(0,0,-1))+rho(-1,-1,1)+rho(1,1,1)+rho(1,-1,-1)+
                rho(-1,0,1)+rho(0,1,1)+rho(0,-1,-1)+rho(-1,-1,0)+rho(1,1,0)+rho(1,0,-1)/12+2*h/K;


            } break;
            case VG_IMJMKM_O: {
                rho(1,1,1)=(rho(-1,0,0)+rho(0,-1,0)+rho(0,0,-1)+rho(-1,1,1)+rho(1,-1,1)+rho(1,1,-1)+
                rho(-1,0,1)+rho(0,-1,1)+rho(0,1,-1)+rho(-1,1,0)+rho(1,-1,0)+rho(1,0,-1))/12+2*h/K;

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
                    rho(-1,-1,0)=(rho(1,0,0)+rho(0,1,0)+rho(1,-1,0)+rho(-1,1,0))/4+2*h/K;

                }
            } break;
            case VG_IPJM_O_KP: {
                if(vg0!=VG_IPKP_O&&vg3!=VG_JMKP_O){
                    rho(0,0,-1)=rho(0,0,1)+2*h/K;

                }
                if(vg4!=VG_IPJMKP_O){
                    rho(-1,1,0)=(rho(1,0,0)+rho(0,-1,0)+rho(-1,-1,0)+rho(1,1,0))/4+2*h/K;

                }
            } break;
            case VG_IMJP_O_KP: {
                if(vg1!=VG_IMKP_O&&vg2!=VG_JPKP_O){
                    rho(0,0,-1)=rho(0,0,1)+2*h/K;

                }
                if(vg4!=VG_IMJPKP_O){
                    rho(1,-1,0)=(rho(-1,0,0)+rho(0,1,0)+rho(-1,-1,0)+rho(1,1,0))/4+2*h/K;

                }
            } break;
            case VG_IMJM_O_KP: {
                if(vg1!=VG_IMKP_O&&vg3!=VG_JMKP_O){
                    rho(0,0,-1)=rho(0,0,1)+2*h/K;

                }
                if(vg4!=VG_IMJMKP_O){
                    rho(1,1,0)=(rho(-1,0,0)+rho(0,-1,0)+rho(1,-1,0)+rho(-1,1,0))/4+2*h/K;
                }
            } break;
            case VG_IPJP_O_KM: {
                if(vg0!=VG_IPKM_O&&vg2!=VG_JPKM_O){
                    rho(0,0,1)=rho(0,0,-1)+2*h/K;

                }
                if(vg5!=VG_IPJPKM_O){
                    rho(-1,-1,0)=(rho(1,0,0)+rho(0,1,0)+rho(1,-1,0)+rho(-1,1,0))/4+2*h/K;

                }
            } break;
            case VG_IPJM_O_KM: {
                if(vg0!=VG_IPKM_O&&vg3!=VG_JMKM_O){
                    rho(0,0,1)=rho(0,0,-1)+2*h/K;

                }
                if(vg5!=VG_IPJMKM_O){
                    rho(-1,1,0)=(rho(1,0,0)+rho(0,-1,0)+rho(-1,-1,0)+rho(1,1,0))/4+2*h/K;

                }
            } break;
            case VG_IMJP_O_KM: {
                if(vg1!=VG_IMKP_O&&vg2!=VG_JPKM_O){
                    rho(0,0,1)=rho(0,0,-1)+2*h/K;

                }
                if(vg5!=VG_IMJPKM_O){
                    rho(1,-1,0)=(rho(-1,0,0)+rho(0,1,0)+rho(-1,-1,0)+rho(1,1,0))/4+2*h/K;

                }
            } break;
            case VG_IMJM_O_KM: {
                if(vg1!=VG_IMKM_O&&vg3!=VG_JMKM_O){
                    rho(0,0,1)=rho(0,0,-1)+2*h/K;

                }
                if(vg5!=VG_IMJMKM_O){
                    rho(1,1,0)=(rho(-1,0,0)+rho(0,-1,0)+rho(1,-1,0)+rho(-1,1,0))/4+2*h/K;
                }
            } break;
            case VG_IPKP_O_JP: {
                if(vg0!=VG_IPJP_O&&vg4!=VG_JPKP_O){
                    rho(0,-1,0)=rho(0,1,0)+2*h/K;

                }
                if(vg2!=VG_IPJPKP_O){
                    rho(-1,0,-1)=(rho(1,0,0)+rho(0,0,1)+rho(1,0,-1)+rho(1,0,-1))/4+2*h/K;
                }
            } break;
            case VG_IPKM_O_JP: {
                if(vg0!=VG_IPJP_O&&vg5!=VG_JPKM_O){
                    rho(0,-1,0)=rho(0,1,0)+2*h/K;

                }
                if(vg2!=VG_IPJPKM_O){
                    rho(-1,0,1)=(rho(1,0,0)+rho(0,0,-1)+rho(-1,0,-1)+rho(1,0,1))/4+2*h/K;
                }
            } break;
            case VG_IMKP_O_JP: {
                if(vg1!=VG_IMJP_O&&vg4!=VG_JPKP_O){
                    rho(0,-1,0)=rho(0,1,0)+2*h/K;

                }
                if(vg2!=VG_IMJPKM_O){
                    rho(1,0,-1)=(rho(-1,0,0)+rho(0,0,1)+rho(-1,0,-1)+rho(1,0,1))/4+2*h/K;
                }
            } break;
            case VG_IMKM_O_JP: {
                if(vg1!=VG_IMJP_O&&vg5!=VG_JPKM_O){
                    rho(0,-1,0)=rho(0,1,0)+2*h/K;

                }
                if(vg2!=VG_IMJPKM_O){
                    rho(1,0,1)=(rho(-1,0,0)+rho(0,0,-1)+rho(1,0,-1)+rho(1,0,-1))/4+2*h/K;
                }
            } break;
            case VG_IPKP_O_JM: {
                if(vg0!=VG_IPJM_O&&vg4!=VG_JMKP_O){
                    rho(0,1,0)=rho(0,-1,0)+2*h/K;

                }
                if(vg3!=VG_IPJMKP_O){
                    rho(-1,0,-1)=(rho(1,0,0)+rho(0,0,1)+rho(1,0,-1)+rho(1,0,-1))/4+2*h/K;
                }
            } break;
            case VG_IPKM_O_JM: {
                if(vg0!=VG_IPJM_O&&vg5!=VG_JMKM_O){
                    rho(0,1,0)=rho(0,-1,0)+2*h/K;

                }
                if(vg3!=VG_IPJMKM_O){
                    rho(-1,0,1)=(rho(1,0,0)+rho(0,0,-1)+rho(-1,0,-1)+rho(1,0,1))/4+2*h/K;
                }
            } break;
            case VG_IMKP_O_JM: {
                if(vg1!=VG_IMJM_O&&vg4!=VG_JMKP_O){
                    rho(0,1,0)=rho(0,-1,0)+2*h/K;

                }
                if(vg3!=VG_IMJMKP_O){
                    rho(1,0,-1)=(rho(-1,0,0)+rho(0,0,1)+rho(-1,0,-1)+rho(1,0,1))/4+2*h/K;
                }
            } break;
            case VG_IMKM_O_JM: {
                if(vg1!=VG_IMJM_O&&vg5!=VG_JMKM_O){
                    rho(0,1,0)=rho(0,-1,0)+2*h/K;

                }
                if(vg3!=VG_IMJMKM_O){
                    rho(1,0,1)=(rho(-1,0,0)+rho(0,0,-1)+rho(1,0,-1)+rho(1,0,-1))/4+2*h/K;
                }
            } break;
            case VG_JPKP_O_IP: {
                if(vg2!=VG_IPJP_O&&vg4!=VG_IPKP_O){
                    rho(-1,0,0)=rho(1,0,0)+2*h/K;
                }
                if(vg0!=VG_IPJPKP_O){
                    rho(0,-1,-1)=(rho(0,1,0)+rho(0,0,1)+rho(0,-1,1)+rho(0,1,-1))/4+2*h/K;
                }
            } break;
            case VG_JPKM_O_IP: {
                if(vg2!=VG_IPJP_O&&vg5!=VG_IPKM_O){
                    rho(-1,0,0)=rho(1,0,0)+2*h/K;
                }
                if(vg0!=VG_IPJPKM_O){
                    rho(0,-1,1)=(rho(0,1,0)+rho(0,0,-1)+rho(0,-1,1)+rho(0,1,-1))/4+2*h/K;
                }
            } break;
            case VG_JMKP_O_IP: {
                if(vg3!=VG_IPJM_O&&vg4!=VG_IPKP_O){
                    rho(-1,0,0)=rho(1,0,0)+2*h/K;
                }
                if(vg0!=VG_IPJMKP_O){
                    rho(0,1,-1)=(rho(0,-1,0)+rho(0,0,1)+rho(0,-1,-1)+rho(0,1,1))/4+2*h/K;
                }
            } break;
            case VG_JMKM_O_IP: {
                if(vg3!=VG_IPJM_O&&vg5!=VG_IPKM_O){
                    rho(-1,0,0)=rho(1,0,0)+2*h/K;
                }
                if(vg0!=VG_IPJMKM_O){
                    rho(0,1,1)=(rho(0,-1,0)+rho(0,0,-1)+rho(0,-1,1)+rho(0,1,-1))/4+2*h/K;
                }
            } break;
            case VG_JPKP_O_IM: {
                if(vg2!=VG_IMJP_O&&vg4!=VG_IMKP_O){
                    rho(1,0,0)=rho(-1,0,0)+2*h/K;
                }
                if(vg1!=VG_IMJPKP_O){
                    rho(0,-1,-1)=(rho(0,1,0)+rho(0,0,1)+rho(0,-1,1)+rho(0,1,-1))/4+2*h/K;
                }
            } break;
            case VG_JPKM_O_IM: {
                if(vg2!=VG_IMJP_O&&vg5!=VG_IMKM_O){
                    rho(1,0,0)=rho(-1,0,0)+2*h/K;
                }
                if(vg1!=VG_IMJPKM_O){
                    rho(0,-1,1)=(rho(0,1,0)+rho(0,0,-1)+rho(0,-1,1)+rho(0,1,-1))/4+2*h/K;
                }
            } break;
            case VG_JMKP_O_IM: {
                if(vg3!=VG_IMJM_O&&vg4!=VG_IMKP_O){
                    rho(1,0,0)=rho(-1,0,0)+2*h/K;
                }
                if(vg1!=VG_IMJMKP_O){
                    rho(0,1,-1)=(rho(0,-1,0)+rho(0,0,1)+rho(0,-1,-1)+rho(0,1,1))/4+2*h/K;
                }
            } break;
            case VG_JMKM_O_IM: {
                if(vg3!=VG_IMJM_O&&vg5!=VG_IMKM_O){
                    rho(1,0,0)=rho(-1,0,0)+2*h/K;
                }
                if(vg1!=VG_IMJMKM_O){
                    rho(0,1,1)=(rho(0,-1,0)+rho(0,0,-1)+rho(0,-1,1)+rho(0,1,-1))/4+2*h/K;
                }
            } break;
        }
    }
}

//Calculate chemical potential within Free Energy Model
void KerCalcMu3D(ACC<Real>&mu,const ACC<Real>&Phi,const ACC<Real>&d2phi, const ACC<int>&nodeType){
    Real phi{Phi(0,0,0)};
    Real phi3 = {phi*phi*phi};
    Real A=0.05;
    Real kappa=0.05;
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (vt != VertexType::ImmersedSolid){
    mu(0,0,0)=-A*phi+A*phi3-kappa*(d2phi(0,0,0,0)+d2phi(1,0,0,0)+d2phi(2,0,0,0));
    }
}

//Function to print a 1D macroscopic variable for debugging purposes
void KerPrintPhi3D(const ACC<Real>& Phi, const ACC<Real>& coordinates){
    if ((coordinates(0,0,0)==60&&coordinates(1,0,0)==10&&coordinates(2,0,0)==16)){
        std::cout<<Phi(-1,0,-1)<<" ";
    }
}

//Update surface tension force between fluid components
void KerUpdateMacroBodyForce3D(ACC<Real>& force, const ACC<Real>& Phi, const ACC<Real>& dmu, const ACC<int>& geometry, const ACC<int>&nodeType, const ACC<Real>& coordinates,
                            const int* idx) {
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
    //std::cout<<"Force: "<<force(0,0,0)<<" "<<force(1,0,0)<<"\n";
}

void KerSetSolid3D(ACC<int>& nodeType, ACC<int>& geometry, const ACC<Real>& coordinates, const int* xMax, const int* yMax, const int* zMax) {
    const Real X{coordinates(0,0,0,0)};
    const Real Y{coordinates(1,0,0,0)};
    const Real Z{coordinates(2,0,0,0)};
    int XMax{*xMax};
    int YMax{*yMax};
    int ZMax{*zMax};
    //std::cout<<XMax<<" ";
    if (X<0.5*Y||X==XMax-1||X==0||Y==YMax-1||Y==0||Z==ZMax-1||Z==0){
        nodeType(0,0,0)=(int)VertexType::ImmersedSolid;
        //std::cout<<nodeType(0,0)<<" ";
        geometry(0,0,0)=(int)VG_ImmersedSolid;
    }
    //std::cout<<nodeType(0,0)<<" ";
}
/*
void KerSetSolid3D(ACC<int>& nodeType, ACC<int>& geometry, const ACC<Real>& coordinates, const int* xMax, const int* yMax, const int* zMax) {
  const Real xk{coordinates(0,0,0,0)};
  const Real yk{coordinates(1,0,0,0)};
  const Real zk{coordinates(2,0,0,0)};
  int NCZ=1;
  int NCX=1;
  int LX{*xMax};
  int LY{*yMax};
  int LZ{*zMax};
  int CZ=(LZ)/NCZ;
  int CX=(LX)/(NCX+1);
  int CY=(LY)-10;
  double bAng=15;
  double oHeight=0.5;
  int pHeight=8;
  int height=oHeight/1.5*CZ;
  int res2=0.44/2*CX;
  int tnss=0.17/1.5*CZ;
  double ang2=bAng*M_PI/180;
  double grad1=0.3/0.4;
  double ang1=atan(1/grad1);
  double grad2=1 / tan (ang2);
  //double grad2=0.8;
  double tnss2=tnss/sin (ang1);
  //Reservoir length in x dir
  //int rsv=80;
  int rightLeft=1;
  int rsv=CX+5; //default
  int ww=10;
  int aPitch=0.11/2*CX;
  int bPitch=CY/2;
  int j=0;
  int i=0;
  if ( (zk>CZ*j+height+tnss	|| zk<CZ*j+height) && (xk==0)&& (yk<CY+5)) {
      nodeType(0,0,0)=(int)VertexType::ImmersedSolid;
      geometry(0,0,0)=(int)VG_ImmersedSolid;
  }
  //Floor/Ceiling of cells
  if ( (zk==0 || zk==LZ-1 ) && (xk>=rightLeft*rsv) && (xk<=rightLeft*rsv+NCX*CX) && (yk<CY+5)) {
      nodeType(0,0,0)=(int)VertexType::ImmersedSolid;
      geometry(0,0,0)=(int)VG_ImmersedSolid;
  }
  
  //Y dir wall
  if ( (yk==0 || yk==LY-1) && (xk>=0) && (xk<LX) )  {
      nodeType(0,0,0)=(int)VertexType::ImmersedSolid;
      geometry(0,0,0)=(int)VG_ImmersedSolid;
  }
  
  //Entrance channel
  if ( (zk>CZ*j+height+tnss	|| zk<CZ*j+height) && (xk>=rightLeft*rsv-5) && (xk<-ww*(rightLeft-1)+rsv*rightLeft) && (zk<(j+1)*CZ) && (zk>=j*CZ) && (yk<CY+5))  {nodeType(0,0,0)=(int)VertexType::ImmersedSolid;geometry(0,0,0)=(int)VG_ImmersedSolid;}
  
  //Entrance slope (upper)
  if ( ((zk-j*CZ)+grad1*(xk-rightLeft*rsv-i*CX)<=height+tnss2) && (zk>CZ*j+height+tnss	) && (xk>=rightLeft*rsv+i*CX) && (xk-rightLeft*rsv<=i*CX+((tnss2-tnss)/grad1)) && (zk<(j+1)*CZ) && (yk<CY+5))  {
      nodeType(0,0,0)=(int)VertexType::ImmersedSolid;
      geometry(0,0,0)=(int)VG_ImmersedSolid;
  }
  
  //Entrance slope (lower)
  if ( (zk<CZ*j+height) && (xk>=rightLeft*rsv+i*CX) && ((zk-j*CZ)+grad1*(xk-rightLeft*rsv-i*CX)<height) && (zk>j*CZ) && (yk<CY+5))  {
      nodeType(0,0,0)=(int)VertexType::ImmersedSolid;
      geometry(0,0,0)=(int)VG_ImmersedSolid;
  }
  
  //Conical left wall
  if ( ((zk-j*CZ)+grad1*(xk-rightLeft*rsv-i*CX)>height+tnss2) && (xk>=rightLeft*rsv+i*CX) && ((zk-2*tnss-j*CZ)+(1 / tan (ang2))*(xk-rightLeft*rsv-i*CX-CX+tnss+res2)<0) && ((xk-rightLeft*rsv-i*CX+res2)>((tnss2-tnss)/grad1)) && ((xk-rightLeft*rsv-i*CX)<CX-tnss-res2) && (zk-CZ*j>tnss) && ((zk-CZ*j-2*tnss)-grad1*0.3*(xk-rightLeft*rsv-i*CX-CX+res2)>0) && (zk<(j+1)*CZ) && (yk<CY+5))  {
      nodeType(0,0,0)=(int)VertexType::ImmersedSolid;
      geometry(0,0,0)=(int)VG_ImmersedSolid;
  }
  if ( ((CZ-2.5*tnss)+(1 / tan (ang2))*(xk-rightLeft*rsv-i*CX-CX+tnss+res2)<0) && ((CZ-2.2*tnss)+(1 / tan (ang2))*(xk-rightLeft*rsv-i*CX-CX+tnss+res2)>=0) && (zk>CZ*j+6*tnss) && (yk<CY+5))  {
      nodeType(0,0,0)=(int)VertexType::ImmersedSolid;
      geometry(0,0,0)=(int)VG_ImmersedSolid;
  }
  
  //Lower channel slope
  if ( ((zk-j*CZ-tnss)-grad1*0.3*(xk-rightLeft*rsv-CX-i*CX+res2)<0) && (xk-rightLeft*rsv-i*CX+res2<CX) && (xk>=rightLeft*rsv+i*CX) && (zk>j*CZ) && (yk<CY+5))  {
      nodeType(0,0,0)=(int)VertexType::ImmersedSolid;
      geometry(0,0,0)=(int)VG_ImmersedSolid;
  }
  
  //Pitch
  if ( (zk>CZ*j+height+tnss	|| zk<CZ*j+height) && (xk<=rightLeft*rsv+i*CX+CX) && (xk>=rightLeft*rsv+i*CX+CX-res2) && (zk<(j+1)*CZ) && (zk>=j*CZ) && (yk<CY+5))  {
      nodeType(0,0,0)=(int)VertexType::ImmersedSolid;
      geometry(0,0,0)=(int)VG_ImmersedSolid;
  }
  if ( (((yk-10+pHeight)*(yk-10+pHeight))+(bPitch*bPitch)/(aPitch*aPitch)*(xk-(rightLeft*rsv+i*CX+CX-res2+aPitch))*(xk-(rightLeft*rsv+i*CX+CX-res2+aPitch))<=bPitch*bPitch) && (yk<CY+5))  {
      nodeType(0,0,0)=(int)VertexType::ImmersedSolid;
      geometry(0,0,0)=(int)VG_ImmersedSolid;
  }
  
  //Exit channel
  if ( (zk>CZ*j+height+tnss	|| zk<CZ*j+height) && (xk<=rightLeft*rsv+NCX*CX+ww) && (xk>=rightLeft*rsv+NCX*CX) && (zk<(j+1)*CZ) && (zk>=j*CZ) && (yk<CY+5))  {
      nodeType(0,0,0)=(int)VertexType::ImmersedSolid;
      geometry(0,0,0)=(int)VG_ImmersedSolid;
  }
  
  
}
*/


void KerSetEmbeddedBodyGeometry3D(ACC<int>& geometryProperty,
                                ACC<int>& nodeType) {
#ifdef OPS_3D
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (VertexType::ImmersedSolid == vt) {
        VertexType neiborVertexType[26];
        /*
                        6*****2*****4
                        *     *     *
                        *     *     *
                        1***********0
                        *     *     *
                        *     *     *
                        5*****3*****7
        */
        neiborVertexType[0] = (VertexType)nodeType(1, 0, 0);
        neiborVertexType[1] = (VertexType)nodeType(-1, 0, 0);
        neiborVertexType[2] = (VertexType)nodeType(0, 1, 0);
        neiborVertexType[3] = (VertexType)nodeType(0, -1, 0);
        neiborVertexType[4] = (VertexType)nodeType(0, 0, 1);
        neiborVertexType[5] = (VertexType)nodeType(0, 0, -1);
        neiborVertexType[6] = (VertexType)nodeType(1, 1, 0);
        neiborVertexType[7] = (VertexType)nodeType(1, -1, 0);
        neiborVertexType[8] = (VertexType)nodeType(-1, 1, 0);
        neiborVertexType[9] = (VertexType)nodeType(-1, -1, 0);
        neiborVertexType[10] = (VertexType)nodeType(1, 0, 1);
        neiborVertexType[11] = (VertexType)nodeType(1, 0, -1);
        neiborVertexType[12] = (VertexType)nodeType(-1, 0, 1);
        neiborVertexType[13] = (VertexType)nodeType(-1, 0, -1);
        neiborVertexType[14] = (VertexType)nodeType(0, 1, 1);
        neiborVertexType[15] = (VertexType)nodeType(0, 1, -1);
        neiborVertexType[16] = (VertexType)nodeType(0, -1, 1);
        neiborVertexType[17] = (VertexType)nodeType(0, -1, -1);
        neiborVertexType[18] = (VertexType)nodeType(1, 1, 1);
        neiborVertexType[19] = (VertexType)nodeType(1, 1, -1);
        neiborVertexType[20] = (VertexType)nodeType(1, -1, 1);
        neiborVertexType[21] = (VertexType)nodeType(1, -1, -1);
        neiborVertexType[22] = (VertexType)nodeType(-1, 1, 1);
        neiborVertexType[23] = (VertexType)nodeType(-1, 1, -1);
        neiborVertexType[24] = (VertexType)nodeType(-1, -1, 1);
        neiborVertexType[25] = (VertexType)nodeType(-1, -1, -1);

        int fluidNeiborNum{0};
        for (int i = 0; i < 26; i++) {
            if (VertexType::ImmersedSolid != neiborVertexType[i]&&VertexType::Wall != neiborVertexType[i]) {
                fluidNeiborNum++;
            }
        }
        int solidNeiborNumatCoord{0};
        for (int i = 0; i < 6; i++) {
            if (VertexType::ImmersedSolid == neiborVertexType[i]||VertexType::Wall == neiborVertexType[i]) {
                solidNeiborNumatCoord++;
            }
        }
        //std::cout<<solidNeiborNumatCoord<<" ";
        if (fluidNeiborNum > 0) {
            
            // outer corner
            if (3 == solidNeiborNumatCoord) {
                if ((VertexType::ImmersedSolid == neiborVertexType[5]||VertexType::Wall == neiborVertexType[5]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[3]||VertexType::Wall == neiborVertexType[3]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]||VertexType::Wall == neiborVertexType[1])) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IPJPKP_O;
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[4]||VertexType::Wall == neiborVertexType[4]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[3]||VertexType::Wall == neiborVertexType[3]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]||VertexType::Wall == neiborVertexType[1])) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IPJPKM_O;
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[5]||VertexType::Wall == neiborVertexType[5]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[2]||VertexType::Wall == neiborVertexType[2]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]||VertexType::Wall == neiborVertexType[1])) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IPJMKP_O;
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[4]||VertexType::Wall == neiborVertexType[4]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[2]||VertexType::Wall == neiborVertexType[2]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]||VertexType::Wall == neiborVertexType[1])) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IPJMKM_O;
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[5]||VertexType::Wall == neiborVertexType[5]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[3]||VertexType::Wall == neiborVertexType[3]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[0]||VertexType::Wall == neiborVertexType[0])) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IMJPKP_O;
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[4]||VertexType::Wall == neiborVertexType[4]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[3]||VertexType::Wall == neiborVertexType[3]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[0]||VertexType::Wall == neiborVertexType[0])) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IMJPKM_O;
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[5]||VertexType::Wall == neiborVertexType[5]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[2]||VertexType::Wall == neiborVertexType[2]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[0]||VertexType::Wall == neiborVertexType[0])) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IMJMKP_O;
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[4]||VertexType::Wall == neiborVertexType[4]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[2]||VertexType::Wall == neiborVertexType[2]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[0]||VertexType::Wall == neiborVertexType[0])) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IMJMKM_O;
                }
            }
            if (4 == solidNeiborNumatCoord) {
                if ((VertexType::ImmersedSolid == neiborVertexType[0] &&
                     VertexType::ImmersedSolid == neiborVertexType[1]) ||
                    (VertexType::ImmersedSolid == neiborVertexType[0] &&
                     VertexType::ImmersedSolid == neiborVertexType[1])) {
                    ops_printf("%s\n",
                               "There appears to be hanged solid points,i.e., "
                               "the solid body may be too thin");
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[3]||VertexType::Wall == neiborVertexType[3]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]||VertexType::Wall == neiborVertexType[1]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[4]||VertexType::Wall == neiborVertexType[4]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[5]||VertexType::Wall == neiborVertexType[5])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[9]||VertexType::Wall == neiborVertexType[9]) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IPJP_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points2");
                    }
                }

                if ((VertexType::ImmersedSolid == neiborVertexType[2]||VertexType::Wall == neiborVertexType[2]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]||VertexType::Wall == neiborVertexType[1]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[4]||VertexType::Wall == neiborVertexType[4]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[5]||VertexType::Wall == neiborVertexType[5])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[8]||VertexType::Wall == neiborVertexType[8]) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IPJM_O;
                        
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points2");
                    }
                }

                if ((VertexType::ImmersedSolid == neiborVertexType[3]||VertexType::Wall == neiborVertexType[3]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[0]||VertexType::Wall == neiborVertexType[0])&&
                    (VertexType::ImmersedSolid == neiborVertexType[4]||VertexType::Wall == neiborVertexType[4]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[5]||VertexType::Wall == neiborVertexType[5])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[7]||VertexType::Wall == neiborVertexType[7]) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IMJP_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points2");
                    }
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[2]||VertexType::Wall == neiborVertexType[2]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[0]||VertexType::Wall == neiborVertexType[0])&&
                    (VertexType::ImmersedSolid == neiborVertexType[4]||VertexType::Wall == neiborVertexType[4]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[5]||VertexType::Wall == neiborVertexType[5])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[6]||VertexType::Wall == neiborVertexType[4]) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IMJM_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points2");
                    }
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[5]||VertexType::Wall == neiborVertexType[5]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]||VertexType::Wall == neiborVertexType[1])&&
                    (VertexType::ImmersedSolid == neiborVertexType[2]||VertexType::Wall == neiborVertexType[2]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[3]||VertexType::Wall == neiborVertexType[3])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[13]||VertexType::Wall == neiborVertexType[13]) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IPKP_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points2");
                    }
                }

                if ((VertexType::ImmersedSolid == neiborVertexType[4]||VertexType::Wall == neiborVertexType[4]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]||VertexType::Wall == neiborVertexType[1])&&
                    (VertexType::ImmersedSolid == neiborVertexType[2]||VertexType::Wall == neiborVertexType[2]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[3]||VertexType::Wall == neiborVertexType[3])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[12]||VertexType::Wall == neiborVertexType[12]) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IPKM_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points2");
                    }
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[5]||VertexType::Wall == neiborVertexType[5]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[0]||VertexType::Wall == neiborVertexType[0])&&
                    (VertexType::ImmersedSolid == neiborVertexType[2]||VertexType::Wall == neiborVertexType[2]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[3]||VertexType::Wall == neiborVertexType[3])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[11]||VertexType::Wall == neiborVertexType[11]) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IMKP_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points2");
                    }
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[4]||VertexType::Wall == neiborVertexType[4]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[0]||VertexType::Wall == neiborVertexType[0])&&
                    (VertexType::ImmersedSolid == neiborVertexType[2]||VertexType::Wall == neiborVertexType[2]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[3]||VertexType::Wall == neiborVertexType[3])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[10]||VertexType::Wall == neiborVertexType[10]) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IMKM_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points2");
                    }
                }

                if ((VertexType::ImmersedSolid == neiborVertexType[5]||VertexType::Wall == neiborVertexType[5]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[3]||VertexType::Wall == neiborVertexType[3])&&
                    (VertexType::ImmersedSolid == neiborVertexType[0]||VertexType::Wall == neiborVertexType[0]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]||VertexType::Wall == neiborVertexType[1])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[17]||VertexType::Wall == neiborVertexType[17]) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_JPKP_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points2");
                    }
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[4]||VertexType::Wall == neiborVertexType[4]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[3]||VertexType::Wall == neiborVertexType[3])&&
                    (VertexType::ImmersedSolid == neiborVertexType[0]||VertexType::Wall == neiborVertexType[0]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]||VertexType::Wall == neiborVertexType[1])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[16]||VertexType::Wall == neiborVertexType[16]) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_JPKM_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points2");
                    }
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[5]||VertexType::Wall == neiborVertexType[5]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[2]||VertexType::Wall == neiborVertexType[2])&&
                    (VertexType::ImmersedSolid == neiborVertexType[0]||VertexType::Wall == neiborVertexType[0]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]||VertexType::Wall == neiborVertexType[1])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[15]||VertexType::Wall == neiborVertexType[15]) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_JMKP_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points2");
                    }
                }

                if ((VertexType::ImmersedSolid == neiborVertexType[4]||VertexType::Wall == neiborVertexType[4]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[2]||VertexType::Wall == neiborVertexType[2])&&
                    (VertexType::ImmersedSolid == neiborVertexType[0]||VertexType::Wall == neiborVertexType[0]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]||VertexType::Wall == neiborVertexType[1])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[14]||VertexType::Wall == neiborVertexType[14]) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_JMKM_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points2");
                    }
                }
                //std::cout<<geometryProperty(0,0,0)<<" ";
            }
            // Planlar corner
            if (5 == solidNeiborNumatCoord) {
                nodeType(0,0,0)=(int)VertexType::Wall;
                if (VertexType::ImmersedSolid != neiborVertexType[0]&&VertexType::Wall != neiborVertexType[0]) {
                    if(VertexType::ImmersedSolid != neiborVertexType[14]&&VertexType::Wall != neiborVertexType[14]){
                        geometryProperty(0,0,0) = (int)VG_JPKP_O_IP;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[15]&&VertexType::Wall != neiborVertexType[15]){
                        geometryProperty(0,0,0) = (int)VG_JPKM_O_IP;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[16]&&VertexType::Wall != neiborVertexType[16]){
                        geometryProperty(0,0,0) = (int)VG_JMKP_O_IP;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[17]&&VertexType::Wall != neiborVertexType[17]) {
                        geometryProperty(0,0,0) = (int)VG_JMKM_O_IP;
                    }   
                    else{
                        geometryProperty(0,0,0) = (int)VG_IP;
                    }
                }
                if (VertexType::ImmersedSolid != neiborVertexType[1]&&VertexType::Wall != neiborVertexType[1]) {
                    if(VertexType::ImmersedSolid != neiborVertexType[14]&&VertexType::Wall != neiborVertexType[14]){
                        geometryProperty(0,0,0) = (int)VG_JPKP_O_IM;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[15]&&VertexType::Wall != neiborVertexType[15]){
                        geometryProperty(0,0,0) = (int)VG_JPKM_O_IM;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[16]&&VertexType::Wall != neiborVertexType[16]){
                        geometryProperty(0,0,0) = (int)VG_JMKP_O_IM;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[17]&&VertexType::Wall != neiborVertexType[17]){
                        geometryProperty(0,0,0) = (int)VG_JMKM_O_IM;
                    } 
                    else{
                        geometryProperty(0,0,0) = (int)VG_IM;
                    }
                }
                if (VertexType::ImmersedSolid != neiborVertexType[2]&&VertexType::Wall != neiborVertexType[2]) {
                    if(VertexType::ImmersedSolid != neiborVertexType[10]&&VertexType::Wall != neiborVertexType[10]){
                        geometryProperty(0,0,0) = (int)VG_IPKP_O_JP;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[11]&&VertexType::Wall != neiborVertexType[11]){
                        geometryProperty(0,0,0) = (int)VG_IPKM_O_JP;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[12]&&VertexType::Wall != neiborVertexType[12]){
                        geometryProperty(0,0,0) = (int)VG_IMKP_O_JP;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[13]&&VertexType::Wall != neiborVertexType[13]) {
                        geometryProperty(0,0,0) = (int)VG_IMKM_O_JP;
                    }
                    else{
                        geometryProperty(0,0,0) = (int)VG_JP;
                    }
                }
                if (VertexType::ImmersedSolid != neiborVertexType[3]&&VertexType::Wall != neiborVertexType[3]) {
                    if(VertexType::ImmersedSolid != neiborVertexType[10]&&VertexType::Wall != neiborVertexType[10]){
                        geometryProperty(0,0,0) = (int)VG_IPKP_O_JP;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[11]&&VertexType::Wall != neiborVertexType[11]){
                        geometryProperty(0,0,0) = (int)VG_IPKM_O_JP;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[12]&&VertexType::Wall != neiborVertexType[12]){
                        geometryProperty(0,0,0) = (int)VG_IMKP_O_JP;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[13]&&VertexType::Wall != neiborVertexType[13]) {
                        geometryProperty(0,0,0) = (int)VG_IMKM_O_JP;
                    }
                    else{
                        geometryProperty(0,0,0) = (int)VG_JM;
                    }
                }
                if (VertexType::ImmersedSolid != neiborVertexType[4]&&VertexType::Wall != neiborVertexType[4]) {
                    if(VertexType::ImmersedSolid != neiborVertexType[6]&&VertexType::Wall != neiborVertexType[6]){
                        geometryProperty(0,0,0) = (int)VG_IPJP_O_KP;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[7]&&VertexType::Wall != neiborVertexType[7]){
                        geometryProperty(0,0,0) = (int)VG_IPJM_O_KP;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[8]&&VertexType::Wall != neiborVertexType[8]){
                        geometryProperty(0,0,0) = (int)VG_IMJP_O_KP;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[9]&&VertexType::Wall != neiborVertexType[9]) {
                        geometryProperty(0,0,0) = (int)VG_IMJM_O_KP;
                    }
                    else{
                        geometryProperty(0,0,0) = (int)VG_KP;
                    }
                }
                if (VertexType::ImmersedSolid != neiborVertexType[5]&&VertexType::Wall != neiborVertexType[5]) {
                    if(VertexType::ImmersedSolid != neiborVertexType[6]&&VertexType::Wall != neiborVertexType[6]){
                        geometryProperty(0,0,0) = (int)VG_IPJP_O_KM;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[7]&&VertexType::Wall != neiborVertexType[7]){
                        geometryProperty(0,0,0) = (int)VG_IPJM_O_KM;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[8]&&VertexType::Wall != neiborVertexType[8]){
                        geometryProperty(0,0,0) = (int)VG_IMJP_O_KM;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[9]&&VertexType::Wall != neiborVertexType[9]) {
                        geometryProperty(0,0,0) = (int)VG_IMJM_O_KM;
                    }
                    else{
                        geometryProperty(0,0,0) = (int)VG_KM;
                    }
                }
            }
            // Inner corner
            if (6 == solidNeiborNumatCoord) {
                nodeType(0,0,0)=(int)VertexType::Wall;
                if (3 == fluidNeiborNum||2==fluidNeiborNum) {
                    if (VertexType::ImmersedSolid != neiborVertexType[6]&&VertexType::Wall != neiborVertexType[6]) {
                        geometryProperty(0,0,0) = (int)VG_IPJP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[7]&&VertexType::Wall != neiborVertexType[7]) {
                        geometryProperty(0,0,0) = (int)VG_IPJM_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[8]&&VertexType::Wall != neiborVertexType[8]) {
                        geometryProperty(0,0,0) = (int)VG_IMJP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[9]&&VertexType::Wall != neiborVertexType[9]) {
                        geometryProperty(0,0,0) = (int)VG_IMJM_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[10]&&VertexType::Wall != neiborVertexType[10]) {
                        geometryProperty(0,0,0) = (int)VG_IPKP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[11]&&VertexType::Wall != neiborVertexType[11]) {
                        geometryProperty(0,0,0) = (int)VG_IPKM_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[12]&&VertexType::Wall != neiborVertexType[12]) {
                        geometryProperty(0,0,0) = (int)VG_IMKP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[13]&&VertexType::Wall != neiborVertexType[13]) {
                        geometryProperty(0,0,0) = (int)VG_IMKM_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[14]&&VertexType::Wall != neiborVertexType[14]) {
                        geometryProperty(0,0,0) = (int)VG_JPKP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[15]&&VertexType::Wall != neiborVertexType[15]) {
                        geometryProperty(0,0,0) = (int)VG_JPKM_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[16]&&VertexType::Wall != neiborVertexType[16]) {
                        geometryProperty(0,0,0) = (int)VG_JMKP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[17]&&VertexType::Wall != neiborVertexType[17]) {
                        geometryProperty(0,0,0) = (int)VG_JMKM_I;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[6]&&VertexType::Wall == neiborVertexType[6])&&(VertexType::ImmersedSolid != neiborVertexType[18]&&VertexType::Wall != neiborVertexType[18])) {
                        geometryProperty(0,0,0) = (int)VG_IPJP_O_KP;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[7]&&VertexType::Wall == neiborVertexType[7])&&(VertexType::ImmersedSolid != neiborVertexType[20]&&VertexType::Wall != neiborVertexType[20])) {
                        geometryProperty(0,0,0) = (int)VG_IPJM_O_KP;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[8]&&VertexType::Wall == neiborVertexType[8])&&(VertexType::ImmersedSolid != neiborVertexType[22]&&VertexType::Wall != neiborVertexType[22])) {
                        geometryProperty(0,0,0) = (int)VG_IMJP_O_KP;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[9]&&VertexType::Wall == neiborVertexType[9])&&(VertexType::ImmersedSolid != neiborVertexType[24]&&VertexType::Wall != neiborVertexType[24])) {
                        geometryProperty(0,0,0) = (int)VG_IMJM_O_KP;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[6]&&VertexType::Wall == neiborVertexType[6])&&(VertexType::ImmersedSolid != neiborVertexType[19]&&VertexType::Wall != neiborVertexType[19])) {
                        geometryProperty(0,0,0) = (int)VG_IPJP_O_KM;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[7]&&VertexType::Wall == neiborVertexType[7])&&(VertexType::ImmersedSolid != neiborVertexType[21]&&VertexType::Wall != neiborVertexType[21])) {
                        geometryProperty(0,0,0) = (int)VG_IPJM_O_KM;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[8]&&VertexType::Wall == neiborVertexType[8])&&(VertexType::ImmersedSolid != neiborVertexType[23]&&VertexType::Wall != neiborVertexType[23])) {
                        geometryProperty(0,0,0) = (int)VG_IMJP_O_KM;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[9]&&VertexType::Wall == neiborVertexType[9])&&(VertexType::ImmersedSolid != neiborVertexType[25]&&VertexType::Wall != neiborVertexType[25])) {
                        geometryProperty(0,0,0) = (int)VG_IMJM_O_KM;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[10]&&VertexType::Wall == neiborVertexType[10])&&(VertexType::ImmersedSolid != neiborVertexType[18]&&VertexType::Wall != neiborVertexType[18])) {
                        geometryProperty(0,0,0) = (int)VG_IPKP_O_JP;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[11]&&VertexType::Wall == neiborVertexType[11])&&(VertexType::ImmersedSolid != neiborVertexType[19]&&VertexType::Wall != neiborVertexType[19])) {
                        geometryProperty(0,0,0) = (int)VG_IPKM_O_JP;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[12]&&VertexType::Wall == neiborVertexType[12])&&(VertexType::ImmersedSolid != neiborVertexType[22]&&VertexType::Wall != neiborVertexType[22])) {
                        geometryProperty(0,0,0) = (int)VG_IMKP_O_JP;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[13]&&VertexType::Wall == neiborVertexType[13])&&(VertexType::ImmersedSolid != neiborVertexType[23]&&VertexType::Wall != neiborVertexType[23])) {
                        geometryProperty(0,0,0) = (int)VG_IMKM_O_JP;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[10]&&VertexType::Wall == neiborVertexType[10])&&(VertexType::ImmersedSolid != neiborVertexType[20]&&VertexType::Wall != neiborVertexType[20])) {
                        geometryProperty(0,0,0) = (int)VG_IPKP_O_JM;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[11]&&VertexType::Wall == neiborVertexType[11])&&(VertexType::ImmersedSolid != neiborVertexType[21]&&VertexType::Wall != neiborVertexType[21])) {
                        geometryProperty(0,0,0) = (int)VG_IPKM_O_JM;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[12]&&VertexType::Wall == neiborVertexType[12])&&(VertexType::ImmersedSolid != neiborVertexType[24]&&VertexType::Wall != neiborVertexType[24])) {
                        geometryProperty(0,0,0) = (int)VG_IMKP_O_JM;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[13]&&VertexType::Wall == neiborVertexType[13])&&(VertexType::ImmersedSolid != neiborVertexType[25]&&VertexType::Wall != neiborVertexType[25])) {
                        geometryProperty(0,0,0) = (int)VG_IMKM_O_JM;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[14]&&VertexType::Wall == neiborVertexType[14])&&(VertexType::ImmersedSolid != neiborVertexType[18]&&VertexType::Wall != neiborVertexType[18])) {
                        geometryProperty(0,0,0) = (int)VG_JPKP_O_IP;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[15]&&VertexType::Wall == neiborVertexType[15])&&(VertexType::ImmersedSolid != neiborVertexType[19]&&VertexType::Wall != neiborVertexType[19])) {
                        geometryProperty(0,0,0) = (int)VG_JPKM_O_IP;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[16]&&VertexType::Wall == neiborVertexType[16])&&(VertexType::ImmersedSolid != neiborVertexType[20]&&VertexType::Wall != neiborVertexType[20])) {
                        geometryProperty(0,0,0) = (int)VG_JMKP_O_IP;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[17]&&VertexType::Wall == neiborVertexType[17])&&(VertexType::ImmersedSolid != neiborVertexType[21]&&VertexType::Wall != neiborVertexType[21])) {
                        geometryProperty(0,0,0) = (int)VG_JMKM_O_IP;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[14]&&VertexType::Wall == neiborVertexType[14])&&(VertexType::ImmersedSolid != neiborVertexType[22]&&VertexType::Wall != neiborVertexType[22])) {
                        geometryProperty(0,0,0) = (int)VG_JPKP_O_IM;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[15]&&VertexType::Wall == neiborVertexType[15])&&(VertexType::ImmersedSolid != neiborVertexType[23]&&VertexType::Wall != neiborVertexType[23])) {
                        geometryProperty(0,0,0) = (int)VG_JPKM_O_IM;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[16]&&VertexType::Wall == neiborVertexType[16])&&(VertexType::ImmersedSolid != neiborVertexType[24]&&VertexType::Wall != neiborVertexType[24])) {
                        geometryProperty(0,0,0) = (int)VG_JMKP_O_IM;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[17]&&VertexType::Wall == neiborVertexType[17])&&(VertexType::ImmersedSolid != neiborVertexType[25]&&VertexType::Wall != neiborVertexType[25])) {
                        geometryProperty(0,0,0) = (int)VG_JMKM_O_IM;
                    }
                }
                else if (1 == fluidNeiborNum) {
                    if (VertexType::ImmersedSolid != neiborVertexType[18]&&VertexType::Wall != neiborVertexType[18]) {
                        geometryProperty(0,0,0) = (int)VG_IPJPKP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[19]&&VertexType::Wall != neiborVertexType[19]) {
                        geometryProperty(0,0,0) = (int)VG_IPJPKM_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[20]&&VertexType::Wall != neiborVertexType[20]) {
                        geometryProperty(0,0,0) = (int)VG_IPJMKP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[21]&&VertexType::Wall != neiborVertexType[21]) {
                        geometryProperty(0,0,0) = (int)VG_IPJMKM_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[22]&&VertexType::Wall != neiborVertexType[22]) {
                        geometryProperty(0,0,0) = (int)VG_IMJPKP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[23]&&VertexType::Wall != neiborVertexType[23]) {
                        geometryProperty(0,0,0) = (int)VG_IMJPKM_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[24]&&VertexType::Wall != neiborVertexType[24]) {
                        geometryProperty(0,0,0) = (int)VG_IMJMKP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[25]&&VertexType::Wall != neiborVertexType[25]) {
                        geometryProperty(0,0,0) = (int)VG_IMJMKM_I;
                    }
                }
                else {

                    //ops_printf("%s\n",
                    //           "There appears to be hanged solid points3");
                }
            }
        }
    }
#endif
}
#endif
#endif // CAVITY2D_KERNEL_INC