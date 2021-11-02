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

//////////////////////////////////////////////////////////////////////////////////
//////////////////////                  2D                  //////////////////////
//////////////////////////////////////////////////////////////////////////////////

#ifndef CAVITY2D_KERNEL_INC
#define CAVITY2D_KERNEL_INC
#ifdef OPS_2D

/**
 * These functions can be adjusted to change the initial setup of the system
*/

//Set initial values of density and velocity
void KerSetInitialMacroVars(ACC<Real>& rho, ACC<Real>& u, ACC<Real>& v,
                            const ACC<Real>& coordinates, const int* idx) {
    rho(0, 0) = 1;
    u(0, 0) = 0;
    v(0, 0) = 0;
}

//Set initial values of the order parameter
void KerCalcConcentration(ACC<Real>& concentration, const ACC<Real>& xiIndex,
                          const Real* TIME, const ACC<Real>& U,
                          const ACC<Real>& V, const Real* Dt) {
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

//Set the solid boundaries in the simulation
void KerSetSolid(ACC<int>& nodeType, ACC<int>& geometry,
                 const ACC<Real>& coordinates, const int* xMax,
                 const int* yMax) {
    const Real X{coordinates(0,0,0)};
    const Real Y{coordinates(1,0,0)};
    int XMax{*xMax};
    int YMax{*yMax};

    if (X==XMax-1||X==0||Y==YMax-1||Y==0){
        nodeType(0,0)=(int)VertexType::ImmersedSolid;
        geometry(0,0)=(int)VG_ImmersedSolid;
    }
}



//////////////////////////////////////////////////////////////////////////////////




/**
 * These functions should be left as is but parameters can be changed if needed
*/



//Function to print a 1D macroscopic variable for debugging purposes
void KerPrintPhi(const ACC<Real>& Phi){
    std::cout << Phi(0,0) << " ";
}

//Labels normal direction for wall nodes
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
        //Count fluid nodes neighboring boundary in 8 nearest nodes
        for (int i = 0; i < 8; i++) {
            if (VertexType::ImmersedSolid != neiborVertexType[i]
                &&VertexType::Wall != neiborVertexType[i]) {
                fluidNeiborNum++;
            }
        }

        //Count solid nodes neighboring boundary in 4 nearest nodes
        //(along lattice axes)
        int solidNeiborNumatCoord{0};
        for (int i = 0; i < 4; i++) {
            if (VertexType::ImmersedSolid == neiborVertexType[i]
                ||VertexType::Wall == neiborVertexType[i]) {
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
                if ((VertexType::ImmersedSolid == neiborVertexType[2]
                     ||VertexType::Wall == neiborVertexType[2]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]
                     ||VertexType::Wall == neiborVertexType[1])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[6]
                        ||VertexType::Wall == neiborVertexType[6]) {
                        nodeType(0,0)=(int)VertexType::Wall;
                        geometryProperty(0, 0) = (int)VG_IPJM_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points");
                    }
                }

                if ((VertexType::ImmersedSolid == neiborVertexType[3]
                     ||VertexType::Wall == neiborVertexType[3]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]
                     ||VertexType::Wall == neiborVertexType[1])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[5]
                        ||VertexType::Wall == neiborVertexType[5]) {
                        nodeType(0,0)=(int)VertexType::Wall;
                        geometryProperty(0, 0) = (int)VG_IPJP_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points");
                    }
                }

                if ((VertexType::ImmersedSolid == neiborVertexType[3]
                     ||VertexType::Wall == neiborVertexType[3]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[0]
                     ||VertexType::Wall == neiborVertexType[0])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[7]
                        ||VertexType::Wall == neiborVertexType[7]) {
                        nodeType(0,0)=(int)VertexType::Wall;
                        geometryProperty(0, 0) = (int)VG_IMJP_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points");
                    }
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[2]
                     ||VertexType::Wall == neiborVertexType[2]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[0]
                     ||VertexType::Wall == neiborVertexType[0])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[4]
                        ||VertexType::Wall == neiborVertexType[4]) {
                        nodeType(0,0)=(int)VertexType::Wall;
                        geometryProperty(0, 0) = (int)VG_IMJM_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points");
                    }
                }
            }
            // Plane
            if (3 == solidNeiborNumatCoord) {
                nodeType(0,0)=(int)VertexType::Wall;
                if (VertexType::ImmersedSolid != neiborVertexType[0]
                    &&VertexType::Wall != neiborVertexType[0]) {
                    geometryProperty(0, 0) = (int)VG_IP;
                }
                if (VertexType::ImmersedSolid != neiborVertexType[1]
                    &&VertexType::Wall != neiborVertexType[1]) {
                    geometryProperty(0, 0) = (int)VG_IM;
                }
                if (VertexType::ImmersedSolid != neiborVertexType[2]
                    &&VertexType::Wall != neiborVertexType[2]) {
                    geometryProperty(0, 0) = (int)VG_JP;
                }
                if (VertexType::ImmersedSolid != neiborVertexType[3]
                    &&VertexType::Wall != neiborVertexType[3]) {
                    geometryProperty(0, 0) = (int)VG_JM;
                }
            }
            // Inner corner
            if (4 == solidNeiborNumatCoord) {
                nodeType(0,0)=(int)VertexType::Wall;
                if (1 == fluidNeiborNum) {
                    if (VertexType::ImmersedSolid != neiborVertexType[4]
                        &&VertexType::Wall != neiborVertexType[4]) {
                        geometryProperty(0, 0) = (int)VG_IPJP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[5]
                        &&VertexType::Wall != neiborVertexType[5]) {
                        geometryProperty(0, 0) = (int)VG_IMJM_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[6]
                        &&VertexType::Wall != neiborVertexType[6]) {
                        geometryProperty(0, 0) = (int)VG_IMJP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[7]
                        &&VertexType::Wall != neiborVertexType[7]) {
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

//////////////////////////////////////////////////////////////////////////////////
//////////////////////                  3D                  //////////////////////
//////////////////////////////////////////////////////////////////////////////////

#ifndef CAVITY3D_KERNEL_INC
#define CAVITY3D_KERNEL_INC
#ifdef OPS_3D

/**
 * These functions can be adjusted to change the initial setup of the system
*/

//Set initial values of density and velocity
void KerSetInitialMacroVars3D(ACC<Real>& rho, ACC<Real>& u, ACC<Real>& v,
                              ACC<Real>& w, const ACC<Real>& coordinates,
                              const int* idx) {
    rho(0, 0, 0) = 1;
    u(0, 0, 0) = 0;
    v(0, 0, 0) = 0;
    w(0, 0, 0) = 0;
}

//Set initial values of the order parameter
void KerCalcConcentration3D(ACC<Real>& concentration, const ACC<Real>& xiIndex,
                            const Real* TIME, const ACC<Real>& U,
                            const ACC<Real>& V, const Real* Dt) {
    Real x0[3] = {40.0,40.0,40.0};
    Real xiIndexVar1 = {xiIndex(0,0,0,0)};
    Real xiIndexVar2 = {xiIndex(1,0,0,0)};
    Real xiIndexVar3 = {xiIndex(2,0,0,0)};
    Real x=(xiIndexVar1-x0[0]);
    Real y=(xiIndexVar2-x0[1]);
    Real z=(xiIndexVar3-x0[2]);
    Real R=sqrt(x*x+y*y+z*z);
    Real r=25.0;
    if (R<=r){
        concentration(0,0,0) = 1;
    }
    else{
        concentration(0,0,0) = -1;
    }
}

//Set the solid boundaries in the simulation
void KerSetSolid3D(ACC<int>& nodeType, ACC<int>& geometry,
                   const ACC<Real>& coordinates, const int* xMax, const int* yMax,
                   const int* zMax) {
    const Real X{coordinates(0,0,0,0)};
    const Real Y{coordinates(1,0,0,0)};
    const Real Z{coordinates(2,0,0,0)};
    Real x0[3] = {40.0,40.0,40.0};
    const Real Xc{X-x0[0]};
    const Real Yc{Y-x0[1]};
    const Real Zc{Z-x0[2]};
    int XMax{*xMax};
    int YMax{*yMax};
    int ZMax{*zMax};
    Real R=sqrt(Xc*Xc+Yc*Yc+Zc*Zc);
    Real r=15.0;
    //if (R<=r){
    //    nodeType(0,0,0)=(int)VertexType::ImmersedSolid;
    //    geometry(0,0,0)=(int)VG_ImmersedSolid;
    //}
    if (Y>=YMax-2||Y<=2){
        nodeType(0,0,0)=(int)VertexType::ImmersedSolid;
        geometry(0,0,0)=(int)VG_ImmersedSolid;
    }
}

void KerUpdateMacroBodyForce(ACC<Real>& force, const ACC<Real>& coordinates,
                            const int* idx) {
    force(0, 0, 0, 0) = 0.001;
    force(1, 0, 0, 0) = 0;
    force(2, 0, 0, 0) = 0;
}

//////////////////////////////////////////////////////////////////////////////////




/**
 * These functions should be left as is but parameters can be changed if needed
*/


/*
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

//Function to print a 1D macroscopic variable for debugging purposes
void KerPrintPhi3D(const ACC<Real>& Phi, const ACC<Real>& coordinates){
    std::cout<<Phi(0,0,0)<<" ";
}

//Labels normal direction for wall nodes
void KerSetEmbeddedBodyGeometry3D(ACC<int>& geometryProperty,
                                ACC<int>& nodeType) {
#ifdef OPS_3D
    VertexType vt = (VertexType)nodeType(0, 0, 0);
    if (VertexType::ImmersedSolid == vt) {
        VertexType neiborVertexType[26];

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
        //Count fluid nodes neighboring boundary in 26 nearest nodes
        for (int i = 0; i < 26; i++) {
            if (VertexType::ImmersedSolid != neiborVertexType[i]
                &&VertexType::Wall != neiborVertexType[i]) {
                fluidNeiborNum++;
            }
        }
        int solidNeiborNumatCoord{0};
        for (int i = 0; i < 6; i++) {
            if (VertexType::ImmersedSolid == neiborVertexType[i]
                ||VertexType::Wall == neiborVertexType[i]) {
                solidNeiborNumatCoord++;
            }
        }
        //std::cout<<solidNeiborNumatCoord<<" ";
        if (fluidNeiborNum > 0) {
            
            //Outer corner
            if (3 == solidNeiborNumatCoord) {
                if ((VertexType::ImmersedSolid == neiborVertexType[5]
                     ||VertexType::Wall == neiborVertexType[5]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[3]
                     ||VertexType::Wall == neiborVertexType[3]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]
                     ||VertexType::Wall == neiborVertexType[1])) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IPJPKP_O;
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[4]
                     ||VertexType::Wall == neiborVertexType[4]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[3]
                     ||VertexType::Wall == neiborVertexType[3]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]
                     ||VertexType::Wall == neiborVertexType[1])) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IPJPKM_O;
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[5]
                     ||VertexType::Wall == neiborVertexType[5]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[2]
                     ||VertexType::Wall == neiborVertexType[2]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]
                     ||VertexType::Wall == neiborVertexType[1])) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IPJMKP_O;
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[4]
                     ||VertexType::Wall == neiborVertexType[4]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[2]
                     ||VertexType::Wall == neiborVertexType[2]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]
                     ||VertexType::Wall == neiborVertexType[1])) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IPJMKM_O;
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[5]
                     ||VertexType::Wall == neiborVertexType[5]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[3]
                     ||VertexType::Wall == neiborVertexType[3]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[0]
                     ||VertexType::Wall == neiborVertexType[0])) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IMJPKP_O;
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[4]
                     ||VertexType::Wall == neiborVertexType[4]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[3]
                     ||VertexType::Wall == neiborVertexType[3]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[0]
                     ||VertexType::Wall == neiborVertexType[0])) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IMJPKM_O;
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[5]
                     ||VertexType::Wall == neiborVertexType[5]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[2]
                     ||VertexType::Wall == neiborVertexType[2]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[0]
                     ||VertexType::Wall == neiborVertexType[0])) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IMJMKP_O;
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[4]
                     ||VertexType::Wall == neiborVertexType[4]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[2]
                     ||VertexType::Wall == neiborVertexType[2]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[0]
                     ||VertexType::Wall == neiborVertexType[0])) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IMJMKM_O;
                }
            }
            //Outer edge
            if (4 == solidNeiborNumatCoord) {
                if ((VertexType::ImmersedSolid == neiborVertexType[0] &&
                     VertexType::ImmersedSolid == neiborVertexType[1]) ||
                    (VertexType::ImmersedSolid == neiborVertexType[0] &&
                     VertexType::ImmersedSolid == neiborVertexType[1])) {
                    ops_printf("%s\n",
                               "There appears to be hanged solid points,i.e., "
                               "the solid body may be too thin");
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[3]
                     ||VertexType::Wall == neiborVertexType[3]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]
                     ||VertexType::Wall == neiborVertexType[1]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[4]
                     ||VertexType::Wall == neiborVertexType[4]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[5]
                     ||VertexType::Wall == neiborVertexType[5])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[9]
                        ||VertexType::Wall == neiborVertexType[9]) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IPJP_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points2");
                    }
                }

                if ((VertexType::ImmersedSolid == neiborVertexType[2]
                     ||VertexType::Wall == neiborVertexType[2]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]
                     ||VertexType::Wall == neiborVertexType[1]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[4]
                     ||VertexType::Wall == neiborVertexType[4]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[5]
                     ||VertexType::Wall == neiborVertexType[5])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[8]
                        ||VertexType::Wall == neiborVertexType[8]) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IPJM_O;
                        
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points2");
                    }
                }

                if ((VertexType::ImmersedSolid == neiborVertexType[3]
                     ||VertexType::Wall == neiborVertexType[3]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[0]
                     ||VertexType::Wall == neiborVertexType[0])&&
                    (VertexType::ImmersedSolid == neiborVertexType[4]
                     ||VertexType::Wall == neiborVertexType[4]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[5]
                     ||VertexType::Wall == neiborVertexType[5])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[7]
                        ||VertexType::Wall == neiborVertexType[7]) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IMJP_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points2");
                    }
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[2]
                     ||VertexType::Wall == neiborVertexType[2]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[0]
                     ||VertexType::Wall == neiborVertexType[0])&&
                    (VertexType::ImmersedSolid == neiborVertexType[4]
                     ||VertexType::Wall == neiborVertexType[4]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[5]
                     ||VertexType::Wall == neiborVertexType[5])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[6]
                        ||VertexType::Wall == neiborVertexType[4]) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IMJM_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points2");
                    }
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[5]
                     ||VertexType::Wall == neiborVertexType[5]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]
                     ||VertexType::Wall == neiborVertexType[1])&&
                    (VertexType::ImmersedSolid == neiborVertexType[2]
                     ||VertexType::Wall == neiborVertexType[2]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[3]
                     ||VertexType::Wall == neiborVertexType[3])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[13]
                        ||VertexType::Wall == neiborVertexType[13]) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IPKP_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points2");
                    }
                }

                if ((VertexType::ImmersedSolid == neiborVertexType[4]
                     ||VertexType::Wall == neiborVertexType[4]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]
                     ||VertexType::Wall == neiborVertexType[1])&&
                    (VertexType::ImmersedSolid == neiborVertexType[2]
                     ||VertexType::Wall == neiborVertexType[2]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[3]
                     ||VertexType::Wall == neiborVertexType[3])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[12]
                        ||VertexType::Wall == neiborVertexType[12]) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IPKM_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points2");
                    }
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[5]
                     ||VertexType::Wall == neiborVertexType[5]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[0]
                     ||VertexType::Wall == neiborVertexType[0])&&
                    (VertexType::ImmersedSolid == neiborVertexType[2]
                     ||VertexType::Wall == neiborVertexType[2]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[3]
                     ||VertexType::Wall == neiborVertexType[3])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[11]
                        ||VertexType::Wall == neiborVertexType[11]) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IMKP_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points2");
                    }
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[4]
                     ||VertexType::Wall == neiborVertexType[4]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[0]
                     ||VertexType::Wall == neiborVertexType[0])&&
                    (VertexType::ImmersedSolid == neiborVertexType[2]
                     ||VertexType::Wall == neiborVertexType[2]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[3]
                     ||VertexType::Wall == neiborVertexType[3])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[10]
                        ||VertexType::Wall == neiborVertexType[10]) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_IMKM_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points2");
                    }
                }

                if ((VertexType::ImmersedSolid == neiborVertexType[5]
                     ||VertexType::Wall == neiborVertexType[5]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[3]
                     ||VertexType::Wall == neiborVertexType[3])&&
                    (VertexType::ImmersedSolid == neiborVertexType[0]
                     ||VertexType::Wall == neiborVertexType[0]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]
                     ||VertexType::Wall == neiborVertexType[1])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[17]
                        ||VertexType::Wall == neiborVertexType[17]) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_JPKP_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points2");
                    }
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[4]
                     ||VertexType::Wall == neiborVertexType[4]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[3]
                     ||VertexType::Wall == neiborVertexType[3])&&
                    (VertexType::ImmersedSolid == neiborVertexType[0]
                     ||VertexType::Wall == neiborVertexType[0]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]
                     ||VertexType::Wall == neiborVertexType[1])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[16]
                        ||VertexType::Wall == neiborVertexType[16]) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_JPKM_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points2");
                    }
                }
                if ((VertexType::ImmersedSolid == neiborVertexType[5]
                     ||VertexType::Wall == neiborVertexType[5]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[2]
                     ||VertexType::Wall == neiborVertexType[2])&&
                    (VertexType::ImmersedSolid == neiborVertexType[0]
                     ||VertexType::Wall == neiborVertexType[0]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]
                     ||VertexType::Wall == neiborVertexType[1])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[15]
                        ||VertexType::Wall == neiborVertexType[15]) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_JMKP_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points2");
                    }
                }

                if ((VertexType::ImmersedSolid == neiborVertexType[4]
                     ||VertexType::Wall == neiborVertexType[4]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[2]
                     ||VertexType::Wall == neiborVertexType[2])&&
                    (VertexType::ImmersedSolid == neiborVertexType[0]
                     ||VertexType::Wall == neiborVertexType[0]) &&
                    (VertexType::ImmersedSolid == neiborVertexType[1]
                     ||VertexType::Wall == neiborVertexType[1])) {
                    if (VertexType::ImmersedSolid == neiborVertexType[14]
                        ||VertexType::Wall == neiborVertexType[14]) {
                        nodeType(0,0,0)=(int)VertexType::Wall;
                        geometryProperty(0,0,0) = (int)VG_JMKM_O;
                    } else {
                        ops_printf("%s\n",
                                   "There appears to be hanged solid points2");
                    }
                }
            }
            //Plane (also inner edge intersection with plane)
            if (5 == solidNeiborNumatCoord) {
                nodeType(0,0,0)=(int)VertexType::Wall;
                if (VertexType::ImmersedSolid != neiborVertexType[0]
                    &&VertexType::Wall != neiborVertexType[0]) {
                    if(VertexType::ImmersedSolid != neiborVertexType[14]
                       &&VertexType::Wall != neiborVertexType[14]){
                        geometryProperty(0,0,0) = (int)VG_JPKP_I_IP;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[15]
                       &&VertexType::Wall != neiborVertexType[15]){
                        geometryProperty(0,0,0) = (int)VG_JPKM_I_IP;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[16]
                       &&VertexType::Wall != neiborVertexType[16]){
                        geometryProperty(0,0,0) = (int)VG_JMKP_I_IP;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[17]
                       &&VertexType::Wall != neiborVertexType[17]) {
                        geometryProperty(0,0,0) = (int)VG_JMKM_I_IP;
                    }   
                    else{
                        geometryProperty(0,0,0) = (int)VG_IP;
                    }
                }
                if (VertexType::ImmersedSolid != neiborVertexType[1]
                    &&VertexType::Wall != neiborVertexType[1]) {
                    if(VertexType::ImmersedSolid != neiborVertexType[14]
                       &&VertexType::Wall != neiborVertexType[14]){
                        geometryProperty(0,0,0) = (int)VG_JPKP_I_IM;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[15]
                       &&VertexType::Wall != neiborVertexType[15]){
                        geometryProperty(0,0,0) = (int)VG_JPKM_I_IM;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[16]
                       &&VertexType::Wall != neiborVertexType[16]){
                        geometryProperty(0,0,0) = (int)VG_JMKP_I_IM;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[17]
                       &&VertexType::Wall != neiborVertexType[17]){
                        geometryProperty(0,0,0) = (int)VG_JMKM_I_IM;
                    } 
                    else{
                        geometryProperty(0,0,0) = (int)VG_IM;
                    }
                }
                if (VertexType::ImmersedSolid != neiborVertexType[2]
                    &&VertexType::Wall != neiborVertexType[2]) {
                    if(VertexType::ImmersedSolid != neiborVertexType[10]
                       &&VertexType::Wall != neiborVertexType[10]){
                        geometryProperty(0,0,0) = (int)VG_IPKP_I_JP;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[11]
                       &&VertexType::Wall != neiborVertexType[11]){
                        geometryProperty(0,0,0) = (int)VG_IPKM_I_JP;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[12]
                       &&VertexType::Wall != neiborVertexType[12]){
                        geometryProperty(0,0,0) = (int)VG_IMKP_I_JP;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[13]
                       &&VertexType::Wall != neiborVertexType[13]) {
                        geometryProperty(0,0,0) = (int)VG_IMKM_I_JP;
                    }
                    else{
                        geometryProperty(0,0,0) = (int)VG_JP;
                    }
                }
                if (VertexType::ImmersedSolid != neiborVertexType[3]
                    &&VertexType::Wall != neiborVertexType[3]) {
                    if(VertexType::ImmersedSolid != neiborVertexType[10]
                       &&VertexType::Wall != neiborVertexType[10]){
                        geometryProperty(0,0,0) = (int)VG_IPKP_I_JP;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[11]
                       &&VertexType::Wall != neiborVertexType[11]){
                        geometryProperty(0,0,0) = (int)VG_IPKM_I_JP;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[12]
                       &&VertexType::Wall != neiborVertexType[12]){
                        geometryProperty(0,0,0) = (int)VG_IMKP_I_JP;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[13]
                       &&VertexType::Wall != neiborVertexType[13]) {
                        geometryProperty(0,0,0) = (int)VG_IMKM_I_JP;
                    }
                    else{
                        geometryProperty(0,0,0) = (int)VG_JM;
                    }
                }
                if (VertexType::ImmersedSolid != neiborVertexType[4]
                    &&VertexType::Wall != neiborVertexType[4]) {
                    if(VertexType::ImmersedSolid != neiborVertexType[6]
                       &&VertexType::Wall != neiborVertexType[6]){
                        geometryProperty(0,0,0) = (int)VG_IPJP_I_KP;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[7]
                       &&VertexType::Wall != neiborVertexType[7]){
                        geometryProperty(0,0,0) = (int)VG_IPJM_I_KP;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[8]
                       &&VertexType::Wall != neiborVertexType[8]){
                        geometryProperty(0,0,0) = (int)VG_IMJP_I_KP;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[9]
                       &&VertexType::Wall != neiborVertexType[9]) {
                        geometryProperty(0,0,0) = (int)VG_IMJM_I_KP;
                    }
                    else{
                        geometryProperty(0,0,0) = (int)VG_KP;
                    }
                }
                if (VertexType::ImmersedSolid != neiborVertexType[5]
                    &&VertexType::Wall != neiborVertexType[5]) {
                    if(VertexType::ImmersedSolid != neiborVertexType[6]
                       &&VertexType::Wall != neiborVertexType[6]){
                        geometryProperty(0,0,0) = (int)VG_IPJP_I_KM;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[7]
                       &&VertexType::Wall != neiborVertexType[7]){
                        geometryProperty(0,0,0) = (int)VG_IPJM_I_KM;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[8]
                       &&VertexType::Wall != neiborVertexType[8]){
                        geometryProperty(0,0,0) = (int)VG_IMJP_I_KM;
                    }
                    if(VertexType::ImmersedSolid != neiborVertexType[9]
                       &&VertexType::Wall != neiborVertexType[9]) {
                        geometryProperty(0,0,0) = (int)VG_IMJM_I_KM;
                    }
                    else{
                        geometryProperty(0,0,0) = (int)VG_KM;
                    }
                }
            }
            //Inner edge
            if (6 == solidNeiborNumatCoord) {
                nodeType(0,0,0)=(int)VertexType::Wall;
                if (3 == fluidNeiborNum||2==fluidNeiborNum) {
                    if (VertexType::ImmersedSolid != neiborVertexType[6]
                        &&VertexType::Wall != neiborVertexType[6]) {
                        geometryProperty(0,0,0) = (int)VG_IPJP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[7]
                        &&VertexType::Wall != neiborVertexType[7]) {
                        geometryProperty(0,0,0) = (int)VG_IPJM_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[8]
                        &&VertexType::Wall != neiborVertexType[8]) {
                        geometryProperty(0,0,0) = (int)VG_IMJP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[9]
                        &&VertexType::Wall != neiborVertexType[9]) {
                        geometryProperty(0,0,0) = (int)VG_IMJM_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[10]
                        &&VertexType::Wall != neiborVertexType[10]) {
                        geometryProperty(0,0,0) = (int)VG_IPKP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[11]
                        &&VertexType::Wall != neiborVertexType[11]) {
                        geometryProperty(0,0,0) = (int)VG_IPKM_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[12]
                        &&VertexType::Wall != neiborVertexType[12]) {
                        geometryProperty(0,0,0) = (int)VG_IMKP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[13]
                        &&VertexType::Wall != neiborVertexType[13]) {
                        geometryProperty(0,0,0) = (int)VG_IMKM_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[14]
                        &&VertexType::Wall != neiborVertexType[14]) {
                        geometryProperty(0,0,0) = (int)VG_JPKP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[15]
                        &&VertexType::Wall != neiborVertexType[15]) {
                        geometryProperty(0,0,0) = (int)VG_JPKM_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[16]
                        &&VertexType::Wall != neiborVertexType[16]) {
                        geometryProperty(0,0,0) = (int)VG_JMKP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[17]
                        &&VertexType::Wall != neiborVertexType[17]) {
                        geometryProperty(0,0,0) = (int)VG_JMKM_I;
                    }
            //Outer edge intersection with plane
                    if ((VertexType::ImmersedSolid == neiborVertexType[6]
                        &&VertexType::Wall == neiborVertexType[6])
                        &&(VertexType::ImmersedSolid != neiborVertexType[18]
                        &&VertexType::Wall != neiborVertexType[18])) {
                        geometryProperty(0,0,0) = (int)VG_IPJP_O_KP;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[7]
                        &&VertexType::Wall == neiborVertexType[7])
                        &&(VertexType::ImmersedSolid != neiborVertexType[20]
                        &&VertexType::Wall != neiborVertexType[20])) {
                        geometryProperty(0,0,0) = (int)VG_IPJM_O_KP;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[8]
                        &&VertexType::Wall == neiborVertexType[8])
                        &&(VertexType::ImmersedSolid != neiborVertexType[22]
                        &&VertexType::Wall != neiborVertexType[22])) {
                        geometryProperty(0,0,0) = (int)VG_IMJP_O_KP;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[9]
                        &&VertexType::Wall == neiborVertexType[9])
                        &&(VertexType::ImmersedSolid != neiborVertexType[24]
                        &&VertexType::Wall != neiborVertexType[24])) {
                        geometryProperty(0,0,0) = (int)VG_IMJM_O_KP;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[6]
                        &&VertexType::Wall == neiborVertexType[6])
                        &&(VertexType::ImmersedSolid != neiborVertexType[19]
                        &&VertexType::Wall != neiborVertexType[19])) {
                        geometryProperty(0,0,0) = (int)VG_IPJP_O_KM;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[7]
                        &&VertexType::Wall == neiborVertexType[7])
                        &&(VertexType::ImmersedSolid != neiborVertexType[21]
                        &&VertexType::Wall != neiborVertexType[21])) {
                        geometryProperty(0,0,0) = (int)VG_IPJM_O_KM;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[8]
                        &&VertexType::Wall == neiborVertexType[8])
                        &&(VertexType::ImmersedSolid != neiborVertexType[23]
                        &&VertexType::Wall != neiborVertexType[23])) {
                        geometryProperty(0,0,0) = (int)VG_IMJP_O_KM;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[9]
                        &&VertexType::Wall == neiborVertexType[9])
                        &&(VertexType::ImmersedSolid != neiborVertexType[25]
                        &&VertexType::Wall != neiborVertexType[25])) {
                        geometryProperty(0,0,0) = (int)VG_IMJM_O_KM;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[10]
                        &&VertexType::Wall == neiborVertexType[10])
                        &&(VertexType::ImmersedSolid != neiborVertexType[18]
                        &&VertexType::Wall != neiborVertexType[18])) {
                        geometryProperty(0,0,0) = (int)VG_IPKP_O_JP;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[11]
                        &&VertexType::Wall == neiborVertexType[11])
                        &&(VertexType::ImmersedSolid != neiborVertexType[19]
                        &&VertexType::Wall != neiborVertexType[19])) {
                        geometryProperty(0,0,0) = (int)VG_IPKM_O_JP;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[12]
                        &&VertexType::Wall == neiborVertexType[12])
                        &&(VertexType::ImmersedSolid != neiborVertexType[22]
                        &&VertexType::Wall != neiborVertexType[22])) {
                        geometryProperty(0,0,0) = (int)VG_IMKP_O_JP;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[13]
                        &&VertexType::Wall == neiborVertexType[13])
                        &&(VertexType::ImmersedSolid != neiborVertexType[23]
                        &&VertexType::Wall != neiborVertexType[23])) {
                        geometryProperty(0,0,0) = (int)VG_IMKM_O_JP;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[10]
                        &&VertexType::Wall == neiborVertexType[10])
                        &&(VertexType::ImmersedSolid != neiborVertexType[20]
                        &&VertexType::Wall != neiborVertexType[20])) {
                        geometryProperty(0,0,0) = (int)VG_IPKP_O_JM;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[11]
                        &&VertexType::Wall == neiborVertexType[11])
                        &&(VertexType::ImmersedSolid != neiborVertexType[21]
                        &&VertexType::Wall != neiborVertexType[21])) {
                        geometryProperty(0,0,0) = (int)VG_IPKM_O_JM;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[12]
                        &&VertexType::Wall == neiborVertexType[12])
                        &&(VertexType::ImmersedSolid != neiborVertexType[24]
                        &&VertexType::Wall != neiborVertexType[24])) {
                        geometryProperty(0,0,0) = (int)VG_IMKP_O_JM;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[13]
                        &&VertexType::Wall == neiborVertexType[13])
                        &&(VertexType::ImmersedSolid != neiborVertexType[25]
                        &&VertexType::Wall != neiborVertexType[25])) {
                        geometryProperty(0,0,0) = (int)VG_IMKM_O_JM;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[14]
                        &&VertexType::Wall == neiborVertexType[14])
                        &&(VertexType::ImmersedSolid != neiborVertexType[18]
                        &&VertexType::Wall != neiborVertexType[18])) {
                        geometryProperty(0,0,0) = (int)VG_JPKP_O_IP;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[15]
                        &&VertexType::Wall == neiborVertexType[15])
                        &&(VertexType::ImmersedSolid != neiborVertexType[19]
                        &&VertexType::Wall != neiborVertexType[19])) {
                        geometryProperty(0,0,0) = (int)VG_JPKM_O_IP;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[16]
                        &&VertexType::Wall == neiborVertexType[16])
                        &&(VertexType::ImmersedSolid != neiborVertexType[20]
                        &&VertexType::Wall != neiborVertexType[20])) {
                        geometryProperty(0,0,0) = (int)VG_JMKP_O_IP;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[17]
                        &&VertexType::Wall == neiborVertexType[17])
                        &&(VertexType::ImmersedSolid != neiborVertexType[21]
                        &&VertexType::Wall != neiborVertexType[21])) {
                        geometryProperty(0,0,0) = (int)VG_JMKM_O_IP;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[14]
                    &&VertexType::Wall == neiborVertexType[14])
                    &&(VertexType::ImmersedSolid != neiborVertexType[22]
                    &&VertexType::Wall != neiborVertexType[22])) {
                        geometryProperty(0,0,0) = (int)VG_JPKP_O_IM;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[15]
                    &&VertexType::Wall == neiborVertexType[15])
                    &&(VertexType::ImmersedSolid != neiborVertexType[23]
                    &&VertexType::Wall != neiborVertexType[23])) {
                        geometryProperty(0,0,0) = (int)VG_JPKM_O_IM;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[16]
                    &&VertexType::Wall == neiborVertexType[16])
                    &&(VertexType::ImmersedSolid != neiborVertexType[24]
                    &&VertexType::Wall != neiborVertexType[24])) {
                        geometryProperty(0,0,0) = (int)VG_JMKP_O_IM;
                    }
                    if ((VertexType::ImmersedSolid == neiborVertexType[17]
                    &&VertexType::Wall == neiborVertexType[17])
                    &&(VertexType::ImmersedSolid != neiborVertexType[25]
                    &&VertexType::Wall != neiborVertexType[25])) {
                        geometryProperty(0,0,0) = (int)VG_JMKM_O_IM;
                    }
                }
            //Inner corner
                else if (1 == fluidNeiborNum) {
                    if (VertexType::ImmersedSolid != neiborVertexType[18]
                    &&VertexType::Wall != neiborVertexType[18]) {
                        geometryProperty(0,0,0) = (int)VG_IPJPKP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[19]
                    &&VertexType::Wall != neiborVertexType[19]) {
                        geometryProperty(0,0,0) = (int)VG_IPJPKM_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[20]
                    &&VertexType::Wall != neiborVertexType[20]) {
                        geometryProperty(0,0,0) = (int)VG_IPJMKP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[21]
                    &&VertexType::Wall != neiborVertexType[21]) {
                        geometryProperty(0,0,0) = (int)VG_IPJMKM_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[22]
                    &&VertexType::Wall != neiborVertexType[22]) {
                        geometryProperty(0,0,0) = (int)VG_IMJPKP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[23]
                    &&VertexType::Wall != neiborVertexType[23]) {
                        geometryProperty(0,0,0) = (int)VG_IMJPKM_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[24]
                    &&VertexType::Wall != neiborVertexType[24]) {
                        geometryProperty(0,0,0) = (int)VG_IMJMKP_I;
                    }
                    if (VertexType::ImmersedSolid != neiborVertexType[25]
                    &&VertexType::Wall != neiborVertexType[25]) {
                        geometryProperty(0,0,0) = (int)VG_IMJMKM_I;
                    }
                }
                else {

                    ops_printf("%s\n",
                               "There appears to be hanged solid points3");
                }
            }
        }
    }
#endif
}
#endif
#endif // CAVITY2D_KERNEL_INC