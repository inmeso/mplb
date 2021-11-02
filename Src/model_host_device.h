#ifndef MODEL_HOST_DEVICE_H
#define MODEL_HOST_DEVICE_H
#ifndef OPS_FUN_PREFIX
#define OPS_FUN_PREFIX
#endif
enum VariableTypes {
    Variable_Rho = 0,
    Variable_U = 1,
    Variable_V = 2,
    Variable_W = 3,
    Variable_T = 4,
    Variable_Qx = 5,
    Variable_Qy = 6,
    Variable_Qz = 7,
    // This is for the force correction needed by using He1998 scheme.
    Variable_U_Force = 8,
    Variable_V_Force = 9,
    Variable_W_Force = 10,
};

/*
 * Calculate the first-order force term
 * Author: Jianping Meng  22-Feb-2019
 */
static inline OPS_FUN_PREFIX Real CalcBodyForce(const int xiIndex, const Real rho,
                   const Real* acceleration) {
    Real cf{0};
    for (int i = 0; i < LATTDIM; i++) {
        cf += CS * XI[xiIndex * LATTDIM + i] * acceleration[i];
    }
    return WEIGHTS[xiIndex] * rho * cf;
}

static inline OPS_FUN_PREFIX Real CalcGuoForce(const int xiIndex, const Real rho,
                   const Real* acceleration, const Real u, const Real v, const Real* dt) {
    Real cs2{1/sqrt(3)};
    Real cs4{cs2*cs2};
    Real cf{0};
    Real uForce=u-((*dt) * acceleration[0] / (2.0*rho));
    Real vForce=v-((*dt) * acceleration[1] / (2.0*rho));
    //Real uForce=u;
    //Real vForce=v;
    Real VEL[]{uForce,vForce};
    //std::cout<<u<<" "<<v<<" "<<*VEL<<"\n";
    for (int i = 0; i < LATTDIM; i++) {
        cf += (3.0*((XI[xiIndex * LATTDIM + i])-VEL[i])+9.0*(XI[xiIndex * LATTDIM + i]*VEL[i]*XI[xiIndex * LATTDIM + i])) * acceleration[i];
    }
    return WEIGHTS[xiIndex] * cf;
}

static inline OPS_FUN_PREFIX Real CalcGuoForce3D(const int xiIndex, const Real rho,
                   const Real* acceleration, const Real u, const Real v, const Real w, const Real* dt) {
    const Real cs2={1/3.0};
    const Real cs4={cs2*cs2};
    Real cf{0};
    Real uForce=u-((*dt) * acceleration[0] / (2.0*rho));
    Real vForce=v-((*dt) * acceleration[1] / (2.0*rho));
    Real wForce=w-((*dt) * acceleration[2] / (2.0*rho));
    Real VEL[]{uForce,vForce,wForce};
    
    //std::cout<<u<<" "<<v<<" "<<*VEL<<"\n";
    for (int i = 0; i < LATTDIM; i++) {
        cf += (3.0*((XI[xiIndex * LATTDIM + i])-VEL[i])+9.0*(XI[xiIndex * LATTDIM + i]*VEL[i]*XI[xiIndex * LATTDIM + i])) * acceleration[i];
	//std::cout<<LATTDIM<<" "<<acceleration[i]<<" "<<rho<<" "<<VEL[i]<<" "<<cf<<" "<<cs2<<" "<<cs4<<" "<<XI[xiIndex * LATTDIM + i]<<"  ";
    }
    return WEIGHTS[xiIndex] * cf;
}


static inline OPS_FUN_PREFIX Real CalcBGKFeq(const int l, const Real rho, const Real u, const Real v,
                const Real T, const int polyOrder) {
    Real cu{(CS * XI[l * LATTDIM] * u + CS * XI[l * LATTDIM + 1] * v)};
    Real c2{(CS * XI[l * LATTDIM] * CS * XI[l * LATTDIM] +
             CS * XI[l * LATTDIM + 1] * CS * XI[l * LATTDIM + 1])};
    Real cu2{cu * cu};
    Real u2{u * u + v * v};
    Real res = 1.0 + cu + 0.5 * (cu2 - u2 + (T - 1.0) * (c2 - LATTDIM));
    if ((polyOrder) >= 3) {
        res = res +
              cu * (cu2 - 3.0 * u2 + 3.0 * (T - 1.0) * (c2 - LATTDIM - 2.0)) /
                  6.0;
    }
    if ((polyOrder) >= 4) {
        res =
            res + (cu2 * cu2 - 6.0 * cu2 * u2 + 3.0 * u2 * u2) / 24.0 +
            (T - 1.0) * ((c2 - (LATTDIM + 2)) * (cu2 - u2) - 2.0 * cu2) / 4.0 +
            (T - 1.0) * (T - 1.0) *
                (c2 * c2 - 2.0 * (LATTDIM + 2) * c2 + LATTDIM * (LATTDIM + 2)) /
                8.0;
    }
    return WEIGHTS[l] * rho * res;
}

static inline OPS_FUN_PREFIX Real CalcBGKFeqAD(const int l, const Real C, const Real U, const Real V,
                const Real T, const int polyOrder) {
    Real u=U;
    Real v=V;
    Real cSound{1/sqrt(3)};
    Real cu{(XI[l * LATTDIM] * u + XI[l * LATTDIM + 1] * v)/(cSound*cSound)};
    //Real c2{(cSound * XI[l * LATTDIM] * cSound * XI[l * LATTDIM] +
    //         cSound * XI[l * LATTDIM + 1] * cSound * XI[l * LATTDIM + 1])};
    Real cu2{(cu * cu)};
    Real u2{(u * u + v * v)/(cSound*cSound)};

    Real res = 1.0 + cu + 0.5 * (cu2 - u2);
    return WEIGHTS[l] * C * res;
}

static inline OPS_FUN_PREFIX Real CalcBGKFeqFE(const int l, const Real C, const Real uForce, const Real vForce, const Real *acceleration, const Real *dt,
                const Real T, const int polyOrder) {
    Real u=uForce-((*dt) * acceleration[0] / (2.0*C));
    Real v=vForce-((*dt) * acceleration[1] / (2.0*C));
    Real cSound{1/sqrt(3)};
    Real cu{(XI[l * LATTDIM] * u + XI[l * LATTDIM + 1] * v)/(cSound*cSound)};
    //Real c2{(cSound * XI[l * LATTDIM] * cSound * XI[l * LATTDIM] +
    //         cSound * XI[l * LATTDIM + 1] * cSound * XI[l * LATTDIM + 1])};
    Real cu2{(cu * cu)};
    Real u2{(u * u + v * v)/(cSound*cSound)};

    Real res = 1.0 + cu + 0.5 * (cu2 - u2);
    return WEIGHTS[l] * C * res;
}


static inline OPS_FUN_PREFIX Real CalcBGKGeqFE(const int l, const Real C, const Real u, const Real v,
                const Real mu, const Real* dt, const int polyOrder) {
    Real cSound{1/sqrt(3)};
    Real cu{(XI[l * LATTDIM] * u + XI[l * LATTDIM + 1] * v)/(cSound*cSound)};
    //Real c2{(cSound * XI[l * LATTDIM] * cSound * XI[l * LATTDIM] +
    //         cSound * XI[l * LATTDIM + 1] * cSound * XI[l * LATTDIM + 1])};
    Real cu2{(cu * cu)};
    Real u2{(u * u + v * v)/(cSound*cSound)};
    Real tau=1.0;
    
    Real css=1.0/sqrt(3);
    Real taul=1.0;
    Real taug=1.0;
    Real Dt{*dt};
    Real vl=css*css*(taul-Dt/2.0);
    Real vg=css*css*(taug-Dt/2.0);
    Real V=vg+(C+1.0)/2.0*(vl-vg);
    Real W=2.0/(6.0*V+1.0);
    Real res = cu + 0.5 * (cu2 - u2);
    Real R=1.0/(W-*dt/2.0);
    /*
    Real res = cu + 0.5 * (cu2 - u2);
    Real R=1.0/(tau-*dt/2.0);
    */
    //std::cout<<mu<<" "<<res<<"   ";
    return WEIGHTS[l] * (C * res+R*mu/(cSound*cSound));
}




static inline OPS_FUN_PREFIX Real CalcBGKFeqAD3D(const int l, const Real C, const Real uForce, const Real vForce, const Real wForce,
                const Real T, const int polyOrder) {
    Real u=uForce;
    Real v=vForce;
    Real w=wForce;
    Real cSound{1/sqrt(3)};
    Real cu{(XI[l * LATTDIM] * u + XI[l * LATTDIM + 1] * v + XI[l * LATTDIM + 2] * w)/(cSound*cSound)};
    //Real c2{(cSound * XI[l * LATTDIM] * cSound * XI[l * LATTDIM] +
    //         cSound * XI[l * LATTDIM + 1] * cSound * XI[l * LATTDIM + 1])};
    Real cu2{(cu * cu)};
    Real u2{(u * u + v * v + w * w)/(cSound*cSound)};
    //std::cout<<u<<" "<<cu<<" "<<u2<<"   ";
    Real res = 1.0 + cu + 0.5 * (cu2 - u2);
    //std::cout<<res<<" "<<C<<"   ";
    return WEIGHTS[l] * C * res;
}

static inline OPS_FUN_PREFIX Real CalcBGKFeqFE3D(const int l, const Real C, const Real uForce, const Real vForce, const Real wForce, const Real *acceleration, const Real *dt,
                const Real T, const int polyOrder) {
    Real u=uForce-((*dt) * acceleration[0] / 2);
    Real v=vForce-((*dt) * acceleration[1] / 2);
    Real w=wForce-((*dt) * acceleration[2] / 2);
    Real cSound{1/sqrt(3)};
    Real cu{(XI[l * LATTDIM] * u + XI[l * LATTDIM + 1] * v + XI[l * LATTDIM + 2] * w)/(cSound*cSound)};
    //Real c2{(cSound * XI[l * LATTDIM] * cSound * XI[l * LATTDIM] +
    //         cSound * XI[l * LATTDIM + 1] * cSound * XI[l * LATTDIM + 1])};
    Real cu2{(cu * cu)};
    Real u2{(u * u + v * v + w * w)/(cSound*cSound)};

    Real res = 1.0 + cu + 0.5 * (cu2 - u2);
    //std::cout<<res<<" "<<C<<"   ";
    return WEIGHTS[l] * C * res;
}


static inline OPS_FUN_PREFIX Real CalcBGKGeqFE3D(const int l, const Real C, const Real u, const Real v,
                const Real w, const Real mu, const Real* dt, const int polyOrder) {
    Real cSound{1/sqrt(3)};
    Real cu{(XI[l * LATTDIM] * u + XI[l * LATTDIM + 1] * v + XI[l * LATTDIM + 2] * w)/(cSound*cSound)};
    //Real c2{(cSound * XI[l * LATTDIM] * cSound * XI[l * LATTDIM] +
    //         cSound * XI[l * LATTDIM + 1] * cSound * XI[l * LATTDIM + 1])};
    Real cu2{(cu * cu)};
    Real u2{(u * u + v * v + w * w)/(cSound*cSound)};
    Real tau=1;
    Real res = cu + 0.5 * (cu2 - u2);
    Real R=1/(tau-*dt/2.0);
    //std::cout<<mu<<" "<<res<<"   ";
    return WEIGHTS[l] * (C * res+R*mu/(cSound*cSound));
}

static inline OPS_FUN_PREFIX Real CalcBGKFeq(const int l, const Real rho, const Real u, const Real v,
                const Real w, const Real T, const int polyOrder) {
    Real cu{(CS * XI[l * LATTDIM] * u + CS * XI[l * LATTDIM + 1] * v +
             CS * XI[l * LATTDIM + 2] * w)};
    Real c2{(CS * XI[l * LATTDIM] * CS * XI[l * LATTDIM] +
             CS * XI[l * LATTDIM + 1] * CS * XI[l * LATTDIM + 1] +
             CS * XI[l * LATTDIM + 2] * CS * XI[l * LATTDIM + 2])};
    Real cu2{cu * cu};
    Real u2{u * u + v * v + w * w};
    Real res = 1.0 + cu + 0.5 * (cu2 - u2 + (T - 1.0) * (c2 - LATTDIM));
    if ((polyOrder) >= 3) {
        res = res +
              cu * (cu2 - 3.0 * u2 + 3.0 * (T - 1.0) * (c2 - LATTDIM - 2.0)) /
                  6.0;
    }
    if ((polyOrder) >= 4) {
        res =
            res + (cu2 * cu2 - 6.0 * cu2 * u2 + 3.0 * u2 * u2) / 24.0 +
            (T - 1.0) * ((c2 - (LATTDIM + 2)) * (cu2 - u2) - 2.0 * cu2) / 4.0 +
            (T - 1.0) * (T - 1.0) *
                (c2 * c2 - 2.0 * (LATTDIM + 2) * c2 + LATTDIM * (LATTDIM + 2)) /
                8.0;
    }
    return WEIGHTS[l] * rho * res;
}

static inline OPS_FUN_PREFIX Real CalcSWEFeq(const int l, const Real h, const Real u, const Real v,
                const int polyOrder) {
    // Implementing the model derived in Please refer to Meng, Gu Emerson, Peng
    // and Zhang, IJMPC 2018(29):1850080

    Real cu{(CS * XI[l * LATTDIM] * u + CS * XI[l * LATTDIM + 1] * v)};
    Real c2{(CS * XI[l * LATTDIM] * CS * XI[l * LATTDIM] +
             CS * XI[l * LATTDIM + 1] * CS * XI[l * LATTDIM + 1])};
    Real cu2{cu * cu};
    Real u2{u * u + v * v};
    Real res = 1.0 + cu + 0.5 * (cu2 - u2 + (h - 1.0) * (c2 - LATTDIM));
    if (polyOrder >= 3) {
        res = res +
              cu * (cu2 - 3.0 * u2 + 3.0 * (h - 1.0) * (c2 - LATTDIM - 2.0)) /
                  6.0;
    }
    if (polyOrder >= 4) {
        res =
            res + (cu2 * cu2 - 6.0 * cu2 * u2 + 3.0 * u2 * u2) / 24.0 +
            (h - 1.0) * ((c2 - (LATTDIM + 2)) * (cu2 - u2) - 2.0 * cu2) / 4.0 +
            (h - 1.0) * (h - 1.0) *
                (c2 * c2 - 2.0 * (LATTDIM + 2) * c2 + LATTDIM * (LATTDIM + 2)) /
                8.0;
    }
    return WEIGHTS[l] * h * res;
}

#endif //MODEL_HOST_DEVICE_H