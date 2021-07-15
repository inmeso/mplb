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

static inline OPS_FUN_PREFIX Real CalcBGKGeqAD(const int l, const Real C, const Real u, const Real v,
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
    return WEIGHTS[l] * C * res;
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