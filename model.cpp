// Copyright 2017 the MPLB team. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

/*! @brief Define discrete model
 *  @author Jianping Meng
 *  @details Define the discrete velocity structure, the macroscopic variables,
 *   and necessary constants
 **/
#include "model.h"

std::string MODELNAME;
int NUMXI = 9;
int FEQORDER = 2;
int LATTDIM = 2;
Real CS = 1;
Real* XI;
Real* WEIGHTS;
int* OPP;
int NUMMACROVAR = 3;
int* VARIABLETYPE;
int* VARIABLECOMPINDEX;
int NUMCOMPONENTS = 1;
int* COMPOINDEX;
Real XIMAXVALUE;
int THERMALPROBLEM = 0;
std::vector<std::string> MACROVARNAME;

// find the particles with opposite directions, for the bounce-back boundary
// Brute-force method, could be slow for large lattice
void FindReverseXi() {
    for (int i = 0; i < NUMXI; i++) {
        for (int j = 0; j < NUMXI; j++) {
            bool isReverse{true};
            for (int k = 0; k < LATTDIM; k++) {
                Real sum{XI[i * LATTDIM + k] + XI[j * LATTDIM + k]};
                isReverse = isReverse && EssentiallyEqual(&sum, &ZERO, EPS);
            }
            if (isReverse) {
                OPP[i] = j;
                break;
            }
        }
    }
}

void SetupD2Q9Latt() {
    MODELNAME = "D2Q9";
    const int nc9 = 9;
    NUMXI = nc9;
    LATTDIM = 2;
    CS = sqrt(3);
    XI = new Real[NUMXI * LATTDIM];
    WEIGHTS = new Real[NUMXI];
    OPP = new int[NUMXI];
    Real t00 = 4.0 / 9.0, t01 = 1.0 / 9.0, t11 = 1.0 / 36.0;
    Real t[nc9] = {t00, t01, t01, t01, t01, t11, t11, t11, t11};
    int cxi[nc9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    int cyi[nc9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
    int op9[nc9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
    XIMAXVALUE = CS;
    for (int l = 0; l < NUMXI; l++) {
        XI[l * LATTDIM] = cxi[l];
        XI[l * LATTDIM + 1] = cyi[l];
        WEIGHTS[l] = t[l];
        OPP[l] = op9[l];
    }
}

void SetupD3Q19Latt() {
    MODELNAME = "D3Q19";
    const int nc19 = 19;
    NUMXI = nc19;
    LATTDIM = 3;
    CS = sqrt(3);
    XI = new Real[NUMXI * LATTDIM];
    WEIGHTS = new Real[NUMXI];
    OPP = new int[NUMXI];
    Real t000{1 / ((Real)3)};
    Real t001{1 / ((Real)18)};
    Real t011{1 / ((Real)36)};
    Real t[nc19] = {t000, t001, t001, t001, t001, t001, t001, t011, t011, t011,
                    t011, t011, t011, t011, t011, t011, t011, t011, t011};
    int cxi[nc19] = {0,  1, -1, 0, 0,  0, 0,  1, -1, 1,
                     -1, 0, 0,  1, -1, 1, -1, 0, 0};
    int cyi[nc19] = {0, 0, 0,  1,  -1, 0, 0, 1, -1, 0,
                     0, 1, -1, -1, 1,  0, 0, 1, -1};
    int czi[nc19] = {0,  0, 0,  0, 0, 1,  -1, 0,  0, 1,
                     -1, 1, -1, 0, 0, -1, 1,  -1, 1};
    XIMAXVALUE = CS;
    for (int l = 0; l < NUMXI; l++) {
        XI[l * LATTDIM] = cxi[l];
        XI[l * LATTDIM + 1] = cyi[l];
        XI[l * LATTDIM + 2] = czi[l];
        WEIGHTS[l] = t[l];
    }
    FindReverseXi();
}

void SetupD3Q15Latt() {
    MODELNAME = "D3Q15";
    const int nc15 = 15;
    NUMXI = nc15;
    LATTDIM = 3;
    CS = sqrt(3);
    XI = new Real[NUMXI * LATTDIM];
    WEIGHTS = new Real[NUMXI];
    OPP = new int[NUMXI];
    Real t000{2 / ((Real)9)};
    Real t001{1 / ((Real)9)};
    Real t111{1 / ((Real)72)};
    Real t[nc15] = {t000, t001, t001, t001, t001, t001, t001, t111,
                    t111, t111, t111, t111, t111, t111, t111};
    int cxi[nc15] = {0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, -1, 1};
    int cyi[nc15] = {0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1};
    int czi[nc15] = {0, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1, 1, -1};
    XIMAXVALUE = CS;
    for (int l = 0; l < NUMXI; l++) {
        XI[l * LATTDIM] = cxi[l];
        XI[l * LATTDIM + 1] = cyi[l];
        XI[l * LATTDIM + 2] = czi[l];
        WEIGHTS[l] = t[l];
    }
    FindReverseXi();
}

/*!
 * Example for 2D isothermal flows
 * Note:
 * This can become automatic by using an input file
 */
void SetupMacroVars() {
    // rho,u,v,w,T m, macroscopic  variables must be stored in a specific order
    NUMMACROVAR = 3;
    VARIABLETYPE = new int[NUMMACROVAR];
    VARIABLETYPE[0] = (int)Variable_Rho;
    VARIABLETYPE[1] = (int)Variable_U;
    VARIABLETYPE[2] = (int)Variable_V;
    //VARIABLETYPE[3] = (int)Variable_W;
    // VARIABLETYPE[3] = (int)Variable_T;
    // VARIABLETYPE[4] = (int)Variable_Qx;
    // VARIABLETYPE[5] = (int)Variable_Qy;
    VARIABLECOMPINDEX = new int[NUMMACROVAR];
    VARIABLECOMPINDEX[0] = 0;
    VARIABLECOMPINDEX[1] = 0;
    VARIABLECOMPINDEX[2] = 0;
    //VARIABLECOMPINDEX[3] = 0;
    // VARIABLECOMPINDEX[3] = 0;
    // VARIABLECOMPINDEX[4] = 0;
    // VARIABLECOMPINDEX[5] = 0;
    MACROVARNAME.reserve(NUMMACROVAR);
    MACROVARNAME.push_back("h");
    MACROVARNAME.push_back("u");
    MACROVARNAME.push_back("v");
    //MACROVARNAME.push_back("w");
    // MACROVARNAME.push_back("T");
    // MACROVARNAME.push_back("qx");
    // MACROVARNAME.push_back("qy");
}

void SetupD2Q16Latt() {
    MODELNAME = "D2Q16";
    const int nc16{16};
    NUMXI = nc16;
    LATTDIM = 2;
    CS = 1;
    XI = new Real[NUMXI * LATTDIM];
    WEIGHTS = new Real[NUMXI];
    OPP = new int[NUMXI];
    const int nc1d{4};
    const Real roots[nc1d] = {-2.3344142183389773, -0.7419637843027259,
                              0.7419637843027258, 2.3344142183389773};
    const Real coeff[nc1d] = {0.045875854768068526, 0.45412414523193156,
                              0.4541241452319317, 0.045875854768068526};
    XIMAXVALUE = 2.3344142183389773;
    // int op16[nc16] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
    int l{0};
    for (int i = 0; i < nc1d; i++) {
        for (int j = 0; j < nc1d; j++) {
            XI[l * LATTDIM] = roots[i];
            XI[l * LATTDIM + 1] = roots[j];
            WEIGHTS[l] = coeff[i] * coeff[j];
            // OPP[l] = op9[l];
            l++;
        }
    }
    FindReverseXi();
}

void SetupD2Q36Latt() {
    MODELNAME = "D2Q36";
    const int nc36{36};
    NUMXI = nc36;
    LATTDIM = 2;
    CS = 1;
    XI = new Real[NUMXI * LATTDIM];
    WEIGHTS = new Real[NUMXI];
    OPP = new int[NUMXI];
    const int nc1d{6};
    const Real roots[nc1d] = {-0.3324257433552119e1, -0.1889175877753711e1,
                              -0.6167065901925942,   0.6167065901925942,
                              0.1889175877753711e1,  0.3324257433552119e1};
    const Real coeff[nc1d] = {0.2555784402056229e-2, 0.8861574604191481e-1,
                              0.4088284695560294,    0.4088284695560294,
                              0.8861574604191481e-1, 0.2555784402056229e-2};
    XIMAXVALUE = 0.3324257433552119e1;
    // int op16[nc16] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
    int l{0};
    for (int i = 0; i < nc1d; i++) {
        for (int j = 0; j < nc1d; j++) {
            XI[l * LATTDIM] = roots[i];
            XI[l * LATTDIM + 1] = roots[j];
            WEIGHTS[l] = coeff[i] * coeff[j];
            // OPP[l] = op9[l];
            l++;
        }
    }
    // find the particles with opposite directions, for the bounce-back boundary
    // Brute-force method, could slow for large lattice
    for (int i = 0; i < NUMXI; i++) {
        for (int j = 0; j < NUMXI; j++) {
            if ((XI[i * LATTDIM] + XI[j * LATTDIM] == 0) &&
                (XI[i * LATTDIM + 1] + XI[j * LATTDIM + 1] == 0)) {
                OPP[i] = j;
            }
        }
    }
}

void DefineModelConstants() {
    // The variable here may cause some confusions for the Python translator
    // while the issues should be solvable.
    ops_decl_const("NUMCOMPONENTS", 1, "int", &NUMCOMPONENTS);
    ops_decl_const("COMPOINDEX", 2 * NUMCOMPONENTS, "int", COMPOINDEX);
    ops_decl_const("THERMALPROBLEM", 1, "int", &THERMALPROBLEM);
    // lattice structure
    ops_decl_const("NUMXI", 1, "int", &NUMXI);
    ops_decl_const("CS", 1, "double", &CS);
    ops_decl_const("LATTDIM", 1, "int", &LATTDIM);
    ops_decl_const("XI", NUMXI * LATTDIM, "double", XI);
    ops_decl_const("WEIGHTS", NUMXI, "double", WEIGHTS);
    ops_decl_const("OPP", NUMXI, "int", OPP);
    ops_decl_const("FEQORDER", 1, "int", &FEQORDER);
    // Macroscopic variables
    ops_decl_const("NUMMACROVAR", 1, "int", &NUMMACROVAR);
    ops_decl_const("VARIABLETYPE", NUMMACROVAR, "int", VARIABLETYPE);
    ops_decl_const("VARIABLECOMPINDEX", NUMMACROVAR, "int", VARIABLECOMPINDEX);
}

void SetupModel() {
    NUMCOMPONENTS = 1;
    FEQORDER = 4;
    THERMALPROBLEM = 0;
    COMPOINDEX = new int[2 * NUMCOMPONENTS];
    COMPOINDEX[0] = 0;
    COMPOINDEX[1] = 15;
    SetupD2Q16Latt();
    //SetupD3Q15Latt();
    SetupMacroVars();
    DefineModelConstants();
}
void DestroyModel() {
    delete[] VARIABLETYPE;
    delete[] VARIABLECOMPINDEX;
    delete[] COMPOINDEX;
    DestroyLatt();
}
void DestroyLatt() {
    delete[] XI;
    delete[] WEIGHTS;
    if (OPP != nullptr) delete[] OPP;
}
/*
 * Two-dimensional version
 */
Real CalcBGKFeq(const int l, const Real rho, const Real u, const Real v,
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
/*
 * three-dimensional version
 */
Real CalcBGKFeq(const int l, const Real rho, const Real u, const Real v,
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
// Equilibrium function for shallow water equations
Real CalcSWEFeq(const int l, const Real h, const Real u, const Real v,
                const int polyOrder) {
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

#include "model_kernel.h"
