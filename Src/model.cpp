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

/*! @brief Define discrete model
 *  @author Jianping Meng
 *  @details Define the discrete velocity structure, the macroscopic variables,
 *   and necessary constants
 **/
#include "model.h"
#include <map>
int NUMXI{9};
int FEQORDER{2};
int LATTDIM{2};
Real CS{1};
Real* XI{nullptr};
Real* WEIGHTS{nullptr};
int* OPP{nullptr};
int NUMMACROVAR{3};
int* VARIABLETYPE{nullptr};
int* VARIABLECOMPINDEX{nullptr};
int NUMCOMPONENTS{1};
int* COMPOINDEX{nullptr};
Real XIMAXVALUE{1};
std::list<std::pair<SizeType, CollisionType>> COLLISIONTERMS;
std::list<std::pair<SizeType, BodyForceType>> FORCETERMS;
std::list<std::pair<SizeType, InitialType>> INITIALTERMS;

const std::list<std::pair<SizeType,CollisionType >> & CollisionTerms(){
    return COLLISIONTERMS;
}

const std::list<std::pair<SizeType,BodyForceType>> & BodyForceTerms(){
    return FORCETERMS;
}

const std::list<std::pair<SizeType,InitialType>> & InitialTerms(){
    return INITIALTERMS;
}

int* VARIABLECOMPPOS{nullptr};
/*!
 *Name of all macroscopic variables
 */
std::vector<std::string> MACROVARNAME;
/*!
 *Name of all Lattices.
 */
std::vector<std::string> LATTICENAME;

struct lattice {
    int lattDim;
    int length;
    Real cs;
};
// Giving the parameters of commonly used lattices.
lattice d2q9{2, 9, sqrt(3)};
lattice d3q19{3, 19, sqrt(3)};
lattice d3q15{3, 15, sqrt(3)};
lattice d2q16{2, 16, 1};
lattice d2q36{2, 36, 1};

std::map<std::string, lattice> latticeSet{
    {"d2q9", d2q9}, {"d3q19", d3q19}, {"d3q15", d3q15}, {"d2q36", d2q36}};

// Find particles with opposite directions, for bounce-back type boundary
// Brute-force method, could be slow for large lattice
void FindReverseXi(const int startPos, const int latticeSize) {
    for (int i = 0; i < latticeSize; i++) {
        for (int j = 0; j < latticeSize; j++) {
            bool isReverse{true};
            for (int k = 0; k < LATTDIM; k++) {
                Real sum{XI[(startPos + i) * LATTDIM + k] +
                         XI[(startPos + j) * LATTDIM + k]};
                isReverse = isReverse && EssentiallyEqual(&sum, &ZERO, EPS);
            }
            if (isReverse) {
                OPP[startPos + i] = startPos + j;
                break;
            }
        }
    }
}

void SetupD2Q9Latt(const int startPos) {
    const int nc9{9};
    Real t00 = 4.0 / 9.0, t01 = 1.0 / 9.0, t11 = 1.0 / 36.0;
    Real t[nc9] = {t00, t01, t01, t01, t01, t11, t11, t11, t11};
    int cxi[nc9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    int cyi[nc9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
    int op9[nc9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
    for (int l = 0; l < nc9; l++) {
        XI[(startPos + l) * LATTDIM] = cxi[l];
        XI[(startPos + l) * LATTDIM + 1] = cyi[l];
        WEIGHTS[startPos + l] = t[l];
        OPP[startPos + l] = op9[l];
    }
}

void SetupD3Q19Latt(const int startPos) {
    const int nc19 = 19;
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
    for (int l = 0; l < nc19; l++) {
        XI[(startPos + l) * LATTDIM] = cxi[l];
        XI[(startPos + l) * LATTDIM + 1] = cyi[l];
        XI[(startPos + l) * LATTDIM + 2] = czi[l];
        WEIGHTS[startPos + l] = t[l];
    }
    FindReverseXi(startPos, nc19);
}

void SetupD3Q15Latt(const int startPos) {
    const int nc15 = 15;
    Real t000{2 / ((Real)9)};
    Real t001{1 / ((Real)9)};
    Real t111{1 / ((Real)72)};
    Real t[nc15] = {t000, t001, t001, t001, t001, t001, t001, t111,
                    t111, t111, t111, t111, t111, t111, t111};
    int cxi[nc15] = {0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, -1, 1};
    int cyi[nc15] = {0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1};
    int czi[nc15] = {0, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1, 1, -1};
    for (int l = 0; l < nc15; l++) {
        XI[(startPos + l) * LATTDIM] = cxi[l];
        XI[(startPos + l) * LATTDIM + 1] = cyi[l];
        XI[(startPos + l) * LATTDIM + 2] = czi[l];
        WEIGHTS[startPos + l] = t[l];
    }
    FindReverseXi(startPos, nc15);
}

void SetupD2Q16Latt(const int startPos) {
    //Gauss-Hermite quadrature from the fourth order polynomial
    const int nc16{16};
    const int nc1d{4};
    const Real roots[nc1d] = {-2.3344142183389773, -0.7419637843027259,
                              0.7419637843027258, 2.3344142183389773};
    const Real coeff[nc1d] = {0.045875854768068526, 0.45412414523193156,
                              0.4541241452319317, 0.045875854768068526};
    int l{0};
    for (int i = 0; i < nc1d; i++) {
        for (int j = 0; j < nc1d; j++) {
            XI[(startPos + l) * LATTDIM] = roots[i];
            XI[(startPos + l) * LATTDIM + 1] = roots[j];
            WEIGHTS[startPos + l] = coeff[i] * coeff[j];
            l++;
        }
    }
    FindReverseXi(startPos, nc16);
}

void SetupD2Q36Latt(const int startPos) {
    //Gauss-Hermite quadrature from the sixth order polynomial
    const int nc36{36};
    const int nc1d{6};
    const Real roots[nc1d] = {-0.3324257433552119e1, -0.1889175877753711e1,
                              -0.6167065901925942,   0.6167065901925942,
                              0.1889175877753711e1,  0.3324257433552119e1};
    const Real coeff[nc1d] = {0.2555784402056229e-2, 0.8861574604191481e-1,
                              0.4088284695560294,    0.4088284695560294,
                              0.8861574604191481e-1, 0.2555784402056229e-2};
    int l{0};
    for (int i = 0; i < nc1d; i++) {
        for (int j = 0; j < nc1d; j++) {
            XI[(startPos + l) * LATTDIM] = roots[i];
            XI[(startPos + l) * LATTDIM + 1] = roots[j];
            WEIGHTS[startPos + l] = coeff[i] * coeff[j];
            l++;
        }
    }
    FindReverseXi(startPos, nc36);
}


void AllocateComponentIndex(const int compoNum) {
    if (compoNum == NUMCOMPONENTS) {
        if (nullptr == COMPOINDEX) {
            COMPOINDEX = new int[2 * compoNum];
        }
        if (nullptr == VARIABLECOMPPOS) {
            VARIABLECOMPPOS = new int[2 * compoNum];
        }
    }
}

void AllocateXi(const int length) {
    if (length == NUMXI) {
        if (nullptr == XI) {
            XI = new Real[length * LATTDIM];
        }
        if (nullptr == WEIGHTS) {
            WEIGHTS = new Real[length];
        }
        if (nullptr == OPP) {
            OPP = new int[length];
        }
    }
}

void AllocateMacroVarProperty(const int macroVarNum) {
    if (macroVarNum == NUMMACROVAR) {
        if (nullptr == VARIABLETYPE) {
            VARIABLETYPE = new int[NUMMACROVAR];
        } else {
            ops_printf("%s\n", "Warning! VARIABLETYPE has been allocated!");
        }
        if (nullptr == VARIABLECOMPINDEX) {
            VARIABLECOMPINDEX = new int[NUMMACROVAR];
        } else {
            ops_printf("%s\n",
                       "Warning! VARIABLECOMPINDEX has been allocated!");
        }
    } else {
        ops_printf("%s\n",
                   "Error! The macroVarNum must be equal to NUMMACROVAR");
        assert(macroVarNum == NUMMACROVAR);
    }
}

void DefineComponents(const std::vector<std::string>& compoNames,
                      const std::vector<SizeType>& compoId,
                      const std::vector<std::string>& lattNames,
                      const SizeType timeStep) {
    NUMCOMPONENTS = compoNames.size();
    if (NUMCOMPONENTS > 0) {
        AllocateComponentIndex(NUMCOMPONENTS);
        ops_printf("There are %i components defined.\n", NUMCOMPONENTS);
    } else {
        ops_printf(
            "Error! There muse be at least one component but we get:%i\n",
            NUMCOMPONENTS);
        assert(NUMCOMPONENTS > 0);
    }
    bool isLattDimSame{true};
    bool isCsSame{true};
    int posCompo{0};
    int totalSize{0};
    int latticeDimension{latticeSet[lattNames[0]].lattDim};
    Real currentCs{latticeSet[lattNames[0]].cs};
    for (int idx = 0; idx < NUMCOMPONENTS; idx++) {
        if (latticeSet.find(lattNames[idx]) != latticeSet.end()) {
            lattice currentLattice{latticeSet[lattNames[idx]]};
            COMPOINDEX[posCompo] = totalSize;
            COMPOINDEX[posCompo + 1] = totalSize + currentLattice.length - 1;
            totalSize += currentLattice.length;
            posCompo += 2;
            isLattDimSame =
                isLattDimSame && (latticeDimension == currentLattice.lattDim);
            isCsSame = isCsSame && (currentCs == currentLattice.cs);
        } else {
            ops_printf("Error! There is no predefined lattice:%s\n",
                       lattNames[idx].c_str());
            assert(latticeSet.find(lattNames[idx]) != latticeSet.end());
        }
    }
    if (!isLattDimSame) {
        ops_printf("%s\n", "Error! The lattice dimension is inconsistent!");
        assert(isLattDimSame);
    }
    if (!isCsSame) {
        ops_printf("%s\n", "Warning: The lattice sound speed is inconsistent!");
    }
    if (isLattDimSame && isCsSame) {
        NUMXI = totalSize;
        CS = currentCs;
        LATTDIM = latticeDimension;
        AllocateXi(totalSize);
        SetLatticeName(lattNames);
        int startPos{0};
        for (int idx = 0; idx < NUMCOMPONENTS; idx++) {
            if ("d3q15" == lattNames[idx]) {
                SetupD3Q15Latt(startPos);
            }
            if ("d3q19" == lattNames[idx]) {
                SetupD3Q19Latt(startPos);
            }
            if ("d2q9" == lattNames[idx]) {
                SetupD2Q9Latt(startPos);
            }
            startPos += latticeSet[lattNames[idx]].length;
            ops_printf("The %s lattice is employed for Component %i.\n",
                       lattNames[idx].c_str(), idx);
        }
        Real maxValue{0};
        for (int l = 0; l < totalSize * LATTDIM; l++) {
            maxValue = maxValue > XI[l] ? maxValue : XI[l];
        }
        XIMAXVALUE = CS * maxValue;
    }
    ops_decl_const("NUMCOMPONENTS", 1, "int", &NUMCOMPONENTS);
    ops_decl_const("COMPOINDEX", 2 * NUMCOMPONENTS, "int", COMPOINDEX);
    ops_decl_const("NUMXI", 1, "int", &NUMXI);
    ops_decl_const("CS", 1, "double", &CS);
    ops_decl_const("LATTDIM", 1, "int", &LATTDIM);
    ops_decl_const("XI", NUMXI * LATTDIM, "double", XI);
    ops_decl_const("WEIGHTS", NUMXI, "double", WEIGHTS);
    ops_decl_const("OPP", NUMXI, "int", OPP);
    g_f.SetDataDim(NUMXI);
    if (timeStep == 0){
    g_f.CreateFieldFromScratch();}
    else {
        g_f.CreateFieldFromCheckPoint(timeStep);
    }
    g_fStage.SetDataDim(NUMXI);
    g_fStage.CreateFieldFromScratch();
}

void DefineMacroVars(std::vector<VariableTypes> types,
                     std::vector<std::string> names,
                     std::vector<SizeType> varId, std::vector<SizeType> compoId,
                     const SizeType timeStep) {
    // It seems varId is not necessary at this moment
    NUMMACROVAR = names.size();
    MACROVARNAME = names;
    if (NUMMACROVAR > 0) {
        AllocateMacroVarProperty(NUMMACROVAR);
        ops_printf("There are %i macroscopic variables defined.\n",
                   NUMMACROVAR);
    } else {
        ops_printf("Warning! There seems no macroscopic variables defined!\n");
    }

    for (int idx = 0; idx < NUMMACROVAR; idx++) {
        VARIABLETYPE[idx] = (int)types[idx];
        VARIABLECOMPINDEX[idx] = (int)compoId[idx];
        ops_printf("The macroscopic variable %s defined for Component %i.\n",
                   names[idx].c_str(), compoId[idx]);

        if (nullptr != VARIABLECOMPPOS) {
            int startPos{0};
            for (int idx = 0; idx < NUMCOMPONENTS; idx++) {
                if (startPos >= NUMMACROVAR) {
                    ops_printf(
                        "Error! There are something wrong in specifying IDs of "
                        "macroscopic variables!");
                    assert(startPos >= NUMMACROVAR);
                }
                VARIABLECOMPPOS[2 * idx] = startPos;
                while (idx == compoId.at(startPos)) {
                    startPos++;
                    if (startPos >= NUMMACROVAR) {
                        break;
                    }
                }
                VARIABLECOMPPOS[2 * idx + 1] = startPos - 1;
            }
        } else {
            ops_printf(
                "Error! It appears that the DefineComponents routine has not "
                "been "
                "called!\n");
            assert(nullptr != VARIABLECOMPPOS);
        }
    }
    ops_decl_const("NUMMACROVAR", 1, "int", &NUMMACROVAR);
    ops_decl_const("VARIABLETYPE", NUMMACROVAR, "int", VARIABLETYPE);
    ops_decl_const("VARIABLECOMPINDEX", NUMMACROVAR, "int", VARIABLECOMPINDEX);
    ops_decl_const("VARIABLECOMPPOS", NUMCOMPONENTS, "int", VARIABLECOMPPOS);
    g_MacroVars.SetDataDim(NUMMACROVAR);
    if (timeStep == 0) {
        g_MacroVars.CreateFieldFromScratch();
    } else {
        g_MacroVars.CreateFieldFromCheckPoint(timeStep);
    }
    if (!IsTransient()) {
        g_MacroVarsCopy.SetDataDim(NUMMACROVAR);
        g_MacroVarsCopy.CreateFieldFromScratch();
        g_ResidualErrorHandle = new ops_reduction[MacroVarsNum()];
        g_ResidualError = new Real[2 * MacroVarsNum()];
        for (int localIdx = 0; localIdx < MacroVarsNum(); localIdx++) {
            g_ResidualErrorHandle[localIdx] = ops_decl_reduction_handle(
                sizeof(Real), "double", MacroVarName()[localIdx].c_str());
        }
    }
}

void DefineCollision(std::vector<CollisionType> types,
                     std::vector<SizeType> compoId) {
    const std::size_t typeSize{types.size()};
    const std::size_t compoSize{compoId.size()};

    if (typeSize != compoSize) {
        ops_printf(
            "Error! There are %i collision types defined for  %i "
            "components\n",
            typeSize, compoSize);
        assert(typeSize == compoSize);
    }

    if (compoSize > NUMCOMPONENTS){
        ops_printf(
            "Error! There are %i collision types defined but only %i "
            "components\n",
            compoSize, NUMCOMPONENTS);
        assert(compoSize < NUMCOMPONENTS);
    }

    for (int idx = 0; idx < typeSize; idx++) {
        ops_printf(
            "The equilibrium function type %i is chosen for "
            "Component "
            "%i\n",
            types.at(idx), compoId.at(idx));
        std::pair<int, CollisionType> pair(compoId.at(idx), types.at(idx));
        COLLISIONTERMS.push_back(pair);
    }
     if (compoSize < NUMCOMPONENTS) {
        ops_printf(
            "Warning! Kernel functions are required for components without "
            "pre-defined collision terms!\n");
    }
}

void DefineBodyForce(std::vector<BodyForceType> types,
                     std::vector<SizeType> compoId) {
    const SizeType typeSize{types.size()};
    const SizeType compoSize{compoId.size()};

    if (typeSize != compoSize) {
        ops_printf(
            "Error! There are %i force types defined for  %i "
            "components\n",
            typeSize, compoSize);
        assert(typeSize == compoSize);
    }

    if (compoSize > NUMCOMPONENTS) {
        ops_printf(
            "Error! There are %i forces types defined but only %i "
            "components\n",
            compoSize, NUMCOMPONENTS);
        assert(compoSize < NUMCOMPONENTS);
    }

    for (int idx = 0; idx < typeSize; idx++) {
        ops_printf(
            "The body force function type %i is chosen for Component "
            "%i\n",
            types.at(idx), compoId.at(idx));
        std::pair<SizeType, BodyForceType> pair(compoId.at(idx), types.at(idx));
        FORCETERMS.push_back(pair);
    }
    if (compoSize < NUMCOMPONENTS) {
        ops_printf(
            "Warning! Kernel functions are required for components without "
            "pre-defined body force terms!\n");
    }
    g_MacroBodyforce.SetDataDim(SPACEDIM * NUMCOMPONENTS);
    g_MacroBodyforce.CreateFieldFromScratch();
}
void DefineInitialCondition(std::vector<InitialType> types,
                            std::vector<SizeType> compoId) {
    const SizeType typeSize{types.size()};
    const SizeType compoSize{compoId.size()};

    if (typeSize != compoSize) {
        ops_printf(
            "Error! There are %i initial types defined for  %i "
            "components\n",
            typeSize, compoSize);
        assert(typeSize == compoSize);
    }

    if (compoSize > NUMCOMPONENTS) {
        ops_printf(
            "Error! There are %i initial types defined but only %i "
            "components\n",
            compoSize, NUMCOMPONENTS);
        assert(compoSize < NUMCOMPONENTS);
    }

    for (int idx = 0; idx < typeSize; idx++) {
        ops_printf(
            "The initial condition type %i is chosen for Component "
            "%i\n",
            types.at(idx), compoId.at(idx));
        std::pair<SizeType, InitialType> pair(compoId.at(idx), types.at(idx));
        INITIALTERMS.push_back(pair);
    }
    if (compoSize < NUMCOMPONENTS) {
        ops_printf(
            "Warning! Kernel functions are required for components without "
            "pre-defined initial contiditions!\n");
    }
}

void DestroyModel() {
    FreeArrayMemory(VARIABLETYPE);
    FreeArrayMemory(VARIABLECOMPINDEX);
    FreeArrayMemory(COMPOINDEX);
    FreeArrayMemory(VARIABLECOMPPOS);
    FreeArrayMemory(XI);
    FreeArrayMemory(WEIGHTS);
    FreeArrayMemory(OPP);
}
/*
* Calculate the first-order force term
* Author: Jianping Meng  22-Feb-2019
*/
Real CalcBodyForce(const int xiIndex, const Real rho,
                   const Real* acceleration) {
    Real cf{0};
    for (int i = 0; i < LATTDIM; i++) {
        cf += CS * XI[xiIndex * LATTDIM + i] * acceleration[i];
    }
    return WEIGHTS[xiIndex] * rho * cf;
}

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

Real CalcSWEFeq(const int l, const Real h, const Real u, const Real v,
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

void SetLatticeName(const std::vector<std::string> &latticeName) {
    LATTICENAME = latticeName;
}
const std::vector<std::string> LatticeName() { return LATTICENAME; }
const std::vector<std::string> MacroVarName() { return MACROVARNAME; }
#include "model_kernel.h"
