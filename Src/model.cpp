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

#include "ops_seq_v2.h"
#include "block.h"
#include "model.h"
#include "model_host_device.h"
#include "flowfield.h"
#include "flowfield_host_device.h"
#include "type.h"

#include <map>
int NUMXI{9};
int FEQORDER{2};
int LATTDIM{2};
Real CS{1};
Real* XI{nullptr};
Real* WEIGHTS{nullptr};
int* OPP{nullptr};
int NUMCOMPONENTS{1};

Real XIMAXVALUE{1};
std::map<int,Component> components;
const std::map<int, Component>& g_Components() { return components; };

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

/**
 * @brief Find particles with opposite directions for bounce-back type boundary
 * using brute-force method. It could be slow for large lattice/discrete velocity set
 * @param startPos the start postion of a set of lattice in XI
 * @param latticeSize the total number of discrete/lattice velocity
 */
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
    // Gauss-Hermite quadrature from the fourth order polynomial
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
    // Gauss-Hermite quadrature from the sixth order polynomial
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

void DefineComponents(const std::vector<std::string>& compoNames,
                      const std::vector<int>& compoId,
                      const std::vector<std::string>& lattNames,
                      const std::vector<Real> tauRef,
                      const SizeType timeStep) {
    if (g_Block().size() < 1) {
        ops_printf("Error:please call DefineBlock first!\n");
        assert(g_Block().size() == 0);
    }
    NUMCOMPONENTS = compoNames.size();
    if (NUMCOMPONENTS > 0) {
        ops_printf("There are %i components defined.\n", NUMCOMPONENTS);
    } else {
        ops_printf(
            "Error! There muse be at least one component but we get:%i\n",
            NUMCOMPONENTS);
        assert(NUMCOMPONENTS > 0);
    }
    bool isLattDimSame{true};
    bool isCsSame{true};
    int totalSize{0};
    int latticeDimension{latticeSet[lattNames[0]].lattDim};
    Real currentCs{latticeSet[lattNames[0]].cs};
    for (int idx = 0; idx < NUMCOMPONENTS; idx++) {
        Component component;
        component.id = compoId.at(idx);
        component.name = compoNames.at(idx);
        component.latticeName = lattNames.at(idx);
        component.tauRef = tauRef.at(idx);
        if (latticeSet.find(lattNames[idx]) != latticeSet.end()) {
            lattice currentLattice{latticeSet[lattNames[idx]]};
            component.index[0] = totalSize;
            component.index[1] = totalSize + currentLattice.length - 1;
            totalSize += currentLattice.length;
            isLattDimSame =
                isLattDimSame && (latticeDimension == currentLattice.lattDim);
            isCsSame = isCsSame && (currentCs == currentLattice.cs);
        } else {
            ops_printf("Error! There is no predefined lattice:%s\n",
                       lattNames[idx].c_str());
            assert(latticeSet.find(lattNames[idx]) != latticeSet.end());
        }
        components.emplace(component.id, component);
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
    ops_decl_const("NUMXI", 1, "int", &NUMXI);
    ops_decl_const("CS", 1, "double", &CS);
    ops_decl_const("LATTDIM", 1, "int", &LATTDIM);
    ops_decl_const("XI", NUMXI * LATTDIM, "double", XI);
    ops_decl_const("WEIGHTS", NUMXI, "double", WEIGHTS);
    ops_decl_const("OPP", NUMXI, "int", OPP);

    for (const auto& pair : components) {
        IntField nodeType{"NodeType_" + pair.second.name};
        g_NodeType().emplace(pair.second.id, nodeType);
    }

    g_f().SetDataDim(NUMXI);
    if (timeStep == 0) {
        g_f().CreateFieldFromScratch(g_Block());
        for (auto& pair : g_NodeType()) {
            pair.second.CreateFieldFromScratch(g_Block());
        }
    } else {
        g_f().CreateFieldFromFile(CaseName(), g_Block(), timeStep);
        for (auto& pair : g_NodeType()) {
            pair.second.CreateFieldFromFile(CaseName(), g_Block(), timeStep);
        }
    }
    g_fStage().SetDataDim(NUMXI);
    g_fStage().CreateFieldFromScratch(g_Block());
    g_fStage().CreateHalos();
}

void DefineMacroVars(std::vector<VariableTypes> types,
                     std::vector<std::string> names, std::vector<int> varId,
                     std::vector<int> compoId, const SizeType timeStep) {
    if (components.size() < 1) {
        ops_printf("Error:please call DefineComponent first!\n");
        assert(components.size() == 0);
    }

    int numMacroVar{static_cast<int>(names.size())};

    if (numMacroVar > 0) {
        ops_printf("There are %i macroscopic variables defined.\n",
                   numMacroVar);
    } else {
        ops_printf("Warning! There seems no macroscopic variables defined!\n");
    }

    for (int idx = 0; idx < numMacroVar; idx++) {
        MacroVariable macroVar;
        macroVar.name = names.at(idx);
        macroVar.id = varId.at(idx);
        macroVar.type = types.at(idx);
        components.at(compoId.at(idx))
            .macroVars.emplace(macroVar.type, macroVar);
        ops_printf("The macroscopic variable %s defined for Component %s.\n",
                   names[idx].c_str(),
                   components.at(compoId.at(idx)).name.c_str());
        RealField macroVarField{macroVar.name};
        g_MacroVars().emplace(macroVar.id, macroVarField);
        if (!IsTransient()) {
            RealField macroVarFieldCopy{macroVar.name + "Copy"};
            g_MacroVarsCopy().emplace(macroVar.id, macroVarFieldCopy);
        }
    }

    for (auto& idCompo : components) {
        Component& compo{idCompo.second};
        const std::map<VariableTypes, MacroVariable>& macroVars{
            compo.macroVars};
        int uId, vId;
        if (macroVars.find(Variable_U) != macroVars.end()) {
            uId = macroVars.at(Variable_U).id;
        }
        if (macroVars.find(Variable_V) != macroVars.end()) {
            vId = macroVars.at(Variable_V).id;
        }

        if (macroVars.find(Variable_U_Force) != macroVars.end()) {
            uId = macroVars.at(Variable_U_Force).id;
        }
        if (macroVars.find(Variable_V_Force) != macroVars.end()) {
            vId = macroVars.at(Variable_V_Force).id;
        }
        compo.uId = uId;
        compo.vId = vId;
#ifdef OPS_3D
        int wId;
        if (macroVars.find(Variable_W) != macroVars.end()) {
            wId = macroVars.at(Variable_W).id;
        }
        if (macroVars.find(Variable_W_Force) != macroVars.end()) {
            wId = macroVars.at(Variable_W_Force).id;
        }
        compo.wId = wId;
#endif
    }

    if (timeStep == 0) {
        for (auto& pair : g_MacroVars()) {
            pair.second.CreateFieldFromScratch(g_Block());
        }

    } else {
        for (auto& pair : g_MacroVars()) {
            pair.second.CreateFieldFromFile(CaseName(), g_Block(), timeStep);
        }
    }

    if (!IsTransient()) {
        for (auto& pair : g_MacroVarsCopy()) {
            pair.second.CreateFieldFromScratch(g_Block());
        }
        for (const auto& compo : components) {
            for (const auto& var : compo.second.macroVars) {
                ops_reduction handle{ops_decl_reduction_handle(
                    sizeof(Real), "double", var.second.name.c_str())};
                g_ResidualErrorHandle().emplace(var.second.id, handle);
                Real error;
                g_ResidualError().emplace(var.second.id, error);
            }
        }
    }
}

void DefineCollision(std::vector<CollisionType> types,
                     std::vector<int> compoId) {
    if (components.size() < 1) {
        ops_printf("Error:please call DefineComponent first!\n");
        assert(components.size() == 0);
    }

    const SizeType typeSize{types.size()};
    const SizeType compoSize{compoId.size()};

    if (typeSize != compoSize) {
        ops_printf(
            "Error! There are %i collision types defined for  %i "
            "components\n",
            typeSize, compoSize);
        assert(typeSize == compoSize);
    }

    if (compoSize > NUMCOMPONENTS) {
        ops_printf(
            "Error! There are %i collision types defined but only %i "
            "components\n",
            compoSize, NUMCOMPONENTS);
        assert(compoSize < NUMCOMPONENTS);
    }

    for (int idx = 0; idx < typeSize; idx++) {
        ops_printf(
            "The collision type %i is chosen for "
            "Component "
            "%i\n",
            types.at(idx), compoId.at(idx));
        components.at(compoId.at(idx)).collisionType = types.at(idx);
    }
    if (compoSize < NUMCOMPONENTS) {
        ops_printf(
            "Warning! User kernel functions are required for components "
            "without pre-defined collision terms!\n");
    }
}

void DefineBodyForce(std::vector<BodyForceType> types,
                     std::vector<SizeType> compoId) {
    if (components.size() < 1) {
        ops_printf("Error:please call DefineComponent first!\n");
        assert(components.size() == 0);
    }

    const SizeType typeSize{types.size()};
    const SizeType compoSize{compoId.size()};

    if (typeSize != compoSize) {
        ops_printf(
            "Error! There are %i force types defined for  %i components\n",
            typeSize, compoSize);
        assert(typeSize == compoSize);
    }

    if (compoSize > NUMCOMPONENTS) {
        ops_printf(
            "Error! There are %i forces types defined but only %i components\n",
            compoSize, NUMCOMPONENTS);
        assert(compoSize < NUMCOMPONENTS);
    }

    for (int idx = 0; idx < typeSize; idx++) {
        ops_printf(
            "The body force function type %i is chosen for Component %i\n",
            types.at(idx), compoId.at(idx));
        components.at(compoId.at(idx)).bodyForceType = types.at(idx);
        RealField force{"Force_" + components.at(compoId.at(idx)).name};
        g_MacroBodyforce().emplace(compoId.at(idx), force);
    }
    if (compoSize < NUMCOMPONENTS) {
        ops_printf(
            "Warning! User kernel functions are required for components "
            "without "
            "pre-defined body force terms!\n");
    }
    for (auto& pair : g_MacroBodyforce()) {
        pair.second.SetDataDim(SpaceDim());
        pair.second.CreateFieldFromScratch(g_Block());
    }
}

void DefineInitialCondition(std::vector<InitialType> types,
                            std::vector<int> compoId) {
    if (components.size() < 1) {
        ops_printf("Error:please call DefineComponent first!\n");
        assert(components.size() == 0);
    }
    const SizeType typeSize{types.size()};
    const SizeType compoSize{compoId.size()};

    if (typeSize != compoSize) {
        ops_printf(
            "Error! There are %i initial types defined for  %i components\n",
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
        ops_printf("The initial condition type %i is chosen for Component %i\n",
                   types.at(idx), compoId.at(idx));
        components.at(compoId.at(idx)).initialType = types.at(idx);
    }
    if (compoSize < NUMCOMPONENTS) {
        ops_printf(
            "Warning! User kernel functions are required for components without "
            "pre-defined initial contiditions!\n");
    }
}

void DestroyModel() {
    FreeArrayMemory(XI);
    FreeArrayMemory(WEIGHTS);
    FreeArrayMemory(OPP);
}

