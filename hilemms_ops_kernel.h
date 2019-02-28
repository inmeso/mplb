#ifndef HILEMMS_OPS_KERNEL
#define HILEMMS_OPS_KERNEL
#include "hilemms.h"

using namespace std;

#ifdef OPS_2D
void KerSetCoordinates(const Real* coordX, const Real* coordY, const int* idx,
                       Real* coordinates) {
    coordinates[OPS_ACC_MD3(0, 0, 0)] = coordX[idx[0]];
    coordinates[OPS_ACC_MD3(1, 0, 0)] = coordY[idx[1]];
}
#endif  // OPS_2D

#ifdef OPS_3D
void KerSetCoordinates3D(const Real* coordX, const Real* coordY,
                         const Real* coordZ, const int* idx,
                         Real* coordinates) {
    coordinates[OPS_ACC_MD4(0, 0, 0, 0)] = coordX[idx[0]];
    coordinates[OPS_ACC_MD4(1, 0, 0, 0)] = coordY[idx[1]];
    coordinates[OPS_ACC_MD4(2, 0, 0, 0)] = coordZ[idx[2]];
}
#endif  // OPS_3D

// Kernel to set initial value for a particlaur component.
// This kernel is multicomponent ready.
void KerSetInitialMacroVarsHilemms(const Real* coordinates, const int* idx,
                                   Real* macroVars, Real* macroVarsInitVal,
                                   const int* componentId) {
#ifdef OPS_2D

    int compoIndex{*componentId};
    for (int m = VARIABLECOMPPOS[2 * compoIndex], i = 0;
         m <= VARIABLECOMPPOS[2 * compoIndex + 1]; m++, i++) {
        macroVars[OPS_ACC_MD2(m, 0, 0)] =
            macroVarsInitVal[OPS_ACC_MD3(i, 0, 0)];
    }
#endif

#ifdef OPS_3D

    int compoIndex{*componentId};

    for (int m = VARIABLECOMPPOS[2 * compoIndex], i = 0;
         m <= VARIABLECOMPPOS[2 * compoIndex + 1]; m++, i++) {
        macroVars[OPS_ACC_MD2(m, 0, 0, 0)] =
            macroVarsInitVal[OPS_ACC_MD3(i, 0, 0, 0)];
    }
#endif
}

#endif  // HILEMMS_OPS_KERNEL