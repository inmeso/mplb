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
#endif // OPS_2D

#ifdef OPS_3D
void KerSetCoordinates3D(const Real* coordX, const Real* coordY,
                         const Real* coordZ, const int* idx,
                         Real* coordinates) {
    //cout<<"Check! Hi I could enter Kerset coordinates 3D"<<endl;
    coordinates[OPS_ACC_MD4(0, 0, 0, 0)] = coordX[idx[0]];
    coordinates[OPS_ACC_MD4(1, 0, 0, 0)] = coordY[idx[1]];
    coordinates[OPS_ACC_MD4(2, 0, 0, 0)] = coordZ[idx[2]];
}
#endif // OPS_3D

// This version differs from the old version in the sense that we are using the value passed from a function
//to initialise the macroscopic variables.
void KerSetInitialMacroVarsHilemms(const Real* coordinates, const int* idx, Real* macroVars, Real *macroVarsInitVal, const int *componentId) 
{
    #ifdef OPS_2D
        
        int compoIndex{*componentId};
        for(int m = VARIABLECOMPPOS[2 * compoIndex], i =0;
                 m <= VARIABLECOMPPOS[2 * compoIndex + 1]; m++, i++)
        {
            macroVars[OPS_ACC_MD2(m,0,0)] = macroVarsInitVal[OPS_ACC_MD3(i,0,0)];
        }
        /*
        macroVars[OPS_ACC_MD2(0,0,0)] = macroVarsInitVal[OPS_ACC_MD3(0,0,0)];//.00005-0.00005*coordinates[OPS_ACC_MD0(0,0,0)]/100;
        macroVars[OPS_ACC_MD2(1,0,0)] = macroVarsInitVal[OPS_ACC_MD3(1,0,0)];
        macroVars[OPS_ACC_MD2(2,0,0)] = macroVarsInitVal[OPS_ACC_MD3(2,0,0)];//-0.2+0.4*coordinates[OPS_ACC_MD0(0,0,0)];
        //macroVars[OPS_ACC_MD2(3,0,0)] =1;
        //macroVars[OPS_ACC_MD2(4,0,0)] =0;
        //macroVars[OPS_ACC_MD2(5,0,0)] =0;
        */
    #endif

    #ifdef OPS_3D

        int compoIndex{*componentId};
        
        //ops_printf("\n Component index for initialisation is %i ",compoIndex);
        //ops_printf("\n Range of m  is from %i to %i +1 ",VARIABLECOMPPOS[2 * compoIndex], VARIABLECOMPPOS[2 * compoIndex+1]);

        for(int m = VARIABLECOMPPOS[2 * compoIndex], i =0;
                 m <= VARIABLECOMPPOS[2 * compoIndex + 1]; m++, i++)
        {
            macroVars[OPS_ACC_MD2(m,0,0,0)] = macroVarsInitVal[OPS_ACC_MD3(i,0,0,0)];

            //ops_printf("\n --------------------");
            //ops_printf("\n Macro vars is = %f", macroVars[OPS_ACC_MD2(m,0,0,0)]);
        }
        
        /*
        macroVars[OPS_ACC_MD2(0, 0, 0, 0)] = macroVarsInitVal[OPS_ACC_MD3(0,0,0,0)]; //1; 
        macroVars[OPS_ACC_MD2(1, 0, 0, 0)] = macroVarsInitVal[OPS_ACC_MD3(1,0,0,0)]; //0.005;//sin(coordinates[OPS_ACC_MD0(0, 0, 0, 0)]); //u
        macroVars[OPS_ACC_MD2(2, 0, 0, 0)] = macroVarsInitVal[OPS_ACC_MD3(2,0,0,0)]; //0;         // v
        macroVars[OPS_ACC_MD2(3, 0, 0, 0)] = macroVarsInitVal[OPS_ACC_MD3(3,0,0,0)]; //0;         // w
        */
    #endif
}

#endif //HILEMMS_OPS_KERNEL