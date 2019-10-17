Real CalcUDFFeqNew (const int XiIdx, const Real* macroVars, const int polyOrder)
{
//Please enter the CPP file with main() function here.
//    CppFileName = lbm3d_hilemms.cpp;

for(int XiIndex =0; XiIndex<=18; XiIndex++ )
{

f[OPS_ACC_MD2(COMPOINDEX[2 *0] + XiIndex,0,0,0)] =

        WEIGHTS[COMPOINDEX[2 *0] + XiIndex] * macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 0,0,0,0)]*
    (
        1.0 +

        pow(3,0.5)* 
                    ( 
                        XI[ (COMPOINDEX[2 *0] + XiIndex) * LATTDIM] * macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 1,0,0,0)] +

                        XI[ (COMPOINDEX[2 *0] + XiIndex) * LATTDIM +1] * macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 2,0,0,0)] +

                        XI[ (COMPOINDEX[2 *0] + XiIndex) * LATTDIM +2] * macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 3,0,0,0)]

                    ) +

        1.5 *
                    (
                        XI[ (COMPOINDEX[2 *0] + XiIndex) * LATTDIM] * macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 1,0,0,0)] +

                        XI[ (COMPOINDEX[2 *0] + XiIndex) * LATTDIM +1] * macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 2,0,0,0)] +

                        XI[ (COMPOINDEX[2 *0] + XiIndex) * LATTDIM +2] * macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 3,0,0,0)]
                    ) *

                    (
                        XI[ (COMPOINDEX[2 *0] + XiIndex) * LATTDIM] * macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 1,0,0,0)] +

                        XI[ (COMPOINDEX[2 *0] + XiIndex) * LATTDIM +1] * macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 2,0,0,0)] +

                        XI[ (COMPOINDEX[2 *0] + XiIndex) * LATTDIM +2] * macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 3,0,0,0)]
                    ) -

        0.5 *
                    (
                        macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 1,0,0,0)] * macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 1,0,0,0)] +

                        macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 2,0,0,0)] * macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 2,0,0,0)] +

                        macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 3,0,0,0)] * macroVars[OPS_ACC_MD1(VARIABLECOMPPOS[2 * 0] + 3,0,0,0)]
                    ) 


    );

} 
//End of for loop. 
return result;

}