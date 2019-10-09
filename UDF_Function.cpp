Real CalcUDFFeqNew (const int XiIdx, const Real* macroVars, const int polyOrder)
{
Real result;
std::vector<int> CompIdsExcluded;
std::vector<int> VelIdsExcluded;
if(XiIdx >= 0 && XiIdx <= 18)
{
	result  =

        WEIGHTS[XiIdx] * macroVars[OPS_ACC_MD1(0,0,0,0)]*
    (
        1.0 +

        pow(3,0.5)* 
                    ( 
                        XI[XiIdx * LATTDIM] * macroVars[OPS_ACC_MD1(1,0,0,0)] +

                        XI[XiIdx * LATTDIM + 1] * macroVars[OPS_ACC_MD1(2,0,0,0)] +

                        XI[XiIdx * LATTDIM + 2] * macroVars[OPS_ACC_MD1(3,0,0,0)]

                    ) +

        1.5 *
                    (
                        XI[XiIdx * LATTDIM] * macroVars[OPS_ACC_MD1(1,0,0,0)] +

                        XI[XiIdx * LATTDIM + 1] * macroVars[OPS_ACC_MD1(2,0,0,0)] +

                        XI[XiIdx * LATTDIM + 2] * macroVars[OPS_ACC_MD1(3,0,0,0)]
                    ) *

                    (
                        XI[XiIdx * LATTDIM] * macroVars[OPS_ACC_MD1(1,0,0,0)] +

                        XI[XiIdx * LATTDIM + 1] * macroVars[OPS_ACC_MD1(2,0,0,0)] +

                        XI[XiIdx * LATTDIM + 2] * macroVars[OPS_ACC_MD1(3,0,0,0)]
                    ) -

        0.5 *
                    (
                        macroVars[OPS_ACC_MD1(1,0,0,0)] * macroVars[OPS_ACC_MD1(1,0,0,0)] +

                        macroVars[OPS_ACC_MD1(2,0,0,0)] * macroVars[OPS_ACC_MD1(2,0,0,0)] +

                        macroVars[OPS_ACC_MD1(3,0,0,0)] * macroVars[OPS_ACC_MD1(3,0,0,0)]
                    ) 


    );
}
return result;

}