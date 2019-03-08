# import relevant functions.
from PostProcess import ReadOPSDataHDF5
from PostProcess import WriteMacroVarsPlainHDF5
import matplotlib.pyplot as  plt
import matplotlib as matplotlib
import numpy as np

#setting up the simulation parameters
#Their meaning are same to those in the code
HaloNum = 1
MacroVarNum = 3
XINUM = 9
SPACEDIM = 2
BlockIndex = 0
MacroVarNames = ['rho', 'u', 'v']
#ratioLD=15
#ratioLH=6
front=4
nx=501
ny=251
# Reading the data at a time step into memory(res)
res = ReadOPSDataHDF5(nx, ny, BlockIndex, HaloNum,SPACEDIM, MacroVarNum, MacroVarNames, XINUM, "Flow_past_2_cylinders_Block_0_12500.h5")
#Convert the data into a layout more friendly to matlab or tecplot
WriteMacroVarsPlainHDF5(res,"testD.h5")
#Density contour plot 
plt.figure(figsize=(10,4))
plt.contour(res['MacroVars']['X'][:,:],res['MacroVars']['Y'][:,:],res['MacroVars']['v'][:,:])#,levels=np.arange(0.996,1.04,0.0001))
plt.show()
#Velocity vector plot 
plt.figure(figsize=(10,4))
plt.quiver(res['MacroVars']['X'][::10,::10].transpose(),res['MacroVars']['Y'][::10,::10].transpose(),res['MacroVars']['u'][::10,::10].transpose(),
res['MacroVars']['v'][::10,::10].transpose(), linewidth=2, cmap='autumn')
plt.show()