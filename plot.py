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
#front=4
nx=51
ny=51

ctr_x = (nx+1)/2
ctr_y = (ny+1)/2

# Reading the data at a time step into memory(res)
res = ReadOPSDataHDF5(nx, ny, BlockIndex, HaloNum,SPACEDIM, MacroVarNum, MacroVarNames, XINUM, "2D_lid_Driven_cavity_Block_0_200.h5")

#Convert the data into a layout more friendly to matlab or tecplot
WriteMacroVarsPlainHDF5(res,"testD.h5")

#Density contour plot 
"""
fig1 = plt.figure(1)
#plt.figure(figsize=(15,5))
plt.contour(res['MacroVars']['X'][:,:],res['MacroVars']['Y'][:,:],res['MacroVars']['v'][:,:])#,levels=np.arange(0.996,1.04,0.0001))
fig1.show()

#Velocity vector plot
fig2 = plt.figure(2) 
#plt.figure(figsize=(15,5))
plt.quiver(res['MacroVars']['X'][::10,::10].transpose(),res['MacroVars']['Y'][::10,::10].transpose(),res['MacroVars']['u'][::10,::10].transpose(),
res['MacroVars']['v'][::10,::10].transpose(), linewidth=2, cmap='autumn')
fig2.show()
input()
"""


fig1 = plt.figure(1)
plt.plot(res['MacroVars']['u'][26,:], res['MacroVars']['Y'][26,:])
plt.xlabel("U")
plt.ylabel("Y")
plt.title("Variation of U on vertical centerline")
fig1.show()

fig2 = plt.figure(2)
plt.plot(res['MacroVars']['X'][:,26], res['MacroVars']['v'][:,26])
plt.xlabel("X")
plt.ylabel("V")
plt.title("Variation of V on horizontal centerline")
fig2.show()
input()
