# Copyright 2018 the MPLB team. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.
"""
    @brief   Post-processing utilities
    @author  Jianping Meng
    @details Providing post-processing utilities including data format transformation, basic visualisation facilities.
    usage: from PostProcess import ReadOPSDataHDF5
           from PostProcess import WriteMacroVarsPlainHDF5
    Specific examples can be found in provide source codes
    Dependency: numpy, h5py, matplotlib for 2D visualisation, and mayavi2 for 3D visualisation.
"""
# python 2 and python 3 compatibility for the print function
from __future__ import print_function

try:
    import numpy as np
    numpyLoaded = True
except ImportError:
    numpyLoaded = False
if (not numpyLoaded):
    print("The numpy module cannot be imported! Without it, all functions cannot work!")

try:
    import h5py as h5
    h5Loaded = True
except ImportError:
    h5Loaded = False
if (not h5Loaded):
    print("The h5py module cannot be imported!  Without it, all functions cannot work!")

try:
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mplLoaded = True
except ImportError:
    mplLoaded = False
if (not mplLoaded):
    print("The matplotlib module cannot be imported! Please install it for two-dimensional visualisation!")

try:
    import math
    mathLoaded = True
except ImportError:
    mathLoaded = False
if (not mathLoaded):
    print("The math module cannot be imported! Please install it for two-dimensional visualisation!")
try:
    from mayavi import mlab
    mlabLoaded = True
except ImportError:
    mlabLoaded = False
if (not mlabLoaded):
    print("The mayavi module cannot be imported! Please install it for three-dimensional visualisation!")

def ChangeShape(data, nx, ny, dataLength, haloNum):
    """Converting the storage order of multidim array in 2D space."""
    data = data.reshape((ny + 2 * haloNum, nx + 2 * haloNum, dataLength))
    return data.transpose((1, 0, 2))

def ChangeShape3D(data, nx, ny, nz, dataLength, haloNum):
    """Converting the storage order of multidim array in 3D space."""
    data = data.reshape((nz + 2 * haloNum, ny + 2 * haloNum,
                         nx + 2 * haloNum, dataLength))
    return data.transpose((2, 1, 0, 3))


def ReadOPSDataHDF5(nx, ny, blockIndex, haloNum, spaceDim, macroVarNum, macroVarNames, xiNum, fileName):
    """Converting a 2D result file into a single dictionary enclosing two sub-dictionaries, MacroVars and Distributions. In particular, all vectors or tensors will be accessed through components. """
    if ((not h5Loaded) or (not numpyLoaded)):
        print("The h5py or numpy is not installed!")
        res = "The h5py or numpy is not installed!"
        return res
    strBlockIdx = str(blockIndex)
    dataFile = h5.File(fileName)
    blockKey = 'Block_' + strBlockIdx
    macroVars = {}
    distributions = {}

    for dataKey in dataFile[blockKey].keys():
        if 'Bodyforce_' + strBlockIdx == dataKey:
            tmpvars = ChangeShape(
                dataFile[blockKey][dataKey][:, :], nx, ny, xiNum, haloNum)
            distributions['Force'] = tmpvars[haloNum:-
                                             haloNum, haloNum:-haloNum, :]
        if 'MacroVars_' + strBlockIdx == dataKey:
            tmpvars = ChangeShape(
                dataFile[blockKey][dataKey][:, :], nx, ny, macroVarNum, haloNum)
            for i in range(len(macroVarNames)):
                macroVars[macroVarNames[i]] = tmpvars[haloNum:-
                                                      haloNum, haloNum:-haloNum, i]
        if 'Tau_' + strBlockIdx == dataKey:
            data = dataFile[blockKey][dataKey][haloNum:-
                                               haloNum, haloNum:-haloNum]
            macroVars['Tau'] = data.transpose()
        if 'Nodetype_' + strBlockIdx == dataKey:
            data = dataFile[blockKey][dataKey][haloNum:-
                                               haloNum, haloNum:-haloNum]
            macroVars['NT'] = data.transpose()
        if 'GeometryProperty_' + strBlockIdx == dataKey:
            data = dataFile[blockKey][dataKey][haloNum:-
                                               haloNum, haloNum:-haloNum]
            macroVars['GP'] = data.transpose()
        if 'CoordinateXYZ_' + strBlockIdx == dataKey:
            tmpvars = ChangeShape(
                dataFile[blockKey][dataKey][:, :], nx, ny, spaceDim, haloNum)
            macroVars['X'] = tmpvars[haloNum:-haloNum, haloNum:-haloNum, 0]
            macroVars['Y'] = tmpvars[haloNum:-haloNum, haloNum:-haloNum, 1]
            distributions['X'] = tmpvars[haloNum:-haloNum, haloNum:-haloNum, 0]
            distributions['Y'] = tmpvars[haloNum:-haloNum, haloNum:-haloNum, 1]
        if 'fStage_' + strBlockIdx == dataKey:
            tmpvars = ChangeShape(
                dataFile[blockKey][dataKey][:, :], nx, ny, xiNum, haloNum)
            distributions['Fstage'] = tmpvars[haloNum:-
                                              haloNum, haloNum:-haloNum, :]
        if 'f_' + strBlockIdx == dataKey:
            tmpvars = ChangeShape(
                dataFile[blockKey][dataKey][:, :], nx, ny, xiNum, haloNum)
            distributions['F'] = tmpvars[haloNum:-
                                         haloNum, haloNum:-haloNum, :]
        if 'feq_' + strBlockIdx == dataKey:
            tmpvars = ChangeShape(
                dataFile[blockKey][dataKey][:, :], nx, ny, xiNum, haloNum)
            distributions['Feq'] = tmpvars[haloNum:-
                                           haloNum, haloNum:-haloNum, :]
    res = {}
    res['MacroVars'] = macroVars
    res['Distributions'] = distributions
    dataFile.close()
    return res


def ReadOPSDataHDF53D(nx, ny, nz, blockIndex, haloNum, spaceDim, macroVarNum, macroVarNames, xiNum, fileName):
    """Converting a 3D result file into a single dictionary enclosing two sub-dictionaries, MacroVars and Distributions. In particular, all vectors or tensors will be accessed through components. """
    if ((not h5Loaded) or (not numpyLoaded)):
        print("The h5py or numpy is not installed!")
        res = "The h5py or numpy is not installed!"
        return res
    strBlockIdx = str(blockIndex)
    dataFile = h5.File(fileName)
    blockKey = 'Block_' + strBlockIdx
    macroVars = {}
    distributions = {}
    for dataKey in dataFile[blockKey].keys():
        if 'Bodyforce_' + strBlockIdx == dataKey:
            tmpvars = ChangeShape3D(
                dataFile[blockKey][dataKey][:, :, :], nx, ny, nz, xiNum, haloNum)
            distributions['Force'] = tmpvars[haloNum:-
                                             haloNum, haloNum:-haloNum, haloNum:-haloNum, :]
        if 'MacroVars_' + strBlockIdx == dataKey:
            tmpvars = ChangeShape3D(
                dataFile[blockKey][dataKey][:, :, :], nx, ny, nz, macroVarNum, haloNum)
            for i in range(len(macroVarNames)):
                macroVars[macroVarNames[i]] = tmpvars[haloNum:-
                                                      haloNum, haloNum:-haloNum, haloNum:-haloNum, i]
        if 'Tau_' + strBlockIdx == dataKey:
            data = dataFile[blockKey][dataKey][haloNum:-
                                               haloNum, haloNum:-haloNum, haloNum:-haloNum]
            macroVars['Tau'] = data.transpose(2, 1, 0)
        if 'Nodetype_' + strBlockIdx == dataKey:
            data = dataFile[blockKey][dataKey][haloNum:-
                                               haloNum, haloNum:-haloNum, haloNum:-haloNum]
            macroVars['NT'] = data.transpose(2, 1, 0)
        if 'GeometryProperty_' + strBlockIdx == dataKey:
            data = dataFile[blockKey][dataKey][haloNum:-
                                               haloNum, haloNum:-haloNum, haloNum:-haloNum]
            macroVars['GP'] = data.transpose(2, 1, 0)
        if 'CoordinateXYZ_' + strBlockIdx == dataKey:
            tmpvars = ChangeShape3D(
                dataFile[blockKey][dataKey][:, :, :], nx, ny, nz, spaceDim, haloNum)
            macroVars['X'] = tmpvars[haloNum:-
                                     haloNum, haloNum:-haloNum, haloNum:-haloNum, 0]
            macroVars['Y'] = tmpvars[haloNum:-
                                     haloNum, haloNum:-haloNum, haloNum:-haloNum, 1]
            macroVars['Z'] = tmpvars[haloNum:-
                                     haloNum, haloNum:-haloNum, haloNum:-haloNum, 2]
            distributions['X'] = tmpvars[haloNum:-
                                         haloNum, haloNum:-haloNum, haloNum:-haloNum, 0]
            distributions['Y'] = tmpvars[haloNum:-
                                         haloNum, haloNum:-haloNum, haloNum:-haloNum, 1]
            distributions['Z'] = tmpvars[haloNum:-
                                         haloNum, haloNum:-haloNum, haloNum:-haloNum, 2]
        if 'fStage_' + strBlockIdx == dataKey:
            tmpvars = ChangeShape3D(
                dataFile[blockKey][dataKey][:, :, :], nx, ny, nz, xiNum, haloNum)
            distributions['Fstage'] = tmpvars[haloNum:-
                                              haloNum, haloNum:-haloNum, haloNum:-haloNum, :]
        if 'f_' + strBlockIdx == dataKey:
            tmpvars = ChangeShape3D(
                dataFile[blockKey][dataKey][:, :, :], nx, ny, nz, xiNum, haloNum)
            distributions['F'] = tmpvars[haloNum:-
                                         haloNum, haloNum:-haloNum, haloNum:-haloNum, :]
        if 'feq_' + strBlockIdx == dataKey:
            tmpvars = ChangeShape3D(
                dataFile[blockKey][dataKey][:, :, :], nx, ny, nz, xiNum, haloNum)
            distributions['Feq'] = tmpvars[haloNum:-
                                           haloNum, haloNum:-haloNum, haloNum:-haloNum, :]
    res = {}
    res['MacroVars'] = macroVars
    res['Distributions'] = distributions
    dataFile.close()
    return res

def WriteMacroVarsPlainHDF5(res, fileName):
    """ Save the data into a plain HDF5 file"""
    if ((not h5Loaded) or (not numpyLoaded)):
        print("The h5py or numpy is not installed!")
        return
    dataFile = h5.File(fileName, "w")
    for key in res['MacroVars'].keys():
        dataFile.create_dataset(key, data=res['MacroVars'][key])
    dataFile.flush()
    dataFile.close()

def WriteMacroVarsTecplotHDF5(res, fileName):
    """ Save the data into a Tecplot HDF5 file"""
    if ((not h5Loaded) or (not numpyLoaded)):
        print("The h5py or numpy is not installed!")
        return
    dataFile = h5.File(fileName, "w")
    spaceDim = len(res['MacroVars']['X'].shape)
    if (2 == spaceDim):
        res['MacroVars']['X'] = res['MacroVars']['X'][:, 0]
        res['MacroVars']['Y'] = res['MacroVars']['Y'][0, :]
    if (3 == spaceDim):
        res['MacroVars']['X'] = res['MacroVars']['X'][:, 0, 0]
        res['MacroVars']['Y'] = res['MacroVars']['Y'][0, :, 0]
        res['MacroVars']['Z'] = res['MacroVars']['Z'][0, 0, :]

    for key in res['MacroVars'].keys():
        dataFile.create_dataset(key, data=res['MacroVars'][key])
    dataFile.flush()
    dataFile.close()

def ContourPlot(res, varName, lineNum, imgSize=1):
    """Plot the Contour of a scalar.
       Note: for better results, the minimal value is taken away from the data, and the data normalised by its maximal difference. """
    if ((not mplLoaded) or (not numpyLoaded)):
        print("The matplotlib or numpy is not installed!")
        return
    x = res['MacroVars']['X']
    y = res['MacroVars']['Y']
    ratio = (np.max(x) - np.min(x)) / (np.max(y) - np.min(y))
    plt.figure(figsize=(imgSize * ratio, imgSize))
    var = res['MacroVars'][varName]
    var = (var - np.min(var)) / (np.max(var) - np.min(var))
    varMin = np.min(var)
    varMax = np.max(var)
    step = (varMax - varMin) / lineNum
    conPlot = plt.contour(x, y, var, levels=np.arange(
        varMin, varMax, step), colors='greenyellow')
    conPlot.ax.set_facecolor('k')
    plt.clabel(conPlot, inline=1, fontsize=10)
    plt.show()


def VectorPlot(res, varName, imgSize=1):
    if ((not mplLoaded) or (not numpyLoaded) or (not mathLoaded)):
        print("The matplotlib, math or numpy is not installed!")
        return
    x = res['MacroVars']['X']
    y = res['MacroVars']['Y']
    ratio = (np.max(x) - np.min(x)) / (np.max(y) - np.min(y))
    plt.figure(figsize=(imgSize * ratio, imgSize))
    varX = res['MacroVars'][varName[0]]
    varY = res['MacroVars'][varName[1]]
    varMax = max(np.max(np.abs(varX)), np.max(np.abs(varY)))
    varX = varX / varMax
    varY = varY / varMax
    nx, ny = x.shape
    skip = math.ceil(min(nx, ny) / 15 / math.sqrt(imgSize))
    vecPlot = plt.quiver(x[::skip, ::skip].transpose(), y[::skip, ::skip].transpose(), varX[::skip, ::skip].transpose(
    ), varY[::skip, ::skip].transpose(), facecolor='greenyellow', edgecolor='greenyellow')
    vecPlot.axes.set_facecolor('k')
    plt.show()

def VectorPlot3D(res, varName):
    if ((not mlabLoaded)):
        print("The mayavi is not installed!")
        return
    x = res['MacroVars']['X'].transpose(1,0,2)
    y = res['MacroVars']['Y'].transpose(1,0,2)
    z = res['MacroVars']['Z'].transpose(1,0,2)
    varX = res['MacroVars'][varName[0]].transpose(1,0,2)
    varY = res['MacroVars'][varName[1]].transpose(1,0,2)
    varZ = res['MacroVars'][varName[2]].transpose(1,0,2)
    mlab.quiver3d(x,y,z,varX,varY,varZ)
    mlab.show()

def ContourPlot3D(res, varName, contours):
    if ((not mlabLoaded)):
        print("The mayavi is not installed!")
        return
    x = res['MacroVars']['X']
    y = res['MacroVars']['Y']
    z = res['MacroVars']['Z']
    var = res['MacroVars'][varName]
    mlab.contour3d(x,y,z,var, contours=contours, transparent=True)
    mlab.show()

def CornerValues(res, varName):
    """get the variable value at corners. """
    var = res['MacroVars'][varName]
    spaceDim = len(var.shape)
    corner = {}
    if (2 == spaceDim):
        nx, ny = var.shape
        corner["left bottom"] = var[0, 0]
        corner["left top"] = var[0, ny-1]
        corner["right bottom"] = var[nx-1, 0]
        corner["right top"] = var[nx-1, 0]
    if (3 == spaceDim):
        nx, ny, nz = var.shape
        corner["left bottom back"] = var[0, 0, 0]
        corner["left bottom front"] = var[0, 0, nz-1]
        corner["left top back"] = var[0, ny-1, 0]
        corner["left top front"] = var[0, ny-1, nz-1]
        corner["right bottom back"] = var[nx-1, 0, 0]
        corner["right bottom front"] = var[nx-1, 0, nz-1]
        corner["right top back"] = var[nx-1, ny-1, 0]
        corner["right top front"] = var[nx-1, ny-1, nz-1]
    return corner


def EdgeValue(res, varName, edge):
    """ get the variable value at a edge : 3D only """
    var = res['MacroVars'][varName]
    if ('left bottom' == edge):
        return var[0, 0, :]
    if ('left top' == edge):
        return var[0, -1, :]
    if ('right bottom' == edge):
        return var[-1, 0, :]
    if ('right top' == edge):
        return var[-1, -1, :]

    if ('left back' == edge):
            return var[0, :, 0]
    if ('left front' == edge):
        return var[0, :, -1]
    if ('right back' == edge):
        return var[-1, :, 0]
    if ('right front' == edge):
        return var[-1, :, -1]

    if ('bottom back' == edge):
            return var[:, 0, 0]
    if ('bottom front' == edge):
        return var[:, 0, -1]
    if ('top back' == edge):
        return var[:, -1, 0]
    if ('top front' == edge):
        return var[:, -1, -1]


def FaceValue(res, varName, face):
    """ get the variable value at a edge : 3D only """
    var = res['MacroVars'][varName]
    spaceDim = len(var.shape)
    if (2 == spaceDim):
        if ('left' == face):
            return var[0, :]
        if ('right' == face):
            return var[-1, :]
        if ('bottom' == face):
            return var[:, 0]
        if ('top' == face):
            return var[:, -1]

    if (3 == spaceDim):
        if ('left' == face):
            return var[0, :, :]
        if ('right' == face):
            return var[-1, :, :]
        if ('bottom' == face):
            return var[:, 0, :]
        if ('top' == face):
            return var[:, -1, :]
        if ('back' == face):
            return var[:, :, 0]
        if ('front' == face):
            return var[:, :, -1]