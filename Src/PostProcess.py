"""
# Copyright 2019 United Kingdom Research and Innovation
 #
 # Authors: See AUTHORS
 #
 # Contact: [jianping.meng@stfc.ac.uk and/or jpmeng@gmail.com]
 #
 # Redistribution and use in source and binary forms, with or without
 # modification, are permitted provided that the #following conditions are met:
 #
 # 1. Redistributions of source code must retain the above copyright notice,
 #    this list of conditions and the following disclaimer.
 # 2. Redistributions in binary form must reproduce the above copyright notice
 #    this list of conditions and the following disclaimer in the documentation
 #    and or other materials provided with the distribution.
 # 3. Neither the name of the copyright holder nor the names of its contributors
 #    may be used to endorse or promote products derived from this software
 #    without specific prior written permission.
 #
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 # ANDANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 # IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 # ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 # LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 # CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 # SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 # INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 # CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE
 # ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 # POSSIBILITY OF SUCH DAMAGE.

 @brief   Post-processing utilities
 @author  Jianping Meng
 @details Providing post-processing utilities including data format transformation, basic visualisation facilities.
 usage: from PostProcess import ReadBlockData
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


def ReadVariableFromHDF5(fileName,varName,varLen=1,haloNum=1,withHalo=False):
    if ((not h5Loaded) or (not numpyLoaded)):
        print("The h5py or numpy is not installed!")
        res = "The h5py or numpy is not installed!"
        return res
    dataFile = h5.File(fileName,"r")
    blockName = list(dataFile.keys())[0]
    dataKey = varName+'_'+blockName
    rawData = np.array(dataFile[blockName][dataKey])
    spaceDim=len(rawData.shape)
    if spaceDim==3:
        nx = int(rawData.shape[2]/varLen)-2*haloNum
        ny = rawData.shape[1]-2*haloNum
        nz = rawData.shape[0]-2*haloNum
        if (varLen == 1):
            if not withHalo:
                data = rawData[haloNum:-haloNum, haloNum:-haloNum, haloNum:-haloNum]
            else:
                data = rawData
            res = data.transpose(2, 1, 0)
        if (varLen > 1):
            data = ChangeShape3D(rawData, nx, ny, nz, varLen, haloNum)
            if not withHalo:
                res = data[haloNum:-haloNum, haloNum:-haloNum, haloNum:-haloNum,:]
            else:
                res = data
    if spaceDim==2:
        nx = int(rawData.shape[1]/varLen)-2*haloNum
        ny = rawData.shape[0]-2*haloNum
        if (varLen == 1):
            if not withHalo:
                data = rawData[haloNum:-haloNum, haloNum:-haloNum]
            else:
                data = rawData
            res = data.transpose()
        if (varLen > 1):
            data = ChangeShape(rawData, nx, ny, varLen, haloNum)
            if not withHalo:
                res = data[haloNum:-haloNum, haloNum:-haloNum,:]
            else:
                res = data
    dataFile.close()
    return res

def ReadBlockData(fileName,variables):
    """Read a series of variables specified by a list of dictionary "variables" on a block from a file specified by "fileName" """
    errorMsg="Please provide a list variables in the format [{'name':'rho','len':1,'haloNum':1,'withHalo':False}"
    if not isinstance(variables,list):
        print(errorMsg)
        return None
    if not all(isinstance(var, dict) for var in variables):
        print(errorMsg)
        return None
    invalidVars=[]
    for var in variables:
        if 'name' not in var.keys():
            invalidVars.append(var)
            print("Please provide the name of varable:",var)
    validVars=[var for var in variables if (var not in invalidVars)]
    res={}
    for var in validVars:
        name = var['name']
        len = 1
        haloNum = 1
        withHalo = False
        if 'len' in var.keys():
            if isinstance(var['len'],int) and var['len']>1:
                len = var['len']
        if 'haloNum' in var.keys():
            if isinstance(var['haloNum'],int) and var['haloNum']>1:
                haloNum = var['haloNum']
        if 'withHalo' in var.keys():
            if isinstance(var['withHalo'],bool):
                withHalo = var['withHalo']
        print("Reading ",var,"...")
        res[name]=ReadVariableFromHDF5(fileName,varName=name,varLen=len,haloNum=haloNum,withHalo=withHalo)
    if "CoordinateXYZ" in res.keys():
        if res['CoordinateXYZ'].shape[-1]==3:
            res['X']=np.copy(res['CoordinateXYZ'][:,:,:,0])
            res['Y']=np.copy(res['CoordinateXYZ'][:,:,:,1])
            res['Z']=np.copy(res['CoordinateXYZ'][:,:,:,2])
        if res['CoordinateXYZ'].shape[-1]==2:
            res['X']=np.copy(res['CoordinateXYZ'][:,:,0])
            res['Y']=np.copy(res['CoordinateXYZ'][:,:,1])
        del res['CoordinateXYZ']
    return res

def WriteVariablesToPlainHDF5(res, fileName):
    """ Save the data into a plain HDF5 file"""
    if ((not h5Loaded) or (not numpyLoaded)):
        print("The h5py or numpy is not installed!")
        return
    dataFile = h5.File(fileName, "w")
    for key in res.keys():
        dataFile.create_dataset(key, data=res[key])
    dataFile.flush()
    dataFile.close()

def WriteMacroVarsTecplotHDF5(res, fileName):
    """
    Save the data into a Tecplot HDF5 file.

    Currently only works for a single block.
    """
    if ((not h5Loaded) or (not numpyLoaded)):
        print("The h5py or numpy is not installed!")
        return
    dataFile = h5.File(fileName, "w")
    spaceDim = len(res['X'].shape)
    if (2 == spaceDim):
        dataFile.create_dataset('X', data=res['X'][:, 0])
        dataFile.create_dataset('Y', data=res['Y'][0, :])
    if (3 == spaceDim):
        dataFile.create_dataset('X', data=res['X'][:, 0, 0])
        dataFile.create_dataset('Y', data=res['Y'][0, :, 0])
        dataFile.create_dataset('Z', data=res['Z'][0, 0, :])
    for key in res.keys():
        if key not in ['X','Y','Z']:
            dataFile.create_dataset(key, data=res[key])
    dataFile.flush()
    dataFile.close()

def contourPlot(x,y,var,lineNum, imgSize=1,labels=('x','y')):
    if ((not mplLoaded) or (not numpyLoaded)):
        print("The matplotlib or numpy is not installed!")
        return
    ratio = (np.max(x) - np.min(x)) / (np.max(y) - np.min(y))
    #var = (var - np.min(var)) / (np.max(var) - np.min(var))
    plt.figure(figsize=(imgSize * ratio, imgSize))
    varMin = np.min(var)
    varMax = np.max(var)
    step = (varMax - varMin) / lineNum
    conPlot = plt.contour(x, y, var, levels=np.arange(
        varMin, varMax, step), colors='black')
    plt.clabel(conPlot, inline=True, fontsize=10)
    plt.imshow(var, extent=[np.min(x), np.max(x), np.min(y), np.max(y)], origin='lower',
           cmap='plasma', alpha=0.5)
    plt.colorbar();
    plt.ylabel(labels[1])
    plt.xlabel(labels[0])
    plt.show()

def ContourPlot(res, varName, lineNum, imgSize=1):
    """ Plot a scalar contour from 2D results at a single block"""
    if ((not mplLoaded) or (not numpyLoaded)):
        print("The matplotlib or numpy is not installed!")
        return
    x = res['X']
    y = res['Y']
    var = res[varName]
    contourPlot(x,y,var,lineNum, imgSize)

def SliceContourPlot(res, varName, slice, lineNum, imgSize=1):
    """
        ContourPlot for a slice [dir,pos] perpendicular to dir='x' (|'y'|'z')
        coordinate at 'x' (|'y'|'z')=pos.

        This routine is for 3D results at a single block.

        Example: SliceContourPlot(right,'rho',['z',16],20,8)
    """
    labels=['x','y']
    if (slice[0] == 'x'):
        x = res['Z'][slice[1],:,:]
        y = res['Y'][slice[1],:,:]
        var = res[varName][slice[1], :, :]
        labels=['z','y']

    if (slice[0] == 'y'):
        x = res['X'][:,slice[1],:]
        y = res['Z'][:,slice[1],:]
        var = res[varName][:, slice[1], :]
        labels=['x','z']

    if (slice[0] == 'z'):
        x = res['X'][:,:,slice[1]]
        y = res['Y'][:,:,slice[1]]
        var = res[varName][:, :, slice[1]]

    contourPlot(x,y,var,lineNum, imgSize,labels)

def vectorPlot(X,Y,U,V,imgSize=1,labels=('x','y')):
    if ((not mplLoaded) or (not numpyLoaded) or (not mathLoaded)):
        print("The matplotlib, math or numpy is not installed!")
        return
    ratio = (np.max(X) - np.min(X)) / (np.max(Y) - np.min(Y))
    plt.figure(figsize=(imgSize * ratio, imgSize))
    varMax = max(np.max(np.abs(U)), np.max(np.abs(V)))
    U = U / varMax
    V = V / varMax
    nx, ny = X.shape
    skip = math.ceil(min(nx, ny) / 15 / math.sqrt(imgSize))
    vecPlot = plt.quiver(X[::skip, ::skip].transpose(), Y[::skip, ::skip].transpose(), U[::skip, ::skip].transpose(
    ), V[::skip, ::skip].transpose(), facecolor='greenyellow', edgecolor='greenyellow')
    vecPlot.axes.set_facecolor('k')
    plt.ylabel(labels[1])
    plt.xlabel(labels[0])
    plt.show()

def VectorPlot(res, varName, imgSize=1):
    """ Plot a vector from 2D results at a single block"""
    if ((not mplLoaded) or (not numpyLoaded) or (not mathLoaded)):
        print("The matplotlib, math or numpy is not installed!")
        return
    x = res['X']
    y = res['Y']
    varX = res[varName[0]]
    varY = res[varName[1]]
    vectorPlot(x,y,varX,varY,imgSize)

def SliceVectorPlot(res, varName, slice, imgSize=1):
    """ VectorPlot for a slice [dir,pos] perpendicular to dir='x' (|'y'|'z')
        coordinate at 'x' (|'y'|'z')=pos.

        This routine is for 3D results at a single block.

        Example: SliceVectorPlot(middle,['u','v','w'],['y',16],8)
    """
    labels=['x','y']
    if (slice[0] == 'x'):
        x = res['Z'][slice[1],:,:]
        y = res['Y'][slice[1],:,:]
        varX = res[varName[2]][slice[1],:,:]
        varY = res[varName[1]][slice[1], :, :]
        labels=['z','y']

    if (slice[0] == 'y'):
        x = res['X'][:,slice[1],:]
        y = res['Z'][:,slice[1],:]
        varX = res[varName[0]][:,slice[1],:]
        varY = res[varName[2]][:, slice[1], :]
        labels=['x','z']

    if (slice[0] == 'z'):
        x = res['X'][:,:,slice[1]]
        y = res['Y'][:,:,slice[1]]
        varX = res[varName[0]][:,:,slice[1]]
        varY = res[varName[1]][:,:,slice[1]]

    vectorPlot(x,y,varX,varY,imgSize,labels)



def VectorPlot3D(res, varName):
    if ((not mlabLoaded)):
        print("The mayavi is not installed!")
        return
    x = res['X'].transpose(1,0,2)
    y = res['Y'].transpose(1,0,2)
    z = res['Z'].transpose(1,0,2)
    varX = res[varName[0]].transpose(1,0,2)
    varY = res[varName[1]].transpose(1,0,2)
    varZ = res[varName[2]].transpose(1,0,2)
    mlab.quiver3d(x,y,z,varX,varY,varZ)
    mlab.show()

def ContourPlot3D(res, varName, contours):
    if ((not mlabLoaded)):
        print("The mayavi is not installed!")
        return
    x = res['X']
    y = res['Y']
    z = res['Z']
    var = res[varName]
    mlab.contour3d(x,y,z,var, contours=contours, transparent=True)
    mlab.show()

def CornerValues(res, varName):
    """Get the macroscopic variable value at corners. """
    var = res[varName]
    corner = {}
    nx, ny = var.shape
    corner["left bottom"] = var[0, 0]
    corner["left top"] = var[0, ny-1]
    corner["right bottom"] = var[nx-1, 0]
    corner["right top"] = var[nx-1, 0]
    return corner

def CornerValues3D(res, varName):
    """Get the macroscopic variable value at corners. """
    var = res[varName]
    corner = {}
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

def EdgeValue3D(res, varName, edge):
    """ Get macroscopic variable value at a edge:3D only """
    var = res[varName]
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
    """ Get macroscopic variable value at a face"""
    var = res[varName]
    if ('left' == face):
        return var[0, :]
    if ('right' == face):
        return var[-1, :]
    if ('bottom' == face):
        return var[:, 0]
    if ('top' == face):
        return var[:, -1]

def FaceValue3D(res, varName, face):
    """ Get macroscopic variable value at a face:3D"""
    var = res[varName]
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
