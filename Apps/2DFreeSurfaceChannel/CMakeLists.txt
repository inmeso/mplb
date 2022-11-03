
cmake_minimum_required(VERSION 3.18)
# Application name
set(AppName Free2D)
# A list of C/C++ source files (.cpp) developed for the application
set(AppSrc lbm2d_freesurface.cpp)
# A list of C/C++ source and head files from the Src direction
# (i.e. provided by MPLB) which are used in the application
set(LibSrc evolution.cpp scheme.cpp scheme_wrapper.cpp configuration.cpp model.cpp model_wrapper.cpp block.cpp flowfield.cpp flowfield_wrapper.cpp boundary.cpp boundary_wrapper.cpp)
set(LibHeadList type.h flowfield_host_device.h boundary_host_device.h model_host_device.h)
# 2D or 3D application
set(SpaceDim 2)
if (NOT OPTIMISE)
    set(LibSrcPath "")
    foreach(Src IN LISTS LibSrc)
        list(APPEND LibSrcPath ${LibDir}/${Src})
    endforeach(Src IN LISTS LibSrc)
    SeqDevTarget("${SpaceDim}" 0)
    MpiDevTarget("${SpaceDim}" 0)
else()
    # set the files needed to be translated by ops.py from the app side
    # source file enclosing wrap function
    set(AppSrcGenList lbm2d_freesurface.cpp)
    # source file for kernel functions
    set(AppKernelGenList freesurface2d_kernel.inc)
    # if any variables are declared in both CPU and GPU memory space.
    set(AppHeadList "")
    CreateTempDir()
    set(TMP_SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/tmp)
    set(HeadList ${LibHeadList} ${AppHeadList})
    WriteJsonConfig(${TMP_SOURCE_DIR} ${AppName} "${LibSrc}" "${AppSrcGenList}" "${AppKernelGenList}" "${HeadList}" ${SpaceDim})
    TranslateSourceCodes(${LibDir} "${LibSrcGenList}" "${AppSrcGenList}" ${TMP_SOURCE_DIR})
    SeqTarget("${SpaceDim}")
    MpiTarget("${SpaceDim}")
    CudaTarget("${SpaceDim}")
endif ()
