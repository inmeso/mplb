
cmake_minimum_required(VERSION 3.18)
set(AppName Cavity3D)
set(AppSrc lbm3d_cavity.cpp)
set(LibSrc evolution3d.cpp scheme.cpp scheme_wrapper.cpp configuration.cpp model.cpp model_wrapper.cpp block.cpp flowfield.cpp flowfield_wrapper.cpp boundary.cpp boundary_wrapper.cpp)

if (NOT OPTIMISE)
    # Set the path for LibSrc
    set(LibSrcPath "")
    foreach(Src IN LISTS LibSrc)
        list(APPEND LibSrcPath ${LibDir}/${Src})
    endforeach(Src IN LISTS LibSrc)
    SeqDevTarget(0)
    MpiDevTarget(0)
else()
    set(SpaceDim 3)
    #set the files needed to be translated by ops.py from the app side
    set(AppSrcGenList lbm3d_cavity.cpp)
    set(AppKernelGenList cavity3d_kernel.inc)
    set(AppHeadGenList "")
    CreateTempDir()
    set(TMP_SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/tmp)

    WriteJsonConfig(${TMP_SOURCE_DIR} ${AppName} "${AppSrcGenList}" "${AppKernelGenList}" "${AppHeadGenList}" ${SpaceDim})

    TranslateSourceCodes(${LibDir} "${LibSrcGenList}" "${AppSrcGenList}" ${TMP_SOURCE_DIR})
    SeqTarget()
    MpiTarget()
endif ()