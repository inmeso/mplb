
cmake_minimum_required(VERSION 3.18)
# Application name
set(AppName Conservation3D)
# A list of C/C++ source files (.cpp) developed for the application
set(AppSrc conservation3d.cpp)
# A list of C/C++ source and head files from the Src direction
# (i.e. provided by MPLB) which are used in the application
set(LibSrc scheme.cpp configuration.cpp model.cpp model_wrapper.cpp block.cpp flowfield.cpp flowfield_wrapper.cpp boundary.cpp boundary_wrapper.cpp)
# 2D or 3D application
set(SpaceDim 3)
if (NOT OPTIMISE)
    set(LibSrcPath "")
    foreach(Src IN LISTS LibSrc)
        list(APPEND LibSrcPath ${LibDir}/${Src})
    endforeach(Src IN LISTS LibSrc)
    SeqDevTarget("${SpaceDim}" 0)
    MpiDevTarget("${SpaceDim}" 0)
endif ()
