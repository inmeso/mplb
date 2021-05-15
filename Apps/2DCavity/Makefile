#OPS_INSTALL_PATH := $(OPS_INSTALL_PATH)/c
OPS_INC   = -I$(OPS_INSTALL_PATH)/include
OPS_LIB   = -L$(OPS_INSTALL_PATH)/lib
CUDA_LIB = -L$(CUDA_INSTALL_PATH)/lib
CUDA_INC = -I$(CUDA_INSTALL_PATH)/include

OPENCL_INC = -I$(OPENCL_INSTALL_PATH)/include
OPENCL_LIB = -L$(OPENCL_INSTALL_PATH)/lib64 -lOpenCL

ifdef MPLB_DIR
	MPLB_INC  = -I$(MPLB_DIR)
	MPLB_LIB  = -L$(MPLB_DIR)
endif


ifdef HDF5_INSTALL_PATH
  HDF5_INC 	  	:= -I$(HDF5_INSTALL_PATH)/include
  HDF5_LIB 	  	:= -L$(HDF5_INSTALL_PATH)/lib
endif
ifdef CHECKPOINTING
  HDF5_INC              := $(HDF5_INC) -DCHECKPOINTING
endif

ifeq ($(OPS_COMPILER),clang)
  CPP   = clang++
ifdef DEBUG
  CPPFLAGS  = -O0 -g -DUNIX -Wall -std=c++11 #-Wall
else
  CPPFLAGS  = -O3  -fPIC -DUNIX -Wall -std=c++11  #-Wall
endif
  OMPFLAGS  = ""#-fopenmp
  MPICPP  = mpic++
  MPIFLAGS  = $(CPPFLAGS)
else
ifeq ($(OPS_COMPILER),gnu)
  CPP   = g++
ifdef DEBUG
  CPPFLAGS  = -ggdb -DUNIX -Wall
else
  CPPFLAGS  = -O3 -std=c++11 -fPIC -DUNIX
endif
  OMPFLAGS  = -fopenmp
  MPICPP  = mpicxx
  MPIFLAGS  = $(CPPFLAGS)
else
ifeq ($(OPS_COMPILER),intel)
  CPP   = icpc
ifdef DEBUG
  CCFLAGS = -O0 -g -no-prec-div -openmp -fp-model strict -fp-model source -prec-div -prec-sqrt -DMPICH_IGNORE_CXX_SEEK #-DOPS_DEBUG
else
  #CCFLAGS  = -O3 -ipo -no-prec-div -fp-model strict -fp-model source -prec-div -prec-sqrt -vec-report2 -xSSE4.2 -parallel #-DCOMM_PERF #-DDEBUG
  CCFLAGS = -O3 -ipo -no-prec-div -restrict -fno-alias -fp-model strict -fp-model source -prec-div -prec-sqrt -DMPICH_IGNORE_CXX_SEEK#-vec-report
endif
  CPPFLAGS  = $(CCFLAGS)
  OMPFLAGS  = -openmp -openmp-report2
  MPICPP  = $(MPI_INSTALL_PATH)/bin/mpicxx
  MPIFLAGS  = $(CCFLAGS)
else
ifeq ($(OPS_COMPILER),mic)
  CPP   = icpc
ifdef DEBUG
  CCFLAGS = -O0 -g -mmic -DMPICH_IGNORE_CXX_SEEK#-vec-report
else
  #CCFLAGS  = -O3 -ipo -no-prec-div -restrict -fno-alias -fp-model strict -fp-model source -prec-div -prec-sqrt -DMPICH_IGNORE_CXX_SEEK#-vec-report
  #CCFLAGS  = -O3 -ipo -no-prec-div -restrict -fno-alias -mmic -vec-report -finline -inline-forceinline -DSIMD_VEC=32 -DMPICH_IGNORE_CXX_SEEK#-vec-report
  CCFLAGS = -O3 -ipo -no-prec-div -restrict -fno-alias -mmic -vec-report -DSIMD_VEC=32 -DMPICH_IGNORE_CXX_SEEK -finline -vec-report
endif
  CPPFLAGS  = $(CCFLAGS)
  OMPFLAGS  = -openmp -openmp-report2
  MPICPP  = $(MPI_INSTALL_PATH)/bin/mpicxx
  MPIFLAGS  = $(CCFLAGS)
  OPENCL_LIB = -L$(OPENCL_INSTALL_PATH)/lib64 -lintelocl -lOpenCL #for offload opencl execution on the phi
else
ifeq ($(OPS_COMPILER),pgi)
  CPP     = pgc++
  CC      = pgcc
  CUDA_ALIGN_FLAG := -D__x86_64 -D__align__\(n\)=__attribute__\(\(aligned\(n\)\)\) -D__location__\(a\)=__annotate__\(a\) -DCUDARTAPI=
ifdef DEBUG
  CCFLAGS   = -O2 -Minline -Kieee -Minform=inform -Minfo=all
else
  CCFLAGS   = -O3 -fast -Kieee -Minline -Minform=severe -Minfo=all -std=c++11
endif
  CPPFLAGS  = $(CCFLAGS)
  OMPFLAGS  = -mp
  MPICPP  = $(MPI_INSTALL_PATH)/bin/mpicxx
  MPICC   = $(MPI_INSTALL_PATH)/bin/mpicc
  MPIFLAGS  = $(CPPFLAGS)
else
ifeq ($(OPS_COMPILER),xl)
  CPP   = xlc++
ifdef DEBUG
  CPPFLAGS  = -O0 -g -qfullpath -Q -qlistopt -qattr=full -qlist -qreport -qxref=full -qsource -qsuppress=1506-224:1500-036 -qsrcmsg -qfloat=nomaf
else
  CPPFLAGS  = -O5 -q64 -qnoeh -qipa=partition=large -qfullpath -Q -qlistopt -qattr=full -qlist -qreport -qxref=full -qsource -qsuppress=1506-224:1500-036 -qsrcmsg -qfloat=nomaf
endif
  OMPFLAGS  = -fopenmp
  MPICPP  = $(MPI_INSTALL_PATH)/bin/mpicxx
  MPIFLAGS  = $(CPPFLAGS)
else
ifeq ($(OPS_COMPILER),cray)
  CPP       = CC
    CUDA_ALIGN_FLAG := -D__x86_64 -D__align__\(n\)=__attribute__\(\(aligned\(n\)\)\) -D__location__\(a\)=__annotate__\(a\) -DCUDARTAPI=
ifdef DEBUG
  CCFLAGS   = -O2
else
  CCFLAGS   = -O3
endif
  CPPFLAGS  = $(CCFLAGS)
  OMPFLAGS  =
  MPICPP    = CC
  MPICC     = cc
  MPIFLAGS  = $(CPPFLAGS)
else
print:
	@echo "unrecognised value for OPS_COMPILER"
endif
endif
endif
endif
endif
endif
endif

NVCC  := $(CUDA_INSTALL_PATH)/bin/nvcc
# flags for nvcc
# set NV_ARCH to select the correct one
ifndef NV_ARCH
  MESSAGE=select an NVIDA device to compile in CUDA, e.g. make NV_ARCH=KEPLER
  NV_ARCH=Kepler
endif
ifeq ($(NV_ARCH),Fermi)
  CODE_GEN_CUDA=-gencode arch=compute_20,code=sm_21
else
ifeq ($(NV_ARCH),Kepler)
  CODE_GEN_CUDA=-gencode arch=compute_35,code=sm_35
endif
endif


ifdef DEBUG
  NVCC_FLAGS = -O0 -g -G
else
  NVCC_FLAGS = -O3 -std=c++11
endif


#
# master to make all versions
#

TARGETS = lbm2d_dev_seq lbm2d_seq lbm2d_dev_mpi lbm2d_mpi lbm2d_openmp lbm2d_mpi_openmp lbm2d_cuda lbm2d_mpi_cuda lbm2d_opencl lbm2d_mpi_opencl

ifeq ($(OPS_COMPILER),pgi)
#TARGETS += lbm2d_openacc lbm2d_mpi_openacc
else
#TARGETS += openacc_error
endif

all: clean $(TARGETS)

env:
	@echo "MPLBDIR = $(MPLB_DIR)"
	@echo "HDF5_INCLUDE_PATH = $(HDF5_INCLUDE_PATH)"
	@echo "HDF5_LIB_PATH = $(HDF5_LIB_PATH)"
	@echo "OPS_INSTALL_PATH = $(OPS_INSTALL_PATH)"
	@echo "Compiler = $(OPS_COMPILER)"

lbm3d_dev_seq :
	$(MPICPP) $(CPPFLAGS)  $(MAINCPP) ${MPLB_DIR}/evolution3d.cpp ${MPLB_DIR}/scheme.cpp ${MPLB_DIR}/scheme_wrapper.cpp ${MPLB_DIR}/configuration.cpp ${MPLB_DIR}/model.cpp ${MPLB_DIR}/model_wrapper.cpp ${MPLB_DIR}/block.cpp ${MPLB_DIR}/flowfield.cpp ${MPLB_DIR}/flowfield_wrapper.cpp ${MPLB_DIR}/boundary.cpp ${MPLB_DIR}/boundary_wrapper.cpp $(OPS_INC) $(HDF5_INC) $(MPLB_INC) $(OPS_LIB) $(HDF5_LIB) -DOPS_3D -DCPU -D$(LEVEL)  -lops_hdf5_seq -lops_seq -lhdf5_hl -lhdf5 -lz  -o lbm3d_dev_seq

lbm3d_dev_seq_lib :
	$(CPP) $(CPPFLAGS) $(MAINCPP) $(OPS_INC) $(HDF5_INC) $(MPLB_INC) $(OPS_LIB) $(MPLB_LIB) -DOPS_3D -DCPU -D$(LEVEL)  -lmplb3d_seq -lops_seq -lops_hdf5_seq  $(HDF5_LIB)  -o lbm3d_dev_seq

lbm3d_dev_mpi :
	$(MPICPP) $(CPPFLAGS)  $(MAINCPP) ${MPLB_DIR}/evolution3d.cpp ${MPLB_DIR}/scheme.cpp ${MPLB_DIR}/scheme_wrapper.cpp ${MPLB_DIR}/configuration.cpp ${MPLB_DIR}/model.cpp ${MPLB_DIR}/model_wrapper.cpp ${MPLB_DIR}/block.cpp ${MPLB_DIR}/flowfield.cpp ${MPLB_DIR}/flowfield_wrapper.cpp ${MPLB_DIR}/boundary.cpp ${MPLB_DIR}/boundary_wrapper.cpp $(OPS_INC) $(HDF5_INC) $(MPLB_INC) $(OPS_LIB) $(MPLB_LIB) $(HDF5_LIB) -DOPS_3D -DOPS_MPI -DCPU -D$(LEVEL)  -lops_hdf5_mpi -lops_mpi -lhdf5_hl -lhdf5 -lz  -o lbm3d_dev_mpi

lbm3d_dev_mpi_lib :
	$(MPICPP) $(MPIFLAGS) $(MAINCPP) $(OPS_INC) $(HDF5_INC) $(MPLB_INC) $(OPS_LIB) $(MPLB_LIB) -DOPS_MPI -DOPS_3D -DCPU -D$(LEVEL)  -lmplb3d_mpi -lops_mpi -lops_hdf5_mpi  $(HDF5_LIB)  -o lbm3d_dev_mpi

lbm2d_dev_seq :
	$(CPP) $(CPPFLAGS) $(MAINCPP) $(shell find $(MPLB_DIR) -type f \( -iname "*.cpp" ! -iname "setup*"  !  -iname "evolution3d.cpp" \))  $(OPS_INC) $(HDF5_INC) $(MPLB_INC) $(OPS_LIB) $(MPLB_LIB) -DOPS_2D -DCPU -D$(LEVEL)  -lops_seq -lops_hdf5_seq  $(HDF5_LIB)  -o lbm2d_dev_seq

lbm2d_dev_seq_lib :
	$(CPP) $(CPPFLAGS) $(MAINCPP) $(OPS_INC) $(HDF5_INC) $(MPLB_INC) $(OPS_LIB) $(MPLB_LIB) -DOPS_2D -DCPU -D$(LEVEL)  -lmplb2d_seq -lops_seq -lops_hdf5_seq  $(HDF5_LIB)  -o lbm2d_dev_seq

lbm2d_dev_mpi :
	$(MPICPP) $(MPIFLAGS) $(MAINCPP) $(shell find $(MPLB_DIR) -type f \( -iname "*.cpp" ! -iname "setup*"  !  -iname "evolution3d.cpp" \))  $(OPS_INC) $(HDF5_INC) $(MPLB_INC) $(OPS_LIB) $(MPLB_LIB) -DOPS_MPI -DOPS_2D -DCPU -D$(LEVEL)  -lops_mpi -lops_hdf5_mpi  $(HDF5_LIB)  -o lbm2d_dev_mpi

lbm2d_dev_mpi_lib :
	$(MPICPP) $(MPIFLAGS) $(MAINCPP) $(OPS_INC) $(HDF5_INC) $(MPLB_INC) $(OPS_LIB) $(MPLB_LIB) -DOPS_MPI -DOPS_2D -DCPU -D$(LEVEL)  -lmplb2d_mpi -lops_mpi -lops_hdf5_mpi  $(HDF5_LIB)  -o lbm2d_dev_mpi

# Old pre-processor 2D
setupdomain: Makefile setup_comput_domain.cpp scheme.cpp model.cpp boundary.cpp Case_Setup.cpp $(OPS_INSTALL_PATH)/lib/libops_seq.a
	$(CPP) $(CPPFLAGS) $(OPS_INC) $(HDF5_INC) $(OPS_LIB) -DOPS_2D type.cpp boundary.cpp setup_comput_domain.cpp scheme.cpp model.cpp -lops_seq -lops_hdf5_seq $(HDF5_LIB) -o setup_comput_domain

# Old pre-processor 3D
setupdomain3D: Makefile setup_comput_domain.cpp scheme.cpp model.cpp boundary.cpp $(OPS_INSTALL_PATH)/lib/libops_seq.a
	$(CPP) $(CPPFLAGS) $(OPS_INC) $(HDF5_INC) $(OPS_LIB) -DOPS_3D type.cpp boundary.cpp setup_comput_domain.cpp scheme.cpp model.cpp -lops_seq -lops_hdf5_seq $(HDF5_LIB) -o setup_comput_domain

#
# mpi version
#

SEQ_KERNELS = lbm2d_seq_kernels.cpp

lbm2d_mpi:  Makefile type_ops.cpp boundary_ops.cpp scheme_ops.cpp flowfield_ops.cpp  model_ops.cpp evolution_ops.cpp lbm2d_ops.cpp $(OPS_INSTALL_PATH)/lib/libops_mpi.a
	$(MPICPP) $(MPIFLAGS) -DOPS_MPI $(OPS_INC) $(HDF5_INC) $(OPS_LIB)  type_ops.cpp boundary_ops.cpp scheme_ops.cpp flowfield_ops.cpp  model_ops.cpp evolution_ops.cpp lbm2d_ops.cpp -I. ./MPI/$(SEQ_KERNELS) -lops_mpi $(HDF5_LIB) -o lbm2d_mpi

#
# OpenMP version
#
OMP_KERNELS = lbm2d_omp_kernels.cpp

lbm2d_mpi_openmp: Makefile type_ops.cpp boundary_ops.cpp scheme_ops.cpp flowfield_ops.cpp  model_ops.cpp lbm2d_ops.cpp $(OPS_INSTALL_PATH)/lib/libops_mpi.a
	$(MPICPP) $(OMPFLAGS) $(MPIFLAGS) -DOPS_MPI $(OPS_INC) $(HDF5_INC) $(OPS_LIB)  type_ops.cpp boundary_ops.cpp scheme_ops.cpp flowfield_ops.cpp  model_ops.cpp lbm2d_ops.cpp -I. ./MPI_OpenMP/$(OMP_KERNELS) -lops_mpi  $(HDF5_LIB) -o lbm2d_mpi_openmp


lbm2d_openmp: Makefile type_ops.cpp boundary_ops.cpp scheme_ops.cpp flowfield_ops.cpp  model_ops.cpp evolution_ops.cpp lbm2d_ops.cpp $(OPS_INSTALL_PATH)/lib/libops_seq.a
	$(MPICPP) $(OMPFLAGS) $(HDF5) $(MPIFLAGS) $(OPS_INC) $(HDF5_INC) $(OPS_LIB)   type_ops.cpp boundary_ops.cpp scheme_ops.cpp flowfield_ops.cpp  model_ops.cpp evolution_ops.cpp lbm2d_ops.cpp -I. ./MPI_OpenMP/$(OMP_KERNELS) -lops_seq $(HDF5_LIB) -o lbm2d_openmp

#
# Sequential version
#
lbm3d_seq :
	$(CPP) $(CPPFLAGS) ./MPI_OpenMP/lbm3d_cavity_cpu_kernels.cpp boundary.cpp evolution3d_ops.cpp scheme.cpp configuration.cpp model.cpp flowfield_ops.cpp lbm3d_cavity_ops.cpp type.cpp  $(OPS_INC) $(HDF5_INC) $(MPLB_INC) $(OPS_LIB) $(MPLB_LIB) -DOPS_3D -DCPU -D$(LEVEL)  -lops_seq -lops_hdf5_seq  $(HDF5_LIB)  -o lbm3d_seq

lbm2d_seq: Makefile type_ops.cpp boundary_ops.cpp scheme_ops.cpp flowfield_ops.cpp  model_ops.cpp evolution_ops.cpp lbm2d_ops.cpp $(OPS_INSTALL_PATH)/lib/libops_seq.a
	$(MPICPP) $(MPIFLAGS) $(OPS_INC) $(HDF5_INC) $(OPS_LIB) -I. ./MPI/$(SEQ_KERNELS) type_ops.cpp boundary_ops.cpp scheme_ops.cpp flowfield_ops.cpp  model_ops.cpp evolution_ops.cpp lbm2d_ops.cpp -lops_seq  $(HDF5_LIB) -o lbm2d_seq

#
# CUDA version
# $(OMPFLAGS)

lbm2d_cuda: Makefile ./CUDA/lbm2d_kernels_cu.o  type_ops.cpp boundary_ops.cpp scheme_ops.cpp flowfield_ops.cpp  model_ops.cpp evolution_ops.cpp lbm2d_ops.cpp $(OPS_INSTALL_PATH)/lib/libops_cuda.a
	$(MPICPP) $(OMPFLAGS) $(MPIFLAGS) $(CUDA_INC) $(HDF5_INC) $(OPS_INC) $(OPS_LIB) $(CUDA_LIB) type_ops.cpp boundary_ops.cpp scheme_ops.cpp flowfield_ops.cpp  model_ops.cpp evolution_ops.cpp lbm2d_ops.cpp ./CUDA/lbm2d_kernels_cu.o -lcudart -lops_cuda $(HDF5_LIB) -o lbm2d_cuda

lbm2d_mpi_cuda: Makefile ./CUDA/lbm2d_kernels_mpi_cu.o type_ops.cpp boundary_ops.cpp scheme_ops.cpp flowfield_ops.cpp  model_ops.cpp lbm2d_ops.cpp $(OPS_INSTALL_PATH)/lib/libops_mpi_cuda.a
	$(MPICPP) $(OMPFLAGS) $(CPPFLAGS) -DOPS_MPI $(CUDA_INC) $(OPS_INC) $(OPS_LIB) $(CUDA_LIB) type_ops.cpp boundary_ops.cpp scheme_ops.cpp flowfield_ops.cpp  model_ops.cpp lbm2d_ops.cpp ./CUDA/lbm2d_kernels_mpi_cu.o -lcudart -lops_mpi_cuda $(HDF5_LIB) -o lbm2d_mpi_cuda

./CUDA/lbm2d_kernels_cu.o: ./CUDA/lbm2d_kernels.cu Makefile
	$(NVCC) $(VAR) $(INC) $(NVCC_FLAGS) $(CODE_GEN_CUDA) $(OPS_INC) -I. -c -o ./CUDA/lbm2d_kernels_cu.o ./CUDA/lbm2d_kernels.cu

./CUDA/lbm2d_kernels_mpi_cu.o: ./CUDA/lbm2d_kernels.cu Makefile
	$(NVCC) $(VAR) $(INC) $(NVCC_FLAGS) -DOPS_MPI $(CODE_GEN_CUDA) $(OPS_INC) -I. -I$(MPI_INSTALL_PATH)/include \
	-DMPICH_IGNORE_CXX_SEEK -c -o ./CUDA/lbm2d_kernels_mpi_cu.o ./CUDA/lbm2d_kernels.cu

openacc_obj_list = $(shell find OpenACC/ -name "*.cpp" | sed s/\\.cpp/\\.o/g | tr "\n" " ")
openacc_c_obj_list = $(shell find OpenACC/ -name "*_c.c" | sed s/\\.c/\\.o/g | tr "\n" " ")

openacc_mpi_obj_list = $(shell find OpenACC/ -name "*.cpp" | sed s/\\.cpp/\\_mpi.o/g | tr "\n" " ")
openacc_mpi_c_obj_list = $(shell find OpenACC/ -name "*_c.c" | sed s/\\.c/\\_mpi.o/g | tr "\n" " ")


./OpenACC/lbm2d_kernels.o: ./OpenACC/lbm2d_kernels.cpp Makefile
	$(CPP) -acc -ta=tesla:cc35 -ta=tesla:maxregcount:64 $(CPPFLAGS) $(CUDA_ALIGN_FLAG) $(OPS_INC) $(CUDA_INC) \
	-I. OpenACC/lbm2d_kernels.cpp -c -o OpenACC/lbm2d_kernels.o


./OpenACC/lbm2d_kernels_mpi.o: ./OpenACC/lbm2d_kernels.cpp Makefile
	$(MPICPP) -acc -ta=tesla:cc35 -ta=tesla:maxregcount:64 -DOPS_MPI $(MPIFLAGS) $(CUDA_ALIGN_FLAG) $(OPS_INC) $(CUDA_INC) \
	-I. OpenACC/lbm2d_kernels.cpp -c -o OpenACC/lbm2d_kernels_mpi.o

./OpenACC/%_mpi.o: ./OpenACC/%.cpp ./OpenACC/%_c.c Makefile
	@echo "Compiling kernel" $*
	$(MPICC) -acc -ta=tesla:cc35 -ta=tesla:maxregcount:64 $(MPIFLAGS) -DOPS_MPI $(CUDA_ALIGN_FLAG) $(OPS_INC) $(CUDA_INC) \
              -I. OpenACC/$*_c.c -c -o OpenACC/$*_c_mpi.o
	$(MPICPP) -acc -ta=tesla:cc35 -ta=tesla:maxregcount:64 $(MPIFLAGS) -DOPS_MPI $(CUDA_ALIGN_FLAG) $(OPS_INC) $(CUDA_INC) \
        -I. OpenACC/$*.cpp -c -o OpenACC/$*_mpi.o

./OpenACC/%.o: ./OpenACC/%.cpp ./OpenACC/%_c.c Makefile
	@echo "Compiling kernel" $*
	$(CC) -acc -ta=tesla:cc35 -ta=tesla:maxregcount:64 $(CCFLAGS) $(CUDA_ALIGN_FLAG) $(OPS_INC) $(CUDA_INC) \
              -I. OpenACC/$*_c.c -c -o OpenACC/$*_c.o
	$(CPP) -acc -ta=tesla:cc35 -ta=tesla:maxregcount:64 $(CPPFLAGS) $(CUDA_ALIGN_FLAG) $(OPS_INC) $(CUDA_INC) \
        -I. OpenACC/$*.cpp -c -o OpenACC/$*.o

lbm2d_openacc: $(openacc_obj_list) ./OpenACC/lbm2d_kernels.o type_ops.cpp boundary_ops.cpp scheme_ops.cpp flowfield_ops.cpp  model_ops.cpp lbm2d_ops.cpp Makefile $(OPS_INSTALL_PATH)/lib/libops_cuda.a
	$(MPICPP) -acc -ta=tesla:cc35 $(MPIFLAGS) $(OPS_INC) $(OPS_LIB) -DOPS_MPI $(CUDA_INC) $(CUDA_LIB) \
    type_ops.cpp boundary_ops.cpp scheme_ops.cpp flowfield_ops.cpp  model_ops.cpp lbm2d_ops.cpp -I. $(openacc_obj_list) $(openacc_c_obj_list) -lcudart -lops_cuda  $(HDF5_LIB -o lbm2d_openacc


lbm2d_mpi_openacc: $(openacc_mpi_obj_list) ./OpenACC/lbm2d_kernels_mpi.o type_ops.cpp boundary_ops.cpp scheme_ops.cpp flowfield_ops.cpp  model_ops.cpp lbm2d_ops.cpp Makefile $(OPS_INSTALL_PATH)/lib/libops_mpi_cuda.a
	$(MPICPP) -acc -ta=tesla:cc35 $(MPIFLAGS) $(OPS_INC) $(OPS_LIB) $(CUDA_INC) $(CUDA_LIB) -DOPS_MPI \
   type_ops.cpp boundary_ops.cpp scheme_ops.cpp flowfield_ops.cpp  model_ops.cpp lbm2d_ops.cpp -I. $(openacc_mpi_obj_list) $(openacc_mpi_c_obj_list) -lcudart -lops_mpi_cuda  $(HDF5_LIB -o lbm2d_mpi_openacc


#
#OpenCL version
#
lbm2d_opencl: ./OpenCL/lbm2d_opencl_kernels.o lbm2d_ops.cpp lbm2d_kernel.h lbm2d_print_kernel.h lbm2d_copy_kernel.h Makefile $(OPS_INSTALL_PATH)/lib/libops_opencl.a
	$(CPP) $(MPIFLAGS) $(OPS_INC) $(OPS_LIB) $(OPENCL_LIB) lbm2d_ops.cpp ./OpenCL/lbm2d_opencl_kernels.o  -lops_opencl -o lbm2d_opencl

lbm2d_mpi_opencl: ./OpenCL/lbm2d_mpi_opencl_kernels.o lbm2d_ops.cpp lbm2d_kernel.h lbm2d_print_kernel.h lbm2d_copy_kernel.h Makefile $(OPS_INSTALL_PATH)/lib/libops_mpi_opencl.a
	$(MPICPP) $(MPIFLAGS) -DOPS_MPI $(OPS_INC) $(OPS_LIB) $(OPENCL_LIB) lbm2d_ops.cpp ./OpenCL/lbm2d_opencl_kernels.o -lops_mpi_opencl -o lbm2d_mpi_opencl

./OpenCL/lbm2d_opencl_kernels.o: ./OpenCL/lbm2d_opencl_kernels.cpp Makefile
	$(CPP) $(CPPFLAGS) $(VAR) $(OPENCL_INC) $(INC) $(OPS_INC) -I. \
                -c -o ./OpenCL/lbm2d_opencl_kernels.o ./OpenCL/lbm2d_opencl_kernels.cpp

./OpenCL/lbm2d_mpi_opencl_kernels.o: ./OpenCL/lbm2d_opencl_kernels.cpp Makefile
	$(MPICPP) $(CPPFLAGS) -DOPS_MPI $(VAR) $(OPENCL_INC) $(INC) $(OPS_INC) -I. \
                -c -o ./OpenCL/lbm2d_opencl_kernels.o ./OpenCL/lbm2d_opencl_kernels.cpp

#
#error
#
openacc_error:
	@echo "lbm2d_mpi_openacc cannot be built with $(OPS_COMPILER) compilers .. requires an OpenACC compiler (e.g. pgi)"


#
# cleanup
#

clean:
	rm -f lbm2d_dev_seq lbm2d_seq lbm2d_dev_seq lbm2d_mpi lbm2d_dev_mpi lbm2d_openmp lbm2d_mpi_openmp lbm2d_cuda lbm2d_mpi_cuda lbm2d_openacc lbm2d_mpi_openacc lbm2d_opencl lbm2d_mpi_opencl ./CUDA/*.o ./OpenACC/*.o *.o lbm2d_opencl lbm3d_dev_seq lbm3d_seq lbm3d_dev_seq lbm3d_mpi lbm3d_dev_mpi./OpenCL/*.o *.o
