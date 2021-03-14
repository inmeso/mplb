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

TARGETS = fieldblock_dev_mpi fieldblock_dev_seq fieldblock_dev_mpi2d fieldblock_dev_seq2d

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

fieldblock_dev_seq :
	$(MPICPP) $(CPPFLAGS)  $(MAINCPP) field_block.cpp $(MPLB_DIR)/block.cpp $(OPS_INC) $(HDF5_INC) $(MPLB_INC) $(OPS_LIB) $(HDF5_LIB)  -DOPS_3D -lops_hdf5_seq -lops_seq -lhdf5_hl -lhdf5 -lz  -o fieldblock_dev_seq

fieldblock_dev_mpi :
	$(MPICPP) $(CPPFLAGS)  $(MAINCPP) field_block.cpp $(MPLB_DIR)/block.cpp $(OPS_INC) $(HDF5_INC) $(MPLB_INC) $(OPS_LIB) $(MPLB_LIB) $(HDF5_LIB) -DOPS_MPI -DOPS_3D -lops_hdf5_mpi -lops_mpi -lhdf5_hl -lhdf5 -lz  -o fieldblock_dev_mpi

fieldblock_dev_seq2d :
	$(MPICPP) $(CPPFLAGS)  $(MAINCPP) field_block.cpp $(MPLB_DIR)/block.cpp $(OPS_INC) $(HDF5_INC) $(MPLB_INC) $(OPS_LIB) $(HDF5_LIB) -DOPS_2D -lops_hdf5_seq -lops_seq -lhdf5_hl -lhdf5 -lz  -o fieldblock_dev_seq2d

fieldblock_dev_mpi2d :
	$(MPICPP) $(CPPFLAGS)  $(MAINCPP) field_block.cpp $(MPLB_DIR)/block.cpp $(OPS_INC) $(HDF5_INC) $(MPLB_INC) $(OPS_LIB) $(MPLB_LIB) $(HDF5_LIB) -DOPS_MPI -DOPS_2D -lops_hdf5_mpi -lops_mpi -lhdf5_hl -lhdf5 -lz  -o fieldblock_dev_mpi2d


#

clean:
	rm -f ${TARGETS}
