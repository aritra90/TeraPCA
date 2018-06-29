MKL_ROOT = /apps/rhel6/intel/mkl
MKL_LIBROOT = $(MKL_ROOT)/lib/intel64
MKL_LIB = -L$(MKL_LIBROOT) -lmkl_blas95_lp64 -lmkl_lapack95_lp64  -Wl,--start-group  $(MKL_LIBROOT)/libmkl_intel_lp64.a $(MKL_LIBROOT)/libmkl_intel_thread.a $(MKL_LIBROOT)/libmkl_core.a $(MKL_LIBROOT)/libmkl_sequential.a -Wl,--end-group -qopenmp -lpthread -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -liomp5 -lpthread -lm

COMP = icpc -mkl -Wall -openmp #-DMKL_ILP64 -I${MKLROOT}/include

SOURCE = terapca.cpp
EXE = TeraPCA.exe
OBJ = terapca.o
COMP_MPI = mpiicpc
SOURCE_MPI = test7_mpi.cpp 
EXE_MPI = test_mpi.exe
OBJ_MPI = test7_mpi.o

test:	$(OBJ)
	$(COMP) $(SOURCE) utilities.c gaussian.c gennorm.c methods.cpp io.c -o $(EXE) $(MKL_LIB)

clean:
	rm *.exe *.o *~
