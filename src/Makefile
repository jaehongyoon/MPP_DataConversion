define colorecho
	@tput setaf 2
	@echo $1
	@tput sgr0
endef

FC     = ifort
MPIFC  = mpiifort
FFLAGS = -O3 -xAVX -openmp
#FFLAGS = -O0 -xAVX -fp-model precise -check bounds -check uninit -check pointers -traceback -g -fpe0
LIB    = -L/opt/intel/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread

all: test_mat_conv.x test_mat_redist.x

test_mat_conv.x: test_mat_conv.o mat_conv.o mat_util.o
	$(call colorecho,"Linking..")
	$(FC) $(FFLAGS) -o $@ $^ 

test_mat_redist.x: test_mat_redist.o mat_redist.o mat_util.o
	$(call colorecho,"Linking..")
	$(MPIFC) $(FFLAGS) -o $@ $^ $(LIB)

test_mat_conv.o: test_mat_conv.f90 mat_conv.o mat_util.o
	$(call colorecho,"Compiling..")
	$(FC) $(FFLAGS) -c $<

test_mat_redist.o: test_mat_redist.f90 mat_redist.o mat_util.o
	$(call colorecho,"Compiling..")
	$(MPIFC) $(FFLAGS) -c $<

mat_redist.o: mat_redist.f90 mat_util.o
	$(call colorecho,"Compiling..")
	$(MPIFC) $(FFLAGS) -c $<

mat_conv.o: mat_conv.f90
	$(call colorecho,"Compiling..")
	$(FC) $(FFLAGS) -c $<

mat_util.o: mat_util.f90
	$(call colorecho,"Compiling..")
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod *.x
