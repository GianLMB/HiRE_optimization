f77files = $(wildcard *.f)
f90files = $(wildcard *.f90)
OBJ= $(f77files:.f=.o) $(f90files:.f90=.o)
MDDRIVER= mddriver
IMD_LINK_GLOB= -I./interactor/ ./interactor/Interactor.o ./interactor/InteractorMDDriver.o ./interactor/OPEP_wrapper.o -I$(MDDRIVER)/include -L$(MDDRIVER)/build/protocol $(MDDRIVER)/build/protocol/libmddriver.a -DIMD_SIMULATION 
CODEVER = test

ifeq ($(arch),gnutest)
  F90 = gfortran -ffree-form
  F77 = gfortran -ffixed-form 
  LINK = $(F90)
  CFLAGS = -Wall -pedantic -O0 -ggdb
  PFLAGS = $(CFLAGS) -x f95-cpp-input
  LIBS = -pthread
  IMD_LINK= $(IMD_LINK_GLOB) -lstdc++
  EXEC = simulator
endif

ifeq ($(arch),gnu2)
  F90 = gfortran -ffree-form -mcmodel=large
  F77 = gfortran -ffixed-form -mcmodel=large 
  LINK = $(F90)
  CFLAGS = -Wall -O3 -flto
  PFLAGS = $(CFLAGS) -x f95-cpp-input
  LIBS = -pthread
  IMD_LINK= $(IMD_LINK_GLOB) -lstdc++
  EXEC = simulator
endif

ifeq ($(arch),nagfor)
  F90 = nagfor -free
  F77 = nagfor -fixed
  LINK = $(F90) -pg
  CFLAGS = -O
  PFLAGS = $(CFLAGS) -fpp
  LIBS = -thread_safe
  IMD_LINK= $(IMD_LINK_GLOB) -lstdc++
  EXEC = simulator
endif

ifeq ($(arch),osx)
  F90 = gfortran -ffree-form 
  F77 = gfortran -ffixed-form 
  LINK = $(F90)# -fuse-linker-plugin
  CFLAGS = -Wall -O3 #-flto
  PFLAGS = $(CFLAGS) -x f95-cpp-input
  LIBS = -pthread
  EXEC = simulator
  IMD_LINK= $(IMD_LINK_GLOB) -lc++
endif

ifeq ($(arch),gnumpi)
  F90 = mpif90 -ffree-form
  F77 = mpif90 -ffixed-form
  LINK = mpif90
  CFLAGS = -cpp -DMPI -O2 #-pg #-ggdb #-fopenmp
  PFLAGS = $(CFLAGS) -ffree-form
  LIBS =
  IMD_LINK= $(IMD_LINK_GLOB) -lstdc++
  EXEC = simulatorMPI
endif

ifeq ($(arch),intelmpi)
  F90 = mpif90 -free
  F77 = mpif90 -fixed
  LINK = mpif90
  CFLAGS = -cpp -DMPI -O2 -g -traceback  #-check bounds -mcmodel=large #-pg #-ggdb #-fopenmp                                               
  PFLAGS = $(CFLAGS) -free 
  LIBS =
  IMD_LINK= $(IMD_LINK_GLOB) -lstdc++
  EXEC = simulatorMPI
endif

ifeq ($(arch),intel)
  F90 = ifort -free
  F77 = ifort -fixed
  LINK = $(F90)
  CFLAGS = -Wall -O3 -g -traceback #-Wall -O3 -flto
  PFLAGS = $(CFLAGS) -free 
  LIBS = 
  IMD_LINK= $(IMD_LINK_GLOB) -lstdc++
  EXEC = simulator
endif

ifeq ($(arch),lucifer)
  F90 = gfortran -ffree-form 
  F77 = gfortran -ffixed-form 
  LINK = $(F90) -L/shared/libraries/lapack/3.4.2/gnu/lib/
  CFLAGS = -Wall -O3 -g -pg -cpp -lpthread
  PFLAGS = $(CFLAGS) 
  LIBS =
  IMD_LINK= $(IMD_LINK_GLOB) -lstdc++
  EXEC = simulator
endif

ifeq ($(arch),luciferMPI)
  F90 = mpif90 -ffree-form
  F77 = mpif90 -ffixed-form
  LINK = mpif90 -L/shared/libraries/lapack/3.4.2/gnu/lib/
  CFLAGS = -cpp -DMPI -O2 #-pg #-ggdb #-fopenmp
  PFLAGS = $(CFLAGS) -ffree-form
  LIBS =
  IMD_LINK= $(IMD_LINK_GLOB) -lstdc++
  EXEC = simulatorMPI
endif


ifeq ($(arch),'')
  error:
        @echo 'ERROR: Please specify the right architecture (see Makefile for choice)'
endif

all:$(EXEC)

gnutest:
	@make arch=gnutest
gnu2:
	@make arch=gnu2
nagfor:
	@make arch=nagfor
osx:
	@make arch=osx
gnumpi:
	@make arch=gnumpi
intelmpi:
	@make arch=intelmpi
lucifer :
	@make arch=lucifer
luciferMPI :
	@make arch=luciferMPI
intel:
	@make arch=intel

$(EXEC): $(OBJ)
#	$(LINK) -Lxdrfile-1.1b/src/ $(CFLAGS) $(OBJ) -o $(EXEC)  $(LIBS) 
#	$(LINK) -L./xdrfile-1.1b/src/ $(CFLAGS) $(OBJ) ./xdrfile-1.1b/src/xdrfile.o ./xdrfile-1.1b/src/xdrfile_xtc.o -o $(EXEC)  $(LIBS) 
	$(LINK) $(IMD_LINK) -L./xdrfile-1.1b/src/ $(CFLAGS) -o $(EXEC)_$(CODEVER)  $(OBJ) $(LIBS) ./xdrfile-1.1b/src/xdrfile.o ./xdrfile-1.1b/src/xdrfile_xtc.o

%.o: %.c
	$(CC) -c  -o $@ $^
%.o: %.f
	$(F77) $(CFLAGS) -c  -o $@ $^
%.o: %.f90
	$(F90) $(PFLAGS) -c  -o   $@ $*.f90 -DIMD_SIMULATION

# Dependencies
#

calcforce.o : calcforce.f90 setup_RNA.o forces_RNA.o
defs.o : defs.f90 
geometric_corrections.o : geometric_corrections.f90 
matrix_mod.o : matrix_mod.f90 defs.o 
md_main_module.o : md_main_module.f90 geometric_corrections.o defs.o 
md_utils.o : md_utils.f90 defs.o geometric_corrections.o md_main_module.o 
md_initialise.o : md_initialise.f90 restart.o md_utils.o md_main_module.o readtop.o
md_statistics.o : md_statistics.f90 geometric_corrections.o md_main_module.o 
write_to_file.o : write_to_file.f90 geometric_corrections.o defs.o 
titration.o : titration.f90
md.o : md.f90 md_utils.o md_statistics.o md_utils.o restart.o md_initialise.o write_to_file.o geometric_corrections.o md_main_module.o random.o defs.o SAXS_scoring.o
md_read_parameters.o : md_read_parameters.f90 md_main_module.o 
SAXS_scoring.o : SAXS_scoring.f90 defs.o hash_int.o geometric_corrections.o
forces_protein.o: forces_protein.f defs.o
intvar.o: intvar.f90 defs.o  setup_RNA.o
forces_RNA.o : forces_RNA.f90 defs.o setup_RNA.o geometric_corrections.o SAXS_scoring.o md_main_module.o
mdmin.o : mdmin.f90 geometric_corrections.o defs.o 
read_parameters.o : read_parameters.f90 md_main_module.o random.o defs.o SAXS_scoring.o
readpar.o : readpar.f defs.o
readtop.o : readtop.f defs.o
restart.o : restart.f90 geometric_corrections.o md_main_module.o fileio.o
setup_RNA.o : setup_RNA.f fileio.o
simulateur.o : simulateur.f90 defs.o md_main_module.o restart.o SAXS_scoring.o geometric_corrections.o md.o
utils_text.o : utils_text.f90 

clean:
	@rm -rf *.o  *.mod  simulateur simulateurMPI

