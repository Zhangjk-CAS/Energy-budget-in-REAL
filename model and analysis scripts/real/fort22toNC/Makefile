#INCLDIR  = -I.  -I/usr/include -I/opt/intel/impi/4.1.0/include -I/usr/local/mpi/include -I/opt/soft/netcdf/c_4.3.3.1_cxx_4.2.1_fortran_4.4.2/2/include
#INCLDIR  = -I.  -I/usr/include -I/THL7/software/intel/compilers_and_libraries_2016.3.210/linux/mpi/intel64/include -I/THL7/home/liuhailong/netcdf-3.6.0-x64/include
#NLIB =   -L/opt/intel/impi/4.1.0/lib64  -L/opt/soft/netcdf/c_4.3.3.1_cxx_4.2.1_fortran_4.4.2/2/lib -lnetcdf -lnetcdff
#NLIB =   -L/THL7/software/intel/compilers_and_libraries_2016.3.210/linux/mpi/intel64/lib  -L/THL7/home/liuhailong/netcdf-3.6.0-x64/lib -lnetcdf
#FFLAGS = -O2 -mcmodel=large -shared-intel -convert big_endian -assume byterecl -fp-model precise -fp-speculation=safe -no-vec 
#FC= /opt/intel/impi/4.1.0/bin64/mpiifort $(FFLAGS)
INCLDIR  = -I.  -I/usr/include -I/usr/local/include -I/opt/hpc/software/compiler/intel/intel-compiler-2017.5.239/compiler/include /public/software/mathlib/libs-intel/netcdf/4.4.1/include
NLIB = -L/opt/hpc/software/compiler/intel/intel-compiler-2017.5.239/compiler/lib/intel64 -L/public/software/mathlib/libs-intel/netcdf/4.4.1/lib -lnetcdf -lnetcdff
FFLAGS = -O2 -mcmodel=large -shared-intel -convert big_endian -assume byterecl
#FFLAGS = -O2 -r8 -i4 -mcmodel=large -shared-intel -convert big_endian -assume byterecl
#NLIB = -L/THL8/software/intel/compilers_and_libraries_2016.3.210/linux/mpi/intel64/lib -L/THL8/software/netcdf/4.4-icc16-mvapich2/lib -lnetcdf -lnetcdff
FC= ifort $(FFLAGS)
OBJS   = mk_data_lhl1.o
#icesnow1.o comqd.o diag_mod.o
###################################################################################
mk_data_lhl: $(OBJS)
	$(FC)  -o mk_data_lhl1  $(OBJS) $(NLIB)

.SUFFIXES : .f90
.f90.o:
	$(FC) -c  $*.f90  
.c.o:
	cc  -c $<
clean:
	rm -f *.o *.i *.il *.mod *.d mk_data_lhl1
