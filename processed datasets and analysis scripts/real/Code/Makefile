INCLDIR  = -I.  -I/usr/include -I/public/software/mathlib/libs-intel/netcdf/4.4.1/include -I/opt/hpc/software/mpi/hpcx/v2.7.4/intel-2017.5.239/include
CPPFLAGS = -P -traditional 
CPP      = /usr/bin/cpp $(CPPFLAGS) $(INCLDIR)
NLIB =   -L/opt/hpc/software/mpi/hpcx/v2.7.4/intel-2017.5.239/lib -L/public/software/mathlib/libs-intel/netcdf/4.4.1/lib -lnetcdff -lnetcdf 
FFLAGS = -O2 -r8 -i4 -fp-model source -convert big_endian -assume byterecl -ftz -traceback -no-vec
#FFLAGS = -O3 -shared-intel -r8 -i4 -convert big_endian -assume byterecl -fp-speculation=safe -no-vec
#FFLAGS = -O3 -shared-intel -r8 -i4 -convert big_endian -assume byterecl -fp-speculation=safe -no-vec 
#FFLAGS = -O3 -mcmodel=large -shared-intel -r8 -i4 -convert big_endian -assume byterecl -fp-model precise -fp-speculation=safe -no-vec  -traceback
#FFLAGS = -O3 -mcmodel=large -shared-intel -r8 -i4 -convert big_endian -assume byterecl -fp-model precise -fp-speculation=safe -no-vec  -traceback
FC= mpif90 -m64 -mcmodel=medium $(FFLAGS)
OBJS   = precision_mod.o  param_mod.o pconst_mod.o  pmix_mod.o isopyc_mod.o cdf_mod.o\
tracer_mod.o dyn_mod.o work_mod.o output_mod.o forc_mod.o msg_mod.o sw_mod.o canuto_mod.o openbc_mod.o tide_mod.o\
ppm_mod.o prm_mod.o accumm.o   diag_mod.o  const.o    filter.o   invtri.o   k1_3.o     rdriver.o   vinteg.o\
addps.o    convadj.o  grids.o    isoadv.o   k2_3.o     readyc.o   smuvh.o    yy00.o  openbc_pre.o openbc_code.o\
openbc_volcons.o openbc_sponge.o openbc_uv_nudge.o openbc_ts_nudge.o smth.o\
barotr.o   dens.o     icesnow.o  isoflux.o  k3_123.o   readyt.o   setbcx.o tide_code.o\
bclinc.o   density.o  inirun.o   isopyc.o   licom.o     tracer.o\
energy.o   intfor.o   isopyi.o   mm00.o     upwell.o  \
ssave-cdf.o smag.o  boundary.o global_to_local.o \
exchange.o local_to_global.o clock0f.o run_time.o sw_absor.o thermalexpansion.o canuto_new.o canuto_ini.o core_daily.o
#icesnow1.o comqd.o diag_mod.o
###################################################################################
licom2: $(OBJS)
	$(FC)  -o licom2  $(OBJS) $(NLIB)

.SUFFIXES : .F90
.F90.o:
	$(CPP) $< > $*.f90  ;  $(FC) -c  $*.f90  ;  mv  $*.f90 $*.i
.c.o:
	cc  -c $<
clean:
	rm -f *.o *.i *.il *.f90 *.f *.mod *.d licom2
