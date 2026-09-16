Data files and their corresponding energy-budget terms
=======================================================


1. Naming convention and common variables
-----------------------------------------

Experiment/product suffixes

  eg       = Reg-EG: ERA5 surface forcing + monthly GLORYS initial/lateral
             boundary forcing.
  jg       = Reg-JG: JRA55-do surface forcing + monthly GLORYS initial/lateral
             boundary forcing.
  jl       = Reg-JL: JRA55-do surface forcing + monthly LICOM3 initial/lateral
             boundary forcing.
  glorys   = the GLORYS global reanalysis result plotted as a dashed line.
  licom3   = the global LICOM3 result plotted as a dashed line.
  dailyBC  = Reg-EG-dailyBC, identical to Reg-EG except that GLORYS lateral
             boundary fields are updated daily rather than monthly. These
             files correspond to the Figure 11 sensitivity experiment, not
             to the five curves in Figures 7 and 8.


2. Eight boundary-flux terms used in Figures 7 and 8
----------------------------------------------------

2.1 edgeumke-*  -> ADV in MKE

2.2 edgeueke-*  -> ADV in EKE

2.3 edgeape-*  -> ADV in MPE

2.4 edgeeape-*  -> ADV in EPE

2.5 edgeuaemke-*  -> CED in KE

2.6 edgerhoaaa-*  -> CED in PE

2.7 edgeup-*  -> PW in MKE

2.8 edgeuapa-*  -> PW in EKE


3. Surface-generation files used for Tables 2 and 3
----------------------------------------------------

Files for each regional experiment

  G-potdens-rhoN0-2014-2018-EG-100_148E-15S_45N.nc.zip
    -> spatial fields for Reg-EG.
  G-potdens-rhoN0-2014-2018-JG-100_148E-15S_45N.nc.zip
    -> spatial fields for Reg-JG.
  G-potdens-rhoN0-2014-2018-JL-100_148E-15S_45N.nc.zip
    -> spatial fields for Reg-JL.

  G-potdens-rhoN0-2014-2018-EG-100_148E-15S_45N.dat.zip
  G-potdens-rhoN0-2014-2018-JG-100_148E-15S_45N.dat.zip
  G-potdens-rhoN0-2014-2018-JL-100_148E-15S_45N.dat.zip
    -> text summaries of the domain-integrated values stored in the matching
       NetCDF files. Values are reported in W and TW.

Variables and article terms in each G-*.nc file

  gkm       -> G_MKE: mean wind-stress work; Table 2, GMKE.
  gke       -> G_EKE: covariance of wind-stress and surface-velocity
               perturbations; Table 2, GEKE.
  gpm_heat  -> heat-flux contribution to G_MPE.
  gpm_fw    -> freshwater-flux contribution to G_MPE.
  gpm       -> total G_MPE = gpm_heat + gpm_fw; Table 3, GMPE.
  gpe_heat  -> heat-flux contribution to G_EPE.
  gpe_fw    -> freshwater-flux contribution to G_EPE.
  gpe       -> total G_EPE = gpe_heat + gpe_fw; Table 3, GEPE.

The remaining variables (rho_ref, n0, alpha0, beta0, rho_star_mean, J_mean,
Gs_mean, cell_area, and diagnostic surface-density fields) are intermediate
or supporting quantities used to calculate the G terms. They are not separate
energy-budget terms in Tables 2 or 3.


4. Reference-density files
--------------------------

Each rhoref-*.txt file contains a 55-level reference density profile rho_r(z)
in kg m-3. The scripts use rho_r to calculate the density anomaly rho_a and
the stratification N^2. These profiles are therefore auxiliary inputs, not
energy-budget terms themselves.

  rhoref-eg.txt       -> reference density for Reg-EG.
  rhoref-jg.txt       -> reference density for Reg-JG.
  rhoref-jl.txt       -> reference density for Reg-JL.
  rhoref-glorys.txt   -> reference density for GLORYS.
  rhoref-dailyBC.txt  -> reference density for Reg-EG-dailyBC.

They support the calculation of:
  - ADV in MPE (edgeape-*),
  - ADV in EPE (edgeeape-*),
  - CED in PE (edgerhoaaa-*), and
  - the pressure anomaly used in PW calculations (especially edgeup-*).


Data are available at https://doi.org/10.5281/zenodo.22753955