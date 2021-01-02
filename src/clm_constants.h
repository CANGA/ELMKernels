#ifndef CLM_CONSTANTS_HH_
#define CLM_CONSTANTS_HH_

// from landunit_varcon.F90
static const int istsoil    = 1;
static const int istcrop    = 2;
static const int istice     = 3;
static const int istice_mec = 4;
static const int istdlak    = 5;
static const int istwet     = 6;
static const int isturb_MIN = 7;
static const int isturb_tbd = 7;
static const int isturb_hd  = 8;
static const int isturb_md  = 9;
static const int isturb_MAX = 9;
static const int max_lunit  = 9;                               
static const int landunit_name_length = 40;
static const int numurbl = isturb_MAX - isturb_MIN + 1;

//from clm_varctl.F90
static const int subgridflag = 0;

//from clm_varpar.F90
static const int nlevsno = 5; // this should be set by driver
static const int nlevgrnd = 15;
static const int nlevurb = 5; // number of urban layers
static const int numrad = 2; // number of solar radiation bands: vis, nir

//from column_varcon.F90
static const int icol_roof        = isturb_MIN*10 + 1;
static const int icol_sunwall     = isturb_MIN*10 + 2;
static const int icol_shadewall   = isturb_MIN*10 + 3;
static const int icol_road_imperv = isturb_MIN*10 + 4;
static const int icol_road_perv   = isturb_MIN*10 + 5;

//from shr_const_mod.F90
static const double SHR_CONST_BOLTZ   = 1.38065e-23;  // Boltzmann's constant ~ J/K/molecule
static const double SHR_CONST_AVOGAD  = 6.02214e26;  // Avogadro's number ~ molecules/kmole
static const double SHR_CONST_MWWV    = 18.016;       // molecular weight water vapor
static const double SHR_CONST_RGAS    = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ;      // Universal gas constant ~ J/K/kmole
static const double SHR_CONST_RWV     = SHR_CONST_RGAS/SHR_CONST_MWWV;          // Water vapor gas constant ~ J/K/kg


//from clm_varcon.F90
static const double snw_rds_min = 54.526;               //minimum allowed snow effective radius (also "fresh snow" value) [microns]
static const double zlnd = 0.01;                        // Roughness length for soil [m]
static const double zsno = 0.0024;                      // Roughness length for snow [m]
static const double spval = 1.0e36;                     // special value for real data
static const int ispval = -9999;                        // special value for int data 
static const double tfrz = 273.15;                      // freezing temperature [K]
static const double denice = 0.917e3;                   // density of ice             ~ kg/m^3
static const double denh2o = 1.000e3;                   // // density of fresh water     ~ kg/m^3
static const double hvap   = 2.501e6;                   // Latent heat of evap for water [J/kg]
static const double hfus   = 3.337e5;                   // Latent heat of fusion for ice [J/kg]
static const double hsub   = hvap + hfus;               // Latent heat of sublimation    [J/kg]
static const double grav   = 9.80616;                   // gravity constant [m/s2]
static const double roverg = SHR_CONST_RWV/grav*1000.0; // Rw/g constant = (8.3144/0.018)/(9.80616)*1000. mm/K
static const double vkc    = 0.4;                       // von Karman constant [-]
static const double cpair  = 1.00464e3;                 // specific heat of dry air   ~ J/kg/K


static const double ELM_PI  = 3.14159265358979323846;   // pi


static const bool perchroot = false;     // calculate root function only in unfrozen soil, based on instantaneous temp
static const bool perchroot_alt = false; // calculate root function only in active layer - estimated from past 2 years

#endif
