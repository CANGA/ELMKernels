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
static const int nlevsno = 5;

//from column_varcon.F90
static const int icol_roof        = isturb_MIN*10 + 1;
static const int icol_sunwall     = isturb_MIN*10 + 2;
static const int icol_shadewall   = isturb_MIN*10 + 3;
static const int icol_road_imperv = isturb_MIN*10 + 4;
static const int icol_road_perv   = isturb_MIN*10 + 5;

#endif
