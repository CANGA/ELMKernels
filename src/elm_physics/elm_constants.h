#pragma once

namespace ELM {

// from SoilStateType
static const double smpmin = -1.e8; // restriction for min of soil potential (mm)

// from CanopyHydrology() - move to CanopyHydrology header when ready
// set shape factor for accumulation of snow
static const double accum_factor = 0.1;

// from landunit_varcon.F90
static const int istsoil = 1;
static const int istcrop = 2;
static const int istice = 3;
static const int istice_mec = 4;
static const int istdlak = 5;
static const int istwet = 6;
static const int isturb_MIN = 7;
static const int isturb_tbd = 7;
static const int isturb_hd = 8;
static const int isturb_md = 9;
static const int isturb_MAX = 9;
static const int max_lunit = 9;
static const int landunit_name_length = 40;
static const int numurbl = isturb_MAX - isturb_MIN + 1;

// from pft_varcon.F90
static const int noveg = 0;
static const int ndllf_evr_tmp_tree = 1;
static const int ndllf_evr_brl_tree = 2;
static const int ndllf_dcd_brl_tree = 3;
static const int nbrdlf_evr_trp_tree = 4;
static const int nbrdlf_evr_tmp_tree = 5;
static const int nbrdlf_dcd_trp_tree = 6;
static const int nbrdlf_dcd_tmp_tree = 7;
static const int nbrdlf_dcd_brl_tree = 8;
static const int nbrdlf_evr_shrub = 9;
static const int nbrdlf_dcd_tmp_shrub = 10;
static const int nbrdlf_dcd_brl_shrub = 11;
static const int nc3_arctic_grass = 12;
static const int nc3_nonarctic_grass = 13;
static const int nc4_grass = 14;
static const int nc3crop = 15;
static const int nc3irrig = 16;
static const int ncorn = 17;
static const int ncornirrig = 18;
static const int nscereal = 19;
static const int nscerealirrig = 20;
static const int nwcereal = 21;
static const int nwcerealirrig = 22;
static const int nsoybean = 23;
static const int nsoybeanirrig = 24;

// from clm_varctl.F90
static const int subgridflag = 1;
static const double co2_ppmv = 355.0; // atmospheric CO2 molar ratio (by volume) (umol/mol) - this is really low

// from clm_varpar.F90
static const int nlevsno = 5; // this should be set by driver
static const int nlevgrnd = 15;
static const int nlevurb = 5;  // number of urban layers
static const int numrad = 2;   // number of solar radiation bands: vis, nir
static const int nlevcan = 1;  // number of leaf layers in canopy layer
static const int nlevsoi = 10; // number of soil layers (hydrologically active)
static const int nlevbed = 15; // number of layers to bedrock (hydrologically inactive below)

static const int mxpft = 25; // maximum number of PFT's for any mode
static const int numveg = 17; // number of veg types (without specific crop)
static const bool use_crop = false;

constexpr int numpft = (use_crop) ? mxpft : numveg; // number of pfts - (use_crop ? 25 : 17)


// from column_varcon.F90
static const int icol_roof = isturb_MIN * 10 + 1;
static const int icol_sunwall = isturb_MIN * 10 + 2;
static const int icol_shadewall = isturb_MIN * 10 + 3;
static const int icol_road_imperv = isturb_MIN * 10 + 4;
static const int icol_road_perv = isturb_MIN * 10 + 5;

static const double ELM_PI = 3.14159265358979323846;   // pi
static const double ELM_BOLTZ = 1.38065e-23;           // Boltzmann's constant ~ J/K/molecule
static const double ELM_AVOGAD = 6.02214e26;           // Avogadro's number ~ molecules/kmole
static const double ELM_MWWV = 18.016;                 // molecular weight water vapor
static const double ELM_RGAS = ELM_AVOGAD * ELM_BOLTZ; // Universal gas constant ~ J/K/kmole
static const double ELM_RWV = ELM_RGAS / ELM_MWWV;     // Water vapor gas constant ~ J/K/kg

static const double ELM_STEBOL = 5.67e-8;              // Stefan-Boltzmann constant ~ W/m^2/K^4
// need to replace 'sb' with 'ELM_STEBOL'
//static const double sb = 5.67e-8;        // Stefan-Boltzmann constant ~ W/m^2/K^4

// from clm_varcon.F90
static const double snw_rds_min = 54.526; // minimum allowed snow effective radius (also "fresh snow" value) [microns]
static const double zlnd = 0.01;          // Roughness length for soil [m]
static const double zsno = 0.0024;        // Roughness length for snow [m]
static const double spval = 1.0e36;       // special value for real data
static const int ispval = -9999;          // special value for int data
inline constexpr double tfrz = 273.15;        // freezing temperature [K]
static const double denice = 0.917e3;     // density of ice             ~ kg/m^3
static const double denh2o = 1.000e3;     // // density of fresh water     ~ kg/m^3
static const double hvap = 2.501e6;       // Latent heat of evap for water [J/kg]
static const double hfus = 3.337e5;       // Latent heat of fusion for ice [J/kg]
static const double hsub = hvap + hfus;   // Latent heat of sublimation    [J/kg]
static const double grav = 9.80616;       // gravity constant [m/s2]
static const double roverg = ELM_RWV / grav * 1000.0; // Rw/g constant = (8.3144/0.018)/(9.80616)*1000. mm/K
static const double vkc = 0.4;                        // von Karman constant [-]
static const double cpair = 1.00464e3;                // specific heat of dry air   ~ J/kg/K
static const double csoilc = 0.004;      // Drag coefficient for soil under canopy [-]
static const double alpha_aero = 1.0; // constant for aerodynamic parameter weighting
static const double tlsai_crit = 2.0;    // critical value of elai+esai for which aerodynamic parameters are maximum
static const double sb = 5.67e-8;        // Stefan-Boltzmann constant ~ W/m^2/K^4
static const double h2osno_max = 1000.0; // max allowed snow thickness (mm H2O)
static const double bdsno = 250.0;       // bulk density snow (kg/m**3)
static const double ELM_MWDAIR = 28.966; // molecular weight dry air ~ kg/kmole
static const double rair = ELM_RGAS / ELM_MWDAIR; // Dry air gas constant     ~ J/K/kg
static const double o2_molar_const = 0.209;       // constant atmospheric O2 molar ratio (mol/mol)

static const bool perchroot = false;     // calculate root function only in unfrozen soil, based on instantaneous temp
static const bool perchroot_alt = false; // calculate root function only in active layer - estimated from past 2 years





// should be read from file - hardwired for now
static const double organic_max = 130.0; // Organic matter content where soil is assumed to act like peat for diffusion ~ kg/m3
static const double secspday = 86400.0; // seconds per day

} // namespace ELM

namespace ELM::ELMconstants {
inline constexpr double TFRZ = 273.15;        // freezing temperature [K]

inline constexpr double PI = 3.14159265358979323846;   // pi
inline constexpr double BOLTZ = 1.38065e-23;           // Boltzmann's constant ~ J/K/molecule
inline constexpr double AVOGAD = 6.02214e26;           // Avogadro's number ~ molecules/kmole
inline constexpr double MWWV = 18.016;                 // molecular weight water vapor
inline constexpr double RGAS = AVOGAD * BOLTZ; // Universal gas constant ~ J/K/kmole
inline constexpr double RWV = RGAS / MWWV;     // Water vapor gas constant ~ J/K/kg

inline constexpr double STEBOL = 5.67e-8;              // Stefan-Boltzmann constant ~ W/m^2/K^4
inline constexpr double MWDAIR = 28.966; // molecular weight dry air ~ kg/kmole
inline constexpr double RAIR = RGAS / MWDAIR; // Dry air gas constant     ~ J/K/kg


inline constexpr double O2_MOLAR_CONST = 0.209;       // constant atmospheric O2 molar ratio (mol/mol)
inline constexpr double CO2_PPMV = 355.0; // atmospheric CO2 molar ratio (by volume) (umol/mol) - this is really low


inline constexpr int numrad = 2;   // number of solar radiation bands: vis, nir

enum class AtmForcType { TBOT, PBOT, QBOT, RH, FLDS, FSDS, PREC, WIND, ZBOT };
} // namespace ELM::constants
