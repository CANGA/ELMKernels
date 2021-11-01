

namespace ELM {
namespace SnowSnicar {

constexpr int numrad_snw = 5; // number of spectral bands used in snow model [nbr]
constexpr int nir_bnd_bgn = 1; // first band index in near-IR spectrum [idx]
constexpr int nir_bnd_end = 4; // ending near-IR band index [idx]
constexpr int sno_nbr_aer = 8; // number of aerosol species in snowpack
constexpr double min_snw = 1.0e-30; // minimum snow mass required for SNICAR RT calculation [kg m-2]
constexpr int ngmax = 8; // gaussian integration index

constexpr int idx_bc_nclrds_min = 1; // minimum index for BC particle size in optics lookup table
constexpr int idx_bc_nclrds_max = 10; // maximum index for BC particle size in optics lookup table
constexpr int idx_bcint_icerds_min = 1; // minimum index for snow grain size in optics lookup table for within-ice BC
constexpr int idx_bcint_icerds_max = 8; // maximum index for snow grain size in optics lookup table for within-ice BC
constexpr int snw_rds_max_tbl = 1500;// maximum effective radius defined in Mie lookup table [microns]
constexpr int snw_rds_min_tbl = 30; // minimium effective radius defined in Mie lookup table [microns]
constexpr double snw_rds_max = 1500.0; // maximum allowed snow effective radius [microns]
constexpr double snw_rds_refrz  = 1000.0; // effective radius of re-frozen snow [microns]

constexpr int idx_Mie_snw_mx = 1471; // number of effective radius indices used in Mie lookup table [idx]
constexpr int idx_T_max      = 11;   // maxiumum temperature index used in aging lookup table [idx]
constexpr int idx_T_min      = 1;    // minimum temperature index used in aging lookup table [idx]
constexpr int idx_Tgrd_max   = 31;   // maxiumum temperature gradient index used in aging lookup table [idx]
constexpr int idx_Tgrd_min   = 1;    // minimum temperature gradient index used in aging lookup table [idx]
constexpr int idx_rhos_max   = 8;    // maxiumum snow density index used in aging lookup table [idx]
constexpr int idx_rhos_min   = 1;    // minimum snow density index used in aging lookup table [idx]


// Gaussian integration angle and coefficients for diffuse radiation
constexpr double difgauspt[8] // gaussian angles (radians)
= { 0.9894009,  0.9445750,
    0.8656312,  0.7554044,
    0.6178762,  0.4580168,
    0.2816036,  0.0950125 };
constexpr double difgauswt[8] //gaussian weights
= { 0.0271525,  0.0622535,
    0.0951585,  0.1246290,
    0.1495960,  0.1691565,
    0.1826034,  0.1894506 };

// constants used in algorithm
constexpr double c0      = 0.0;
constexpr double c1      = 1.0;
constexpr double c3      = 3.0;
constexpr double c4      = 4.0;
constexpr double c6      = 6.0;
constexpr double cp01    = 0.01;
constexpr double cp5     = 0.5;
constexpr double cp75    = 0.75;
constexpr double c1p5    = 1.5;
constexpr double trmin   = 0.001;
constexpr double argmax  = 10.0; // maximum argument of exponential

// cconstant coefficients used for SZA parameterization
constexpr double sza_a0 =  0.085730;
constexpr double sza_a1 = -0.630883;
constexpr double sza_a2 =  1.303723;
constexpr double sza_b0 =  1.467291;
constexpr double sza_b1 = -3.338043;
constexpr double sza_b2 =  6.807489;
constexpr double puny   =  1.0e-11;
constexpr double mu_75  =  0.2588; // cosine of 75 degree

constexpr exp_min = exp(-argmax);

// shorter name for pi
constexpr double pi = ELM_PI;
// always use Delta approximation for snow
constexpr int DELTA = 1;




} // namespace ELM
} // namespace SnowSnicar
