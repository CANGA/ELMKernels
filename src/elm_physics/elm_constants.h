#pragma once

namespace ELM::constants {
    
static constexpr double TFRZ{273.15}; // freezing temperature [K]
static constexpr double ELM_PI{3.14159265358979323846}; // pi
static constexpr double BOLTZ{1.38065e-23};         // Boltzmann's constant ~ J/K/molecule
static constexpr double AVOGAD{6.02214e26};         // Avogadro's number ~ molecules/kmole
static constexpr double MWWV{18.016};               // molecular weight water vapor
static constexpr double RGAS{AVOGAD * BOLTZ};       // Universal gas constant ~ J/K/kmole
static constexpr double RWV{RGAS / MWWV};           // Water vapor gas constant ~ J/K/kg
static constexpr double STEBOL{5.67e-8};     // Stefan-Boltzmann constant ~ W/m^2/K^4
static constexpr double MWDAIR{28.966};      // molecular weight dry air ~ kg/kmole
static constexpr double RAIR{RGAS / MWDAIR}; // Dry air gas constant     ~ J/K/kg
static constexpr double O2_MOLAR_CONST{0.209}; // constant atmospheric O2 molar ratio (mol/mol)
static constexpr double CO2_PPMV{355.0}; // atmospheric CO2 molar ratio (by volume) (umol/mol) - this is really low

} // namespace ELM::ELMconstants

namespace ELM {

enum class AtmForcType { TBOT, PBOT, QBOT, RH, FLDS, FSDS, PREC, WIND, ZBOT };

// from snow_snicar - needed here for array sizing
static constexpr int sno_nbr_aer{8};      // number of aerosol species in snowpack
static constexpr int numrad_snw{5};            // number of spectral bands used in snow model [nbr]


// from landunit_varcon.F90
static constexpr int istsoil{1};
static constexpr int istcrop{2};
static constexpr int istice{3};
static constexpr int istice_mec{4};
static constexpr int istdlak{5};
static constexpr int istwet{6};
static constexpr int isturb_MIN{7};
static constexpr int isturb_tbd{7};
static constexpr int isturb_hd{8};
static constexpr int isturb_md{9};
static constexpr int isturb_MAX{9};
static constexpr int max_lunit{9};
static constexpr int landunit_name_length{40};
static constexpr int numurbl{isturb_MAX - isturb_MIN + 1};

// from SoilStateType
static const double smpmin{-1.e8}; // restriction for min of soil potential (mm)


// from pft_varcon.F90
static constexpr int noveg{0};
static constexpr int ndllf_evr_tmp_tree{1};
static constexpr int ndllf_evr_brl_tree{2};
static constexpr int ndllf_dcd_brl_tree{3};
static constexpr int nbrdlf_evr_trp_tree{4};
static constexpr int nbrdlf_evr_tmp_tree{5};
static constexpr int nbrdlf_dcd_trp_tree{6};
static constexpr int nbrdlf_dcd_tmp_tree{7};
static constexpr int nbrdlf_dcd_brl_tree{8};
static constexpr int nbrdlf_evr_shrub{9};
static constexpr int nbrdlf_dcd_tmp_shrub{10};
static constexpr int nbrdlf_dcd_brl_shrub{11};
static constexpr int nc3_arctic_grass{12};
static constexpr int nc3_nonarctic_grass{13};
static constexpr int nc4_grass{14};
static constexpr int nc3crop{15};
static constexpr int nc3irrig{16};
static constexpr int ncorn{17};
static constexpr int ncornirrig{18};
static constexpr int nscereal{19};
static constexpr int nscerealirrig{20};
static constexpr int nwcereal{21};
static constexpr int nwcerealirrig{22};
static constexpr int nsoybean{23};
static constexpr int nsoybeanirrig{24};

// from clm_varctl.F90
static constexpr int subgridflag{1};
static constexpr double co2_ppmv{355.0}; // atmospheric CO2 molar ratio (by volume) (umol/mol) - this is really low

// from clm_varpar.F90
static constexpr int nlevsno{5}; // this should be set by driver
static constexpr int nlevgrnd{15};
static constexpr int nlevurb{5};  // number of urban layers
static constexpr int numrad{2};   // number of solar radiation bands: vis, nir
static constexpr int nlevcan{1};  // number of leaf layers in canopy layer
static constexpr int nlevsoi{10}; // number of soil layers (hydrologically active)
static constexpr int nlevbed{15}; // number of layers to bedrock (hydrologically inactive below)
static constexpr int mxpft{25};  // maximum number of PFT's for any mode
static constexpr int numveg{17}; // number of veg types (without specific crop)
static constexpr bool use_crop{false};
static constexpr int numpft{(use_crop) ? mxpft : numveg}; // number of pfts - (use_crop ? 25 : 17)

// from column_varcon.F90
static constexpr int icol_roof{isturb_MIN * 10 + 1};
static constexpr int icol_sunwall{isturb_MIN * 10 + 2};
static constexpr int icol_shadewall{isturb_MIN * 10 + 3};
static constexpr int icol_road_imperv{isturb_MIN * 10 + 4};
static constexpr int icol_road_perv{isturb_MIN * 10 + 5};

// from clm_varcon.F90
static constexpr double snw_rds_min{54.526}; // minimum allowed snow effective radius (also "fresh snow" value) [microns]
static constexpr double snw_rds_max{1500.0}; // maximum allowed snow effective radius [microns]
static constexpr double zlnd{0.01};          // Roughness length for soil [m]
static constexpr double zsno{0.0024};        // Roughness length for snow [m]
static constexpr double spval{1.0e36};       // special value for real data
static constexpr int ispval{-9999};          // special value for int data
static constexpr double tfrz{273.15};    // freezing temperature [K]
static constexpr double denice{0.917e3};     // density of ice             ~ kg/m^3
static constexpr double denh2o{1.000e3};     // // density of fresh water     ~ kg/m^3
static constexpr double hvap{2.501e6};       // Latent heat of evap for water [J/kg]
static constexpr double hfus{3.337e5};       // Latent heat of fusion for ice [J/kg]
static constexpr double hsub{hvap + hfus};   // Latent heat of sublimation    [J/kg]
static constexpr double grav{9.80616};       // gravity constant [m/s2]
static constexpr double roverg{constants::RWV / grav * 1000.0}; // Rw/g constant = (8.3144/0.018)/(9.80616)*1000. mm/K
static constexpr double vkc{0.4};                        // von Karman constant [-]
static constexpr double cpair{1.00464e3};                // specific heat of dry air   ~ J/kg/K
static constexpr double csoilc{0.004};                   // Drag coefficient for soil under canopy [-]
static constexpr double alpha_aero{1.0};                 // constant for aerodynamic parameter weighting
static constexpr double tlsai_crit{2.0};    // critical value of elai+esai for which aerodynamic parameters are maximum
static constexpr double h2osno_max{1000.0}; // max allowed snow thickness (mm H2O)
static constexpr double bdsno{250.0};       // bulk density snow (kg/m**3)

static constexpr bool perchroot{false};     // calculate root function only in unfrozen soil, based on instantaneous temp
static constexpr bool perchroot_alt{false}; // calculate root function only in active layer - estimated from past 2 years

// should be read from file - hardwired for now
static constexpr double organic_max{130.0}; // Organic matter content where soil is assumed to act like peat for diffusion ~ kg/m3
static constexpr double secspday{86400.0}; // seconds per day

} // namespace ELM
