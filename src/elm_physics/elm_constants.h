#pragma once

namespace ELM::constants {
    
inline constexpr double TFRZ{273.15}; // freezing temperature [K]
inline constexpr double ELM_PI{3.14159265358979323846}; // pi
inline constexpr double BOLTZ{1.38065e-23};         // Boltzmann's constant ~ J/K/molecule
inline constexpr double AVOGAD{6.02214e26};         // Avogadro's number ~ molecules/kmole
inline constexpr double MWWV{18.016};               // molecular weight water vapor
inline constexpr double RGAS{AVOGAD * BOLTZ};       // Universal gas constant ~ J/K/kmole
inline constexpr double RWV{RGAS / MWWV};           // Water vapor gas constant ~ J/K/kg
inline constexpr double STEBOL{5.67e-8};     // Stefan-Boltzmann constant ~ W/m^2/K^4
inline constexpr double MWDAIR{28.966};      // molecular weight dry air ~ kg/kmole
inline constexpr double RAIR{RGAS / MWDAIR}; // Dry air gas constant     ~ J/K/kg
inline constexpr double O2_MOLAR_CONST{0.209}; // constant atmospheric O2 molar ratio (mol/mol)
inline constexpr double CO2_PPMV{355.0}; // atmospheric CO2 molar ratio (by volume) (umol/mol) - this is really low

} // namespace ELM::ELMconstants

namespace ELM {

enum class AtmForcType { TBOT, PBOT, QBOT, RH, FLDS, FSDS, PREC, WIND, ZBOT };

// from landunit_varcon.F90
inline constexpr int istsoil{1};
inline constexpr int istcrop{2};
inline constexpr int istice{3};
inline constexpr int istice_mec{4};
inline constexpr int istdlak{5};
inline constexpr int istwet{6};
inline constexpr int isturb_MIN{7};
inline constexpr int isturb_tbd{7};
inline constexpr int isturb_hd{8};
inline constexpr int isturb_md{9};
inline constexpr int isturb_MAX{9};
inline constexpr int max_lunit{9};
inline constexpr int landunit_name_length{40};
inline constexpr int numurbl{isturb_MAX - isturb_MIN + 1};

// from SoilStateType
inline const double smpmin{-1.e8}; // restriction for min of soil potential (mm)

// from CanopyHydrology() - move to CanopyHydrology header when ready
inline const double accum_factor{0.1}; // set shape factor for accumulation of snow

// from pft_varcon.F90
inline constexpr int noveg{0};
inline constexpr int ndllf_evr_tmp_tree{1};
inline constexpr int ndllf_evr_brl_tree{2};
inline constexpr int ndllf_dcd_brl_tree{3};
inline constexpr int nbrdlf_evr_trp_tree{4};
inline constexpr int nbrdlf_evr_tmp_tree{5};
inline constexpr int nbrdlf_dcd_trp_tree{6};
inline constexpr int nbrdlf_dcd_tmp_tree{7};
inline constexpr int nbrdlf_dcd_brl_tree{8};
inline constexpr int nbrdlf_evr_shrub{9};
inline constexpr int nbrdlf_dcd_tmp_shrub{10};
inline constexpr int nbrdlf_dcd_brl_shrub{11};
inline constexpr int nc3_arctic_grass{12};
inline constexpr int nc3_nonarctic_grass{13};
inline constexpr int nc4_grass{14};
inline constexpr int nc3crop{15};
inline constexpr int nc3irrig{16};
inline constexpr int ncorn{17};
inline constexpr int ncornirrig{18};
inline constexpr int nscereal{19};
inline constexpr int nscerealirrig{20};
inline constexpr int nwcereal{21};
inline constexpr int nwcerealirrig{22};
inline constexpr int nsoybean{23};
inline constexpr int nsoybeanirrig{24};

// from clm_varctl.F90
inline constexpr int subgridflag{1};
inline constexpr double co2_ppmv{355.0}; // atmospheric CO2 molar ratio (by volume) (umol/mol) - this is really low

// from clm_varpar.F90
inline constexpr int nlevsno{5}; // this should be set by driver
inline constexpr int nlevgrnd{15};
inline constexpr int nlevurb{5};  // number of urban layers
inline constexpr int numrad{2};   // number of solar radiation bands: vis, nir
inline constexpr int nlevcan{1};  // number of leaf layers in canopy layer
inline constexpr int nlevsoi{10}; // number of soil layers (hydrologically active)
inline constexpr int nlevbed{15}; // number of layers to bedrock (hydrologically inactive below)
inline constexpr int mxpft{25};  // maximum number of PFT's for any mode
inline constexpr int numveg{17}; // number of veg types (without specific crop)
inline constexpr bool use_crop{false};
inline constexpr int numpft{(use_crop) ? mxpft : numveg}; // number of pfts - (use_crop ? 25 : 17)

// from column_varcon.F90
inline constexpr int icol_roof{isturb_MIN * 10 + 1};
inline constexpr int icol_sunwall{isturb_MIN * 10 + 2};
inline constexpr int icol_shadewall{isturb_MIN * 10 + 3};
inline constexpr int icol_road_imperv{isturb_MIN * 10 + 4};
inline constexpr int icol_road_perv{isturb_MIN * 10 + 5};

// from clm_varcon.F90
inline constexpr double snw_rds_min{54.526}; // minimum allowed snow effective radius (also "fresh snow" value) [microns]
inline constexpr double zlnd{0.01};          // Roughness length for soil [m]
inline constexpr double zsno{0.0024};        // Roughness length for snow [m]
inline constexpr double spval{1.0e36};       // special value for real data
inline constexpr int ispval{-9999};          // special value for int data
inline constexpr double tfrz{273.15};    // freezing temperature [K]
inline constexpr double denice{0.917e3};     // density of ice             ~ kg/m^3
inline constexpr double denh2o{1.000e3};     // // density of fresh water     ~ kg/m^3
inline constexpr double hvap{2.501e6};       // Latent heat of evap for water [J/kg]
inline constexpr double hfus{3.337e5};       // Latent heat of fusion for ice [J/kg]
inline constexpr double hsub{hvap + hfus};   // Latent heat of sublimation    [J/kg]
inline constexpr double grav{9.80616};       // gravity constant [m/s2]
inline constexpr double roverg{constants::RWV / grav * 1000.0}; // Rw/g constant = (8.3144/0.018)/(9.80616)*1000. mm/K
inline constexpr double vkc{0.4};                        // von Karman constant [-]
inline constexpr double cpair{1.00464e3};                // specific heat of dry air   ~ J/kg/K
inline constexpr double csoilc{0.004};                   // Drag coefficient for soil under canopy [-]
inline constexpr double alpha_aero{1.0};                 // constant for aerodynamic parameter weighting
inline constexpr double tlsai_crit{2.0};    // critical value of elai+esai for which aerodynamic parameters are maximum
inline constexpr double h2osno_max{1000.0}; // max allowed snow thickness (mm H2O)
inline constexpr double bdsno{250.0};       // bulk density snow (kg/m**3)

inline constexpr bool perchroot{false};     // calculate root function only in unfrozen soil, based on instantaneous temp
inline constexpr bool perchroot_alt{false}; // calculate root function only in active layer - estimated from past 2 years

// should be read from file - hardwired for now
inline constexpr double organic_max{130.0}; // Organic matter content where soil is assumed to act like peat for diffusion ~ kg/m3
inline constexpr double secspday{86400.0}; // seconds per day

} // namespace ELM
