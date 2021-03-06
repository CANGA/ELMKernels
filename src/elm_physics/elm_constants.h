#pragma once


namespace ELM {
  // enum for using forcing data utility
  enum class AtmForcType { TBOT, PBOT, QBOT, RH, FLDS, FSDS, PREC, WIND, ZBOT };
} // namespace ELM

// configuration options
namespace ELM::ELMconfig {
  static constexpr int subgridflag{1};        
  static constexpr bool use_crop{false};      // perform crop calculations
  static constexpr bool perchroot{false};     // calculate root function only in unfrozen soil, based on instantaneous temp
  static constexpr bool perchroot_alt{false}; // calculate root function only in active layer - estimated from past 2 years
} // namespace ELM::ELMconfig

// physical constants
namespace ELM::ELMconst {
  static constexpr double TFRZ{273.15};                   // freezing temperature [K]
  static constexpr double ELM_PI{3.14159265358979323846}; // pi
  static constexpr double BOLTZ{1.38065e-23};             // Boltzmann's constant ~ J/K/molecule
  static constexpr double AVOGAD{6.02214e26};             // Avogadro's number ~ molecules/kmole
  static constexpr double MWWV{18.016};                   // molecular weight water vapor
  static constexpr double RGAS{AVOGAD * BOLTZ};           // Universal gas constant ~ J/K/kmole
  static constexpr double RWV{RGAS / MWWV};               // Water vapor gas constant ~ J/K/kg
  static constexpr double STEBOL{5.67e-8};                // Stefan-Boltzmann constant ~ W/m^2/K^4
  static constexpr double MWDAIR{28.966};                 // molecular weight dry air ~ kg/kmole
  static constexpr double RAIR{RGAS / MWDAIR};            // Dry air gas constant     ~ J/K/kg
  static constexpr double GRAV{9.80616};                  // gravity constant [m/s2]
  static constexpr double ROVERG{RWV / GRAV * 1000.0};    // Rw/g constant = (8.3144/0.018)/(9.80616)*1000. mm/K
  static constexpr double O2_MOLAR_CONST{0.209};          // constant atmospheric O2 molar ratio (mol/mol)
  static constexpr double CO2_PPMV{355.0};                // atmospheric CO2 molar ratio (by volume) (umol/mol) - this is really low
  static constexpr double DENICE{0.917e3};                // density of ice             ~ kg/m^3
  static constexpr double DENH2O{1.000e3};                // // density of fresh water     ~ kg/m^3
  static constexpr double HVAP{2.501e6};                  // Latent heat of evap for water [J/kg]
  static constexpr double HFUS{3.337e5};                  // Latent heat of fusion for ice [J/kg]
  static constexpr double HSUB{HVAP + HFUS};              // Latent heat of sublimation    [J/kg]
  static constexpr double VKC{0.4};                       // von Karman constant [-]
  static constexpr double CPAIR{1.00464e3};               // specific heat of dry air   ~ J/kg/K
  static constexpr double CPICE{2.11727e3};               // specific heat of fresh ice ~ J/kg/K
  static constexpr double CPWAT{4.188e3};                 // specific heat of fresh h2o ~ J/kg/K
  static constexpr double CSOILC{0.004};                  // Drag coefficient for soil under canopy [-]
  static constexpr double ZLND{0.01};                     // Roughness length for soil [m]
  static constexpr double ZSNO{0.0024};                   // Roughness length for snow [m]
  static constexpr double SNW_RDS_MIN{54.526};            // minimum allowed snow effective radius (also "fresh snow" value) [microns]
  static constexpr double SNW_RDS_MAX{1500.0};            // maximum allowed snow effective radius [microns]
  static constexpr double H2OSNO_MAX{1000.0};             // max allowed snow thickness (mm H2O)
  static constexpr double BDSNO{250.0};                   // bulk density snow (kg/m**3)
  static constexpr double SECSPDAY{86400.0};              // seconds per day
  static constexpr double SPVAL{1.0e36};                  // special value for real data
  static constexpr int ISPVAL{-9999};                     // special value for int data
} // namespace ELM::ELMconst

// pft type data
namespace ELM::PFT {
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
} // namespace ELM::PFT

// dimensions of the problem
namespace ELM::ELMdims {
  static constexpr int nlevsno{5};     // maximum number of snow layers
  static constexpr int nlevgrnd{15};   // number of total subsurface layers
  static constexpr int nlevurb{5};     // number of urban layers
  static constexpr int numrad{2};      // number of solar radiation bands: vis, nir
  static constexpr int nlevcan{1};     // number of leaf layers in canopy layer
  static constexpr int nlevsoi{10};    // number of soil layers (hydrologically active)
  static constexpr int nlevbed{15};    // number of layers to bedrock (hydrologically inactive below)
  static constexpr int mxpft{25};      // maximum number of PFT's for any mode
  static constexpr int numveg{17};     // number of veg types (without specific crop)
  static constexpr int numpft{(ELMconfig::use_crop) ? mxpft : numveg}; // number of pfts - (use_crop ? 25 : 17)
  static constexpr int sno_nbr_aer{8}; // number of aerosol species in snowpack
  static constexpr int numrad_snw{5};  // number of spectral bands used in snow model [nbr]
  static constexpr int nband{5};       // number of bands of the tridigonal matrix
} // namespace ELM::ELMdims


