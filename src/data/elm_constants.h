#pragma once


namespace ELM {
  // enum for using forcing data utility
  enum class AtmForcType { TBOT, PBOT, QBOT, RH, FLDS, FSDS, PREC, WIND, ZBOT };
} // namespace ELM

// configuration options
namespace ELM::ELMconfig {
  inline constexpr int subgridflag() { return 1; }        
  inline constexpr int use_crop() { return 0; }      // perform crop calculations
  inline constexpr int perchroot() { return 0; }     // calculate root function only in unfrozen soil, based on instantaneous temp
  inline constexpr int perchroot_alt() { return 0; } // calculate root function only in active layer - estimated from past 2 years
} // namespace ELM::ELMconfig

// physical constants
namespace ELM::ELMconst {
  inline constexpr double TFRZ() { return 273.15; }                   // freezing temperature [K]
  inline constexpr double ELM_PI() { return 3.14159265358979323846; } // pi
  inline constexpr double BOLTZ() { return 1.38065e-23; }             // Boltzmann's constant ~ J/K/molecule
  inline constexpr double AVOGAD() { return 6.02214e26; }             // Avogadro's number ~ molecules/kmole
  inline constexpr double MWWV() { return 18.016; }                   // molecular weight water vapor
  inline constexpr double RGAS() { return AVOGAD() * BOLTZ(); }       // Universal gas constant ~ J/K/kmole
  inline constexpr double RWV() { return RGAS() / MWWV(); }           // Water vapor gas constant ~ J/K/kg
  inline constexpr double STEBOL() { return 5.67e-8; }                // Stefan-Boltzmann constant ~ W/m^2/K^4
  inline constexpr double MWDAIR() { return 28.966; }                 // molecular weight dry air ~ kg/kmole
  inline constexpr double RAIR() { return RGAS() / MWDAIR(); }        // Dry air gas constant     ~ J/K/kg
  inline constexpr double GRAV() { return 9.80616; }                  // gravity constant [m/s2]
  inline constexpr double ROVERG() { return RWV() / GRAV() * 1000.; } // Rw/g constant = (8.3144/0.018)/(9.80616)*1000. mm/K
  inline constexpr double O2_MOLAR_CONST() { return 0.209; }          // constant atmospheric O2 molar ratio (mol/mol)
  inline constexpr double CO2_PPMV() { return 355.0; }                // atmospheric CO2 molar ratio (by volume) (umol/mol) - designed for long spinup, really low - TODO consider changing to a modern value 
  inline constexpr double DENICE() { return 0.917e3; }                // density of ice             ~ kg/m^3
  inline constexpr double DENH2O() { return 1.000e3; }                // // density of fresh water     ~ kg/m^3
  inline constexpr double HVAP() { return 2.501e6; }                  // Latent heat of evap for water [J/kg]
  inline constexpr double HFUS() { return 3.337e5; }                  // Latent heat of fusion for ice [J/kg]
  inline constexpr double HSUB() { return HVAP() + HFUS(); }          // Latent heat of sublimation    [J/kg]
  inline constexpr double VKC() { return 0.4; }                       // von Karman constant [-]
  inline constexpr double CPAIR() { return 1.00464e3; }               // specific heat of dry air   ~ J/kg/K
  inline constexpr double CPICE() { return 2.11727e3; }               // specific heat of fresh ice ~ J/kg/K
  inline constexpr double CPWAT() { return 4.188e3; }                 // specific heat of fresh h2o ~ J/kg/K
  inline constexpr double CSOILC() { return 0.004; }                  // Drag coefficient for soil under canopy [-]
  inline constexpr double ZLND() { return 0.01; }                     // Roughness length for soil [m]
  inline constexpr double ZSNO() { return 0.0024; }                   // Roughness length for snow [m]
  inline constexpr double SNW_RDS_MIN() { return 54.526; }            // minimum allowed snow effective radius (also "fresh snow" value) [microns]
  inline constexpr double SNW_RDS_MAX() { return 1500.0; }            // maximum allowed snow effective radius [microns]
  inline constexpr double H2OSNO_MAX() { return 1000.0; }             // max allowed snow thickness (mm H2O)
  inline constexpr double BDSNO() { return 250.0; }                   // bulk density snow (kg/m**3)
  inline constexpr double SECSPDAY() { return 86400.0; }              // seconds per day
  inline constexpr double SPVAL() { return 1.0e36; }                  // special value for real data
  inline constexpr int ISPVAL() { return -9999; }                     // special value for int data
} // namespace ELM::ELMconst

//  plant functional trait keys
namespace ELM::PFT {
  inline constexpr int noveg() { return 0; }
  inline constexpr int ndllf_evr_tmp_tree() { return 1; }
  inline constexpr int ndllf_evr_brl_tree() { return 2; }
  inline constexpr int ndllf_dcd_brl_tree() { return 3; }
  inline constexpr int nbrdlf_evr_trp_tree() { return 4; }
  inline constexpr int nbrdlf_evr_tmp_tree() { return 5; }
  inline constexpr int nbrdlf_dcd_trp_tree() { return 6; }
  inline constexpr int nbrdlf_dcd_tmp_tree() { return 7; }
  inline constexpr int nbrdlf_dcd_brl_tree() { return 8; }
  inline constexpr int nbrdlf_evr_shrub() { return 9; }
  inline constexpr int nbrdlf_dcd_tmp_shrub() { return 10; }
  inline constexpr int nbrdlf_dcd_brl_shrub() { return 11; }
  inline constexpr int nc3_arctic_grass() { return 12; }
  inline constexpr int nc3_nonarctic_grass() { return 13; }
  inline constexpr int nc4_grass() { return 14; }
  inline constexpr int nc3crop() { return 15; }
  inline constexpr int nc3irrig() { return 16; }
  inline constexpr int ncorn() { return 17; }
  inline constexpr int ncornirrig() { return 18; }
  inline constexpr int nscereal() { return 19; }
  inline constexpr int nscerealirrig() { return 20; }
  inline constexpr int nwcereal() { return 21; }
  inline constexpr int nwcerealirrig() { return 22; }
  inline constexpr int nsoybean() { return 23; }
  inline constexpr int nsoybeanirrig() { return 24; }
} // namespace ELM::PFT

// dimensions of the problem
namespace ELM::ELMdims {
  inline constexpr int nlevsno() { return 5; }     // maximum number of snow layers
  inline constexpr int nlevgrnd() { return 15; }   // number of total subsurface layers
  inline constexpr int nlevurb() { return 5; }     // number of urban layers
  inline constexpr int numrad() { return 2; }      // number of solar radiation bands: vis, nir
  inline constexpr int nlevcan() { return 1; }     // number of leaf layers in canopy layer
  inline constexpr int nlevsoi() { return 10; }    // number of soil layers (hydrologically active)
  inline constexpr int nlevbed() { return 15; }    // number of layers to bedrock (hydrologically inactive below)
  inline constexpr int mxpft() { return 25; }      // maximum number of PFT's for any mode
  inline constexpr int numveg() { return 17; }     // number of veg types (without specific crop)
  inline constexpr int numpft() { return (ELMconfig::use_crop()) ? mxpft() : numveg(); } // number of pfts - (use_crop ? 25 : 17)
  inline constexpr int sno_nbr_aer() { return 8; } // number of aerosol species in snowpack
  inline constexpr int numrad_snw() { return 5; }  // number of spectral bands used in snow model [nbr]
  inline constexpr int nband() { return 5; }       // number of bands of the tridigonal matrix
} // namespace ELM::ELMdims


