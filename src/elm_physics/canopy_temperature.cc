// functions derived from CanopyTemperatureMod.F90

#include "elm_constants.h"
#include "landtype.h"
#include <algorithm>
#include <cmath>

namespace ELM {
namespace canopy_temperature {

void forcing_height(const LandType &Land, const bool &veg_active, const int &frac_veg_nosno,
                            const double &forc_hgt_u, const double &forc_hgt_t, const double &forc_hgt_q,
                            const double &z0m, const double &z0mg, const double &z_0_town, const double &z_d_town,
                            const double &forc_t, const double &displa, double &forc_hgt_u_patch,
                            double &forc_hgt_t_patch, double &forc_hgt_q_patch, double &thm) {
  // Make forcing height a pft-level quantity that is the atmospheric forcing
  // height plus each pft's z0m+displa
  if (veg_active) {
    if (Land.ltype == istsoil || Land.ltype == istcrop) {
      if (frac_veg_nosno == 0) {
        forc_hgt_u_patch = forc_hgt_u + z0mg + displa;
        forc_hgt_t_patch = forc_hgt_t + z0mg + displa;
        forc_hgt_q_patch = forc_hgt_q + z0mg + displa;
      } else {
        forc_hgt_u_patch = forc_hgt_u + z0m + displa;
        forc_hgt_t_patch = forc_hgt_t + z0m + displa;
        forc_hgt_q_patch = forc_hgt_q + z0m + displa;
      }
    } else if (Land.ltype == istwet || Land.ltype == istice || Land.ltype == istice_mec) {
      forc_hgt_u_patch = forc_hgt_u + z0mg;
      forc_hgt_t_patch = forc_hgt_t + z0mg;
      forc_hgt_q_patch = forc_hgt_q + z0mg;
    } else if (Land.ltype == istdlak) {
      forc_hgt_u_patch = forc_hgt_u;
      forc_hgt_t_patch = forc_hgt_t;
      forc_hgt_q_patch = forc_hgt_q;
    } else if (Land.urbpoi) {
      forc_hgt_u_patch = forc_hgt_u + z_0_town + z_d_town;
      forc_hgt_t_patch = forc_hgt_t + z_0_town + z_d_town;
      forc_hgt_q_patch = forc_hgt_q + z_0_town + z_d_town;
    }
  }

  thm = forc_t + 0.0098 * forc_hgt_t_patch;
} // forcing_height

void init_energy_fluxes(const LandType &Land, double &eflx_sh_tot, double &eflx_sh_tot_u, double &eflx_sh_tot_r,
                            double &eflx_lh_tot, double &eflx_lh_tot_u, double &eflx_lh_tot_r, double &eflx_sh_veg,
                            double &qflx_evap_tot, double &qflx_evap_veg, double &qflx_tran_veg) {
  // Initial set (needed for history tape fields)
  eflx_sh_tot = 0.0;
  if (Land.urbpoi) {
    eflx_sh_tot_u = 0.0;
  } else if (Land.ltype == istsoil || Land.ltype == istcrop) {
    eflx_sh_tot_r = 0.0;
  }
  eflx_lh_tot = 0.0;
  if (Land.urbpoi) {
    eflx_lh_tot_u = 0.0;
  } else if (Land.ltype == istsoil || Land.ltype == istcrop) {
    eflx_lh_tot_r = 0.0;
  }
  eflx_sh_veg = 0.0;
  qflx_evap_tot = 0.0;
  qflx_evap_veg = 0.0;
  qflx_tran_veg = 0.0;
} // init_energy_fluxes

} // namespace canopy_temperature
} // namespace ELM
