/* functions derived from CanopyFluxesMod.F90
DESCRIPTION:
Calculates leaf temperature leaf fluxes, transpiration, photosynthesis and updates the dew accumulation due to
evaporation. Also calculates the irrigation rate.

The member functions should be called in the following order:
                     InitializeFlux()
       StabilityIteration()   Irrigation()
ComputeFlux()
*/
#include "Photosynthesis.h"
#include "QSat.h"
#include "SoilMoistStress.hh"
#include "clm_constants.h"
#include "frictionvelocity.h"
#include "landtype.h"
#include "vegproperties.h"
#include <algorithm>
#include <assert.h>
#include <cmath>

namespace ELM {

class CanopyFluxes {

private:
  double btran0 = 0.0; // hardwired??
  double rssun;        // leaf sunlit stomatal resistance (s/m) (output from Photosynthesis)
  double rssha;        // leaf shaded stomatal resistance (s/m) (output from Photosynthesis)
  double lbl_rsc_h2o;  // laminar boundary layer resistance for h2o
  // double rresis[nlevgrnd]; // root soil water stress (resistance) by layer (0-1)
  double del;         // absolute change in leaf temp in current iteration [K]
  double efeb;        // latent heat flux from leaf (previous iter) [mm/s]
  double efe;         // water flux from leaf [mm/s]
  double wtg;         // heat conductance for ground [m/s]
  double wtl0;        // normalized heat conductance for leaf [-]
  double wta0;        // normalized heat conductance for air [-]
  double wtlq0;       // normalized latent heat conductance for leaf [-]
  double wtgq;        // latent heat conductance for ground [m/s]
  double wtal;        // normalized heat conductance for air and leaf [-]
  double wtalq;       // normalized latent heat cond. for air and leaf [-]
  double wtaq0;       // normalized latent heat conductance for air [-]
  double wtgaq;       // normalized latent heat cond. for air and ground [-]
  double obuold;      // monin-obukhov length from previous iteration
  double dayl_factor; // scalar (0-1) for daylength effect on Vcmax
  double air;         // atmos. radiation temporay set
  double bir;         // atmos. radiation temporay set
  double cir;         // atmos. radiation temporay set
  double el;          // vapor pressure on leaf surface [pa]
  double deldT;       // derivative of "el" on "t_veg" [pa/K]
  double qsatl;       // leaf specific humidity [kg/kg]
  double qsatldT;     // derivative of "qsatl" on "t_veg"
  double co2;         // atmospheric co2 partial pressure (pa)
  double o2;          // atmospheric o2 partial pressure (pa)
  int nmozsgn;        // number of times stability changes sign
  double taf;         // air temperature within canopy space [K]
  double qaf;         // humidity of canopy air [kg/kg]
  double um;          // wind speed including the stablity effect [m/s]
  double ur;          // wind speed at reference height [m/s]
  double dth;         // diff of virtual temp. between ref. height and surface
  double dthv;        // diff of vir. poten. temp. between ref. height and surface
  double dqh;         // diff of humidity between ref. height and surface
  double obu;         // Monin-Obukhov length (m)
  double zldis;       // reference height "minus" zero displacement height [m]
  double temp1;       // relation for potential temperature profile
  double temp2;       // relation for specific humidity profile
  double temp12m;     // relation for potential temperature profile applied at 2-m
  double temp22m;     // relation for specific humidity profile applied at 2-m
  double tlbef;       // leaf temperature from previous iteration [K]
  double delq;        // temporary
  double dt_veg;      // change in t_veg, last iteration (Kelvin)

  int frac_veg_nosno;   // local copy of fraction of vegetation not covered by snow (0 OR 1) [-]
  double frac_sno;      // local copy of fraction of ground covered by snow (0 to 1)
  double displa;        // local copy of displa
  double t_veg;         // local copy of vegetation temperature (Kelvin)
  double elai;          // local copy of one-sided leaf area index with burying by snow
  double esai;          // local copy of one-sided stem area index with burying by snow
  double t_grnd;        // local copy of ground temperature (Kelvin)
  double forc_pbot;     // local copy of atmospheric pressure (Pa)
  double qflx_tran_veg; // local copy of vegetation transpiration (mm H2O/s) (+ = to atm)
  double qflx_evap_veg; // local copy of vegetation evaporation (mm H2O/s) (+ = to atm)
  double forc_u;        // local copy of atmospheric wind speed in east direction (m/s)
  double forc_v;        // local copy of atmospheric wind speed in north direction (m/s)
  double forc_q;        // local copy of atmospheric specific humidity (kg/kg)
  double forc_th;       // local copy of atmospheric potential temperature (K)
  double forc_rho;      // local copy of air density (kg/m**3)
  double z0mg;          // local copy of roughness length over ground, momentum [m]
  double z0mv;          // local copy of roughness length over vegetation, momentum [m]
  double z0hv;          // local copy of roughness length over vegetation, sensible heat [m]
  double z0qv;          // local copy of roughness length over vegetation, latent heat [m]
  int snl;              // local copy of number of active snow layers
  double dtime;         // local copy of timestep size (sec)
  double thm;           // local copy of intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
  double qg;            // local copy of ground specific humidity [kg/kg]
  double eflx_sh_veg;   // local copy of sensible heat flux from leaves (W/m**2) [+ to atm]
  double thv;           // local copy of virtual potential temperature (kelvin)
  double emv;           // local copy of vegetation emissivity
  double emg;           // local copy of ground emissivity
  double forc_lwrad;    // downward infrared (longwave) radiation (W/m**2)

public:
  /*
  We start applying the irrigation in the time step FOLLOWING this time,
  since we won't begin irrigating until the next call to CanopyHydrology



  */
  // void Irrigation(
  //  const LandType& Land,
  //  const double elai,
  //  const double *irrigated, // doesn't need to be in array form if 1 pft/cell -- fix
  //
  //  int& n_irrig_steps_left,
  //  double& irrig_rate)
  //{
  //  if (!Land.lakpoi && !Land.urbpoi && frac_veg_nosno != 0) {
  //    // Determines target soil moisture level for irrigation. If h2osoi_liq_so is the soil moisture level at
  //    // which stomata are fully open and h2osoi_liq_sat is the soil moisture level at saturation (eff_porosity),
  //    // then the target soil moisture level is (h2osoi_liq_so + irrig_factor*(h2osoi_liq_sat - h2osoi_liq_so)).
  //    // A value of 0 means that the target soil moisture level is h2osoi_liq_so.
  //    // A value of 1 means that the target soil moisture level is h2osoi_liq_sat
  //    double irrig_factor = 0.7;
  //
  //    double irrig_min_lai = 0.0;  // Minimum LAI for irrigation
  //    double irrig_btran_thresh = 0.999999; // Irrigate when btran falls below 0.999999 rather than 1 to allow for
  //    round-off error int irrig_start_time = isecspday/4   // Time of day to check whether we need irrigation, seconds
  //    (0 = midnight). int irrig_length = isecspday/6; // Desired amount of time to irrigate per day (sec). Actual time
  //    may differ if this is not a multiple of dtime. Irrigation won't work properly if dtime > secsperday int
  //    irrig_nsteps_per_day = ((irrig_length + (dtime - 1))/dtime);  // number of time steps per day in which we
  //    irrigate
  //
  //    // Determine if irrigation is needed (over irrigated soil columns)
  //    // First, determine in what grid cells we need to bother 'measuring' soil water, to see if we need irrigation
  //    // Also set n_irrig_steps_left for these grid cells
  //    // n_irrig_steps_left(p) > 0 is ok even if irrig_rate(p) ends up = 0
  //    // in this case, we'll irrigate by 0 for the given number of time steps
  //
  //    // get_prev_date(yr, mon, day, time)  ! get time as of beginning of time step --- figure out!! -- need variable
  //    'time' if (irrigated[Land.vtype] == 1.0 && elai > irrig_min_lai && btran < irrig_btran_thresh) {
  //      //see if it's the right time of day to start irrigating:
  //      int local_time = (((time + round(londeg/degpsec)) % isecspday) + isecspday) % isecspday;
  //      int seconds_since_irrig_start_time = (((local_time - irrig_start_time) % isecspday) + isecspday) % isecspday;
  //      if (seconds_since_irrig_start_time < dtime) { // it's time to start irrigating
  //        bool check_for_irrig = true;
  //        n_irrig_steps_left = irrig_nsteps_per_day;
  //        irrig_rate = 0.0  // reset; we'll add to this later
  //      } else {
  //        bool check_for_irrig = false;
  //      }
  //    } else { // non-irrig pft or elai<=irrig_min_lai or btran>irrig_btran_thresh
  //      bool check_for_irrig = false;
  //    }
  //
  //    // Now 'measure' soil water for the grid cells identified above and see if the
  //    // soil is dry enough to warrant irrigation
  //    // (Note: frozen_soil could probably be a column-level variable, but that would be
  //    // slightly less robust to potential future modifications)
  //    // This should not be operating on FATES patches (see is_fates filter above, pushes
  //    // check_for_irrig = false
  //    bool frozen_soil = false;
  //    if (check_for_irrig && !frozen_soil) {
  //      for (int i = 0; i < nlevgrnd; i++) {
  //        // if level i was frozen, then we don't look at any levels below L
  //        if (t_soisno[nlevsno+i] <= tfrz) {
  //          frozen_soil = true;
  //        } else if (rootfr > 0.0) {
  //          // determine soil water deficit in this layer:
  //          // Calculate vol_liq_so - i.e., vol_liq at which smp_node = smpso - by inverting the above equations
  //          // for the root resistance factors
  //          vol_liq_so = eff_porosity[i] * pow((-smpso[Land.vtype]/sucsat[i]), (-1.0/bsw[i]));
  //          // Translate vol_liq_so and eff_porosity into h2osoi_liq_so and h2osoi_liq_sat and calculate deficit
  //          h2osoi_liq_so = vol_liq_so * denh2o * dz[i];
  //          h2osoi_liq_sat = eff_porosity[i] * denh2o * dz[i];
  //          deficit = std::max((h2osoi_liq_so + irrig_factor * (h2osoi_liq_sat - h2osoi_liq_so)) -
  //          h2osoi_liq[nlevsno+i], 0.0);
  //          // Add deficit to irrig_rate, converting units from mm to mm/sec
  //          irrig_rate += deficit / (dtime*irrig_nsteps_per_day);
  //        }
  //      }
  //    }
  //  } // if (!Land.lakpoi && !Land.urbpoi && frac_veg_nosno != 0)
  //} // Irrigation

  /*
  INPUTS:
  Land             [LandType] struct containing information about landtype
  dtime_in                   [double] timestep size (sec)
  snl_in                        [int] number of snow layers
  frac_veg_nosno_in [int]  fraction of vegetation not covered by snow (0 OR 1) [-]
  frac_sno_in       [double] fraction of ground covered by snow (0 to 1)
  forc_hgt_u_patch [double] observational height of wind at pft level [m]
  thm_in            [double]  intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
  thv_in                          [double] virtual potential temperature (kelvin)
  max_dayl  [double] maximum daylength for this grid cell (s)
  dayl      [double] daylength (seconds)
  altmax_indx                     [int] index corresponding to maximum active layer depth from current year
  altmax_lastyear_indx            [int] index corresponding to maximum active layer depth from prior year
  t_soisno[nlevgrnd+nlevsno] [double] col soil temperature (Kelvin)
  h2osoi_ice[nlevgrnd+nlevsno] [double] ice lens (kg/m2)
  h2osoi_liq[nlevgrnd+nlevsno] [double] liquid water (kg/m2)
  dz[nlevgrnd+nlevsno]         [double] layer thickness (m)
  rootfr[nlevgrnd]             [double] fraction of roots in each soil layer
  tc_stress                       [double] critical soil temperature for soil water stress (C)
  sucsat[nlevgrnd]                [double] minimum soil suction (mm)
  watsat[nlevgrnd]                [double] volumetric soil water at saturation (porosity)
  bsw[nlevgrnd]                   [double] Clapp and Hornberger "b
  smpso[numpft]                   [double] soil water potential at full stomatal opening (mm)
  smpsc[numpft]                   [double] soil water potential at full stomatal closure (mm)
  elai_in                         [double] one-sided leaf area index with burying by snow
  esai_in                         [double] one-sided stem area index with burying by snow
  emv_in                          [double] vegetation emissivity
  emg_in                          [double] ground emissivity
  qg_in                         [double] ground specific humidity [kg/kg]
  t_grnd_in                     [double] ground temperature (Kelvin
  forc_t         [double]  atmospheric temperature (Kelvin)
  forc_pbot_in      [double]  atmospheric pressure (Pa)
  forc_lwrad_in           [double]  downward infrared (longwave) radiation (W/m**2)
  forc_pco2            [double]  partial pressure co2 (Pa)
  forc_po2             [double]  partial pressure o2 (Pa))
  forc_u_in                     [double] atmospheric wind speed in east direction (m/s)
  forc_v_in                     [double] atmospheric wind speed in north direction (m/s)
  forc_q_in                     [double] atmospheric specific humidity (kg/kg)
  forc_th_in                    [double] atmospheric potential temperature (Kelvin)
  z0mg_in       [double] roughness length over ground, momentum [m]
  t_veg_in          [double]  vegetation temperature (Kelvin)


  OUTPUTS:
  btran          [double]  transpiration wetness factor (0 to 1)
  z0mv_out       [double] roughness length over vegetation, momentum [m]
  z0hv_out       [double] roughness length over vegetation, sensible heat [m]
  z0qv_out       [double] roughness length over vegetation, latent heat [m]
  displa_out              [double]  displacement height (m)
  rootr[nlevgrnd]  [double] effective fraction of roots in each soil layer
  eff_porosity[nlevgrnd]       [double] effective soil porosity
  h2osoi_liqvol[nlevgrnd+nlevsno] [double] liquid volumetric moisture



  */
  void InitializeFlux(const LandType &Land, const double &dtime_in, const int &snl_in, const int &frac_veg_nosno_in,
                      const double &frac_sno_in, const double &forc_hgt_u_patch, const double &thm_in,
                      const double &thv_in, const double &max_dayl, const double &dayl, const int &altmax_indx,
                      const int &altmax_lastyear_indx, const double t_soisno[nlevgrnd + nlevsno],
                      const double h2osoi_ice[nlevgrnd + nlevsno], const double h2osoi_liq[nlevgrnd + nlevsno],
                      const double dz[nlevgrnd + nlevsno], const double rootfr[nlevgrnd], const double &tc_stress,
                      const double *sucsat, const double *watsat, const double *bsw, const double *smpso,
                      const double *smpsc, const double &elai_in, const double &esai_in, const double &emv_in,
                      const double &emg_in, const double &qg_in, const double &t_grnd_in, const double &forc_t,
                      const double &forc_pbot_in, const double &forc_lwrad_in, const double &forc_pco2,
                      const double &forc_po2, const double &forc_u_in, const double &forc_v_in, const double &forc_q_in,
                      const double &forc_th_in, const double &z0mg_in, const double &t_veg_in,

                      double &btran, double &displa_out, double &z0mv_out, double &z0hv_out, double &z0qv_out,

                      double *rootr, double eff_porosity[nlevgrnd], double h2osoi_liqvol[nlevgrnd + nlevsno])
  // need to decide where to place photosyns_vars%TimeStepInit() calc_effective_soilporosity() calc_volumetric_h2oliq()
  // calc_root_moist_stress()
  {

    // -----------------------------------------------------------------
    // Time step initialization of photosynthesis variables
    // -----------------------------------------------------------------
    // call photosyns_vars%TimeStepInit(bounds)

    forc_u = forc_u_in;
    forc_v = forc_v_in;
    forc_q = forc_q_in;
    t_veg = t_veg_in;
    elai = elai_in;
    esai = esai_in;
    t_grnd = t_grnd_in;
    forc_pbot = forc_pbot_in;
    frac_veg_nosno = frac_veg_nosno_in;
    frac_sno = frac_sno_in;
    snl = snl_in;
    forc_th = forc_th_in;
    dtime = dtime_in;
    thm = thm_in;
    qg = qg_in;
    thv = thv_in;
    emv = emv_in;
    emg = emg_in;
    forc_lwrad = forc_lwrad_in;

    if (!Land.lakpoi && !Land.urbpoi) {
      double lt; // elai+esai
      if (frac_veg_nosno == 0) {
        btran = 0.0;
        t_veg = forc_t;
        double cf_bare = forc_pbot / (ELM_RGAS * 0.001 * thm) * 1.e06; // heat transfer coefficient from bare ground [-]
        rssun = 1.0 / 1.e15 * cf_bare;
        rssha = 1.0 / 1.e15 * cf_bare;
        lbl_rsc_h2o = 0.0;
        for (int i = 0; i < nlevgrnd; i++) {
          rootr[i] = 0.0;
          // rresis[i] = 0.0;
        }
      } else {
        del = 0.0;  // change in leaf temperature from previous iteration
        efeb = 0.0; // latent head flux from leaf for previous iteration
        wtlq0 = 0.0;
        wtalq = 0.0;
        wtgq = 0.0;
        wtaq0 = 0.0;
        obuold = 0.0;
        btran = btran0;

        // calculate dayl_factor as the ratio of (current:max dayl)^2
        // set a minimum of 0.01 (1%) for the dayl_factor
        dayl_factor = std::min(1.0, std::max(0.01, (dayl * dayl) / (max_dayl * max_dayl)));

        // compute effective soil porosity
        calc_effective_soilporosity(watsat, h2osoi_ice, dz, eff_porosity);
        // compute volumetric liquid water content
        calc_volumetric_h2oliq(eff_porosity, h2osoi_liq, dz, h2osoi_liqvol);
        // calculate root moisture stress
        calc_root_moist_stress(Land.vtype, h2osoi_liqvol, rootfr, t_soisno, tc_stress, sucsat, watsat, bsw, smpso,
                               smpsc, eff_porosity, altmax_indx, altmax_lastyear_indx, rootr, btran);

        // Modify aerodynamic parameters for sparse/dense canopy (X. Zeng)
        lt = std::min(elai + esai, tlsai_crit);
        double egvf = (1.0 - exp(-lt)) / (1.0 - exp(-tlsai_crit));
        displa_out *= egvf;
        displa = displa_out;
        z0mv_out = exp(egvf * std::log(z0mv_out) + (1.0 - egvf) * std::log(z0mg_in));
        z0hv_out = z0mv_out;
        z0qv_out = z0mv_out;

        z0mg = z0mg_in;
        z0mv = z0mv_out;
        z0hv = z0hv_out;
        z0qv = z0qv_out;

        // Net absorbed longwave radiation by canopy and ground
        air = emv * (1.0 + (1.0 - emv) * (1.0 - emg)) * forc_lwrad;
        bir = -(2.0 - emv * (1.0 - emg)) * emv * sb;
        cir = emv * emg * sb;

        // Saturated vapor pressure, specific humidity, and their derivatives at the leaf surface
        QSat(t_veg, forc_pbot, el, deldT, qsatl, qsatldT);
        // Determine atmospheric co2 and o2
        co2 = forc_pco2;
        o2 = forc_po2;

        // Initialize flux profile
        nmozsgn = 0;
        taf = (t_grnd + thm) / 2.0;
        qaf = (forc_q + qg) / 2.0;
        ur = std::max(1.0, std::sqrt(forc_u * forc_u + forc_v * forc_v));
        dth = thm - taf;
        dqh = forc_q - qaf;
        delq = qg - qaf;
        dthv = dth * (1.0 + 0.61 * forc_q) + 0.61 * forc_th * dqh;
        zldis = forc_hgt_u_patch - displa;

        // Check to see if the forcing height is below the canopy height
        assert(zldis >= 0.0);
        // Initialize Monin-Obukhov length and wind speed
        MoninObukIni(ur, thv, dthv, zldis, z0mv, um, obu);
      }
    }
  }

  /* CanopyFluxes::StabilityIteration()
  DESCRIPTION:
  calculate Monin-Obukhov length and wind speed, call photosynthesis, calculate ET & SH flux
  Iterates until convergence, up to 40 iterations, calling friction velocity functions, then
  photosynthesis for sun & shade.

  INPUTS:
  Land             [LandType] struct containing information about landtype
  forc_hgt_u_patch [double] observational height of wind at pft level [m]
  forc_hgt_t_patch [double] observational height of temperature at pft level [m]
  forc_hgt_q_patch [double] observational height of specific humidity at pft level [m]
  dleaf[numpft]            [double] characteristic leaf dimension (m)
  fwet             [double] fraction of canopy that is wet (0 to 1)
  fdry             [double] fraction of foliage that is green and dry [-]
  laisun           [double] sunlit leaf area
  laisha           [double] shaded leaf area
  forc_rho_in         [double] air density (kg/m**3)
  snow_depth       [double] snow height (m)
  soilbeta                   [double] soil wetness relative to field capacity
  frac_h2osfc                [double] fraction of ground covered by surface water (0 to 1)
  t_h2osfc                   [double] surface water temperature
  sabv               [double] solar radiation absorbed by vegetation (W/m**2)
  h2ocan                     [double] canopy water (mm H2O)
  htop               [double] canopy top(m)
  t_soisno[nlevgrnd+nlevsno] [double] col soil temperature (Kelvin)

  added for photosynthesis:
  veg                [VegProperties] struct containing vegetation constant parameters
  nrad                [int]  number of canopy layers above snow for radiative transfer
  t10                  [double] 10-day running mean of the 2 m temperature (K)
  tlai_z[nlevcan]       [double] pft total leaf area index for canopy layer
  vcmaxcintsha          [double] leaf to canopy scaling coefficient - shade
  vcmaxcintsun          [double] leaf to canopy scaling coefficient - sun
  parsha_z[nlevcan]     [double] par absorbed per unit lai for canopy layer (w/m**2) - shade
  parsun_z[nlevcan]     [double] par absorbed per unit lai for canopy layer (w/m**2) - sun
  laisha_z[nlevcan]     [double] leaf area index for canopy layer, sunlit or shaded - shade
  laisun_z[nlevcan]     [double] leaf area index for canopy layer, sunlit or shaded - sun


  OUTPUTS:
  btran            [double]  transpiration wetness factor (0 to 1)
  qflx_tran_veg_out    [double] vegetation transpiration (mm H2O/s) (+ = to atm)
  qflx_evap_veg_out    [double] vegetation evaporation (mm H2O/s) (+ = to atm)
  eflx_sh_veg_out      [double] sensible heat flux from leaves (W/m**2) [+ to atm]
  */
  void
  StabilityIteration(const LandType &Land, const double &forc_hgt_u_patch, const double &forc_hgt_t_patch,
                     const double &forc_hgt_q_patch, const double *dleaf, const double &fwet, const double &fdry,
                     const double &laisun, const double &laisha, const double &forc_rho_in, const double &snow_depth,
                     const double &soilbeta, const double &frac_h2osfc, const double &t_h2osfc, const double &sabv,
                     const double &h2ocan, const double &htop, const double t_soisno[nlevgrnd + nlevsno],

                     const VegProperties &veg, const int &nrad, const double &t10, const double tlai_z[nlevcan],
                     const double &vcmaxcintsha, const double &vcmaxcintsun, const double parsha_z[nlevcan],
                     const double parsun_z[nlevcan], const double laisha_z[nlevcan], const double laisun_z[nlevcan],

                     double &btran, double &qflx_tran_veg_out, double &qflx_evap_veg_out, double &eflx_sh_veg_out)

  // clm uses a decreasing index filter to act as a stop criteria
  // we operate on single cells, so we will replace the filter criteria with something else
  {
    forc_rho = forc_rho_in;
    if (!Land.lakpoi && !Land.urbpoi && frac_veg_nosno != 0) {
      bool stop = false;
      int itmax = 40; // maximum number of iteration [-]
      int itmin = 2;  // minimum number of iteration [-]
      int itlef = 0;
      double csoilc = 0.004; // Drag coefficient for soil under canopy [-]
      double ustar, del2, uaf, cf, rb, ram, rah[2], raw[2];
      double csoilcn, csoilb, ri, ricsoilc, w, svpts, eah;
      double wta, wtl, wtshi, wtg0, wtga, wtsqi, wtgq0;
      double snow_depth_c, fsno_dl, elai_dl, rdl, rppdry, efpot, rpp;
      double dc1, dc2, efsh, erre, err, efeold, lw_grnd, dels, ecidif;
      double tstar, qstar, thvstar, wc, zeta, wtaq, wtlq, dele, det;

      while (itlef <= itmax && !stop) {
        // Determine friction velocity, and potential temperature and humidity profiles of the surface boundary layer
        FrictionVelocityWind(forc_hgt_u_patch, displa, um, obu, z0mv, ustar);
        FrictionVelocityTemperature(forc_hgt_t_patch, displa, obu, z0hv, temp1);
        FrictionVelocityHumidity(forc_hgt_q_patch, forc_hgt_t_patch, displa, obu, z0hv, z0qv, temp1, temp2);
        FrictionVelocityTemperature2m(obu, z0hv, temp12m);
        FrictionVelocityHumidity2m(obu, z0hv, z0qv, temp12m, temp22m);

        // save leaf temp and leaf temp delta from previous iteration
        tlbef = t_veg;
        del2 = del;
        // Determine aerodynamic resistances
        ram = 1.0 / (ustar * ustar / um);
        rah[0] = 1.0 / (temp1 * ustar);
        raw[0] = 1.0 / (temp2 * ustar);
        // Bulk boundary layer resistance of leaves
        uaf = um * std::sqrt(1.0 / (ram * um));
        // Use pft parameter for leaf characteristic width dleaf
        cf = 0.01 / (std::sqrt(uaf) * std::sqrt(dleaf[Land.vtype]));
        rb = 1.0 / (cf * uaf);

        // Parameterization for variation of csoilc with canopy density from X. Zeng, University of Arizona
        w = exp(-(elai + esai));
        // changed by K.Sakaguchi from here
        // transfer coefficient over bare soil is changed to a local variable
        // just for readability of the code (from line 680)
        csoilb = (vkc / (0.13 * pow((z0mg * uaf / 1.5e-5), 0.45)));
        // compute the stability parameter for ricsoilc  ("S" in Sakaguchi&Zeng,2008)
        ri = (grav * htop * (taf - t_grnd)) / pow((taf * uaf), 2.0);

        // modify csoilc value (0.004) if the under-canopy is in stable condition
        if ((taf - t_grnd) > 0.0) {
          // decrease the value of csoilc by dividing it with (1+gamma*min(S, 10.0))
          // ria ("gmanna" in Sakaguchi&Zeng, 2008) is a constant (=0.5)
          ricsoilc = csoilc / (1.0 + 0.5 * std::min(ri, 10.0));
          csoilcn = csoilb * w + ricsoilc * (1.0 - w);
        } else {
          csoilcn = csoilb * w + csoilc * (1.0 - w);
        }

        // Sakaguchi changes for stability formulation ends here
        rah[1] = 1.0 / (csoilcn * uaf);
        raw[1] = rah[1];
        // Stomatal resistances for sunlit and shaded fractions of canopy.
        // Done each iteration to account for differences in eah, tv.
        svpts = el;                    // pa
        eah = forc_pbot * qaf / 0.622; // pa

        if (Land.vtype == nsoybean || Land.vtype == nsoybeanirrig) {
          btran = std::min(1.0, btran * 1.25);
        }

        // call photosynthesis (phase=sun)
        Photosynthesis(veg, Land.vtype, nrad, forc_pbot, t_veg, t10, svpts, eah, o2, co2, rb, btran, dayl_factor, thm,
                       tlai_z, vcmaxcintsun, parsun_z, laisun_z, rssun);

        if (Land.vtype == nsoybean || Land.vtype == nsoybeanirrig) {
          btran = std::min(1.0, btran * 1.25);
        }

        // call photosynthesis (phase=shade)
        Photosynthesis(veg, Land.vtype, nrad, forc_pbot, t_veg, t10, svpts, eah, o2, co2, rb, btran, dayl_factor, thm,
                       tlai_z, vcmaxcintsha, parsha_z, laisha_z, rssha);

        // Sensible heat conductance for air, leaf and ground
        wta = 1.0 / rah[0];       // air
        wtl = (elai + esai) / rb; // leaf
        wtg = 1.0 / rah[1];       // ground
        wtshi = 1.0 / (wta + wtl + wtg);
        wtl0 = wtl * wtshi; // leaf
        wtg0 = wtg * wtshi; // ground
        wta0 = wta * wtshi; // air
        wtga = wta0 + wtg0; // ground + air
        wtal = wta0 + wtl0; // air + leaf

        // Fraction of potential evaporation from leaf
        if (fdry > 0.0) {
          rppdry = fdry * rb * (laisun / (rb + rssun) + laisha / (rb + rssha)) / elai;
        } else {
          rppdry = 0.0;
        }

        efpot = forc_rho * wtl * (qsatl - qaf);
        if (efpot > 0.0) {
          if (btran > btran0) {
            qflx_tran_veg_out = efpot * rppdry;
            rpp = rppdry + fwet;
          } else {
            // No transpiration if btran below 1.e-10
            rpp = fwet;
            qflx_tran_veg_out = 0.0;
          }
          // Check total evapotranspiration from leaves
          rpp = std::min(rpp, (qflx_tran_veg_out + h2ocan / dtime) / efpot);
        } else {
          // No transpiration if potential evaporation less than zero
          rpp = 1.0;
          qflx_tran_veg_out = 0.0;
        }

        // Update conductances for changes in rpp
        // Latent heat conductances for ground and leaf.
        // Air has same conductance for both sensible and latent heat.
        wtaq = frac_veg_nosno / raw[0];                   // air
        wtlq = frac_veg_nosno * (elai + esai) / rb * rpp; // leaf
        // Litter layer resistance. Added by K.Sakaguchi
        snow_depth_c = 0.05; // critical depth for 100% litter burial by snow (=litter thickness) -- litter thickness
                             // hardwired as 0.05 m
        fsno_dl = snow_depth / snow_depth_c; // effective snow cover for (dry)plant litter
        elai_dl = 0.5 * (1.0 - std::min(fsno_dl,
                                        1.0)); // exposed (dry)litter area index -- liter area hardwired as 0.5 m^2 m^-2
        rdl = (1.0 - exp(-elai_dl)) / (0.004 * uaf); // dry litter layer resistance
        // add litter resistance and Lee and Pielke 1992 beta
        if (delq < 0.0) { // dew. Do not apply beta for negative flux (follow old rsoil)
          wtgq = frac_veg_nosno / (raw[1] + rdl);
        } else {
          wtgq = soilbeta * frac_veg_nosno / (raw[1] + rdl);
        }

        wtsqi = 1.0 / (wtaq + wtlq + wtgq);
        wtgq0 = wtgq * wtsqi;  // ground
        wtlq0 = wtlq * wtsqi;  // leaf
        wtaq0 = wtaq * wtsqi;  // air
        wtgaq = wtaq0 + wtgq0; // air + ground
        wtalq = wtaq0 + wtlq0; // air + leaf
        dc1 = forc_rho * cpair * wtl;
        dc2 = hvap * forc_rho * wtlq;
        efsh = dc1 * (wtga * t_veg - wtg0 * t_grnd - wta0 * thm);
        efe = dc2 * (wtgaq * qsatl - wtgq0 * qg - wtaq0 * forc_q);

        // Evaporation flux from foliage
        erre = 0.0;
        if (efe * efeb < 0.0) {
          efeold = efe;
          efe = 0.1 * efeold;
          erre = efe - efeold;
        }

        // fractionate ground emitted longwave
        lw_grnd = (frac_sno * pow(t_soisno[nlevsno - snl], 4.0) +
                   (1.0 - frac_sno - frac_h2osfc) * pow(t_soisno[nlevsno], 4.0) + frac_h2osfc * pow(t_h2osfc, 4.0));
        dt_veg = (sabv + air + bir * pow(t_veg, 4.0) + cir * lw_grnd - efsh - efe) /
                 (-4.0 * bir * pow(t_veg, 3.0) + dc1 * wtga + dc2 * wtgaq * qsatldT);
        t_veg = tlbef + dt_veg;
        dels = dt_veg;
        del = std::abs(dels);
        err = 0.0;
        if (del > 1.0) { // maxchange in  leaf temperature [K] --  hardwired as 1.0
          dt_veg = dels / del;
          t_veg = tlbef + dt_veg;
          err = sabv + air + bir * pow(tlbef, 3.0) * (tlbef + 4.0 * dt_veg) + cir * lw_grnd -
                (efsh + dc1 * wtga * dt_veg) - (efe + dc2 * wtgaq * qsatldT * dt_veg);
        }

        // Fluxes from leaves to canopy space
        // "efe" was limited as its sign changes frequently.  This limit may
        // result in an imbalance in "hvap*qflx_evap_veg" and
        // "efe + dc2*wtgaq*qsatdt_veg"
        efpot = forc_rho * wtl * (wtgaq * (qsatl + qsatldT * dt_veg) - wtgq0 * qg - wtaq0 * forc_q);
        qflx_evap_veg_out = rpp * efpot;
        // Calculation of evaporative potentials (efpot) and interception losses; flux in kg m**-2 s-1.
        // ecidif holds the excess energy if all intercepted water is evaporated during the timestep.
        // This energy is later added to the sensible heat flux.
        ecidif = 0.0;
        if (efpot > 0.0 && btran > btran0) {
          qflx_tran_veg_out = efpot * rppdry;
        } else {
          qflx_tran_veg = 0.0;
        }
        ecidif = std::max(0.0, qflx_evap_veg_out - qflx_tran_veg_out - h2ocan / dtime);
        qflx_evap_veg_out = std::min(qflx_evap_veg_out, qflx_tran_veg_out + h2ocan / dtime);
        qflx_evap_veg = qflx_evap_veg_out;
        qflx_tran_veg = qflx_tran_veg_out;
        // The energy loss due to above two limits is added to the sensible heat flux.
        eflx_sh_veg_out = efsh + dc1 * wtga * dt_veg + err + erre + hvap * ecidif;
        eflx_sh_veg = eflx_sh_veg_out;
        // Re-calculate saturated vapor pressure, specific humidity, and their derivatives at the leaf surface
        QSat(t_veg, forc_pbot, el, deldT, qsatl, qsatldT);

        // Update vegetation/ground surface temperature, canopy air
        // temperature, canopy vapor pressure, aerodynamic temperature, and
        // Monin-Obukhov stability parameter for next iteration.
        taf = wtg0 * t_grnd + wta0 * thm + wtl0 * t_veg;
        qaf = wtlq0 * qsatl + wtgq0 * qg + forc_q * wtaq0;
        // Update Monin-Obukhov length and wind speed including the stability effect
        dth = thm - taf;
        dqh = forc_q - qaf;
        delq = wtalq * qg - wtlq0 * qsatl - wtaq0 * forc_q;
        tstar = temp1 * dth;
        qstar = temp2 * dqh;
        thvstar = tstar * (1.0 + 0.61 * forc_q) + 0.61 * forc_th * qstar;
        zeta = zldis * vkc * grav * thvstar / (pow(ustar, 2.0) * thv);
        if (zeta >= 0.0) { // stable
          zeta = std::min(2.0, std::max(zeta, 0.01));
          um = std::max(ur, 0.1);
        } else {
          zeta = std::max(-100.0, std::min(zeta, -0.01));
          wc = pow((-grav * ustar * thvstar * 1000.0 / thv),
                   0.333); // convective boundary layer height hardwired as 1000.0
          um = std::sqrt(ur * ur + wc * wc);
        }
        obu = zldis / zeta;
        if (obuold * obu < 0.0) {
          nmozsgn += 1;
        }
        if (nmozsgn >= 4) {
          obu = zldis / (-0.01);
        }
        obuold = obu;

        // Test for convergence
        itlef += 1;
        if (itlef > itmin) {
          dele = std::abs(efe - efeb);
          efeb = efe;
          det = std::max(del, del2);
          if ((det < 0.01) && (dele < 0.1)) {
            stop = true;
          }
        }
      } // stability iteration
    }   // land type
  }     // CanopyFluxes::StabilityIteration()

  /*
  INPUTS:
  Land             [LandType] struct containing information about landtype
  t_soisno[nlevgrnd+nlevsno] [double] col soil temperature (Kelvin)
  frac_h2osfc                [double] fraction of ground covered by surface water (0 to 1)
  t_h2osfc                   [double] surface water temperature
  sabv               [double] solar radiation absorbed by vegetation (W/m**2)
  qg_snow             [double] specific humidity at snow surface [kg/kg]
  qg_soil             [double] specific humidity at soil surface [kg/kg]
  qg_h2osfc           [double] specific humidity at h2osfc [kg/kg]
  dqgdT                      [double] d(qg)/dT
  htvp                       [double] latent heat of vapor of water (or sublimation) [j/kg]


  OUTPUTS:
  h2ocan                     [double] canopy water (mm H2O)
  t_veg_out          [double]  vegetation temperature (Kelvin)
  eflx_sh_grnd        [double] sensible heat flux from ground (W/m**2) [+ to atm]
  eflx_sh_snow        [double] sensible heat flux from snow (W/m**2) [+ to atm]
  eflx_sh_soil        [double] sensible heat flux from soil (W/m**2) [+ to atm]
  eflx_sh_h2osfc      [double] sensible heat flux from h2osfc (W/m**2) [+ to atm]

  qflx_evap_soi       [double] soil evaporation (mm H2O/s) (+ = to atm)
  qflx_ev_snow        [double] evaporation flux from snow (W/m**2) [+ to atm]
  qflx_ev_soil        [double] evaporation flux from soil (W/m**2) [+ to atm]
  qflx_ev_h2osfc      [double] evaporation flux from h2osfc (W/m**2) [+ to atm]

  dlrad               [double] downward longwave radiation below the canopy [W/m2]
  ulrad               [double] upward longwave radiation above the canopy [W/m2]
  cgrnds                     [double] deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
  cgrndl                     [double] deriv of soil latent heat flux wrt soil temp [w/m**2/k]
  cgrnd                      [double] deriv. of soil energy flux wrt to soil temp [w/m2/k]
  t_ref2m                    [double]  2 m height surface air temperature (Kelvin)
  t_ref2m_r                  [double]  Rural 2 m height surface air temperature (Kelvin)
  q_ref2m                    [double]  2 m height surface specific humidity (kg/kg)
  rh_ref2m_r                 [double]  Rural 2 m height surface relative humidity (%)
  rh_ref2m                   [double]  2 m height surface relative humidity (%)
  */
  void ComputeFlux(const LandType &Land, const double t_soisno[nlevgrnd + nlevsno], const double &frac_h2osfc,
                   const double &t_h2osfc, const double &sabv, const double &qg_snow, const double &qg_soil,
                   const double &qg_h2osfc, const double &dqgdT, const double &htvp,

                   double &h2ocan, double &t_veg_out, double &eflx_sh_grnd, double &eflx_sh_snow, double &eflx_sh_soil,
                   double &eflx_sh_h2osfc,

                   double &qflx_evap_soi, double &qflx_ev_snow, double &qflx_ev_soil, double &qflx_ev_h2osfc,

                   double &dlrad, double &ulrad, double &cgrnds, double &cgrndl, double &cgrnd, double &t_ref2m,
                   double &t_ref2m_r, double &q_ref2m, double &rh_ref2m, double &rh_ref2m_r) {

    if (!Land.lakpoi && !Land.urbpoi && frac_veg_nosno != 0) {

      t_veg_out = t_veg;
      double e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT;

      // Energy balance check in canopy
      double lw_grnd =
          (frac_sno * pow(t_soisno[nlevsno - snl], 4.0) + (1.0 - frac_sno - frac_h2osfc) * pow(t_soisno[nlevsno], 4.0) +
           frac_h2osfc * pow(t_h2osfc, 4.0));
      double err = sabv + air + bir * pow(tlbef, 3.0) * (tlbef + 4.0 * dt_veg) + cir * lw_grnd - eflx_sh_veg -
                   hvap * qflx_evap_veg;

      // Fluxes from ground to canopy space
      double delt = wtal * t_grnd - wtl0 * t_veg - wta0 * thm;
      eflx_sh_grnd = cpair * forc_rho * wtg * delt;

      // compute individual sensible heat fluxes
      double delt_snow = wtal * t_soisno[nlevsno - snl] - wtl0 * t_veg - wta0 * thm;
      eflx_sh_snow = cpair * forc_rho * wtg * delt_snow;
      double delt_soil = wtal * t_soisno[nlevsno] - wtl0 * t_veg - wta0 * thm;
      eflx_sh_soil = cpair * forc_rho * wtg * delt_soil;
      double delt_h2osfc = wtal * t_h2osfc - wtl0 * t_veg - wta0 * thm;
      eflx_sh_h2osfc = cpair * forc_rho * wtg * delt_h2osfc;
      qflx_evap_soi = forc_rho * wtgq * delq;

      // compute individual latent heat fluxes
      double delq_snow = wtalq * qg_snow - wtlq0 * qsatl - wtaq0 * forc_q;
      qflx_ev_snow = forc_rho * wtgq * delq_snow;
      double delq_soil = wtalq * qg_soil - wtlq0 * qsatl - wtaq0 * forc_q;
      qflx_ev_soil = forc_rho * wtgq * delq_soil;
      double delq_h2osfc = wtalq * qg_h2osfc - wtlq0 * qsatl - wtaq0 * forc_q;
      qflx_ev_h2osfc = forc_rho * wtgq * delq_h2osfc;

      // 2 m height air temperature
      t_ref2m = thm + temp1 * dth * (1.0 / temp12m - 1.0 / temp1);
      t_ref2m_r = t_ref2m;
      // 2 m height specific humidity
      q_ref2m = forc_q + temp2 * dqh * (1.0 / temp22m - 1.0 / temp2);
      // 2 m height relative humidity
      QSat(t_ref2m, forc_pbot, e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT);
      rh_ref2m = std::min(100.0, (q_ref2m / qsat_ref2m) * 100.0);
      rh_ref2m_r = rh_ref2m;

      // Downward longwave radiation below the canopy
      dlrad = (1.0 - emv) * emg * forc_lwrad + emv * emg * sb * pow(tlbef, 3.0) * (tlbef + 4.0 * dt_veg);
      // Upward longwave radiation above the canopy
      ulrad = ((1.0 - emg) * (1.0 - emv) * (1.0 - emv) * forc_lwrad +
               emv * (1.0 + (1.0 - emg) * (1.0 - emv)) * sb * pow(tlbef, 3.0) * (tlbef + 4.0 * dt_veg) +
               emg * (1.0 - emv) * sb * lw_grnd);
      // Derivative of soil energy flux with respect to soil temperature
      cgrnds = cgrnds + cpair * forc_rho * wtg * wtal;
      cgrndl = cgrndl + forc_rho * wtgq * wtalq * dqgdT;
      cgrnd = cgrnds + cgrndl * htvp;
      // Update dew accumulation (kg/m2)
      h2ocan = std::max(0.0, h2ocan + (qflx_tran_veg - qflx_evap_veg) * dtime);

      // Determine total photosynthesis -- need to implement -- vars don't get used - diagnostics, maybe??
      // PhotosynthesisTotal(fn, filterp, atm2lnd_vars, cnstate_vars, canopystate_vars, photosyns_vars)

      // evaluate error and write?? -- best way to write?
      //  if (abs(err(p)) > 0.1_r8)
    }
  }
};

} // namespace ELM
