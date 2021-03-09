struct SnowData {
  double dz[5] = {0.0};
  double z[5] = {0.0};
  double zi[6] = {0.0};

  int snl;
  double snow_depth;
  double swe_old[5], double h2osoi_liq[5], double h2osoi_ice[5], double t_soisno[5], double frac_iceold[5], Array_d dz,
      Array_d z, Array_d zi,
};

struct SurfRadData {
  double sabg_lyr[6]; // sabg_lyr has values for all snowpack layers (sabg_lyr[1,..snl]) AND the first soil cell
                      // (sabg_lyr[0]);
};