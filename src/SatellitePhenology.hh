#pragma once

#include "read_input.hh"
#include "utils.hh"
namespace ELM {

/*
Read monthly vegetation data for two consec. months.
*/
template <typename Array_t>
void readMonthlyVegetation(const std::string &data_dir, const std::string &basename_phen, const Utils::Date &time,
                           const Utils::DomainDecomposition<2> &dd, Array_t &mlai, Array_t &msai, Array_t &mhgtt,
                           Array_t &mhgtb) {

  const int n_months = 2;

  // -- read lai
  ELM::IO::read_and_reshape_phenology(data_dir, basename_phen, "MONTHLY_LAI", time, n_months, dd, mlai);
  // -- read sai
  ELM::IO::read_and_reshape_phenology(data_dir, basename_phen, "MONTHLY_SAI", time, n_months, dd, msai);
  // -- read max veg height
  ELM::IO::read_and_reshape_phenology(data_dir, basename_phen, "MONTHLY_HEIGHT_TOP", time, n_months, dd, mhgtt);
  // -- read min veg height
  ELM::IO::read_and_reshape_phenology(data_dir, basename_phen, "MONTHLY_HEIGHT_BOT", time, n_months, dd, mhgtb);
}

/*
Determine if 2 new months of data are to be read.
*/
template <typename ArrayI1, typename ArrayD1, typename ArrayD2>
void interpMonthlyVeg(const Utils::Date &time_start, int dtime, int n_pfts, const std::string &dir,
                      const std::string &basename, const Utils::DomainDecomposition<2> &dd, const ArrayI1 &vtype,
                      ArrayD1 &timwt, ArrayD2 &mlai2t, ArrayD2 &msai2t, ArrayD2 &mhvt2t, ArrayD2 &mhvb2t) {
  int it[2];     //  month 1 and month 2 (step 1)
  int months[2]; //  months to be interpolated (1 to 12)
  int veg_type;
  static int InterpMonths1 = -999;                                                         // saved month index
  const auto ndaypm = std::array<int, 12>{31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}; // days per month
  auto n_grid_cells = dd.n_local[0] * dd.n_local[1];

  Utils::Date time_end(time_start);
  time_end.increment_seconds(dtime);
  auto date = time_end.date();
  int kmo = std::get<1>(date); // month
  int kda = std::get<2>(date); // day

  double t = (kda - 0.5) / ndaypm[kmo - 1];
  it[0] = t + 0.5;
  it[1] = it[0] + 1;
  months[0] = kmo + it[0] - 1;
  months[1] = kmo + it[1] - 1;

  if (months[0] < 1) {
    months[0] = 12;
  }
  if (months[1] > 12) {
    months[1] = 1;
  }

  timwt[0] = (it[0] + 0.5) - t;
  timwt[1] = 1.0 - timwt[0];

  if (InterpMonths1 != months[0]) {
    // allocate 2 months of LAI, SAI, HTOP, HBOT
    Array<double, 3> mlai(2, n_grid_cells, n_pfts);
    Array<double, 3> msai(2, n_grid_cells, n_pfts);
    Array<double, 3> mhgtt(2, n_grid_cells, n_pfts);
    Array<double, 3> mhgtb(2, n_grid_cells, n_pfts);

    readMonthlyVegetation(dir, basename, time_end, dd, mlai, msai, mhgtt, mhgtb);
    InterpMonths1 = months[0];

    for (int i = 0; i < n_grid_cells; ++i) {
      veg_type = vtype[i];
      for (int month = 0; month < 2; ++month) {
        mlai2t(i, month) = mlai(month, i, veg_type);
        msai2t(i, month) = msai(month, i, veg_type);
        mhvt2t(i, month) = mhgtt(month, i, veg_type);
        mhvb2t(i, month) = mhgtb(month, i, veg_type);
      }
    }
  }
}

template <typename ArrayD1>
void SatellitePhenology(const ArrayD1 &mlai, const ArrayD1 &msai, const ArrayD1 &mhvt, const ArrayD1 &mhvb,
                        const ArrayD1 &timwt, const int &vtype, const double &snow_depth, const double &frac_sno,

                        double &tlai, double &tsai, double &htop, double &hbot, double &elai, double &esai,
                        int &frac_veg_nosno_alb) {

  double ol, fb;

  // need to update elai and esai only every albedo time step so do not
  // have any inconsistency in lai and sai between SurfaceAlbedo calls (i.e.,
  // if albedos are not done every time step).
  // leaf phenology
  // Set leaf and stem areas based on day of year
  // Interpolate leaf area index, stem area index, and vegetation heights
  // between two monthly
  // The weights below (timwt(1) and timwt(2)) were obtained by a call to
  // routine InterpMonthlyVeg in subroutine NCARlsm.
  //                 Field   Monthly Values
  //                -------------------------
  // leaf area index LAI  <- mlai1 and mlai2
  // leaf area index SAI  <- msai1 and msai2
  // top height      HTOP <- mhvt1 and mhvt2
  // bottom height   HBOT <- mhvb1 and mhvb2

  tlai = timwt[0] * mlai[0] + timwt[1] * mlai[1];
  tsai = timwt[0] * msai[0] + timwt[1] * msai[1];
  htop = timwt[0] * mhvt[0] + timwt[1] * mhvt[1];
  hbot = timwt[0] * mhvb[0] + timwt[1] * mhvb[1];

  // adjust lai and sai for burying by snow. if exposed lai and sai
  // are less than 0.05, set equal to zero to prevent numerical
  // problems associated with very small lai and sai.
  // snow burial fraction for short vegetation (e.g. grasses) as in
  // Wang and Zeng, 2007.

  if (vtype > noveg && vtype <= nbrdlf_dcd_brl_shrub) {
    ol = std::min(std::max(snow_depth - hbot, 0.0), htop - hbot);
    fb = 1.0 - ol / std::max(1.e-06, htop - hbot);
  } else {
    fb = 1.0 - std::max(std::min(snow_depth, 0.2), 0.0) / 0.2; // 0.2m is assumed
    // depth of snow required for complete burial of grasses
  }

  // area weight by snow covered fraction
  elai = std::max(tlai * (1.0 - frac_sno) + tlai * fb * frac_sno, 0.0);
  esai = std::max(tsai * (1.0 - frac_sno) + tsai * fb * frac_sno, 0.0);
  if (elai < 0.05) {
    elai = 0.0;
  }
  if (esai < 0.05) {
    esai = 0.0;
  }

  // Fraction of vegetation free of snow
  if ((elai + esai) >= 0.05) {
    frac_veg_nosno_alb = 1;
  } else {
    frac_veg_nosno_alb = 0;
  }
}

} // namespace ELM
