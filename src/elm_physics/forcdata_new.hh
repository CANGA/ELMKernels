#include <type_traits>
#include <cmath>
#include <string>

#include <iostream>

#include "elm_constants.h"
#include "array.hh"
#include "date_time.hh"

namespace forc_data_manager {

enum class forcDataType { TBOT, PBOT, QBOT, FLDS, FSDS, PREC, WIND, ZBOT };

template<forcDataType type>
using forc_tag = std::integral_constant<forcDataType, type>;





constexpr double tdc(const double &t) {
  return std::min(50.0, std::max(-50.0, (t - ELM::tfrz)));
}

constexpr double esatw(const double &t) {
  constexpr double a[7] = {6.107799961,     4.436518521e-01, 1.428945805e-02, 2.650648471e-04,
                       3.031240396e-06, 2.034080948e-08, 6.136820929e-11};
  return 100.0 * (a[0] + t * (a[1] + t * (a[2] + t * (a[3] + t * (a[4] + t * (a[5] + t * a[6]))))));
}

constexpr double esati(const double &t) {
  constexpr double b[7] = {6.109177956,     5.034698970e-01, 1.886013408e-02, 4.176223716e-04,
                       5.824720280e-06, 4.838803174e-08, 1.838826904e-10};
  return 100.0 * (b[0] + t * (b[1] + t * (b[2] + t * (b[3] + t * (b[4] + t * (b[5] + t * b[6]))))));
}

void forc_rh_to_qbot(const double& forc_tbot, const double& forc_pbot, const double& forc_rh, double& forc_qbot) {
  double e = (forc_tbot > ELM::tfrz) ? esatw(tdc(forc_tbot)) : esati(tdc(forc_tbot));
  double qsat = 0.622 * e / (forc_pbot - 0.378 * e);
  forc_qbot = qsat * forc_rh / 100.0;
}

void forc_qbot_to_rh(const double& forc_tbot, const double& forc_pbot, const double& forc_qbot, double& forc_rh) {
  double e = (forc_tbot > ELM::tfrz) ? esatw(tdc(forc_tbot)) : esati(tdc(forc_tbot));
  double qsat = 0.622 * e / (forc_pbot - 0.378 * e);
  forc_rh = 100.0 * (forc_qbot / qsat);
}

 // rho, pO2, pCO2

  void derive_forc_vp (const double& forc_qbot, const double& forc_pbot, double& forc_vp) 
  { forc_vp = forc_qbot * forc_pbot / (0.622 + 0.378 * forc_qbot); }

  void derive_forc_rho (const double& forc_pbot, const double& forc_vp, const double& forc_tbot, double& forc_rho) 
  { forc_rho = (forc_pbot - 0.378 * forc_vp) / (ELM::rair * forc_tbot); }

  void derive_forc_po2 (const double& forc_pbot, double& forc_po2) 
  { forc_po2 = ELM::o2_molar_const * forc_pbot; }

  void derive_forc_pco2 (const double& forc_pbot, double& forc_pco2) 
  { forc_pco2 = ELM::co2_ppmv * 1.0e-6 * forc_pbot; }





constexpr double interp_forcing(double wt1, double wt2, double forc1, double forc2) {
  return forc1 * wt1 + forc2 * wt2;
}


template<typename ArrayD1, typename ArrayD2>
struct ProcessTBOT {
  ProcessTBOT(const int& t_idx, const double& wt1, const double& wt2, const ArrayD2& atm_tbot, ArrayD1& forc_tbot) : 
  t_idx_(t_idx), wt1_(wt1), wt2_(wt2), atm_tbot_(atm_tbot), forc_tbot_(forc_tbot) {}

  constexpr void operator()(const int i) const {
    forc_tbot_(i) = std::min(interp_forcing(wt1_, wt2_, atm_tbot_(t_idx_, i), atm_tbot_(t_idx_+1, i)), 323.0);
    std::cout << "inside tbot()" << std::endl;
  }
private:
  int t_idx_;
  double wt1_, wt2_;
  ArrayD2 atm_tbot_;
  ArrayD1 forc_tbot_;
};


template<typename ArrayD1, typename ArrayD2>
struct ProcessPBOT {
  ProcessPBOT(const int& t_idx, const double& wt1, const double& wt2, const ArrayD2& atm_pbot, ArrayD1& forc_pbot) : 
  t_idx_(t_idx), wt1_(wt1), wt2_(wt2), atm_pbot_(atm_pbot), forc_pbot_(forc_pbot) {}

  constexpr void operator()(const int i) const {
    forc_pbot_(i) = std::max(interp_forcing(wt1_, wt2_, atm_pbot_(t_idx_, i), atm_pbot_(t_idx_+1, i)), 4.0e4);
    std::cout << "inside pbot()" << std::endl;
  }
private:
  int t_idx_;
  double wt1_, wt2_;
  ArrayD2 atm_pbot_;
  ArrayD1 forc_pbot_;
};


template<typename ArrayD1, typename ArrayD2>
struct ProcessQBOT {
  ProcessQBOT(const int& t_idx, const double& wt1, const double& wt2, const ArrayD2& atm_qbot, ArrayD1& forc_qbot) : 
  t_idx_(t_idx), wt1_(wt1), wt2_(wt2), atm_qbot_(atm_qbot), forc_qbot_(forc_qbot) {}

  constexpr void operator()(const int i) const {
    forc_qbot_(i) = std::max(interp_forcing(wt1_, wt2_, atm_qbot_(t_idx_, i), atm_qbot_(t_idx_+1, i)), 1.0e-9);
    std::cout << "inside qbot()" << std::endl;
  }
private:
  int t_idx_;
  double wt1_, wt2_;
  ArrayD2 atm_qbot_;
  ArrayD1 forc_qbot_;
};


template<typename ArrayD1, typename ArrayD2>
struct ProcessFLDS {
  ProcessFLDS(
    const int& t_idx, const double& wt1, const double& wt2, const ArrayD2& atm_flds, const ArrayD1& forc_pbot, const ArrayD1& forc_qbot, 
    const ArrayD1& forc_tbot, ArrayD1& forc_lwrad) 
  : t_idx_(t_idx), wt1_(wt1), wt2_(wt2), atm_flds_(atm_flds), forc_pbot_(forc_pbot),
    forc_qbot_(forc_qbot), forc_tbot_(forc_tbot), forc_lwrad_(forc_lwrad) {}

  constexpr void operator()(const int i) const {
    const double flds = interp_forcing(wt1_, wt2_, atm_flds_(t_idx_, i), atm_flds_(t_idx_+1, i));
    if (flds <= 50.0 || flds >= 600.0) {
      const double e = forc_pbot_(i) * forc_qbot_(i) / (0.622 + 0.378 * forc_qbot_(i));
      const double ea = 0.70 + 5.95e-5 * 0.01 * e * exp(1500.0 / forc_tbot_(i));
      forc_lwrad_(i) = ea * ELM::ELM_STEBOL * pow(forc_tbot_(i), 4.0);
    } else {
      forc_lwrad_(i) = flds;
    }
    std::cout << "inside flds()" << std::endl;
  }
private:
  int t_idx_;
  double wt1_, wt2_;
  ArrayD2 atm_flds_;
  ArrayD1 forc_pbot_, forc_qbot_, forc_tbot_;
  ArrayD1 forc_lwrad_;
};


template<typename ArrayD1, typename ArrayD2>
struct ProcessFSDS {
  ProcessFSDS(const ArrayD1& atm_fsds, const ArrayD1& coszen, ArrayD2& forc_solai, ArrayD2& forc_solad) : 
  atm_fsds_(atm_fsds), coszen_(coszen), forc_solai_(forc_solai), forc_solad_(forc_solad) {}

  constexpr void operator()(const int i) const {
    // need to impement model for coszen factor
    // ELM uses fac = (cosz > 0.001) ? min(cosz/avg_forc_cosz, 10) : 0.0
    // ATS uses a slope based factor
    // ELM's method could probably be calculated outside the parallel region
    // ATS's method should probably be calculated inside parallel region
    const double swndr = std::max(atm_fsds_(i) * coszen_(i) * 0.5, 0.0);
    const double& swndf = swndr;
    const double& swvdr = swndr; // these vars are only used with a specific forcing data stream
    const double& swvdf = swndr; // maybe implement later? placeholders for now
    const double ratio_rvrf_vis = std::min(
        0.99, std::max(0.17639 + 0.00380 * swvdr - 9.0039e-06 * pow(swvdr, 2.0) + 8.1351e-09 * pow(swvdr, 3.0), 0.01));
    const double ratio_rvrf_nir = std::min(
        0.99, std::max(0.29548 + 0.00504 * swndr - 1.4957e-05 * pow(swndr, 2.0) + 1.4881e-08 * pow(swndr, 3.0), 0.01));
    forc_solad_(i,0) = ratio_rvrf_vis * swvdr;
    forc_solad_(i,1) = ratio_rvrf_nir * swndr;
    forc_solai_(i,0) = (1.0 - ratio_rvrf_vis) * swvdf;
    forc_solai_(i,1) = (1.0 - ratio_rvrf_nir) * swndf;
    std::cout << "inside fsds()" << std::endl;
  }
private:
  ArrayD1 atm_fsds_, coszen_;
  ArrayD2 forc_solai_, forc_solad_;
};


template<typename ArrayD1>
struct ProcessPREC {
  ProcessPREC(const ArrayD1& atm_prec, const ArrayD1& forc_tbot, 
    ArrayD1& forc_rain, ArrayD1& forc_snow) 
  : atm_prec_(atm_prec), forc_tbot_(forc_tbot), 
    forc_rain_(forc_rain), forc_snow_(forc_snow) {}

  constexpr void operator()(const int i) const {
    double frac = (forc_tbot_(i) - ELM::tfrz) * 0.5; // ramp near freezing
    frac = std::min(1.0, std::max(0.0, frac)); // bound in [0,1]
    forc_rain_(i) = frac * std::max(atm_prec_(i), 0.0);
    forc_snow_(i) = (1.0 - frac) * std::max(atm_prec_(i), 0.0);
  }
private:
  ArrayD1 atm_prec_, forc_tbot_;
  ArrayD1 forc_rain_, forc_snow_;
};


template<typename ArrayD1, typename ArrayD2>
struct ProcessWIND {
  ProcessWIND(const int& t_idx, const double& wt1, const double& wt2, const ArrayD2& atm_wind, ArrayD1& forc_u, ArrayD1& forc_v) : 
  t_idx_(t_idx), wt1_(wt1), wt2_(wt2), atm_wind_(atm_wind), forc_u_(forc_u), forc_v_(forc_v) {}

  constexpr void operator()(const int i) const {
    forc_u_(i) = interp_forcing(wt1_, wt2_, atm_wind_(t_idx_, i), atm_wind_(t_idx_+1, i));
    forc_v_(i) = 0.0;
  }
private:
  int t_idx_;
  double wt1_, wt2_;
  ArrayD2 atm_wind_;
  ArrayD1 forc_u_, forc_v_;
};


template<typename ArrayD1>
struct ProcessZBOT {
  //ProcessZBOT(const int& t_idx, const double& wt1, const double& wt2, const ArrayD2& atm_hgt, ArrayD1& forc_hgt, ArrayD1& forc_hgt_u, ArrayD1& forc_hgt_t, ArrayD1& forc_hgt_q) : 
  ProcessZBOT(ArrayD1& forc_hgt, ArrayD1& forc_hgt_u, ArrayD1& forc_hgt_t, ArrayD1& forc_hgt_q) 
  : forc_hgt_(forc_hgt), forc_hgt_u_(forc_hgt_u), forc_hgt_t_(forc_hgt_t), forc_hgt_q_(forc_hgt_q) {}

  constexpr void operator()(const int i) const {
    forc_hgt_(i) = 30.0; // hardwired? what about zbot from forcing file?
    forc_hgt_u_(i) = forc_hgt_(i);  // observational height of wind [m]
    forc_hgt_t_(i) = forc_hgt_(i);  // observational height of temperature [m]
    forc_hgt_q_(i) = forc_hgt_(i);  // observational height of humidity [m]
  }
private:
  ArrayD1 forc_hgt_, forc_hgt_u_, forc_hgt_t_, forc_hgt_q_;
};








template<typename ArrayD1, typename ArrayD2, forcDataType type>
struct ForcDataBase {


  //ForcDataBase(const std::string &filename, const std::string &varname, const Utils::Date &file_start_time, int ntimes, int ncells) :
  //    fname_(filename), file_start_time_(file_start_time), ntimes_(ntimes),
  //    ncells_(ncells), data_(varname, ntimes_, ncells_)  { }


ForcDataBase(const std::string& varname, int ntimes, int ncells) :
    ntimes_(ntimes), ncells_(ncells), data_(varname, ntimes, ncells, 9915.0)  { }



void update_chunk_info(const ELM::Utils::Date& new_chunk_start_time, int ntimes) {
    chunk_start_time_ = new_chunk_start_time;
    ntimes_ = ntimes;
  }



// need to copy into forc_thbot
void process_data(forc_tag<forcDataType::TBOT>, const double& atm_tbot, double& forc_tbot) const {
  forc_tbot = std::min(atm_tbot, 323.0);
}

void process_data(forc_tag<forcDataType::PBOT>, const double& atm_pbot, double& forc_pbot) const {
  forc_pbot = std::max(atm_pbot, 4.0e4);
}


void process_data(forc_tag<forcDataType::QBOT>, const double& atm_qbot, double& forc_qbot) const {
  forc_qbot = std::max(atm_qbot, 1.0e-9);
}

constexpr void process_data(forc_tag<forcDataType::FLDS>, const double& atm_flds, double& forc_lwrad, const double& forc_pbot, 
  const double& forc_qbot, const double& forc_tbot) const {
  if (atm_flds <= 50.0 || atm_flds >= 600.0) {
    const double e = forc_pbot * forc_qbot / (0.622 + 0.378 * forc_qbot);
    const double ea = 0.70 + 5.95e-5 * 0.01 * e * exp(1500.0 / forc_tbot);
    forc_lwrad = ea * ELM::ELM_STEBOL * pow(forc_tbot, 4.0);
  } else {
    forc_lwrad = atm_flds;
  }
}


void process_data(forc_tag<forcDataType::FSDS>, const double& atm_fsds, ArrayD1 forc_solad, ArrayD1 forc_solai, 
  const double& coszen) const {
  // need to impement model for coszen factor
  // ELM uses fac = (cosz > 0.001) ? min(cosz/avg_forc_cosz, 10) : 0.0
  // ATS uses a slope based factor
  // ELM's method could probably be calculated outside the parallel region
  // ATS's method should probably be calculated inside parallel region
  const double swndr = std::max(atm_fsds * coszen * 0.5, 0.0);
  const double& swndf = swndr;
  const double& swvdr = swndr; // these vars are only used with a specific forcing data stream
  const double& swvdf = swndr; // maybe implement later? placeholders for now
  const double ratio_rvrf_vis = std::min(
      0.99, std::max(0.17639 + 0.00380 * swvdr - 9.0039e-06 * pow(swvdr, 2.0) + 8.1351e-09 * pow(swvdr, 3.0), 0.01));
  const double ratio_rvrf_nir = std::min(
      0.99, std::max(0.29548 + 0.00504 * swndr - 1.4957e-05 * pow(swndr, 2.0) + 1.4881e-08 * pow(swndr, 3.0), 0.01));
  forc_solad[0] = ratio_rvrf_vis * swvdr;
  forc_solad[1] = ratio_rvrf_nir * swndr;
  forc_solai[0] = (1.0 - ratio_rvrf_vis) * swvdf;
  forc_solai[1] = (1.0 - ratio_rvrf_nir) * swndf;
}



void process_data(forc_tag<forcDataType::PREC>, const double& atm_prec, double& forc_rain, double& forc_snow, 
  const double& forc_tbot) const {
  double frac = (forc_tbot - ELM::tfrz) * 0.5;       // ramp near freezing
  frac = std::min(1.0, std::max(0.0, frac)); // bound in [0,1]
  forc_rain = frac * std::max(atm_prec, 0.0);
  forc_snow = (1.0 - frac) * std::max(atm_prec, 0.0);
}

void process_data(forc_tag<forcDataType::WIND>, const double& atm_wind, double& forc_u, double& forc_v) const {
  forc_u = atm_wind;
  forc_v = 0.0;
}

void process_data(forc_tag<forcDataType::ZBOT>, const double& atm_zbot, double& forc_hgt_u, double& forc_hgt_t, 
  double& forc_hgt_q) const {
  // all hardwired at 30 m as ELM default
  // it seems like many forcing files provide a different value, so build it into reader
  // the value of atm_zbot is likely constant throughout the domain, so this could probably just be called for a single variable and shared
  // need to trace the forc_hgt_* params and see if they ever change
  forc_hgt_u = 30.0;  // observational height of wind [m]
  forc_hgt_t = 30.0;  // observational height of temperature [m]
  forc_hgt_q = 30.0;  // observational height of humidity [m]
}







// need to copy into forc_thbot
//void process_data(forc_tag<forcDataType::TBOT>, const int& i, const int& t_idx, 
//  const double& wt1, const double& wt2, double& forc_tbot) const {
//  forc_tbot = std::min(wt1 * data_(t_idx,i) + wt2 * data_(t_idx+1,i), 323.0);
//}










double elapsed_forc_dt(const ELM::Utils::Date& date, const ELM::Utils::Date& prior_date) const {
  double delta_days = ELM::Utils::days_since(date, prior_date);
  assert(delta_days >= 0.0 && "difference in dates is negative");
  return delta_days / forc_dt_;
}



std::pair<double,double> forcing_time_weights(const int& t_idx, const ELM::Utils::Date& model_time) const {
  ELM::Utils::Date forc_start(chunk_start_time_);
  forc_start.increment_seconds(static_cast<int>(86400 * forc_dt_) * t_idx);
  double elapsed_dt = elapsed_forc_dt(model_time, forc_start);
  assert(elapsed_dt <= 1.0 && "model_time is greater than forc_t_start + forc_dt");
  double wt1 = 1.0 - elapsed_dt;
  double wt2 = 1.0 - wt1;
  return std::make_pair(wt1, wt2);
}


template < typename T, typename Index >
using subscript_t = decltype(std::declval<T>()[std::declval<Index>()]);

template < typename, typename, typename = void >
struct has_subscript : std::false_type {};

template < typename T, typename Index >
struct has_subscript< T, Index, std::void_t< subscript_t<T,Index> > > : std::true_type {};

//template<typename Args>
//constexpr bool has_op

//template<typename... Args>
//void get_current_forcing_data(Args&&...args) { process_data(forc_tag<type>{}, std::forward<Args>(args)...); }

// Variable template that checks if a type has begin() and end() member functions
//template <typename, typename = void>
//constexpr bool is_iterable{};
// 
//template <typename T>
//constexpr bool is_iterable<T, std::void_t< decltype(std::declval<T>().begin()), decltype(std::declval<T>().end()) >> = true;

template<typename Args>
constexpr decltype(auto) unpack_elems(const int& i, Args&& arg) const {
  if constexpr (std::is_same_v<std::decay_t<Args>, ArrayD2>) {
    // or kokkos subview
    return arg[i];
  } else if constexpr (std::is_same_v<std::decay_t<Args>, ArrayD1>) {
    return std::forward<decltype(std::forward<Args>(arg)[i])>(std::forward<Args>(arg)[i]);
  } else {
    return std::forward<Args>(arg);
  }
}

//template<typename Args>
//constexpr decltype(auto) unpack_elems(const int& i, Args&& arg) const {
//  if constexpr (has_subscript<Args>::value) {
//    return std::forward<decltype(arg[i])>(arg[i]);
//  } else {
//    return std::forward<Args>(arg);
//  }
//}


//const double data_interp = data_(t_idx, i) * wts.first + data_(t_idx+1, i) * wts.second;


template<typename... Args>
void get_forcing(const ELM::Utils::Date& model_time, Args&&...args) const {
  const int t_idx = static_cast<int>(elapsed_forc_dt(model_time, chunk_start_time_));

  if constexpr (type == forcDataType::FSDS || type == forcDataType::PREC) {
    for (int i = 0; i != ncells_; ++i) {
      process_data(forc_tag<type>{}, data_(t_idx, i), unpack_elems(i, std::forward<Args>(args))...);
    }
  } else {
    const auto& [wt1, wt2] = forcing_time_weights(t_idx, model_time);
    for (int i = 0; i != ncells_; ++i) {
      const double data_interp = data_(t_idx, i) * wt1 + data_(t_idx+1, i) * wt2;
      process_data(forc_tag<type>{}, data_interp, unpack_elems(i, std::forward<Args>(args))...);
    }
  }
}


//template<typename... Args>
//void get_forcing_test(const ELM::Utils::Date& model_time, Args&&...args) const {
//  const int t_idx = static_cast<int>(elapsed_forc_dt(model_time, chunk_start_time_));
//  const auto& wts = forcing_time_weights(t_idx, model_time);
//  if constexpr(type == forcDataType::TBOT) {
//    ProcessTBOT process_tbot(t_idx, wts.first, wts.second, data_, std::forward<Args>(args)...);
//    for (int i = 0; i != ncells_; ++i) {
//      process_tbot<forcDataType::TBOT>(i);
//      std::cout << "tbot" << std::endl;
//    //invoke1(t_idx, wts.first, wts.second, data_, test);
//    }
//  }
//}


 



template<typename... Args>
constexpr void get_forcing_test(const ELM::Utils::Date& model_time, Args&&...args) {
  const int t_idx = static_cast<int>(elapsed_forc_dt(model_time, chunk_start_time_));
  const auto& wts = forcing_time_weights(t_idx, model_time);

  auto myObj = [&] {
    if constexpr(type == forcDataType::TBOT) {
      return ProcessTBOT(t_idx, wts.first, wts.second, data_, std::forward<Args>(args)...);
    } else if constexpr (type == forcDataType::PBOT) {
      return ProcessPBOT(t_idx, wts.first, wts.second, data_, std::forward<Args>(args)...);
    } else if constexpr (type == forcDataType::QBOT) {
      return ProcessQBOT(t_idx, wts.first, wts.second, data_, std::forward<Args>(args)...);
    } else if constexpr (type == forcDataType::FLDS) {
      return ProcessFLDS(t_idx, wts.first, wts.second, data_, std::forward<Args>(args)...);
    } else if constexpr (type == forcDataType::FSDS) {
      return ProcessFSDS(data_[t_idx], std::forward<Args>(args)...);
    } else if constexpr (type == forcDataType::PREC) {
      return ProcessPREC(data_[t_idx], std::forward<Args>(args)...);
    } else if constexpr (type == forcDataType::WIND) {
      return ProcessWIND(t_idx, wts.first, wts.second, data_, std::forward<Args>(args)...);
    } else if constexpr (type == forcDataType::ZBOT) {
      return ProcessZBOT(std::forward<Args>(args)...);
    }
  }();

  for (int i = 0; i != ncells_; ++i) {
    //myObj(i);
    std::invoke(myObj, i);
  }
}

protected:
ArrayD2 data_;
int ntimes_{0}, ncells_{0};
double forc_dt_{0.125};
ELM::Utils::Date chunk_start_time_;
ELM::Utils::Date chunk_end_time_;

};






} // namespace forc_data_manager
