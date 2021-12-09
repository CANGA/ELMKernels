

#include "mpi_types.hh"
#include "read_input.hh"
#include "utils.hh"
#include <string.h>

namespace ELM {

template<class ArrayD1, class ArrayD2>
struct ForcData {

const ArrayD2 atm_tbot; // air temp
const ArrayD2 atm_pbot; // pressure
const ArrayD2 atm_qbot; // specific humidity
const ArrayD2 atm_fsds; // shortwave radiation
const ArrayD2 atm_prec; // precip
const ArrayD2 atm_wind; // wind speed
const ArrayD2 atm_flds; // longwave radiation
//const ArrayD2 atm_zbot;


//const ArrayD2 atm_rh;







//ArrayD1 forc_t
//ArrayD1 forc_th
//ArrayD1 forc_pbot
//ArrayD1 forc_q
//ArrayD1 forc_lwrad
//ArrayD1 forc_u
//ArrayD1 forc_v
//ArrayD1 forc_hgt
//ArrayD2 forc_solai
//ArrayD2 forc_solad

//void resize_arrays() {
//atm_zbot
//atm_tbot
//atm_rh
//atm_qbot
//atm_wind
//atm_fsds
//atm_flds
//atm_psrf
//atm_prec 
//}

  ForcData(int ncells, int ntimes) : ncells_(ncells), ntimes_(ntimes), 
  atm_zbot("atm_zbot", ntimes_, ncells_), atm_tbot("atm_tbot", ntimes_, ncells_), atm_rh("atm_rh", ntimes_, ncells_), 
  atm_qbot("atm_qbot", ntimes_, ncells_), atm_wind("atm_wind", ntimes_, ncells_), atm_fsds("atm_fsds", ntimes_, ncells_), 
  atm_flds("atm_flds", ntimes_, ncells_), atm_pbot("atm_pbot", ntimes_, ncells_), 
  atm_prec("atm_prec", ntimes_, ncells_)  {}

  // default constructor
  // shouldn't be used
  ForcData() {};
  
  // default destructor
  ~ForcData() {};

  // serial I/O function
  //template <class ArrayD2>
  void read_atm_forcing(const std::string &data_dir, const std::string &basename_atm, const Utils::Date &time,
                      const Utils::DomainDecomposition<2> &dd, const int n_months) {
  
    
    IO::read_and_reshape_forcing(data_dir, basename_atm, "TBOT", time, n_months, dd, this->atm_tbot);
    IO::read_and_reshape_forcing(data_dir, basename_atm, "PSRF", time, n_months, dd, this->atm_pbot);
    IO::read_and_reshape_forcing(data_dir, basename_atm, "QBOT", time, n_months, dd, this->atm_qbot);
    IO::read_and_reshape_forcing(data_dir, basename_atm, "FSDS", time, n_months, dd, this->atm_fsds);
    IO::read_and_reshape_forcing(data_dir, basename_atm, "PRECTmms", time, n_months, dd, this->atm_prec);
    IO::read_and_reshape_forcing(data_dir, basename_atm, "WIND", time, n_months, dd, this->atm_wind);
    IO::read_and_reshape_forcing(data_dir, basename_atm, "FLDS", time, n_months, dd, this->atm_flds);

    //IO::read_and_reshape_forcing(data_dir, basename_atm, "RH", time, n_months, dd, atm_rh);
    //IO::read_and_reshape_forcing(data_dir, basename_atm, "ZBOT", time, n_months, dd, atm_zbot);
  }


  void interp_forc_tbot(const double& model_time, ArrayD1 lsm_forc) {
    const double& offset = this->offset_tbot_;
    const double& forc_dt = this->dt_tbot_;
    auto atm_forc = this->atm_tbot;
    interp_forcing(model_time, offset, forc_dt, atm_forc, lsm_forc);
  }


  void interp_forc_pbot(const double& model_time, ArrayD1 lsm_forc) {
    const double& offset = this->offset_pbot_;
    const double& forc_dt = this->dt_pbot_;
    auto atm_forc = this->atm_pbot;
    interp_forcing(model_time, offset, forc_dt, atm_forc, lsm_forc);
  }


  void interp_forc_t(const double& model_time, ArrayD1 lsm_forc) {
    const double& offset = this->offset_t_;
    const double& forc_dt = this->dt_t_;
    auto atm_forc = this->atm_tbot;
    interp_forcing(model_time, offset, forc_dt, atm_forc, lsm_forc);
  }


  void interp_forc_t(const double& model_time, ArrayD1 lsm_forc) {
    const double& offset = this->offset_t_;
    const double& forc_dt = this->dt_t_;
    auto atm_forc = this->atm_tbot;
    interp_forcing(model_time, offset, forc_dt, atm_forc, lsm_forc);
  }


  void interp_forc_t(const double& model_time, ArrayD1 lsm_forc) {
    const double& offset = this->offset_t_;
    const double& forc_dt = this->dt_t_;
    auto atm_forc = this->atm_tbot;
    interp_forcing(model_time, offset, forc_dt, atm_forc, lsm_forc);
  }

int forc_start_idx(const double& model_time, const double& offset, const double& forc_dt) const {
  return (int)((model_time - offset)/forc_dt);
}

double forc_start_time(const int& t_idx, const double& forc_dt) const {
  return forc_dt * t_idx;
}

std::pair<double,double> forcing_time_weights(const double& model_time, const double& forc_t_start, const double& forc_dt) {
  assert(model_time >= forc_t_start && "(model_time >= forc_t_start) is false");
  assert(model_time <= forc_t_start + forc_dt && "(model_time <= forc_t_start + forc_dt) is false");
  double wt1 = 1.0 - (model_time - forc_t_start) / forc_dt;
  double wt2 = 1.0 - wt1;
  return std::make_pair(wt1, wt2);
}

// interpolate forcing data at t + dt/2
typename ArrayD2::value_type interp_forcing_array (const int& i, const int& t_idx, const double& wt1, const double& wt2, const ArrayD2& atm_forc) {
  return wt1 * atm_forc(t_idx, i) + wt2 * atm_forc(t_idx+1, i);
}

void interp_forcing(const double& model_time, const double& offset, const double& forc_dt, const ArrayD2 atm_forc, ArrayD1 lsm_forc) {
  int t_idx = forc_start_idx(model_time, offset, forc_dt);
  auto wt = forcing_time_weights(model_time, forc_start_time(t_idx, forc_dt), forc_dt);
  // this should be parallelized
  for (int i = 0; i < lsm_forc.extent(0); ++i) {
    lsm_forc(i) = interp_forcing_array(i, t_idx, wt.first, wt.second, atm_forc);
  }
}

private:
  double offset_tbot_, dt_tbot_;
  double offset_pbot_, dt_pbot_;
  double offset_qbot_, dt_qbot_;
  double offset_lwrad_, dt_lwrad_;
  double offset_u_, dt_u_;
  double offset_v_, dt_v_;
  double offset_hgt_, dt_hgt_;
  double offset_solai_, dt_solai_;
  double offset_solad_, dt_solad_;

int ntimes_, ncells_;



};

} // namespace ELM

