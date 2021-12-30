
const int forcdim = 3;
struct ForcDataTimeManager {

ForcDataTimeManager(const Date& )

}

enum class forcDataType { TBOT, PBOT, QBOT, FSDS, FLDS, PRECTmms, WIND, ZBOT };

template<forcDataType type>
using forc_tag = std::integral_constant<forcDataType, type>;

double tdc(const double &t) {
  return std::min(50.0, std::max(-50.0, (t - tfrz)));
}

double esatw(const double &t) {
  const double a[7] = {6.107799961,     4.436518521e-01, 1.428945805e-02, 2.650648471e-04,
                       3.031240396e-06, 2.034080948e-08, 6.136820929e-11};
  return 100.0 * (a[0] + t * (a[1] + t * (a[2] + t * (a[3] + t * (a[4] + t * (a[5] + t * a[6]))))));
}

double esati(const double &t) {
  const double b[7] = {6.109177956,     5.034698970e-01, 1.886013408e-02, 4.176223716e-04,
                       5.824720280e-06, 4.838803174e-08, 1.838826904e-10};
  return 100.0 * (b[0] + t * (b[1] + t * (b[2] + t * (b[3] + t * (b[4] + t * (b[5] + t * b[6]))))));
}

void forc_rh_to_qbot(const double& forc_t, const double& forc_pbot, const double& forc_rh, double& forc_qbot) {
  double e = (forc_t > tfrz) ? esatw(tdc(forc_t)) : esati(tdc(forc_t));
  double qsat = 0.622 * e / (forc_pbot - 0.378 * e);
  forc_qbot = qsat * forc_rh / 100.0;
}


void forc_qbot_to_rh() {


}


template<>
void ForcDataBase<ArrayD2, ForcType::TBOT>::process_data(const int& i, const int& t_idx, 
  const double& wt1, const double& wt2, double& forc_tbot) const {
  forc_tbot = std::min(wt1 * data_(t_idx,i) + wt2 * data_(t_idx+1,i), 323.0);
}

template<>
void ForcDataBase<ArrayD2, ForcType::PBOT>::process_data(const int& i, const int& t_idx, 
  const double& wt1, const double& wt2, double& forc_pbot) const {
  forc_pbot = std::max(wt1 * data_(t_idx,i) + wt2 * data_(t_idx+1,i), 4.0e4);
}

template<>
void ForcDataBase<ArrayD2, ForcType::QBOT>::process_data(const int& i, const int& t_idx, 
  const double& wt1, const double& wt2, double& forc_qbot) const {
  forc_qbot = std::max(wt1 * data_(t_idx,i) + wt2 * data_(t_idx+1,i), 1.0e-9);
}

template<>
void ForcDataBase<ArrayD2, ForcType::RH>::process_data(const int& i, const int& t_idx, 
  const double& wt1, const double& wt2, double& forc_rh) const {
  forc_rh = std::max(wt1 * data_(t_idx,i) + wt2 * data_(t_idx+1,i), 1.0e-9);
}

template<>
void ForcDataBase<ArrayD2, ForcType::FLDS>::process_data(const int& i, const int& t_idx, 
  const double& wt1, const double& wt2, double& forc_lwrad) const {
  forc_lwrad = wt1 * data_(t_idx,i) + wt2 * data_(t_idx+1,i);
}


void process_forc_t()

template<typename ArrayD2>
struct ForcDataBase {


  ForcDataBase(const std::string &filename, const std::string &varname, const Utils::Date &file_start_time, int ntimes, int ncells) :
      fname_(filename), file_start_time_(file_start_time), ntimes_(ntimes),
      ncells_(ncells), data_(varname, ntimes_, ncells_)  { }




//std::pair<double,double> forcing_time_weights(const double& model_time, const double& forc_t_start) {
//  assert(model_time >= forc_t_start && "(model_time >= forc_t_start) is false");
//  assert(model_time <= forc_t_start + forc_dt_ && "(model_time <= forc_t_start + forc_dt) is false");
//  double wt1 = 1.0 - (model_time - forc_t_start) / forc_dt_;
//  double wt2 = 1.0 - wt1;
//  return std::make_pair(wt1, wt2);
//}

std::pair<double,double> forcing_time_weights(const int& forc_idx, const Utils::Date& model_time) {
  ELM::Utils::Date start(chunk_start_time_);
  start.increment_seconds((int)(86400 * forc_dt_) * forc_idx );
  double elapsed_ratio = ELM::Utils::days_since(model_time, start) / forc_dt_;
  assert(elapsed_ratio >= 0.0 && "model_time is less than forc_t_start");
  assert(elapsed_ratio <= 1.0 && "model_time is greater than forc_t_start + forc_dt");
  double wt1 = 1.0 - elapsed_ratio;
  double wt2 = 1.0 - wt1;
  return std::make_pair(wt1, wt2);
}

  void process_timestep(const Utils::Date& model_time, ArrayD1& forc_processed) {
    int forc_idx = current_time_idx(model_time);
    //ELM::Utils::Date start(chunk_start_time_);
    //double forc_start_time =
    auto wt = forcing_time_weights(forc_idx, model_time);
    process_data()

  }
  void process_data();
  void set_varname();


  void update_file_info(const Utils::Date& new_file_start_time, const std::string& new_filename) {
    file_start_time_ = new_file_start_time;
    fname_ = new_filename;
  }

  void update_chunk_info(const Utils::Date& new_chunk_start_time, int ntimes) {
    chunk_start_time_ = new_chunk_start_time;
    chunk_end_time_ = ntimes * forc_dt_ + new_chunk_start_time;
  }


  // serial I/O function
  // make generic interface for arbitrary timesteps - either resize array time dimension and use that for read:
  // read_atm_forcing(const std::string &filename, const std::string &varname, const Utils::Date &time,
  // const Utils::DomainDecomposition<2> &dd)    - or -
  // use a parameter n_times
  // read_atm_forcing(const std::string &filename, const std::string &varname, const Utils::Date &time,
  // const Utils::DomainDecomposition<2> &dd, const int& n_times)

  // no, just use this and calc start[] and count[] outside
  //read_atm_forcing(const std::string &filename, const std::string &varname,
  //const std::array<size_t, D> &start, const std::array<size_t, D> &count);

  void read_atm_forcing(const std::string& filename, const std::string& varname,
  const Utils::DomainDecomposition<2>& dd, const Utils::Date& model_time, int nsteps) {
    if (model_time >= chunk_end_time_) { // check here for now
      assert(nsteps <= ntimes_ && "forcing input: nsteps > array size")
      int offset = get_file_offset(model_time);
      std::array<GO, 3> start = {offset, dd.start[0], dd.start[1]};
      std::array<GO, 3> count = {offset+nsteps, dd.n_local[0], dd.n_local[1]};


    }



  }

  // get forc_dt from this read?
  void read_atm_forcing(const std::string& data_dir, const std::string& basename_atm, const Utils::Date& time,
                        const Utils::DomainDecomposition<2>& dd, const int n_months) {
    IO::read_and_reshape_forcing(data_dir, basename_atm, "TBOT", time, n_months, dd, this->atm_data_);
    // need new read_forc method - template off size_t D from incoming start/count arrays
  }


//int get_file_offset (const Date& model_time) const {
//  double decimal_day_offset = ELM::Utils::days_since(model_time, file_start_time_);
//  assert(decimal_day_offset >= 0.0 && "decimal_day_offset < 0.0")
//  return (int)(decimal_day_offset / forc_dt);
//}


int current_time_idx (const Date& model_time) const {
  double decimal_day_offset = ELM::Utils::days_since(model_time, chunk_start_time_);
  assert(decimal_day_offset >= 0.0 && "decimal_day_offset < 0.0")
  return (int)(decimal_day_offset / forc_dt_);
}
// these can be combined/generalized later
int file_start_idx(const Date& model_time) const {
  //return (int)((model_time - offset)/forc_dt);
  double decimal_day_offset = ELM::Utils::days_since(model_time, file_start_time_);
  assert(decimal_day_offset >= 0.0 && "decimal_day_offset < 0.0")
  return (int)(decimal_day_offset / forc_dt_);
}


  template <size_t D>
  std::array<int, D> get_data_dims(const Date& model_time, int ntimes) {
    int offset = get_file_offset(const Date& model_time);
    start = {offset, dd.start[0], dd.start[1]};
    count = {offset + ntimes, dd.n_local[0], dd.n_local[1]};




    , const std::array<size_t, D> &count,
  }


protected:
ArrayD2 data_;
Date file_start_time_, file_end_time_, chunk_start_time_, chunk_end_time_;
std::string fname_;
double forc_dt_{0.0};
int ntimes_{0}, ncells_{0};//, file_start_idx{0};


};
