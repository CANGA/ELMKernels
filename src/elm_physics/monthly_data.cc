
#include "monthly_data.h"

namespace ELM {

double MonthlyDataManager::month_frac(const Utils::Date& model_time) {
  static constexpr auto ndaypm = std::array<int, 12>{31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}; // days per month
  const int kmo = std::get<1>(model_time.date()); // month
  const int kda = std::get<2>(model_time.date()) - 1; // day -
  double t = (kda + model_time.sec/86400.0) / ndaypm[kmo-1]; // elapsed fraction of current month
  return t;
}
  
std::pair<int, int> MonthlyDataManager::month_indices(const Utils::Date& model_time) {
  const double t = month_frac(model_time);
  const int t1 = t + 0.5; // implicit rounding - t1 = 0 if t < 0.5 and t1 = 1 if t >= 0.5
  const int t2 = t1 + 1;
  const int kmo = std::get<1>(model_time.date()); // month
  int m1 = kmo + t1 - 1;
  int m2 = kmo + t2 - 1;
  if (m1 < 1) {
    m1 = 12;
  }
  if (m2 > 12) {
    m2 = 1;
  }
  m1 -= 1;
  m2 -= 1;
  return std::make_pair(m1, m2);
}
  
std::pair<double, double> MonthlyDataManager::monthly_data_weights(const Utils::Date& model_time) {
  const double t = month_frac(model_time);
  const int t1 = t + 0.5; // implicit rounding - t1 = 0 if t < 0.5 and t1 = 1 if t >= 0.5
  double wt1 = (t1 + 0.5) - t; // 
  return std::make_pair(wt1, 1.0 - wt1);
}

} // namespace ELM
