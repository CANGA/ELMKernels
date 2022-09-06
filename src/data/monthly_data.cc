
#include "monthly_data.h"
#include <array>
#include <string>

namespace ns = ELM::monthly_data;

double ns::month_frac(const Utils::Date& model_time)
{
  static constexpr std::array<int, 12> ndaypm{31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}; // days per month
  const int kmo{std::get<1>(model_time.date())};                                              // month
  const int kda{std::get<2>(model_time.date()) - 1};                                          // day -
  double t = (kda + model_time.sec / 86400.0) / ndaypm[kmo - 1]; // elapsed fraction of current month
  return t;
}

int ns::first_month_idx(const Utils::Date& model_time)
{
  const double t{month_frac(model_time)};
  const int t1 = t < 0.5 ? 0 : 1; // t1 = 0 if t < 0.5 and t1 = 1 if t >= 0.5
  const int kmo{std::get<1>(model_time.date())}; // month
  int m1 = kmo + t1 - 2;
  if (m1 < 0) {
    m1 = 11;
  }
  return m1;
}

std::pair<int, int> ns::month_indices(const Utils::Date& model_time)
{
  const int m1{first_month_idx(model_time)};
  int m2 = m1 + 1;
  if (m2 > 11) {
    m2 = 0;
  }
  return std::make_pair(m1, m2);
}

int ns::third_month_idx(const Utils::Date& model_time)
{
  const auto [m1, m2] = month_indices(model_time);
  int m3 = m2 + 1;
  if (m3 > 11) {
    m3 = 0;
  }
  return m3;
}

std::tuple<int, int, int> ns::triple_month_indices(const Utils::Date& model_time)
{
  const auto [m1, m2] = month_indices(model_time);
  const auto m3{third_month_idx(model_time)};
  return std::make_tuple(m1, m2, m3);
}

std::pair<double, double> ns::monthly_data_weights(const Utils::Date& model_time)
{
  const double t{month_frac(model_time)};
  const int t1 = t < 0.5 ? 0 : 1; // t1 = 0 if t < 0.5 and t1 = 1 if t >= 0.5
  double wt1 = (t1 + 0.5) - t;
  return std::make_pair(wt1, 1.0 - wt1);
}
