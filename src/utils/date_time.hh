#ifndef ELM_UTILS_DATE_TIME_HH_
#define ELM_UTILS_DATE_TIME_HH_

#include <iomanip>

namespace ELM {
namespace Utils {

inline int to_doy(int month, int day) {
  const auto dy_per_mo = std::array<int, 12>{31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

  int doy = 0;
  for (int i = 0; i != month - 1; ++i) {
    doy += dy_per_mo[i];
  }
  doy += day - 1;
  return std::move(doy);
}

inline std::tuple<int, int, int> to_date(int year, int doy) {
  const auto dy_per_mo = std::array<int, 12>{31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

  assert(doy < 365 && "ELM::Utils::to_date() expects day_of_year < 365");
  assert(doy >= 0 && "ELM::Utils::to_date() expects day_of_year >= 0");
  int month = 0;
  while (doy >= 0) {
    doy -= dy_per_mo[month];
    ++month;
  }

  doy += dy_per_mo[month - 1];
  int day = doy + 1;
  return std::make_tuple(year, month, day);
}

//
// A no-leap calendar.
struct Date {
  int year;
  int doy;
  int sec;
  static const int sec_per_day = 86400;

  Date(int year_, int doy_) : year(year_), doy(doy_), sec(0) {}

  Date(int year_, int month, int day) : Date(year_, to_doy(month, day)) {}

  Date(int year_, int month, int day, int seconds) : Date(year_, to_doy(month, day)) { this->increment_seconds(seconds); }

  Date(int year_) : Date(year_, 0) {}
  Date() : Date(0, 0) {}

  Date(const Date &other) = default;
  Date &operator=(const Date &other) = default;

  // year,month,day
  std::tuple<int, int, int> date() const { return to_date(year, doy); }

  //
  // increment
  //
  Date &increment_day(size_t days = 1) { return operator+=((int)days); }
  Date &increment_year(size_t years = 1) {
    year += (int)years;
    return *this;
  }
  Date &increment_month(size_t months = 1) {
    const auto dy_per_mo = std::array<int, 12>{31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    for (size_t m = 0; m != months; ++m) {
      auto d = date();
      operator+=(dy_per_mo[std::get<1>(d) - 1]);
    }
    return *this;
  }

  Date &increment_seconds(size_t seconds) {
    sec += (int)seconds;
    while (sec >= sec_per_day) {
      operator+=(1);
      sec -= sec_per_day;
    }
    return *this;
  }

  Date &decrement_seconds(size_t seconds) {
    sec -= (int)seconds;
    while (sec < 0) {
      operator-=(1);
      sec += sec_per_day;
    }
    return *this;
  }

  // pre-increment
  Date &operator++() { return this->operator+=(1); }

  // post-increment
  Date operator++(int) {
    Date old(*this);
    operator++();
    return old;
  }

  // pre-decrement
  Date &operator--() { return this->operator+=(-1); }

  // post-decrement
  Date operator--(int) {
    Date old(*this);
    operator--();
    return old;
  }

  //
  // add
  //
  Date &operator+=(int days) {
    doy += days;
    while (doy >= 365) {
      ++year;
      doy -= 365;
    }

    while (doy < 0) {
      --year;
      doy += 365;
    }
    return *this;
  }

  //
  // subtract
  //
  Date &operator-=(int days) { return operator+=(-days); }
};

//
// Date operators
//

//
// stream
//
inline std::ostream &operator<<(std::ostream &out, const Date &d) {
  auto date = d.date();
  out << std::setfill('0') << std::setw(4) << std::get<0>(date) << "-" << std::get<1>(date) << "-" << std::get<2>(date);
  return out;
}

//
// math
//
inline Date operator+(Date lhs, int rhs) {
  lhs += rhs;
  return lhs;
}
inline Date operator-(Date lhs, int rhs) {
  lhs -= rhs;
  return lhs;
}

//
// decimal Julian date
//
inline double decimal_doy(int month, int day, int seconds) {
  //int doy = to_doy(month, day);
  //double d_doy = doy + seconds/86400.0;
  return to_doy(month, day) + seconds/86400.0;
}

inline double decimal_doy(const Date& date) {
  //double doy = (double)date.doy;
  //double d_ doy += date.sec/86400.0;
  //return std::move(doy);
  return static_cast<double>(date.doy) + static_cast<double>(date.sec)/86400.0;
}

//
// Subtract
//
// NOTE: adding two dates isn't supported, but diffing two dates is ok!
//
inline double days_since(const Date &lhs, const Date &rhs) {
  int d_year = lhs.year - rhs.year;
  double d_doy = decimal_doy(lhs) - decimal_doy(rhs);
  return d_doy + 365.0 * d_year;
}

inline int months_since(const Date &lhs, const Date &rhs) {
  auto lhs_d = lhs.date();
  auto rhs_d = rhs.date();
  int d_year = std::get<0>(lhs_d) - std::get<0>(rhs_d);
  int d_mo = std::get<1>(lhs_d) - std::get<1>(rhs_d);
  return d_mo + 12 * d_year;
}

inline int years_since(const Date &lhs, const Date &rhs) { return lhs.year - rhs.year; }

inline double operator-(const Date &lhs, const Date &rhs) { return days_since(lhs, rhs); }

//
// relational
//
inline bool operator<(const Date &lhs, const Date &rhs) {
  // orders by year first, doy second, sec last
  return std::tie(lhs.year, lhs.doy, lhs.sec) < std::tie(rhs.year, rhs.doy, rhs.sec);
}
inline bool operator>(const Date &lhs, const Date &rhs) { return rhs < lhs; }
inline bool operator<=(const Date &lhs, const Date &rhs) { return !(lhs > rhs); }
inline bool operator>=(const Date &lhs, const Date &rhs) { return !(lhs < rhs); }

inline bool operator==(const Date &lhs, const Date &rhs) { return lhs.year == rhs.year && lhs.doy == rhs.doy && lhs.sec == rhs.sec; }
inline bool operator!=(const Date &lhs, const Date &rhs) { return !(lhs == rhs); }

//
// Ticker class handles timesteps.
//
// NOTE: this implies all timesteps are at most a day!
//
struct Ticker {
  Date start;
  int days;
  int ticks;
  unsigned short ticks_per_day;

  Ticker(const Date &start_, unsigned short ticks_per_day_)
      : start(start_), days(0), ticks(0), ticks_per_day(ticks_per_day_) {}

  Ticker(unsigned short ticks_per_day_) : start(), days(0), ticks(0), ticks_per_day(ticks_per_day_) {}

  Ticker(const Ticker &other) = default;
  Ticker &operator=(const Ticker &other) = default;

  Date now() const { return start + days; }

  int ticks_since() const { return ticks + ticks_per_day * days; }

  int days_since() const { return days; }

  //
  // increment
  //
  // pre-increment
  Ticker &operator++() { return this->operator+=(1); }

  // post-increment
  Ticker operator++(int) {
    Ticker old(*this);
    operator++();
    return old;
  }

  // pre-decrement
  Ticker &operator--() { return this->operator+=(-1); }

  // post-decrement
  Ticker operator--(int) {
    Ticker old(*this);
    operator--();
    return old;
  }

  //
  // add
  //
  Ticker &operator+=(int d_ticks) {
    ticks += d_ticks;

    while (ticks >= ticks_per_day) {
      ++days;
      ticks -= ticks_per_day;
    }

    while (ticks < 0) {
      --days;
      ticks += ticks_per_day;
    }
    return *this;
  }

  //
  // subtract
  //
  Ticker &operator-=(int d_ticks) { return operator+=(-d_ticks); }
};

//
// stream
//
inline std::ostream &operator<<(std::ostream &out, const Ticker &d) {
  out << d.start << "_" << d.ticks;
  return out;
}

} // namespace Utils
} // namespace ELM

#endif
