

#pragma once

#include "date_time.hh"

#include <utility>
#include <tuple>

/*

monthly weights and indexing is reliant on implicit type conversions - could be replaced by round()
assumes data is defined at the center of the month:  t_data = thismonth_start + thismonth_length/2
when t < t_data (t < midpoint of current month):  t_indices = prev_month, curr_month
when t >= t_data (t >= midpoint of current month):  t_indices = curr_month, next_month

the time weights change in concert with the indices
with t = elapsed fraction of current month
weights start at (0.5,0.5) for t=0 and approach (0,1) as t approaches 0.5
weights switch to (1,0) when t=0.5 and approach (0.5,0.5) as t approaches 1


INDICES:
           t<0.5    t>=0.5
kmo        m1 m2    m1 m2
 1 (0)     11 0     0  1
 2 (1)     0  1     1  2
 3 (2)     1  2     2  3
 4 (3)     2  3     3  4
 5 (4)     3  4     4  5
 6 (5)     4  5     5  6
 7 (6)     5  6     6  7
 8 (7)     6  7     7  8
 9 (8)     7  8     8  9
 10 (9)    8  9     9  10
 11 (10)   9  10    10 11
 12 (11)   10 11    11 0


                switch to next
 WEIGHTS:         index pair
(0.5, 0.5) -> (0,1)   |  (1,0) -> (0.5, 0.5)  
    t=0      t=0.25       t=0.5       t=0.75
 wt1  wt2   wt1   wt2    wt1   wt2   wt1   wt2
 0.5  0.5   0.25  0.75   1.0   0.0   0.75  0.25
*/

namespace ELM::monthly_data {

// utility namespace for monthly data (phenology, aerosols)
// month_frac() - returns elapsed month fraction
// month_indices() - return monthly data timeseries indices bracketing model_time input parameter
// monthly_data_weights() - return weights for monthly data 

double month_frac(const Utils::Date& model_time);
int first_month_idx(const Utils::Date& model_time);
std::pair<int, int> month_indices(const Utils::Date& model_time);
int third_month_idx(const Utils::Date& model_time);
std::tuple<int,int,int> triple_month_indices(const Utils::Date& model_time);
std::pair<double, double> monthly_data_weights(const Utils::Date& model_time);


} // namespace ELM::monthly_data
