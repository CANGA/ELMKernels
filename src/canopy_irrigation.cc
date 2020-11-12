/*
DESCRIPTION:
Determine whether we're irrigating here; set qflx_irrig appropriately.

INPUTS:
irrig_rate         [double] current irrigation rate (applied if n_irrig_steps_left > 0) [mm/s]

INOUT:
n_irrig_steps_left [int] number of time steps for which we still need to irrigate today

OUTPUTS:
qflx_irrig         [double] irrigation amount (mm/s)
*/

void CanopyIrrigation(const double& irrig_rate,
    int& n_irrig_steps_left,
    double& qflx_irrig)
{
  if (n_irrig_steps_left > 0) {
     qflx_irrig = irrig_rate;
     n_irrig_steps_left -= 1;
   } else {
     qflx_irrig = 0.0;
   }
}
