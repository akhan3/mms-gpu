#ifndef _ODE_FUNCTIONS_H_
#define _ODE_FUNCTIONS_H_

#include "Vector3.hpp"


int time_marching(  Vector3 *M, // initial state. This will be overwritten in each time step
                    const float finaltime,
                    const int xdim, const int ydim, const int zdim, const float meshwidth,
                    const float mu_0, const float Ms, const float Aexch, const float alfa, const float gamma,
                    const int verbose_level );

#endif // #ifndef  _ODE_FUNCTIONS_H_
