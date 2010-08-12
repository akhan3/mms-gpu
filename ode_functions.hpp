#ifndef _ODE_FUNCTIONS_H_
#define _ODE_FUNCTIONS_H_

#include "mydefs.hpp"
#include "Vector3.hpp"


int time_marching(  Vector3 *M, // initial state. This will be overwritten in each time step
                    const fptype finaltime,
                    const int xdim, const int ydim, const int zdim, const fptype meshwidth,
                    const fptype mu_0, const fptype Ms, const fptype Aexch, const fptype alfa, const fptype gamma,
                    const int verbose_level );

#endif // #ifndef  _ODE_FUNCTIONS_H_
