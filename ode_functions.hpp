#ifndef _ODE_FUNCTIONS_H_
#define _ODE_FUNCTIONS_H_

#include "mydefs.hpp"
#include "Vector3.hpp"


int time_marching(  byte *material, Vector3 *M, // initial state. This will be overwritten in each time step
                    const fptype finaltime, const fptype timestep,
                    const int xdim, const int ydim, const int zdim, const fptype meshwidth, const int P,
                    const fptype mu_0, const fptype Ms, const fptype Aexch, const fptype alfa, const fptype gamma,
                    const int coupling, const int exchange, const int external, const int use_fmm,
                    const int verbose_level );

#endif // #ifndef  _ODE_FUNCTIONS_H_
