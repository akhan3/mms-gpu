#ifndef _ODE_FUNCTIONS_H_
#define _ODE_FUNCTIONS_H_

#include "mydefs.hpp"
#include "Vector3.hpp"


int Hfield (    const Vector3 *M, Vector3 *H, fptype *charge, fptype *potential,
                const int xdim, const int ydim, const int zdim, const fptype meshwidth,
                const fptype mu_0, const fptype Ms, const fptype Aexch,
                const int demag, const int exchange, const int external, const int use_fmm, const int P,
                const int use_gpu, const int verbosity );

int time_marching(  byte *material, Vector3 *M, // initial state. This will be overwritten in each time step
                    const fptype finaltime, const fptype timestep,
                    const int xdim, const int ydim, const int zdim, const fptype meshwidth, const int P,
                    const fptype mu_0, const fptype Ms, const fptype Aexch, const fptype alfa, const fptype gamma,
                    const int coupling, const int exchange, const int external, const int use_fmm,
                    const int use_gpu, const int verbosity );

#endif // #ifndef  _ODE_FUNCTIONS_H_
