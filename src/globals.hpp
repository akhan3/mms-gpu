#ifndef _GLOBALS_H_
#define _GLOBALS_H_

#include <cmath>
#include <gflags/gflags.h>

// physical constants
const double mu_0 = 4 * M_PI * 1e-7; // magnetic permeability of vacuum
const double hbar = 1.05457148e-34; // Reduced Planck's constant
const double e = 1.60217646e-19; // electron charge


DECLARE_string   (simName);

DECLARE_int32    (Nx);
DECLARE_int32    (Ny);
DECLARE_int32    (Nlayers);
DECLARE_double   (cellSize);
DECLARE_double   (timestep);
DECLARE_double   (finaltime);
DECLARE_double   (terminatingTorque);
DECLARE_bool     (adjust_step);

DECLARE_string   (MinitFile);
DECLARE_double   (Ms);
DECLARE_double   (alpha);
DECLARE_double   (gamma);
DECLARE_double   (Aexch);

DECLARE_string   (Bext);
DECLARE_double   (STO_I);
DECLARE_string   (STO_Pdir);
DECLARE_double   (STO_A);
DECLARE_double   (STO_P);
DECLARE_double   (STO_Lambda);
DECLARE_double   (STO_t0);
// DECLARE_int32    (STO_radius, 0, "Radius of STO contact in cells");

DECLARE_bool     (demag);
DECLARE_int32    (subsample_demag);
DECLARE_bool     (exchange);
DECLARE_bool     (external);
DECLARE_bool     (useGPU);
DECLARE_bool     (useFMM);
DECLARE_int32    (fmmP);

DECLARE_bool     (log_Mfield);
DECLARE_int32    (subsample);

DECLARE_int32    (seed);
DECLARE_int32    (verbosity);
DECLARE_bool     (printArgsAndExit);

#endif // #ifndef  _GLOBALS_H_
