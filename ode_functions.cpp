#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdlib.h>
#include "ode_functions.hpp"
#include "helper_functions.hpp"
#include "vector_functions.hpp"


Vector3 LLG_Mprime( const Vector3 &M,
                    const Vector3 &H,
                    const fptype alfa, const fptype gamma, const fptype Ms )
{
    Vector3 McrossH = M.cross(H);
    Vector3 Mprime = -gamma * McrossH - (alfa * gamma / Ms) * M.cross(McrossH);
    return Mprime;
}


int Hfield(const Vector3 *M, Vector3 *H,
                const int xdim, const int ydim, const int zdim, const fptype meshwidth,
                const fptype mu_0, const fptype Ms, const fptype Aexch,
                const int coupling, const int exchange, const int external, const int use_fmm, const int P,
                const int verbose_level )
{
    int status = 0;
    if(verbose_level) {}
// reset H before beginning
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i = 0; i < xdim*ydim*zdim; i++)
        H[i] = Vector3(0,0,0);

    if(coupling) {
    // allocate memory for charge and potential field
        fptype *charge = new fptype[zdim*ydim*xdim]();
        fptype *potential = new fptype[zdim*ydim*xdim]();
        if(charge == NULL || potential == NULL) {
            fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
            return EXIT_FAILURE;
        }
    // magnetic volume charge density
        divergence_3d(M, xdim, ydim, zdim, charge);
    // calculate potential
        // calc_H_nearest_neighbor(M, H, xdim, ydim, zdim); // nearest neighbor coupling only
        if(use_fmm)
            status |= fmm_calc(charge, potential, xdim, ydim, zdim, P, verbose_level);
        else
            calc_potential_exact(charge, xdim, ydim, zdim, potential); // Exact O(N^2) calculation
    // magnetostatic field from potential = H_demag + H_coupl
        gradient_3d(potential, xdim, ydim, zdim, 1/(4.0*M_PI), H);
    // clean-up
        delete []charge;
        delete []potential;
    }

    if(exchange) {
    // Exchange field from magnetization = H_exch
        add_exchange_field(M, Ms, Aexch, mu_0, xdim, ydim, zdim, meshwidth, H);
    }

    if(external) {
    // add external field = H_ext
        const Vector3 Hext = 1*Ms * Vector3(0,0,1);
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int i = 0; i < zdim*ydim*xdim; i++)
            H[i] += Hext;
    }

    return status ? EXIT_FAILURE : EXIT_SUCCESS;
}


int postprocess_M(  const Vector3 *M,
                    const int tindex, const fptype t, const fptype dt,
                    const fptype energy, fptype *energy_new,
                    const int xdim, const int ydim, const int zdim, const fptype meshwidth,
                    const fptype mu_0, const fptype Ms, const fptype Aexch,
                    const int coupling, const int exchange, const int external, const int use_fmm, const int P,
                    FILE *fh, FILE *fhM,
                    fptype *torque_max,
                    int verbose_level )
{
    int status = 0;
// append M to Mdynamics file
    status |= append_vector3d(M, zdim, ydim, xdim, tindex, t, fhM, verbose_level);
// update H for current M
    Vector3 *H = new Vector3[zdim*ydim*xdim]();
    if(H == NULL) {
        fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
    status |= Hfield(M, H, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, coupling, exchange, external, use_fmm, P, verbose_level);
// calculate energy and average magnetization
    *torque_max = 0;
    *energy_new = 0;
    fptype Mmag_avg = 0;
    int count = 0;
    Vector3 M_avg(0,0,0);
    #ifdef _OPENMP
    // #pragma omp parallel for
    #endif
    for(int i = 0; i < zdim*ydim*xdim; i++) {
        if(M[i].magnitude()) {
            fptype torque =  M[i].cross(H[i]).magnitude() * (1/Ms/Ms);
            *torque_max = (*torque_max > torque) ? *torque_max : torque;
            *energy_new += M[i].dot(H[i]);
            Mmag_avg += M[i].magnitude();
            M_avg = M_avg + M[i];
            count++;
        }
    }
    *energy_new *= -0.5*mu_0 * (meshwidth*meshwidth*meshwidth) / 1.6e-19;
    fptype dE = *energy_new - energy;
    Mmag_avg /= count;
    M_avg = M_avg * (1.0 / count);
    fprintf(fh, "%d, %g, %g, %g, %g, %g, %g, %g, %g, %g \n",
            tindex, t, dt, energy, dE, M_avg.x/Ms, M_avg.y/Ms, M_avg.z/Ms, Mmag_avg/Ms, *torque_max);
    fprintf(stdout, "%d, %g, %g, %g, %g, %g, %g, %g, %g, %g \n",
            tindex, t, dt, energy, dE, M_avg.x/Ms, M_avg.y/Ms, M_avg.z/Ms, Mmag_avg/Ms, *torque_max);
    if(dE > 0 && tindex != 0) {
        status |= save_vector3d(M, zdim, ydim, xdim, "M.dat", verbose_level);
        if(status) return EXIT_FAILURE;
        printf("PANIC!! dE = %g eV,  dE/E = %g \n", dE, dE/energy);
        // status = EXIT_FAILURE;
        // break;
    }
    delete []H;
    if(!(tindex % 100))
    {
        status |= save_vector3d(M, zdim, ydim, xdim, "M.dat", verbose_level);
        if(status) return EXIT_FAILURE;
    }
    return status;
}

// update M from current value of M (Runge-Kutta 4th order)
// =================================================================
int rk4_step(   const fptype t, const Vector3 *M, const fptype dt,
                const int xdim, const int ydim, const int zdim, const fptype meshwidth,
                const fptype mu_0, const fptype Ms, const fptype Aexch, const fptype alfa, const fptype gamma,
                const int coupling, const int exchange, const int external, const int use_fmm, const int P,
                const int normalize,
                fptype *t_new, Vector3 *M_new,    // output
                const int verbose_level)
{
    const int xyzdim = zdim*ydim*xdim;
    Vector3 *H = new Vector3[xyzdim]();
    Vector3 *slope = new Vector3[xyzdim]();
    Vector3 *this_slope = new Vector3[xyzdim]();
    if(H == NULL || slope == NULL || this_slope == NULL) {
        fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
// k1 @ t1
    Hfield(M, H, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, coupling, exchange, external, use_fmm, P, verbose_level);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i = 0; i < xyzdim; i++) {
        if(M[i].magnitude()) {
            this_slope[i] = LLG_Mprime(M[i], H[i], alfa, gamma, Ms);
            M_new[i] = M[i] + (this_slope[i] * 0.5 * dt);
            slope[i] += this_slope[i];
        }
    }
// k2 @ t1 + dt/2
    Hfield(M_new, H, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, coupling, exchange, external, use_fmm, P, verbose_level);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i = 0; i < xyzdim; i++) {
        if(M[i].magnitude()) {
            this_slope[i] = LLG_Mprime(M_new[i], H[i], alfa, gamma, Ms);
            M_new[i] = M[i] + (this_slope[i] * 0.5 * dt);
            slope[i] += 2 * this_slope[i];
        }
    }
// k3 @ t1 + dt/2
    Hfield(M_new, H, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, coupling, exchange, external, use_fmm, P, verbose_level);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i = 0; i < xyzdim; i++) {
        if(M[i].magnitude()) {
            this_slope[i] = LLG_Mprime(M_new[i], H[i], alfa, gamma, Ms);
            M_new[i] = M[i] + (this_slope[i] * dt);
            slope[i] += 2 * this_slope[i];
        }
    }
// k4 @ t1 + dt
    Hfield(M_new, H, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, coupling, exchange, external, use_fmm, P, verbose_level);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i = 0; i < xyzdim; i++) {
        if(M[i].magnitude()) {
            this_slope[i] = LLG_Mprime(M_new[i], H[i], alfa, gamma, Ms);
            slope[i] += this_slope[i];
            M_new[i] = M[i] + dt * slope[i] * (1/6.0);
            if(normalize)
                M_new[i] = M_new[i] * (Ms / M_new[i].magnitude());  // re-normalize with Ms
        }
    }
    *t_new = t + dt;
    delete []H;
    delete []slope;
    delete []this_slope;
    return EXIT_SUCCESS;
}




int rk4_step_adaptive(   const fptype t, const Vector3 *M, const fptype dt,
                const fptype dt_min, const fptype dt_max, const fptype tolerance, const fptype tolerance_hyst, const fptype safety_factor,
                const int xdim, const int ydim, const int zdim, const fptype meshwidth,
                const fptype mu_0, const fptype Ms, const fptype Aexch, const fptype alfa, const fptype gamma,
                const int coupling, const int exchange, const int external, const int use_fmm, const int P,
                const int normalize,
                const int adjust_step,
                fptype *t_new, Vector3 *M_new, fptype *dt_new,    // output
                const int verbose_level)
{
    int status = 0;
    const int xyzdim = zdim*ydim*xdim;
    fptype energy = 0;
    Vector3 *H = new Vector3[xyzdim]();
    status |= Hfield(M, H, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, coupling, exchange, external, use_fmm, P, verbose_level);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i = 0; i < zdim*ydim*xdim; i++)
        if(M[i].magnitude())
            energy += M[i].dot(H[i]);
    energy *= -0.5*mu_0 * (meshwidth*meshwidth*meshwidth) / 1.6e-19;

// provisional RK step
    // printf("rk4_step_adaptive: %g, %g, %g, %g\n", t, dt, *t_new, *dt_new);
    status |= rk4_step( t, M, dt,
                        xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, alfa, gamma,
                        coupling, exchange, external, use_fmm, P,
                        normalize,
                        t_new, M_new,   // output from RK step
                        verbose_level);
    if(!adjust_step) {
        *dt_new = dt;
    }
    else {
    // calculate energy and maximum torque
        fptype torque_max = 0;
        fptype energy_new = 0;
        status |= Hfield(M_new, H, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, coupling, exchange, external, use_fmm, P, verbose_level);
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int i = 0; i < zdim*ydim*xdim; i++) {
            if(M_new[i].magnitude()) {
                fptype torque =  M_new[i].cross(H[i]).magnitude() * (1/Ms/Ms);
                torque_max = (torque_max > torque) ? torque_max : torque;
                energy_new += M_new[i].dot(H[i]);
            }
        }
        energy_new *= -0.5*mu_0 * (meshwidth*meshwidth*meshwidth) / 1.6e-19;
        if(energy_new > energy) {
            *t_new = t;  // reject and invalidate this step
            *dt_new = dt / 2.5;
            printf("PANIC!! energy is increasing, so halving the step\n");
        }
        else {
            // if(torque_max * dt < 1e-13)
                *dt_new = 1e-13 / torque_max/2;
                // *dt_new = dt * 1.1;
            // else
                // *dt_new = 1e-13 / torque_max/2;
            // *dt_new = (*dt_new < dt_max) ? *dt_new : dt_max;
        }
        if(torque_max < 1e-6) {
            printf("Simulation is practically finished!\n");
            return EXIT_FAILURE;
        }
    }
    return status;
}


// =============================================
// main time marching function
// =============================================
int time_marching(  Vector3 *M, // initial state. This will be overwritten in each time step
                    const fptype finaltime, const fptype timestep,
                    const int xdim, const int ydim, const int zdim, const fptype meshwidth, const int P,
                    const fptype mu_0, const fptype Ms, const fptype Aexch, const fptype alfa, const fptype gamma,
                    const int coupling, const int exchange, const int external, const int use_fmm,
                    const int verbose_level )
{
    int status = 0;
// open file for logging runtime statistics
    char filename[] = "dynamics.dat";
    FILE *fh = fopen(filename, "w");
    if(fh == NULL) {
        printf("FATAL ERROR: Error opening file %s\n", filename);
        return EXIT_FAILURE;
    }
// open file for logging time evolution of M
    FILE *fhM = fopen("Mdynamics.dat", "w");
    if(fhM == NULL) {
        printf("FATAL ERROR: Error opening file %s\n", filename);
        return EXIT_FAILURE;
    }
// adaptive rk4 parameters
    // const fptype dt_min = 1e-17;   // OOMMF has a default of zero
    const fptype dt_min = 0;   // OOMMF has a default of zero
    const fptype dt_max = 1e-10;
    fptype dt = timestep; // start very optimistically at dt_max
    // const fptype tolerance = Ms/1e3;
    const fptype tolerance = 1 * M_PI/180;  // in radians
    const fptype tolerance_hyst = .5 * tolerance;  // in radians
    const fptype safety_factor = 0.5; // 1 is no safety at all, while 0 is infinite safety
    const int normalize = true;
    const int adjust_step = true;
// H field parameters
    // const int coupling = true;
    // const int exchange = true;
    // const int external = false;
// starting point
    int tindex = 0;
    fptype t = 0;
    fptype energy = 0;
    fptype energy_new = 0;
    fptype torque_max = 0;
// post-process initial state
    status |= postprocess_M( M, tindex, t, dt,
                             energy, &energy_new,
                             xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch,
                             coupling, exchange, external, use_fmm, P,
                             fh, fhM, &torque_max, verbose_level );
    printf("If you want to see the initial state of M, now is the time! \n"); fflush(NULL);
    // getchar();
    printf("\n");

// Time-marching loop
// ======================================
    while(t < finaltime)
    {
        fptype t_new = -99;
        fptype dt_new = -99;
        Vector3 *M_new = new Vector3[zdim*ydim*xdim]();
        if(M_new == NULL) {
            fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
            return EXIT_FAILURE;
        }
    // execute one step of RK
        // printf("\nexecuting RK step...\n");
        status |= rk4_step_adaptive( t, M, dt,
                            dt_min, dt_max, tolerance, tolerance_hyst, safety_factor,
                            xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, alfa, gamma,
                            coupling, exchange, external, use_fmm, P,
                            normalize,
                            adjust_step,
                            &t_new, M_new, &dt_new,  // output from RK step
                            verbose_level);
        if(status) return EXIT_FAILURE;
        dt = dt_new; // update the step-size anyway
        if(t_new > t) { // if(step_valid) do additional book-keeping
            tindex++;
            t = t_new;
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for(int i = 0; i < zdim*ydim*xdim; i++)
                M[i] = M_new[i];
            fptype energy_new = 0;
            // if(!(tindex % 50))
            status |= postprocess_M( M, tindex, t, dt,
                                     energy, &energy_new,
                                     xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch,
                                     coupling, exchange, external, use_fmm, P,
                                     fh, fhM, &torque_max, verbose_level );
            energy = energy_new;
        } // if(step_valid)
        delete []M_new;
        fflush(NULL);
        // if(tindex >= 5) break;
    } // time marching while loop

    status |= postprocess_M( M, tindex, t, dt,
                             energy, &energy_new,
                             xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch,
                             coupling, exchange, external, use_fmm, P,
                             fh, fhM, &torque_max, verbose_level );

// closing
    fclose(fh);
    fclose(fhM);
    return status ? EXIT_FAILURE : EXIT_SUCCESS;
}
