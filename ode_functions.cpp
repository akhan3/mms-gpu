#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdlib.h>
#include <sys/time.h>
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


int Hfield (    const Vector3 *M, Vector3 *H, fptype *charge, fptype *potential,
                const int xdim, const int ydim, const int zdim, const fptype meshwidth,
                const fptype mu_0, const fptype Ms, const fptype Aexch,
                const int demag, const int exchange, const int external, const int use_fmm, const int P,
                const int use_gpu, const int verbosity )
{
    int status = 0;
    const int xyzdim = zdim*ydim*xdim;
    if(verbosity) {}
// reset H before beginning
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i = 0; i < xyzdim; i++)
        H[i] = Vector3(0,0,0);

    if(demag) {
        // magnetic volume charge density
        divergence_3d(M, xdim, ydim, zdim, charge);

        // calculate potential from charge
        if(use_fmm)
            status |= fmm_calc(charge, potential, xdim, ydim, zdim, P, use_gpu, verbosity);
        else
            status |= calc_potential_exact(charge, xdim, ydim, zdim, potential, use_gpu, verbosity); // Exact O(N^2) calculation

        // char filename_pot[200];
        // if     (use_gpu == 0) sprintf(filename_pot, "potential_cpu.dat");
        // else if(use_gpu == 1) sprintf(filename_pot, "potential_gpu.dat");
        // else if(use_gpu == 2) sprintf(filename_pot, "potential_gpuemu.dat");
        // else                  sprintf(filename_pot, "potential_gpugpu.dat");
        // status |= save_scalar3d(potential, zdim, ydim, xdim, filename_pot, 100);
        // status |= matrix2file(charge, ydim, xdim, "charge.dat", 100);

        // magnetostatic field from potential = H_demag
        gradient_3d(potential, xdim, ydim, zdim, 1/(4.0*M_PI), H);
    }

    if(exchange) {
        // add exchange field from magnetization = H_exch
        add_exchange_field(M, Ms, Aexch, mu_0, xdim, ydim, zdim, meshwidth, H);
    }

    if(external) {
        // add external field = H_ext
        const Vector3 Hext = 1*Ms * Vector3(0,0,1);
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int i = 0; i < xyzdim; i++)
            H[i] += Hext;
    }

    return status ? EXIT_FAILURE : EXIT_SUCCESS;
}



// update M from current value of M (Runge-Kutta 4th order)
// =================================================================
int rk4_step(   const fptype t, const fptype dt, fptype *t2,
                const byte *material, const Vector3 *M, Vector3 *M2,
                Vector3 *H, Vector3 *slope, Vector3 *this_slope,
                fptype *charge, fptype *potential,
                const int xdim, const int ydim, const int zdim, const fptype meshwidth,
                const fptype mu_0, const fptype Ms, const fptype Aexch, const fptype alfa, const fptype gamma,
                const int demag, const int exchange, const int external, const int use_fmm, const int P,
                const int normalize,
                const int use_gpu,
                const int verbosity)
{
    const int xyzdim = zdim*ydim*xdim;
// reset slope before beginning
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i = 0; i < xyzdim; i++)
        slope[i] = Vector3(0,0,0);
// k1 @ t1
    // Hfield has already been calculated, so save time
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i = 0; i < xyzdim; i++) {
        if(material[i]) {
            this_slope[i] = LLG_Mprime(M[i], H[i], alfa, gamma, Ms);
            M2[i] = M[i] + (this_slope[i] * 0.5 * dt);
            slope[i] += this_slope[i];
        }
    }
// k2 @ t1 + dt/2
    Hfield(M2, H, charge, potential, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, demag, exchange, external, use_fmm, P, use_gpu, verbosity);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i = 0; i < xyzdim; i++) {
        if(material[i]) {
            this_slope[i] = LLG_Mprime(M2[i], H[i], alfa, gamma, Ms);
            M2[i] = M[i] + (this_slope[i] * 0.5 * dt);
            slope[i] += 2 * this_slope[i];
        }
    }
// k3 @ t1 + dt/2
    Hfield(M2, H, charge, potential, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, demag, exchange, external, use_fmm, P, use_gpu, verbosity);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i = 0; i < xyzdim; i++) {
        if(material[i]) {
            this_slope[i] = LLG_Mprime(M2[i], H[i], alfa, gamma, Ms);
            M2[i] = M[i] + (this_slope[i] * dt);
            slope[i] += 2 * this_slope[i];
        }
    }
// k4 @ t1 + dt
    Hfield(M2, H, charge, potential, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, demag, exchange, external, use_fmm, P, use_gpu, verbosity);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i = 0; i < xyzdim; i++) {
        if(material[i]) {
            this_slope[i] = LLG_Mprime(M2[i], H[i], alfa, gamma, Ms);
            slope[i] += this_slope[i];
            M2[i] = M[i] + dt * slope[i] * (1/6.0);
            if(normalize)
                M2[i] = M2[i] * (Ms / M2[i].magnitude());  // re-normalize with Ms
        }
    }
    *t2 = t + dt;
    return EXIT_SUCCESS;
}






// =============================================
// main time marching function
// =============================================
int time_marching(  byte *material, Vector3 *M, // initial state. This will be overwritten in each time step
                    const fptype finaltime, const fptype timestep,
                    const int xdim, const int ydim, const int zdim, const fptype meshwidth, const int P,
                    const fptype mu_0, const fptype Ms, const fptype Aexch, const fptype alfa, const fptype gamma,
                    const int demag, const int exchange, const int external, const int use_fmm,
                    const int use_gpu, const char *sim_name, const int verbosity )
{
    int status = 0;
    const int xyzdim = zdim*ydim*xdim;

// allocate memory
    Vector3 *M2         = new Vector3[xyzdim]();
    Vector3 *H          = new Vector3[xyzdim]();
    Vector3 *H2         = new Vector3[xyzdim]();
    Vector3 *H3         = new Vector3[xyzdim]();
    Vector3 *slope      = new Vector3[xyzdim]();
    Vector3 *this_slope = new Vector3[xyzdim]();
    fptype *charge      = new fptype[xyzdim]();
    fptype *potential   = new fptype[xyzdim]();
    if(M2 == NULL && H == NULL && H2 == NULL && slope == NULL && this_slope == NULL && charge == NULL && potential == NULL) {
        fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }

    char filename[1000];

// open file for logging runtime statistics
    sprintf(filename, "%s/%s", sim_name, "dynamics.dat");
    FILE *fh = fopen(filename, "w");
    if(fh == NULL) {
        printf("FATAL ERROR: Error opening file %s\n", filename);
        return EXIT_FAILURE;
    }
// open file for logging time evolution of M
    sprintf(filename, "%s/%s", sim_name, "Mdynamics.dat");
    FILE *fhM = fopen(filename, "w");
    if(fhM == NULL) {
        printf("FATAL ERROR: Error opening file %s\n", filename);
        return EXIT_FAILURE;
    }
// open file for logging panics
    sprintf(filename, "%s/%s", sim_name, "panic.dat");
    FILE *fhp = fopen(filename, "w");
    if(fhp == NULL) {
        printf("FATAL ERROR: Error opening file %s\n", filename);
        return EXIT_FAILURE;
    }

// adaptive rk4 parameters
    fptype dt = timestep; // start very optimistically at dt_max
    // const fptype dt_min = 0;   // OOMMF has a default of zero
    // const fptype dt_max = 1e-10;
    // const fptype tolerance = 1 * M_PI/180;  // in radians
    // const fptype tolerance_hyst = .5 * tolerance;  // in radians
    // const fptype safety_factor = 0.5; // 1 is no safety at all, while 0 is infinite safety
    const int normalize = true;
    const int adjust_step = true;
// starting point
    int tindex = 0;
    fptype t = 0;
    fptype t2 = t;

    printf("If you want to see the initial state of M, now is the time! \n"); fflush(NULL);
    // getchar();
    printf("\n");



// Time-marching loop
// ======================================
    status |= Hfield(M, H, charge, potential, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, demag, exchange, external, use_fmm, P, use_gpu, verbosity);
    while(t <= finaltime)
    {
    // energies before
        fptype E = 0;
        #ifdef _OPENMP
        // #pragma omp parallel for
        #endif
        for(int i = 0; i < xyzdim; i++) {
            if(material[i]) {
                E += M[i].dot(H[i]);
            }
        }
        E *= -0.5*mu_0 * (meshwidth*meshwidth*meshwidth) / 1.6e-19;

    // RK4 step
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int i = 0; i < xyzdim; i++)
            H3[i] = H[i];
        timeval time1, time2;
        status |= gettimeofday(&time1, NULL);
        status |= rk4_step( t, dt, &t2,
                            material, M, M2,
                            H3, slope, this_slope,
                            charge, potential,
                            xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, alfa, gamma,
                            demag, exchange, external, use_fmm, P,
                            normalize,
                            use_gpu, verbosity );
        if(status) return EXIT_FAILURE;
        status |= gettimeofday(&time2, NULL);
        double deltatime = (time2.tv_sec + time2.tv_usec/1e6) - (time1.tv_sec + time1.tv_usec/1e6);
        if(verbosity >= 5)
            printf("RK4: took %f seconds\n", deltatime);

    // energies and torque after
        status |= Hfield(M2, H2, charge, potential, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, demag, exchange, external, use_fmm, P, use_gpu, verbosity);
        fptype E2 = 0;
        fptype torque = 0;
        #ifdef _OPENMP
        // #pragma omp parallel for
        #endif
        for(int i = 0; i < xyzdim; i++) {
            if(material[i]) {
                E2 += M2[i].dot(H2[i]);
                fptype this_torque =  M2[i].cross(H2[i]).magnitude();
                torque = (torque > this_torque) ? torque : this_torque;
            }
        }
        E2 *= -0.5*mu_0 * (meshwidth*meshwidth*meshwidth) / 1.6e-19;
        torque *= (1/Ms/Ms);
        // printf("dt=%g, E=%f , E2=%f\n", dt, E, E2);

    // check energies and adjust stepsize
        if(E2 > E) { // reject and invalidate this step
            dt = dt / 2.5; // reduce the stepsize
            printf("PANIC!! energy is increasing, so halving the step\n");
            fprintf(fhp, "%d, %g, %g, %g, %g PANIC!! energy is increasing, so halving the step\n",
                        tindex, t, dt, E2, torque);
            // continue;
        }
        else { // accept the step
            if(adjust_step)
                dt = 1e-13 / torque/2;
            t = t2;
            tindex++;
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for(int i = 0; i < xyzdim; i++) {
                M[i] = M2[i];
                H[i] = H2[i];
            }

        // log all the important variables
            if(1)
            // if(E2 <= E) // if step was validated
            {
                // append M to Mdynamics file
                status |= append_vector3d(M, zdim, ydim, xdim, fhM, verbosity);
                // calculate average magnetization
                fptype Mmag_avg = 0;
                int count = 0;
                Vector3 M_avg(0,0,0);
                #ifdef _OPENMP
                // #pragma omp parallel for
                #endif
                for(int i = 0; i < xyzdim; i++) {
                    if(material[i]) {
                        Mmag_avg += M[i].magnitude();
                        M_avg = M_avg + M[i];
                        count++;
                    }
                }
                Mmag_avg /= count;
                M_avg = M_avg * (1.0 / count);

                fprintf(fh, "%d, %g, %g, %g, %g, %g, %g, %g, %g \n",
                        tindex, t, dt, E2, M_avg.x/Ms, M_avg.y/Ms, M_avg.z/Ms, Mmag_avg/Ms, torque);
                fprintf(stdout, "%d, %g, %g, %g, %g, %g, %g, %g, %g \n",
                        tindex, t, dt, E2, M_avg.x/Ms, M_avg.y/Ms, M_avg.z/Ms, Mmag_avg/Ms, torque);
            }

        // check stopping torque criteria
            if(torque < 1e-4) {
                printf("Torque is too low. Breaking simulation...\n");
                break;
            }
        } // valid step


        fflush(NULL);
        // if(tindex >= 10) break;
    } // time marching while loop

// closing
    delete []M2;
    delete []H;
    delete []H2;
    delete []H3;
    delete []slope;
    delete []this_slope;
    delete []charge;
    delete []potential;
    fclose(fh);
    fclose(fhM);
    return status ? EXIT_FAILURE : EXIT_SUCCESS;
}
