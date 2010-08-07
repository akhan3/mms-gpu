// #include <stdio.h>
// #include <assert.h>
// #include <cmath>
#include <stdlib.h>
#include "ode_functions.hpp"
#include "helper_functions.hpp"
#include "vector_functions.hpp"


Vector3 LLG_Mprime( const Vector3 &M,
                    const Vector3 &H,
                    const float alfa, const float gamma, const float Ms )
{
    Vector3 McrossH = M.cross(H);
    Vector3 Mprime = -gamma * McrossH - (alfa * gamma / Ms) * M.cross(McrossH);
    return Mprime;
}



int calc_Hfield(const Vector3 *M, Vector3 *H,
                const int xdim, const int ydim, const int zdim, const float meshwidth,
                const float mu_0, const float Ms, const float Aexch,
                const int verbose_level )
{
    int status = 0;

// reset H before beginning
    for(int i = 0; i < xdim*ydim*zdim; i++)
        H[i] = Vector3(0,0,0);

// nearest neighbor coupling only
    // calc_H_nearest_neighbor(M, H, xdim, ydim, zdim);

// allocate memory for charge and potential field
    float *charge = new float[zdim*ydim*xdim]();
    float *potential = new float[zdim*ydim*xdim]();
    if(charge == NULL || potential == NULL) {
        fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }

// // magnetic volume charge density
    divergence_3d(M, xdim, ydim, zdim, charge);

// // calculate potential in each layer through FMM
    status |= fmm_calc(charge, potential, xdim, ydim, zdim, 3, verbose_level);
    // calc_potential_exact(charge, xdim, ydim, zdim, potential); // Exact O(N^2) calculation

// magnetostatic field from potential = H_demag + H_coupl
    // gradient_3d(potential, xdim, ydim, zdim, 1/(4*M_PI), H);

// Exchange field from magnetization = H_exch
    add_exchange_field(M, Ms, Aexch, mu_0, xdim, ydim, zdim, meshwidth, H);

// // add external field = H_ext
    // const Vector3 Hext = 0.1*Ms * Vector3(0,1,0);
    // for(unsigned int i = 0; i < zdim*ydim*xdim; i++)
        // H[i] += Hext;

    delete []charge;
    delete []potential;
    return status ? EXIT_FAILURE : EXIT_SUCCESS;
}

// update M from current value of M (Runge-Kutta 4th order)
// =================================================================
int rk4_step(   Vector3 *M, // input and output state
                const int normalize, const float dt,
                const int xdim, const int ydim, const int zdim, const float meshwidth,
                const float mu_0, const float Ms, const float Aexch, const float alfa, const float gamma,
                const int verbose_level )
{
    const int xyzdim = zdim*ydim*xdim;
    Vector3 *M_ = new Vector3[xyzdim]();
    Vector3 *H_ = new Vector3[xyzdim]();
    Vector3 *average_slope = new Vector3[xyzdim]();
    Vector3 *this_slope = new Vector3[xyzdim]();
    if(M_ == NULL || H_ == NULL || average_slope == NULL || this_slope == NULL) {
        fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
// k1 @ t1
    calc_Hfield(M, H_, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, verbose_level);
    for(int i = 0; i < xyzdim; i++) {
        this_slope[i] = LLG_Mprime(M[i], H_[i], alfa, gamma, Ms);
        M_[i] = M[i] + (this_slope[i] * 0.5 * dt);
        average_slope[i] += this_slope[i];
    }
// k2 @ t1 + dt/2
    calc_Hfield(M_, H_, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, verbose_level);
    for(int i = 0; i < xyzdim; i++) {
        this_slope[i] = LLG_Mprime(M_[i], H_[i], alfa, gamma, Ms);
        M_[i] = M[i] + (this_slope[i] * 0.5 * dt);
        average_slope[i] += 2 * this_slope[i];
    }
// k3 @ t1 + dt/2
    calc_Hfield(M_, H_, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, verbose_level);
    for(int i = 0; i < xyzdim; i++) {
        this_slope[i] = LLG_Mprime(M_[i], H_[i], alfa, gamma, Ms);
        M_[i] = M[i] + (this_slope[i] * dt);
        average_slope[i] += 2 * this_slope[i];
    }
// k4 @ t1 + dt
    calc_Hfield(M_, H_, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, verbose_level);
    for(int i = 0; i < xyzdim; i++) {
        this_slope[i] = LLG_Mprime(M_[i], H_[i], alfa, gamma, Ms);
        average_slope[i] += this_slope[i];
        M[i] = M[i] + dt * average_slope[i] * (1/6.0);
        if(normalize)
            M[i] = M[i] * (Ms / M[i].magnitude());  // re-normalize with Ms
    }
    return EXIT_SUCCESS;
}

// main time marching loop
// =============================================
int time_marching(  Vector3 *M, // initial state. This will be overwritten in each time step
                    const int tdim, const float dt, const float finaltime,
                    const int xdim, const int ydim, const int zdim, const float meshwidth,
                    const float mu_0, const float Ms, const float Aexch, const float alfa, const float gamma,
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

// allocate memory for H
    Vector3 *H = new Vector3[zdim*ydim*xdim]();
    if(H == NULL) {
        fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }

// Time-marching loop
// ======================================
    for(int t = 0; t < tdim; t++)
    {
        printf("**t = %d, time = %g\n", t, t*dt);

    // update H for current M
        status |= calc_Hfield(M, H, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, verbose_level);

    // // write H vectorfield to file
        // status |= save_vector3d(H, zdim, ydim, xdim, "H", verbose_level);
        // if(status) return EXIT_FAILURE;

    // calculate energy and average magnetization
        float energy = 0;
        float Ms_avg = 0;
        int count = 0;
        Vector3 M_avg;
        for(int i = 0; i < zdim*ydim*xdim; i++) {
            if(M[i].magnitude()) {
                energy += M[i].dot(H[i]);
                Ms_avg += M[i].magnitude();
                M_avg = M_avg + M[i];
                count++;
            }
        }
        energy *= -0.5*mu_0 * (meshwidth*meshwidth*meshwidth);
        Ms_avg /= count;
        M_avg = M_avg * (1.0 / count);
        fprintf(fh, "%d, %g, %g, %g, %g, %g, %g\n", t, t*dt, energy, M_avg.x, M_avg.y, M_avg.z, Ms_avg);

    // update M from current value of M by executing one step of RK
        if(t == tdim-1) break;
        printf("  executing RK step...\n");
        status |= rk4_step(M, false, dt, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, alfa, gamma, verbose_level);

        if(!(t % 500)) {
            printf("step completed! ");
            // write M vectorfield to file
            status |= save_vector3d(M, zdim, ydim, xdim, "M", verbose_level);
            if(status) return EXIT_FAILURE;
            fflush(NULL);
            getchar();
        }

    } // time loop

// freeing up memory
    fclose(fh);
    delete []H;
    return status ? EXIT_FAILURE : EXIT_SUCCESS;
}
