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

// magnetic volume charge density
    divergence_3d(M, xdim, ydim, zdim, charge);

// calculate potential in each layer through FMM
    // status |= fmm_calc(charge, potential, xdim, ydim, zdim, 3, verbose_level);
    // calc_potential_exact(charge, xdim, ydim, zdim, potential); // Exact O(N^2) calculation

// magnetostatic field from potential = H_demag + H_coupl
    // gradient_3d(potential, xdim, ydim, zdim, 1/(4*M_PI), H);

// Exchange field from magnetization = H_exch
    add_exchange_field(M, Ms, Aexch, mu_0, xdim, ydim, zdim, meshwidth, H);

// // add external field = H_ext
    // const Vector3 Hext = 0.1*Ms * Vector3(0,0,1);
    // for(int i = 0; i < zdim*ydim*xdim; i++)
        // H[i] += Hext;

    delete []charge;
    delete []potential;
    return status ? EXIT_FAILURE : EXIT_SUCCESS;
}


// adjust step size if necesssary
float adjust_step(  const Vector3 *slope,
                    const float dt, const int xyzdim, const float Ms,
                    const float control_point, const float control_point_hysteresis,
                    const float step_correction )
{
    Vector3 Mmax;
    for(int i = 0; i < xyzdim; i++) {
        Vector3 dM = slope[i] * dt * (1/Ms);
        Mmax.x = (fabs(dM.x) > Mmax.x) ? fabs(dM.x) : Mmax.x;
        Mmax.y = (fabs(dM.y) > Mmax.y) ? fabs(dM.y) : Mmax.y;
        Mmax.z = (fabs(dM.z) > Mmax.z) ? fabs(dM.z) : Mmax.z;
    }
    float maxchange =  (Mmax.x > Mmax.y) ? Mmax.x : ((Mmax.y > Mmax.z) ? Mmax.y : Mmax.z);
    const float llimit = control_point - control_point_hysteresis/2;
    const float ulimit = control_point + control_point_hysteresis/2;
    float dt_new = dt;
    if(maxchange < llimit) { // increase step size
        dt_new += dt * step_correction;
        printf("%g < %g ", maxchange, llimit);
        printf("step too small, increasing... ");
        printf("%g, %g\n", dt, dt_new);
        dt_new = adjust_step(slope, dt_new, xyzdim, Ms, control_point, control_point_hysteresis, step_correction);
    }
    else if(maxchange < ulimit)
        return dt_new;
    if(maxchange > ulimit) { // decrease step size
        dt_new -= dt * step_correction;
        printf("%g > %g ", maxchange, ulimit);
        printf("step too large, decreasing... ");
        printf("%g, %g\n", dt, dt_new);
        dt_new = adjust_step(slope, dt_new, xyzdim, Ms, control_point, control_point_hysteresis, step_correction);
    }
    else if(maxchange > llimit)
        return dt_new;
    return dt_new;
}


// update M from current value of M (Runge-Kutta 4th order)
// =================================================================
int rk4_step(   Vector3 *M, // input and output state
                float *dt, const int normalize,
                const int xdim, const int ydim, const int zdim, const float meshwidth,
                const float mu_0, const float Ms, const float Aexch, const float alfa, const float gamma,
                const int verbose_level )
{
    const int xyzdim = zdim*ydim*xdim;
    Vector3 *M_ = new Vector3[xyzdim]();
    Vector3 *H_ = new Vector3[xyzdim]();
    Vector3 *slope = new Vector3[xyzdim]();
    Vector3 *this_slope = new Vector3[xyzdim]();
    if(M_ == NULL || H_ == NULL || slope == NULL || this_slope == NULL) {
        fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
// initial slope @ t1 with current step size
    calc_Hfield(M, H_, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, verbose_level);
    for(int i = 0; i < xyzdim; i++) {
        this_slope[i] = LLG_Mprime(M[i], H_[i], alfa, gamma, Ms);
    }
// adjust step size
    *dt = adjust_step(this_slope, *dt, xyzdim, Ms, 0.01, 0.001, 0.5);
// k1 @ t1
    for(int i = 0; i < xyzdim; i++) {
        M_[i] = M[i] + (this_slope[i] * 0.5 * *dt);
        slope[i] += this_slope[i];
    }
// k2 @ t1 + dt/2
    calc_Hfield(M_, H_, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, verbose_level);
    for(int i = 0; i < xyzdim; i++) {
        this_slope[i] = LLG_Mprime(M_[i], H_[i], alfa, gamma, Ms);
        M_[i] = M[i] + (this_slope[i] * 0.5 * *dt);
        slope[i] += 2 * this_slope[i];
    }
// k3 @ t1 + dt/2
    calc_Hfield(M_, H_, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, verbose_level);
    for(int i = 0; i < xyzdim; i++) {
        this_slope[i] = LLG_Mprime(M_[i], H_[i], alfa, gamma, Ms);
        M_[i] = M[i] + (this_slope[i] * *dt);
        slope[i] += 2 * this_slope[i];
    }
// k4 @ t1 + dt
    calc_Hfield(M_, H_, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, verbose_level);
    for(int i = 0; i < xyzdim; i++) {
        this_slope[i] = LLG_Mprime(M_[i], H_[i], alfa, gamma, Ms);
        slope[i] += this_slope[i];
        M[i] = M[i] + *dt * slope[i] * (1/6.0);
        if(normalize)
            M[i] = M[i] * (Ms / M[i].magnitude());  // re-normalize with Ms
    }
    return EXIT_SUCCESS;
}

// main time marching function
// =============================================
int time_marching(  Vector3 *M, // initial state. This will be overwritten in each time step
                    float dt, const float finaltime,
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
    int tindex = 0;
    float time = 0;
    float energy_prev = 0;
    while(time <= finaltime)
    {
        printf("**tindex = %d, time = %g, dt = %g\n", tindex, time, dt);

    // update H for current M
        status |= calc_Hfield(M, H, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, verbose_level);

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
        fprintf(fh, "%d, %g, %g, %g, %g, %g, %g, %g, %g\n",
                tindex, time, dt, energy, energy-energy_prev, M_avg.x, M_avg.y, M_avg.z, Ms_avg);
        energy_prev = energy;

    // update M from current value of M by executing one step of RK
        // printf("executing RK step...\n");
        status |= rk4_step(M, &dt, false, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, alfa, gamma, verbose_level);
        if(status) return EXIT_FAILURE;

        if(!(tindex % 500))
        {
            // write M vectorfield to file
            status |= save_vector3d(M, zdim, ydim, xdim, "M", verbose_level);
            if(status) return EXIT_FAILURE;
            printf("step completed! ");
            fflush(NULL);
            // getchar();
        }

    // increment time
        tindex++;
        time += dt;
    } // time loop

// freeing up memory
    fclose(fh);
    delete []H;
    return status ? EXIT_FAILURE : EXIT_SUCCESS;
}
