#!/bin/csh

# set OpenMP threads
setenv OMP_NUM_THREADS=4

# launch the simulation
./main                      \
                dummy.png   \
                3           \
                5e-10       \
                5e-14       \
                1e-9        \
                32          \
                32          \
                3           \
                1           \
                1           \
                0           \
                0           \
                sim_coupling_exchange_cel1nm_nofmm     \
                1985        \
    | tee log

echo "program exit staus = $?" >> log


# default command line args
# =====================================
# imagefile = dummy.png
# P = 3
# finaltime = 1e-9
# timestep = 1e-14
# meshwidth = 1e-9
# xdim = 16
# ydim = 16
# zdim = 3
# coupling = 0
# exchange = 1
# external = 0
# use_fmm = 0
# sim_name = sim_untitled
# SEED = time(NULL)
