#!/bin/bash

# set OpenMP threads
export OMP_NUM_THREADS=8

# launch the simulation
#valgrind --leak-check=full
./main                      \
                dummy.png   \
                3           \
                10e-9      \
                1e-13       \
                20e-9       \
                32         \
                32          \
                3           \
                1           \
                1           \
                0           \
                0           \
                1           \
                sim_cell20nm_dots_16x16     \
                1985        \
| tee -a log

exit_status=$?
echo "program exit staus = $exit_status"
exit $exit_status

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
# demag = 1
# exchange = 1
# external = 0
# use_fmm = 0
# use_gpu = 0
# sim_name = sim_untitled
# SEED = time(NULL)
