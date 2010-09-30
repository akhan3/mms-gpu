#!/bin/bash

# set OpenMP threads
export OMP_NUM_THREADS=8

# launch the simulation
#valgrind --leak-check=full \

./main                      \
                disc_40x40.png   \
                3           \
                            \
                10e-9       \
                1e-13       \
                            \
                320e-9       \
                320e-9       \
                30e-9        \
                            \
                32          \
                32          \
                3           \
                            \
                1           \
                1           \
                0           \
                            \
                0           \
                1           \
                            \
                sim_square_32x32_10nm     \
                2010        \

#| tee -a log.dat

exit_status=$?
echo "program exit staus = $exit_status"
exit $exit_status

# default command line args
# =====================================
# imagefile = dummy.png
# P = 3

# finaltime = 1e-9
# timestep = 1e-14

# sample_width  = 10e-9
# sample_height = 10e-9
# sample_depth  = 1e-9

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
