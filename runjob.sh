#!/bin/bash

# set OpenMP threads
export OMP_NUM_THREADS=8

# launch the simulation
#valgrind --leak-check=full \

./main                      \
                dummy_disc_40x40.png   \
                3           \
                            \
                10e-9       \
                1e-13       \
                            \
                250e-9       \
                250e-9       \
                10e-9        \
                            \
                100          \
                100          \
                4           \
                            \
                1           \
                1           \
                0           \
                            \
                0           \
                1           \
                            \
                sim_square_250x10nm_randomIC     \
                2010        \
                2       \

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
# IC = 0; // 0 = SD, 1 = Vortex, 2 = Random
