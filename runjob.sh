#!/bin/bash

# set OpenMP threads
export OMP_NUM_THREADS=16


SIM_NAME="results/sim_spin_50x50x3_testing"
FIGURE="shapes/disc.png"

mkdir -p $SIM_NAME
cp $0 $SIM_NAME/$SIM_NAME.runjob.sh
cp $FIGURE $SIM_NAME

# launch the simulation
#valgrind --leak-check=full \

src/main                      \
                $FIGURE   \
                3           \
                            \
                10e-9       \
                2e-13       \
                            \
                65e-9       \
                65e-9       \
                15e-9        \
                            \
                13          \
                13          \
                3           \
                            \
                1           \
                1           \
                0           \
                            \
                0           \
                0           \
                            \
                $SIM_NAME     \
                2011        \
                0       \

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
