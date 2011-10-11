#!/bin/bash

# set OpenMP threads
export OMP_NUM_THREADS=8


SIM_NAME="sim_spin_150x150x3_disc_10GHz_constant_demag"
FIGURE="disc_150x150.png"

mkdir $SIM_NAME
cp $0 $SIM_NAME/$SIM_NAME.runjob.sh
cp $FIGURE $SIM_NAME

# launch the simulation
#valgrind --leak-check=full \

./main                      \
                $FIGURE   \
                3           \
                            \
                0.5e-9       \
                2.5e-13       \
                            \
                750e-9       \
                750e-9       \
                15e-9        \
                            \
                150          \
                150          \
                3           \
                            \
                0           \
                1           \
                1           \
                            \
                0           \
                1           \
                            \
                $SIM_NAME     \
                2010        \
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
