#!/bin/bash

# set OpenMP threads
export OMP_NUM_THREADS=8

SIM_NAME="sim_Test"
FIGURE="shapes/disc_50x50.png"

#mkdir $SIM_NAME
#cp $0 $SIM_NAME/runjob.sh.$SIM_NAME
#cp $FIGURE $SIM_NAME

# launch the simulation
#valgrind --leak-check=full \

src/main    --maskTxt shapes/disc_13x13.txt         \
            --P 3                       \
            --finaltime 0.5e-12         \
            --timestep 2.5e-13          \
            --cellSize 5e-9       \
            --Nlayers 1        \
            --demag                     \
            --exchange                  \
            --external                  \
            --nouse_fmm                 \
            --use_gpu                   \
            --sim_name $SIM_NAME        \
            --seed 2010                 \
            --IC 0                      \

#| tee -a log.dat

exit_status=$?
echo "program exit staus = $exit_status"
exit $exit_status
