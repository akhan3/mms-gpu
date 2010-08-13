#!/bin/bash

# clean-up
rm -rf log*
rm -rf *.dat
make clean
# fresh compile
make

# launch the simulation
time ./main                 \
                dummy.png   \
                3           \
                5e-10       \
                2e-14       \
                1e-9        \
                1           \
                1           \
                0           \
                0           \
                sim_coupling_exchange_cell1nm_nofmm     \
                1985        \
    &> log
    echo "program exit staus = $?" >> log


# default command line args
# =====================================
# imagefile = dummy.png
# P = 3
# finaltime = 1e-9
# timestep = 1e-14
# meshwidth = 1e-9
# coupling = 0
# exchange = 1
# external = 0
# use_fmm = 0
# sim_name = sim_untitled
# SEED = time(NULL)
