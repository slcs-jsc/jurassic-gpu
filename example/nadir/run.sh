#! /bin/bash

# Setup...
jurassic=../../src
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../../lib/build/lib:../../lib/build/lib64

# Create atmospheric data file...
OMP_NUM_THREADS=1 $jurassic/climatology nadir.ctl atm.tab

# Create observation geometry...
OMP_NUM_THREADS=1 $jurassic/nadir nadir.ctl obs.tab T1 10

rm -f rad.tab
# Call forward model...
# reference=~/Codes/JURASSIC/reference-jurassic/jurassic/src
# OMP_NUM_THREADS=1 $reference/formod nadir.ctl obs.tab atm.tab rad.tab TASK time
OMP_NUM_THREADS=1 $jurassic/formod nadir.ctl obs.tab atm.tab rad.tab TASK time CHECKMODE $1

# Plot results... brightness temperature vs. tangent point latitude
gnuplot <<EOF
set term png enh truecolor font "Helvetica,28" size 1600,1200 crop lw 2
set out "plot.png"
set xla "latitude [deg]"
set yla "brightness temperature [K]"
set mxtics
set mytics
set key spac 1.5
set key outside
plot "rad.org" u 10:11 w lp pt 1 t  "ref (667.8 cm^{-1})", \
     "rad.tab" u 10:11 w lp pt 2 t "test (667.8 cm^{-1})", \
     "rad.org" u 10:12 w lp pt 1 t  "ref (668.5 cm^{-1})", \
     "rad.tab" u 10:12 w lp pt 2 t "test (668.5 cm^{-1})", \
     "rad.org" u 10:13 w lp pt 1 t  "ref (669.8 cm^{-1})", \
     "rad.tab" u 10:13 w lp pt 2 t "test (669.8 cm^{-1})"
EOF

# Plot results... transmittance vs. tangent point latitude
gnuplot <<EOF
set term png enh truecolor font "Helvetica,28" size 1600,1200 crop lw 2
set out "plot2.png"
set xla "latitude [deg]"
set yla "transmittance [%]"
set mxtics
set mytics
set key spac 1.5
set key outside
plot "rad.org" u 10:(\$14*100) w lp pt 1 t  "ref (667.8 cm^{-1})", \
     "rad.tab" u 10:(\$14*100) w lp pt 2 t "test (667.8 cm^{-1})", \
     "rad.org" u 10:(\$15*100) w lp pt 1 t  "ref (668.5 cm^{-1})", \
     "rad.tab" u 10:(\$15*100) w lp pt 2 t "test (668.5 cm^{-1})", \
     "rad.org" u 10:(\$16*100) w lp pt 1 t  "ref (669.8 cm^{-1})", \
     "rad.tab" u 10:(\$16*100) w lp pt 2 t "test (669.8 cm^{-1})"
EOF

# Get differences...
echo -e "\nCheck for differences..."
diff -sq rad.tab rad.org
