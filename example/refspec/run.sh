#! /bin/bash

# Setup...
jurassic=jurassic/src

# Loop over windows...
for nu in $(seq 650 100 2350) ; do
    
    # Modify control file...
    cp template.ctl limb_$nu.ctl
    echo $nu | awk '{
      for(i=0; i<100; i++)
        print "NU["i"] = "$1+i
    }' >> limb_$nu.ctl
    
    # Create atmospheric data file...
    $jurassic/climatology limb_$nu.ctl atm.tab
    
    # Create observation geometry...
    $jurassic/limb limb_$nu.ctl obs.tab Z0 3 Z1 68 DZ 1.0
    
    # Call forward model...
    $jurassic/formod limb_$nu.ctl obs.tab atm.tab rad_$nu.tab TASK contrib

    # Convert spectra...
    for f in $(ls rad_$nu*) ; do
	$jurassic/obs2spec limb_$nu.ctl $f spec.$f
	rm $f
    done
    
done
