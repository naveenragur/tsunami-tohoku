16/11/2021
Folders:
gis - contains all the raw and formatted datasets - topo dtopo fault historical observations etc
_input and _input/B0 - contains geoclaw setup files for simulation created forced with the displacement files dtopo ( B0 needs no dtopo ie without and coseismic dispalcement
_output - contains simulation run files like tohoku, synthetic events
_plots/runname/fgmax... - contains postprocessed plots for each simulation and gis outputs for hazard
_tsunami - input event forcing for each tsunami - tohoku etc
others: test (initial run for tohoku), scripts(initial backup)

Scripts:
Makefile - contains settings to use make and links with makefile.common in clawpack utils
process_dat.py - used to format raw elevation files from JP tsunami project to usable ascii files for geoclaw
maketopo.py - used to create topo plots and dtopo file for tohoku and synthetic events
setrun.py and setrunB0 - used to setup geoclaw simulaiton - will create the input folders
make_B0.py - to process the bathy from B0 simulation and save in inputB0 folder as npy file
setplot.py - to postprocess simulations into usable plots
startbatch.py - to start parallised batch of events to simulate 
process_fg_max.py - to postprocess a simulation result ie fgmax file to usable point and raster file(geotiff)
xgeoclaw - compiled fortran exe for running geoclaw

Others import scripts modified to work with batch setup etc
file:///home/nrr/opt/clawpack/clawpack-v5.8.0/clawutil/src/Makefile.common
file:///home/nrr/opt/clawpack/clawpack-v5.8.0/clawutil/src/python/clawutil/runclaw.py
file:///home/nrr/opt/clawpack/clawpack-v5.8.0/visclaw/src/python/visclaw/plotclaw.py



 
 
