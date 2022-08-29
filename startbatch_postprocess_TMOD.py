"""
Module to set up batch of geoclaw runs using model from input folder RUNDIR and dtopo files in folder DTOPODIR.
"""

import os, threading, queue,glob
from runclaw import runclaw
from plotclaw import plotclaw
from process_fg_max import process_fg_max

root_dir = os.environ['PTHA']  # should point to main repo directory
print('root_dir = ', root_dir)

#make fresh inputs if changes are made to setrun like refinement or change in dataset
#os.system('make data')

# Default run parameters passed to runclaw and parallel runs
exe = 'xgeoclaw'
rundir = '_input'
outdir = '_results/_output_DART'
plotdir = '_results/_plots_DART2'
dtopodir = 'gis/dtopo_TMOD'
proc = 1

# List of Tsunami Events to run defined by the dtopo file in the DTOPO folder:
# dtopo_batch = glob.glob(dtopodir + '/*.tt3')
# dtopo_batch = [os.path.basename(dtopo_batch[i]) for i in range(len(dtopo_batch))]
# dtopo_done = [(item + '.tt3') for item in os.listdir(plotdir)]
# dtopo_batch = [item for item in dtopo_batch if item not in dtopo_done]

#dtopo_batch = ['Gusman.txydz', 'Ammon.txydz', 'GUSMAN2012.tt3', 'LAY2011.tt3', 'YUE2013.tt3', 'UCSB3.txydz',
               # 'GCMT.txydz', 'YAMAZAKI2011.tt3', 'YAMAZAKI2018_TMOD_Static.tt3', 'YAMAZAKI2018_TMOD.tt3',
               # 'Hayes.txydz', 'Fujii.txydz', 'HAYES2011.tt3', 'YAMAZAKI2018_TPMOD_Static.tt3',
               # 'YAMAZAKI2018_TPMOD.tt3', 'SATAKE2013.tt3', 'Caltech.txydz', 'GusmanAU.txydz', 'HAYES2017.tt3',
               # 'SATAKE2013_Static.tt3']

dtopo_batch = ['FUJI2011_42']

# Queue containing input dtopo file for parallel  geoclaw jobs:
q = queue.Queue()
for dtopo in dtopo_batch:
    q.put_nowait(dtopo)

# Runclaw job function - run <proc> times at once
def job1(threadid):
    while not q.empty():
        inpfile = q.get_nowait()
        print("Thread ", threadid, " Starts job to plot ", inpfile)
        p2 = plotclaw(outdir=outdir, plotdir=plotdir, dtopo=inpfile)
        # print("Thread ", threadid, " Starts job to process fgmax ", inpfile)
        # p3 = process_fg_max(outdir=outdir, plotdir=plotdir, dtopo=inpfile, save_figs=True, testdata=False,dtopodir=dtopodir)

# Main Code:
threads = []
for tid in range(proc):
    t = threading.Thread(target=job1, args=(tid,))
    threads.append(t)
    t.start()

print('\n All jobs pushed!')
