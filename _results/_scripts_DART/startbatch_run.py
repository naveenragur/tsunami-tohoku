"""
Module to set up batch of geoclaw runs using model from input folder RUNDIR and dtopo files in folder DTOPODIR.
"""

import os, threading, queue
from runclaw import runclaw
from plotclaw import plotclaw
from process_fg_max import process_fg_max

root_dir = os.environ['PTHA']  # should point to main repo directory
print('root_dir = ', root_dir)


#make fresh inputs if changes are made to setrun like refinement or change in dataset
os.system('make data')

# Default run parameters passed to runclaw and parallel runs
exe = 'xgeoclaw'
rundir = '_input'
outdir = '_output'
plotdir = '_plots'
dtopodir = '_tsunami'
proc = 10

# List of Tsunami Events to run defined by the dtopo file in the DTOPO folder:
dtopo_batch = os.listdir(dtopodir)
# dtopo_batch = ['FUJI2011_46.tt3']
print(dtopo_batch)
# Queue containing input dtopo file for parallel  geoclaw jobs:
q = queue.Queue()
for dtopo in dtopo_batch:
    q.put_nowait(dtopo)

# Runclaw job function - run <proc> times at once
def job1(threadid):
    while not q.empty():
        inpfile = q.get_nowait()
        print("Thread ", threadid, " Starts job to run ", inpfile)
        p1 = runclaw(xclawcmd=exe, outdir=outdir, rundir=rundir, dtopo=inpfile)
        # print("Thread ", threadid, " Starts job to plot ", inpfile)
        # p2 = plotclaw(outdir=outdir, plotdir=plotdir, dtopo=inpfile)
        # print("Thread ", threadid, " Starts job to process fgmax ", inpfile)
        # p3 = process_fg_max(outdir=outdir, plotdir=plotdir, dtopo=inpfile, save_figs=True, testdata=False)

# Main Code:
threads = []
for tid in range(proc):
    t = threading.Thread(target=job1, args=(tid,))
    threads.append(t)
    t.start()

print('\n All jobs pushed!')
