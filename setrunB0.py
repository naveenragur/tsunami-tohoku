"""
Module to set up run time parameters for Clawpack -- AMRClaw code.

The values set in the function setrun are then written out to data files in the rundir location mentioned below
that will be read in by the Fortran code.

"""

import os
import numpy as np
from clawpack.geoclaw import fgmax_tools
#from clawpack.amrclaw.data import FlagRegion

try:
    CLAW = os.environ['CLAW']
    HOME = os.environ['HOME']
except:
    raise Exception("*** Must first set CLAW environment variable")

topodir = os.path.join(HOME, 'Tsunami', 'Tohoku','gis','topo','ASC')
dtopodir = os.path.join(HOME, 'Tsunami', 'Tohoku','gis','dtopo')
if not os.path.isdir(topodir):
    raise Exception("*** Missing topo directory: %s" % topodir)

#initial run without a tsunami forcing to derive bathymetry
making_B0 = True
print("its the B0 set run file here!!!!!!!!!!")
#------------------------------
def setrun(claw_pkg='geoclaw'):
#------------------------------

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData

    """

    from clawpack.clawutil import data
    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)

    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------
    #probdata = rundata.new_UserData(name='probdata',fname='setprob.data')


    #------------------------------------------------------------------
    # GeoClaw specific parameters:
    #------------------------------------------------------------------
    rundata = setgeo(rundata)

    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #------------------------------------------------------------------

    clawdata = rundata.clawdata  # initialized when rundata instantiated


    # Set single grid parameters first.
    # See below for AMR parameters.


    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim

    # Lower and upper edge of computational domain:
    clawdata.lower[0] = 136.0          # xlower 10degrees
    clawdata.upper[0] = 146.0          # xupper
    clawdata.lower[1] = 32.0          # ylower 8degrees
    clawdata.upper[1] = 46.0          # yupper

    # Number of grid cells: Res of raw data {"1350": 0.01215, "0450" : 0.00405, "0150" : 0.00135, "0050": 0.00045}
    # region0:dx & dy = ~ 0.06075 ~6.68km
    #region1:dx & dy = ~ 0.06075 / 5 =.01215 ~ 1.336km
    clawdata.num_cells[0] = 164 # mx
    clawdata.num_cells[1] = 132 # my

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 3

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 3

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 2

    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0


    # Restart from checkpoint file of a previous run?
    # Note: If restarting, you must also change the Makefile to set:
    #    RESTART = True
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False               # True to restart from prior results
    clawdata.restart_file = 'fort.chk00006'  # File to use for restart data


    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.

    clawdata.output_style = 1

    if clawdata.output_style == 1:
        # Output ntimes frames at equally spaced times up to tfinal:
        # Can specify num_output_times = 0 for no output
        if making_B0:
            clawdata.num_output_times = 1  # for making B0
            clawdata.tfinal = 1.  # for making B0
            clawdata.output_t0 = True  # for making B0
        else:
            clawdata.num_output_times = 4
            clawdata.tfinal = 1*3600.
            clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list or numpy array of output times:
        # Include t0 if you want output at the initial time.
        clawdata.output_times =  np.linspace(7,13,4)*3600.

    elif clawdata.output_style == 3:
        # Output every step_interval timesteps over total_steps timesteps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 1
        clawdata.output_t0 = False  # output at initial (or restart) time?

    clawdata.output_format == 'binary'      # 'ascii', 'binary', 'netcdf'

    clawdata.output_q_components = 'all'   # could be list such as [True,True]
    clawdata.output_aux_components = 'none'  # could be list
    clawdata.output_aux_onlyonce = True    # output aux arrays only at t0


    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 0



    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==True:  variable time steps used based on cfl_desired,
    # if dt_variable==Falseixed time steps dt = dt_initial always used.
    clawdata.dt_variable = True

    # Initial time step for variable dt.
    # (If dt_variable==0 then dt=dt_initial for all steps)
    clawdata.dt_initial = 0.016

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99

    # Desired Courant number if variable dt used
    clawdata.cfl_desired = 0.75
    # max Courant number to allow without retaking step with a smaller dt:
    clawdata.cfl_max = 1.0

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 50000


    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2

    # Use dimensional splitting? (not yet available for AMR)
    clawdata.dimensional_split = 'unsplit'

    # For unsplit method, transverse_waves can be
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 2


    # Number of waves in the Riemann solution:
    clawdata.num_waves = 3

    # List of limiters to use for each wave family:
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'vanleer'  ==> van Leer
    #   4 or 'mc'       ==> MC limiter
    clawdata.limiter = ['vanleer', 'vanleer', 'vanleer']

    clawdata.use_fwaves = True    # True ==> use f-wave version of algorithms

    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used,
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 1


    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2

    # Choice of BCs at xlower and xupper:
    #   0 or 'user'     => user specified (must modify bcNamr.f to use this option)
    #   1 or 'extrap'   => extrapolation (non-reflecting outflow)
    #   2 or 'periodic' => periodic (must specify this at both boundaries)
    #   3 or 'wall'     => solid wall for systems where q(2) is normal velocity

    clawdata.bc_lower[0] = 'extrap'   # at xlower
    clawdata.bc_upper[0] = 'extrap'   # at xupper

    clawdata.bc_lower[1] = 'extrap'   # at ylower
    clawdata.bc_upper[1] = 'extrap'   # at yupper


    # ---------------
    # Gauges:
    # ---------------
    gauges = rundata.gaugedata.gauges

    # for gauges append lines of the form [gaugeno, x, y, t1, t2]
    gauges.append([1001, 141.45, 38.08, 0, 1.e9])
    gauges.append([1002, 141.40, 37.90, 0, 1.e9])
    gauges.append([301, 141.255, 38.298, 0, 1.e9])
    gauges.append([302, 141.011, 38.005, 0, 1.e9])
    gauges.append([51, 141.2244, 38.3912, 0, 1.e9])
    gauges.append([52, 140.9216, 38.0026, 0, 1.e9])

    # --------------
    # Checkpointing:
    # --------------

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = 0

    if clawdata.checkpt_style == 0:
      # Do not checkpoint at all
      pass

    elif clawdata.checkpt_style == 1:
      # Checkpoint only at tfinal.
      pass

    elif clawdata.checkpt_style == 2:
      # Specify a list of checkpoint times.
      clawdata.checkpt_times = [0.1,0.15]

    elif clawdata.checkpt_style == 3:
      # Checkpoint every checkpt_interval timesteps (on Level 1)
      # and at the final time.
      clawdata.checkpt_interval = 5



    # ---------------
    # AMR parameters:   (written to amr.data)
    # ---------------
    amrdata = rundata.amrdata

    # max number of refinement levels:
    amrdata.amr_levels_max = 4   # Set to 5 originally

    # List of refinement ratios at each level (length at least amr_level_max-1)
    #Res of raw data {"1350": 0.01215, "0450" : 0.00405, "0150" : 0.00135, "0050": 0.00045}
    # gebco = 0.06075/5 level 1
    # 1350  = 0.01215/3 level 2
    # 0450  = 0.00405/3 level 3
    # 0150  = 0.00135/3 level 4
    # 0050  = 0.00045   level 5

    amrdata.refinement_ratios_x = [5,3,3]#,3]
    amrdata.refinement_ratios_y = [5,3,3]#,3]
    amrdata.refinement_ratios_t = [5,3,3]#,3]

    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length num_aux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).
    amrdata.aux_type = ['center', 'capacity', 'yleft']

    # Flag for refinement based on Richardson error estimater:
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag_richardson_tol = 1.0  # Richardson tolerance

    # Flag for refinement using routine flag2refine:
    amrdata.flag2refine = True      # use this?
    amrdata.flag2refine_tol = 0.5  # tolerance used in this routine
    # Note: in geoclaw the refinement tolerance is set as wave_tolerance below
    # and flag2refine_tol is unused!

    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 3

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width  = 2

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    amrdata.clustering_cutoff = 0.7

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 0


    # ---------------
    # Regions:
    # ---------------
    regions = rundata.regiondata.regions
    # to specify regions of refinement append lines of the form
    # [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]

    inf = 1e9

    # Region 0 : Global region.  This assures a maximum refinement in any
    # regions not covered by other regions listed below.(gebco).
    regions.append([1, 1, 0., inf, (clawdata.lower[0] - 0.1), (clawdata.upper[0] + 0.1),
                (clawdata.lower[1] - 0.1), (clawdata.upper[1] + 0.1)])

    # Region 1 : (from the dtopo file, below)(gebco).
    # Time interval taken from dtopo file : [0,1]
    regions.append([2, 2, 0, 1, 136., 146., 32., 45.])

    # Region 2 : Region encompassing Sendai coast(150m-09) .
    # Time interval :  (0,18000)
    regions.append([4, 4, 0., inf, 140.84, 141.10, 37.75, 38.34])

    # Region 3 : Region encompassing N.Sendai coast(50m-21).
    # Time interval :  (0,18000)
    #regions.append([5, 5, 0., inf, 140.902, 141.486, 38.123, 38.451])

    # Region 4 : Region encompassing S.Sendai coast(50m-20) .
    # Time interval :  (0,18000)
    #regions.append([5, 5, 0., inf, 140.845, 141.152, 37.737, 38.115])

    # Region 5 : Region encompassing S.Sendai coast(150m) for fgmax readings.
    # Time interval :  (0,18000)
    regions.append([4, 4, 0., inf, 140.82, 141.03, 37.725, 38.325])

    #  ----- For developers -----
    # Toggle debugging print statements:
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = False      # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting

    return rundata

    # end of function setrun
    # ----------------------


#-------------------
def setgeo(rundata):
#-------------------
    """
    Set GeoClaw specific runtime parameters.
    """

    try:
        geo_data = rundata.geo_data
    except:
        print("*** Error, this rundata has no geo_data attribute")
        raise AttributeError("Missing geo_data attribute")

    # == Physics ==
    geo_data.gravity = 9.81
    geo_data.coordinate_system =  2
    geo_data.earth_radius = 6367500.0

    # == Forcing Options
    geo_data.coriolis_forcing = False

    # == Algorithm and Initial Conditions ==
    geo_data.sea_level = 0.0
    geo_data.dry_tolerance = 0.001
    geo_data.friction_forcing = True
    geo_data.manning_coefficient = 0.035
    geo_data.friction_depth = 500.0

    # Refinement settings
    refinement_data = rundata.refinement_data
    refinement_data.variable_dt_refinement_ratios = True
    refinement_data.wave_tolerance = 0.1
    refinement_data.deep_depth = 200.0
    refinement_data.max_level_deep = 3

    # == settopo.data values ==
    #[topotype, fname]
    topofiles = rundata.topo_data.topofiles
    topofiles.append([3, os.path.join(topodir, 'depth_wgs_gebco.asc')])
    topofiles.append([3, os.path.join(topodir, 'depth_wgs_0150-09.asc')])
    # topofiles.append([3, os.path.join(topodir, 'depth_wgs_0050-21.asc')])
    # topofiles.append([3, os.path.join(topodir, 'depth_wgs_0050-22.asc')])

    # == setdtopo.data values ==
    # [dtype, minlevel, maxlevel, 'dtopofname']
    # for moving topography, append lines of the form :   (<= 1 allowed for now!)

    if making_B0:
        rundata.dtopo_data.dtopofiles = []  # no dtopo for making B0
        rundata.dtopo_data.dt_max_dtopo = 1.0
    else:
        # Region 1, above handles this region.
        rundata.dtopo_data.dtopofiles = [[3, os.path.join(dtopodir, 'tohoku_dtopo_2011.tt3')]]

    # == setqinit.data values ==
    rundata.qinit_data.qinit_type =  0
    rundata.qinit_data.qinitfiles = []

    # == fixedgrids.data values ==
    # fixedgrids = []
    # rundata.fixed_grid_data.fixedgrids = []
    # fixedgrids = rundata.fixed_grid_data.fixedgrids

    # == fgmax_grids.data values ==
    # NEW STYLE STARTING IN v5.7.0

    # set num_fgmax_val = 1 to save only max depth,
    #                     2 to also save max speed,
    #                     5 to also save max hs,hss,hmin

    if making_B0:
        rundata.fgmax_data.num_fgmax_val = 1  # for making B0 only depth

    else:
        rundata.fgmax_data.num_fgmax_val = 2  # Save depth and speed

    fgmax_grids = rundata.fgmax_data.fgmax_grids  # empty list to start

    # Now append to this list objects of class fgmax_tools.FGmaxGrid
    # specifying any fgmax grids.

    tstart_max = 0.  # for making max depth
    dx150 = 0.00135  # 9/3600.
    ###### Sendai Bay
    # Points on a uniform 2d grid:
    fg = fgmax_tools.FGmaxGrid()
    fg.point_style = 2  # uniform rectangular x-y grid
    fg.x1 = 140.85 + dx150 / 2.0
    fg.x2 = 141.09 - dx150 / 2.0
    fg.y1 = 37.76 + dx150 / 2.0
    fg.y2 = 38.33 - dx150 / 2.0
    fg.dx = dx150  # 150m level ~4.86 sec
    fg.min_level_check = 4
    fg.tstart_max = tstart_max
    fg.tend_max = 1.e10  # when to stop monitoring max values
    fg.dt_check = 20.  # how often to update max values
    fg.interp_method = 0  # 0 ==> pw const in cells, recommended
    fgmax_grids.append(fg)  # written to fgmax_grids.data

     # Now append to this list objects of class fgmax_tools.FGmaxGrid
    # specifying any fgmax grids.
    if making_inundation:
        tstart_max = 0.  # for making max depth
        dx150 = 0.00135
        dx050 = 0.00045

        ###### Sendai Bay
        # Points on a uniform 2d grid:
        fg = fgmax_tools.FGmaxGrid()
        fg.point_style = 2  # uniform rectangular x-y grid
        fg.x1 = 140.85 + dx050 / 2.0
        fg.x2 = 141.09 - dx050 / 2.0
        fg.y1 = 37.76 + dx050 / 2.0
        fg.y2 = 38.33 - dx050 / 2.0
        fg.dx = dx050  # 50m level
        fg.min_level_check = 4
        fg.tstart_max = tstart_max
        fg.tend_max = 1.e10  # when to stop monitoring max values
        fg.dt_check = 30.  # how often to update max values
        fg.interp_method = 0  # 0 ==> pw const in cells, recommended
        fgmax_grids.append(fg)  # written to fgmax_grids.data

        # # ###### Ishinomaki Bay
        # Points on a uniform 2d grid:
        fg = fgmax_tools.FGmaxGrid()
        fg.point_style = 2  # uniform rectangular x-y grid
        fg.x1 = 141.14 + dx050 / 2.0
        fg.x2 = 141.39 - dx050 / 2.0
        fg.y1 = 38.35 + dx050 / 2.0
        fg.y2 = 38.49 - dx050 / 2.0
        fg.dx = dx050  # 50m level
        fg.min_level_check = 4
        fg.tstart_max = tstart_max
        fg.tend_max = 1.e10  # when to stop monitoring max values
        fg.dt_check = 30.  # how often to update max values
        fg.interp_method = 0  # 0 ==> pw const in cells, recommended
        fgmax_grids.append(fg)  # written to fgmax_grids.data

        # ####### Rikuzentakata Bay
        # # Points on a uniform 2d grid:
        fg = fgmax_tools.FGmaxGrid()
        fg.point_style = 2  # uniform rectangular x-y grid
        fg.x1 = 141.59 + dx050 / 2.0
        fg.x2 = 141.69 - dx050 / 2.0
        fg.y1 = 38.96 + dx050 / 2.0
        fg.y2 = 39.04 - dx050 / 2.0
        fg.dx = dx050  # 50m level
        fg.min_level_check = 4
        fg.tstart_max = tstart_max
        fg.tend_max = 1.e10  # when to stop monitoring max values
        fg.dt_check = 30.  # how often to update max values
        fg.interp_method = 0  # 0 ==> pw const in cells, recommended
        fgmax_grids.append(fg)  # written to fgmax_grids.data

    return rundata
    # end of function setgeo
    # ----------------------

if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    from clawpack.geoclaw import kmltools

    if making_B0:
        rundir='_input/B0'
    else:
        rundir='_input'

    os.makedirs(rundir, exist_ok = True)
    os.chdir(rundir)

    rundata = setrun(*sys.argv[1:])
    rundata.write()

    kmltools.make_input_data_kmls(rundata)
