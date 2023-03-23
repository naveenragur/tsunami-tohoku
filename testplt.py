import matplotlib.pyplot as plt
import numpy as np
from clawpack.geoclaw import dtopotools, topotools
import os
try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("*** Must first set CLAW environment variable")

try:
    root_dir = os.environ['PTHA']
except:
    raise Exception("*** Must first set PTHA enviornment variable")

# xpoints = np.array([1, 8])
# ypoints = np.array([3, 10])

# plt.plot(xpoints, ypoints)
# plt.show()
# plt.savefig("dummy_name.png")

subfault_fname = root_dir + '/gis/dtopo_sift/dtopofil/SL_0641.csv'
input_units = {"length":"km", "width":"km", "depth":"km", "slip":"m"}
fault = dtopotools.CSVFault()
fault.read(subfault_fname, input_units=input_units, coordinate_specification="bottom center")

print("The seismic moment is %g N-m" % fault.Mo())
print("The Moment magnitude is %g" % fault.Mw())
print("  (Assuming the rigidity mu of all subfaults is the default value %g Pa)"\
      % fault.subfaults[0].mu)
shorelines_file = root_dir + '/gis/coastlines/JPshoreline_xy.npy'
shore = np.load(shorelines_file)

fault.plot_subfaults(plot_rake=True)
fault.plot_subfaults(slip_color=True)
plt.plot(shore[:, 0], shore[:, 1], 'g')
plt.show() #to print plot on screen

