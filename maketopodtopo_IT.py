"""
Create topo and dtopo needed for this tohoku test:
Call functions with makeplots==True to create plots of topo, slip, and dtopo.
"""

from __future__ import absolute_import
from __future__ import print_function
from clawpack.geoclaw import dtopotools, topotools
from matplotlib import pyplot as plt
import math, pandas, os, numpy, glob

try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("*** Must first set CLAW environment variable")

try:
    root_dir = os.environ['PTHA']
except:
    raise Exception("*** Must first set PTHA enviornment variable")


def Mw_2_M0(Mw, constant=9.05):
    """ Function that computes the Seismic Moment(M0) from Moment Magnitude(Mw). """
    M0 = 10 ** (Mw * 3 / 2 + constant)
    return M0

def M0_2_Mw(M0, constant=9.05):
    """ Function that computes the Seismic Moment(M0) from Moment Magnitude(Mw). """
    Mw = (2 / 3) * (math.log10(M0) - constant)
    return Mw

def Mw_2_RuptureDim(Mw):
#    Strasser, F. O., Arango, M. C., & Bommer, J. J. (2010).
    area_absigma = [-3.476, 0.952, 0.304]
    width_absigma = [-0.882, 0.351, 0.173]
    length_absigma = [-2.477, 0.585, 0.180]
    area = 10 ** (area_absigma[0] + Mw * area_absigma[1])
    width = 10 ** (width_absigma[0] + Mw * width_absigma[1])
    length = 10 ** (length_absigma[0] + Mw * length_absigma[1])
    output = [area, width, length]
    return output

def slip_from_Mw(Mw, mu=3.5e10, constant=9.05):
    M0 = Mw_2_M0(Mw, constant=constant)
    area = Mw_2_RuptureDim(Mw)[0] * 1e+06  # m^2
    slip = M0 / (area * mu)
    return slip

def make_dtopo(EveID, EveMw, EveLat, EveLon, EveDep, EveRak, EveStr ,EveDip,Length,Area,Slip, makeplots=False):
    
    GRID_SIZ = 0.01215
    X0, X1 = 10, 20
    Y0, Y1 = 30, 44

    XUNITS = round((X1 - X0) / GRID_SIZ)
    YUNITS = round((Y1 - Y0) / GRID_SIZ)

    dtopo_fname = os.path.join(DtopoDir,EveID + '.tt3')
    png_fname = os.path.join(PngDir, EveID + '.png')
    asc_fname = os.path.join(ASCDir, EveID + '.asc')
    prj_fname = os.path.join(ASCDir, EveID + 'prj')

    # For tohoku source defined by [Hayes (USGS, Tohoku 2011) ] , DIP= 10.21 degree, STRIKE = 194 degree and RAKE = 87.52 degree
    #create a subfault
    tohoku_fault = dtopotools.SubFault()
    tohoku_fault.strike = EveStr
    tohoku_fault.rake = EveRak
    tohoku_fault.dip = EveDip
    tohoku_fault.length = Length * 1000 #Mw_2_RuptureDim(EveMw)[2] 
    tohoku_fault.width = Area/Length * 1000 #Mw_2_RuptureDim(EveMw)[1]
    tohoku_fault.depth = EveDep * 1000
    tohoku_fault.slip = Slip #slip_from_Mw(EveMw)
    tohoku_fault.longitude = EveLon
    tohoku_fault.latitude = EveLat
    tohoku_fault.coordinate_specification = "top center"
    # create a fault and add the above single subfault to it
    fault = dtopotools.Fault()
    fault.subfaults = [tohoku_fault]

    x = numpy.linspace(X0, X1, XUNITS)
    y = numpy.linspace(Y0, Y1, YUNITS)
    times = [1.]

    fault.create_dtopography(x, y, times)
    dtopo = fault.dtopo
    dtopo.write(dtopo_fname, dtopo_type=3)
    print("Created ", dtopo_fname)

    tt3_file = open(dtopo_fname, "r")
    lines = tt3_file.readlines()
    tt3_file.close()
    asc_file = open(asc_fname, "w")
    asc_file.write('ncols ' + str.split(lines[0])[0] + '\n')
    asc_file.write('nrows ' + str.split(lines[1])[0] + '\n')
    asc_file.write('xllcorner ' + str.split(lines[3])[0] + '\n')
    asc_file.write('yllcorner ' + str.split(lines[4])[0] + '\n')
    asc_file.write('cellsize ' + str(GRID_SIZ) + '\n')
    asc_file.write('nodata_value -32768' + '\n')

    for line in lines[9:]:
        asc_file.write(line)
    asc_file.close()
    print("Created ", asc_fname)

    prjline = 'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],' \
              'PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]'
    with open(prj_fname, 'w') as f:
        f.write(prjline)

    if makeplots:
        plt.figure(figsize=(12, 7))
        ax1 = plt.subplot(121)
        ax2 = plt.subplot(122)
        fault.plot_subfaults(axes=ax1, slip_color=True)
        ax1.plot(EveLon, EveLat, 'ro')
        ax1.text(142, 32, ('ID:' + str(EveID)), fontsize=11)
        ax1.text(142, 31.2, ('Loc:' + str("{:.3f}".format(EveLat)) + ',' + str("{:.3f}".format(EveLon))), fontsize=11)
        ax1.text(142, 30.4, ('Mw,Depth :' + str(EveMw) + ',' + str(EveDep) + 'km'), fontsize=11)
        ax1.set_xlim(x.min(), x.max())
        ax1.set_ylim(y.min(), y.max())
        dtopo.plot_dZ_colors(1., axes=ax2)
        ax2.plot(EveLon, EveLat, 'ro')
        ax2.text(142, 32, ('Dip:' + str(EveDip) + ' degree'), fontsize=11)
        ax2.text(142, 31.2, ('Rak:' + str(EveRak) + ' degree'), fontsize=11)
        ax2.text(142, 30.4, ('Str:' + str(EveStr) + ' degree'), fontsize=11)
        ax2.set_xlim(x.min(), x.max())
        ax2.set_ylim(y.min(), y.max())
        plt.savefig(png_fname)
        # plt.clf()
        # plt.close('all')
        print("Created ", png_fname)

if __name__ == "__main__":
    # Directory for storing topo and dtopo files:
    DtopoDir = os.path.join(root_dir, 'gis', 'dtopo_IT')
    TopoDir = os.path.join(root_dir, 'gis', 'topo', 'ASC')
    PngDir = os.path.join(DtopoDir, 'dtopopng')
    ASCDir = os.path.join(DtopoDir, 'dtopoasc')
    os.makedirs(DtopoDir, exist_ok=True)
    os.makedirs(PngDir, exist_ok=True)
    os.makedirs(ASCDir, exist_ok=True)

     # Sythetic Events in Tohoku Region based on SLAB.2 data and custom source region with points at .25 degree spacing
    subfaultdf = pandas.read_csv('/mnt/data/nragu/Tsunami/Tohoku/gis/dtopo_IT/slab_dtopo_list.csv')
    for i in range(len(subfaultdf[:-2])):
        print("Creating dtopo file for ", subfaultdf.EveID[i])
        #create dtopo files
        make_dtopo(subfaultdf.EveID[i], subfaultdf.EveMw[i], subfaultdf.lat[i], subfaultdf.lon[i], abs(subfaultdf.dep[i]), subfaultdf.rak[i],
                    subfaultdf.str[i], subfaultdf.dip[i],subfaultdf.length[i],subfaultdf.area[i],subfaultdf.slip[i], makeplots=True)
    print("Completed Creating dtopo files")
