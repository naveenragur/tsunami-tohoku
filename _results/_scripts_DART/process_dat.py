"""
To preprocess dat file from JP cabinet office and gebco
Needs python-gis env in conda for geopandas and osgeo
"""

import os
import pandas as pd
import geopandas as gp
from osgeo import gdal

if __name__ == "__main__":
    print(os.getcwd())
    TsunamiDir = '/home/nrr/Tsunami/Tohoku'
    tileinfofile = os.path.join(TsunamiDir, 'gis', 'topo', 'topotileinfo.csv')
    tileinfo = pd.read_csv(tileinfofile)
    exfile = gp.read_file('gis/topo/Extents.shp')

    for i in range(len(tileinfo)):
        print(tileinfo.tile[i].split("-")[0] + '-' + tileinfo.tile[i].split("-")[1])
        filepath_in = 'gis/topo/depth' + tileinfo.tile[i].split("-")[0] + '/' + tileinfo.tile[i].split("-")[0] \
                      + 'm/depth_' + tileinfo.tile[i] + '.dat'
        raster = 'depth_wgs_' + tileinfo.tile[i] + '.asc'
        filepath_out = 'gis/topo/ASC/depth_' + tileinfo.tile[i] + '.asc'
        filepath_prjout = 'gis/topo/ASC/depth_' + tileinfo.tile[i] + '.prj'
        filepath_wgsout = 'gis/topo/ASC/depth_wgs_' + tileinfo.tile[i] + '.asc'

        NCOLS = tileinfo.nx[i]
        NROWS = tileinfo.ny[i]
        XLLCORNER = tileinfo.xtl[i] - 300
        YLLCORNER = tileinfo.ytl[i] - (tileinfo.ny[i]* tileinfo.dx[i]) + 350
        CELLSIZE = tileinfo.dx[i]
        NODATA_VALUE = -32768
        RES_VALUE = {"1350": 0.01215, "0450": 0.00405, "0150": 0.00135, "0050": 0.00045}
        print(NCOLS,NROWS,XLLCORNER,YLLCORNER,CELLSIZE)

        prjline = 'PROJCS["unknown",GEOGCS["GCS_unknown",DATUM["D_Unknown_based_on_Bessel_1841_ellipsoid",' \
                  'SPHEROID["Bessel_1841",6377397.155,299.1528128]],PRIMEM["Greenwich",0.0],' \
                  'UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],' \
                  'PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",143.0],PARAMETER["Scale_Factor",0.9996],' \
                  'PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]'
        with open(filepath_prjout, 'w') as f:
            f.write(prjline)

        topo = pd.read_fwf(filepath_in, header=None, widths=[8, 8, 8, 8, 8, 8, 8, 8, 8, 8], delimiter="\n\t")

        lines = ['NCOLS ' + str(NCOLS), 'NROWS ' + str(NROWS), 'XLLCORNER ' + str(XLLCORNER),
                 'YLLCORNER ' + str(YLLCORNER),
                 'CELLSIZE ' + str(CELLSIZE), 'NODATA_VALUE ' + str(NODATA_VALUE)]
        with open(filepath_out, 'w') as f:
            for line in lines:
                f.write(line)
                f.write('\n')

        rtopo = pd.DataFrame(topo.values.reshape(NROWS, NCOLS)) * -1
        rtopo.round(decimals=3)
        rtopo.to_csv(filepath_out, sep=' ', header=False, index=None, mode='a')

        # Find the extent to clip with from the shapefile
        onetile = exfile.loc[exfile['Rastername'] == raster]
        corners = onetile.iloc[0].Extent[16:-3].split(',')
        xs = []
        ys = []

        for row in range(4):
            x = corners[row].split()[0]
            y = corners[row].split()[1]
            xs.append(float(x))
            ys.append(float(y))
        gdal.Warp(filepath_wgsout, filepath_out, dstSRS='EPSG:4326', xRes=RES_VALUE[tileinfo.tile[i].split("-")[0]],
                  yRes=RES_VALUE[tileinfo.tile[i].split("-")[0]],
                  outputBounds=[round(min(xs), 3) + .005, round(min(ys), 3) + .005,
                                round(max(xs), 3) - .005, round(max(ys), 3) - .005])

    # Gebco main tile resample resolution
    # filepath_gebco_in = 'gis/topo/gebco/gebco_2021_n50.0_s30.0_w130.0_e150.0.asc'
    # filepath_gebco_out = 'gis/topo/ASC/depth_wgs_gebco.asc'
    # gdal.Warp(filepath_gebco_out, filepath_gebco_in, dstSRS='EPSG:4326', xRes=0.00405, yRes=0.00405)

