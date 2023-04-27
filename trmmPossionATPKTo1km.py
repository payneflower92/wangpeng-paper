import calendar
import os.path
import re, pathlib
import sys

import numpy
import pyinterpolate
'''

reference Page：https://pyinterpolate.readthedocs.io/en/latest/usage/tutorials/Poisson%20Kriging%20-%20Area%20to%20Point%20%28Advanced%29.html

'''
try:
    from osgeo import gdal
    from osgeo import ogr
except Exception:
    import gdal
    import ogr
import pandas, geopandas
from geopandas import points_from_xy
from pyinterpolate import build_experimental_variogram, Blocks, PointSupport, Deconvolution, TheoreticalVariogram, \
    smooth_blocks
from parrallelSmoothBlocks import *

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'


def atpkProcess(srcFile, tifShpFile, newpath, monthHours):
    newRoot = os.path.dirname(newpath)
    newCurrentYearMonthTempDir = os.path.join(newRoot, yearmonth + "_tempDir", "areaIdPickles")
    pathlib.Path(newCurrentYearMonthTempDir).mkdir(parents=True, exist_ok=True)

    modelFileName = re.sub(r'TRMM1KmAfterQHshpMaskAndATPK.tif',
                           'TRMM1KmAfterQHshpMaskAndATPK.regularized_model.json', os.path.basename(newpath))
    regularized_model_json = os.path.join(newRoot, modelFileName)

    dr = gdal.Open(srcFile)
    transform = dr.GetGeoTransform()
    trmm_x_origin = transform[0]
    trmm_y_origin = transform[3]
    pixel_width = transform[1]
    pixel_height = transform[5]
    prj = dr.GetProjection()
    tiffDriver = gdal.GetDriverByName("Gtiff")
    srcArr = None

    try:
        srcArr = dr.ReadAsArray()
    except Exception:
        pass

    arr = []
    for j in range(dr.RasterYSize):
        for i in range(dr.RasterXSize):
            lat = trmm_y_origin + j * pixel_height
            lon = trmm_x_origin + i * pixel_width
            if srcArr[j, i] >= 0:
                arr.append([lat, lon, srcArr[j, i] * monthHours])

    srcDf = pandas.DataFrame(arr, columns=["latitude", "longtitude", "TrmmValue"])
    tdf = geopandas.GeoDataFrame(srcDf.values[:, 0:3], columns=srcDf.columns[0:3], )
    tdf["geometry"] = points_from_xy(tdf.longtitude, tdf.latitude, crs="epsg:4326")
    tdf.geometry = tdf["geometry"]
    # 获取block dataframe
    gdf = geopandas.GeoDataFrame.from_file(tifShpFile, encoding='gb18030')
    gdf.columns = ["TrmmValue", "geometry"]
    filterGdf = gdf[gdf["TrmmValue"] >= 0]
    filterGdf["TrmmValue"] = filterGdf["TrmmValue"] * monthHours

    blocks = pyinterpolate.Blocks()
    blocks.from_geodataframe(filterGdf, "TrmmValue", geometry_col='geometry', index_col=None)
    point_support = pyinterpolate.PointSupport()
    point_support.from_geodataframes(tdf, blocks.data, "geometry", "TrmmValue", "geometry", "index",
                                     use_point_support_crs=True)

    if not os.path.isfile(regularized_model_json):
        # experimental_variogram = build_experimental_variogram(input_array=dem, step_size=step_radius,
        #                                                       max_range=max_range)
        # semivariogram_model = TheoreticalVariogram()
        # fitted = semivariogram_model.autofit(
        #     experimental_variogram=experimental_variogram,
        #     model_types='all',
        #     nugget=0,
        #     rang=var_range,
        #     sill=sill)

        maximum_range = 6
        step_size = 0.25
        dt = blocks.data[[blocks.cx, blocks.cy, blocks.value_column_name]]  # x, y, val
        exp_semivar = build_experimental_variogram(input_array=dt, step_size=step_size, max_range=maximum_range)

        pt = point_support.point_support[[point_support.x_col, point_support.y_col, point_support.value_column]].values
        exp_semivarPoint = build_experimental_variogram(input_array=pt, step_size=step_size, max_range=maximum_range)

        reg_mod = Deconvolution(verbose=True)
        reg_mod.fit(agg_dataset=blocks,
                    point_support_dataset=point_support,
                    agg_step_size=step_size,
                    agg_max_range=maximum_range,
                    variogram_weighting_method='closest',
                    model_types='basic')

        reg_mod.transform(max_iters=5)
        # exp_semivar.plot(plot_semivariance=True, plot_covariance=True, plot_variance=True)
        # exp_semivarPoint.plot(plot_semivariance=True, plot_covariance=True, plot_variance=True)
        # reg_mod.plot_deviations()
        # reg_mod.plot_weights()
        # reg_mod.plot_variograms()

        reg_mod.export_model(regularized_model_json)

    smoothedFile = re.sub(r'.TRMM1KmAfterQHshpMaskAndATPK.tif', '.TRMM1KmAfterQHshpMaskAndATPK.smoothedPointMatrix.csv',
                          newpath)

    if not os.path.isfile(smoothedFile):

        semivariogram = TheoreticalVariogram()
        semivariogram.from_json(regularized_model_json)
        smoothed = smooth_area_to_point_pk(semivariogram_model=semivariogram,
                                           blocks=blocks,
                                           point_support=point_support,
                                           number_of_neighbors=8,
                                           max_range=6,
                                           crs=blocks.data.crs,
                                           raise_when_negative_error=False,
                                           raise_when_negative_prediction=True)



        smoothed.to_csv(smoothedFile,sep="\t")

    # if not os.path.isfile(newpath):

    # base = blocks.data.plot(figsize=(14, 14), color='white', edgecolor='black')
    # smooth_plot_data = smoothed.copy()
    # smooth_plot_data[smooth_plot_data['pred'] == 0] = numpy.NAN
    # smooth_plot_data.plot(ax=base, column='pred', cmap='plasma', legend=True, markersize=30, alpha=0.7)


if __name__ == '__main__':

    dataRoot = r"C:\Users\crawl\Desktop\workForWP\wangpeng-master-thesis\reorganizedRawData"
    if sys.platform == "linux":
        dataRoot = "/mnt/sata3_storage/jinlinfang/wpThesis"
    outDir = os.path.join(dataRoot, "TRMMtifATPK")
    tifSrcDir = os.path.join(dataRoot, "TRMMtifMaskedByQH")
    tifShpDir = os.path.join(dataRoot, "TRMMtif2ShpPolygon")


    # processAreaScripts = r'C:\Users\crawl\Desktop\workForWP\wangpeng-master-thesis\jinlinfangScript\smooth_area_to_point_pk_per_areaID.py'
    pythonBin=r"C:\Users\crawl\anaconda3\envs\geoProcessing\python.exe"
    if sys.platform == "linux":
        # processAreaScripts = os.path.join(dataRoot,
        #                                   "wangpengPaperCode/jinlinfangScript/smooth_area_to_point_pk_per_areaID.py")
        pythonBin="/mnt/sata3_storage/pipe_app/anaconda3/envs/geoProcessing310/bin/python"
    srcFileDic=dict()
    for d in os.listdir(tifSrcDir):
        if not re.search(r'^\d+$', d):
            continue
        for f in os.listdir(os.path.join(tifSrcDir, d)):
            if not re.search(r'HDF.tif.clip.tif$', f):
                continue
            srcFile = os.path.join(tifSrcDir, d, f)
            yearmonth = re.search(r'(\d+).HDF.tif.clip.tif', f).group(1)
            srcFileDic[int(yearmonth)] = (srcFile,d,f)

    sortedSrcFileDic = dict(sorted(srcFileDic.items()))
    for yearmonth in sortedSrcFileDic.keys():
        (srcFile, d, f) = sortedSrcFileDic[yearmonth]
        yearmonth = str(yearmonth)
        year = int(str(yearmonth)[0:4])
        month = int(str(yearmonth)[4:])


        newName = re.sub(r'.HDF.tif.clip.tif', '.TRMM1KmAfterQHshpMaskAndATPK.tif', f)
        newSubDir = os.path.join(outDir, d)
        newpath = os.path.join(newSubDir, newName)
        pathlib.Path(newSubDir).mkdir(parents=True, exist_ok=True)


        shpName = yearmonth + "TRMMPolygonShp.shp"
        tifShpFile = os.path.join(tifShpDir, d, shpName)
        allHours = 24 * calendar.monthrange(year, month)[1]

        print(re.sub(r'\\', "/", "srcFile=\"" + srcFile + "\""))
        print(re.sub(r'\\', "/", "tifShpFile=\"" + tifShpFile + "\""))
        print(re.sub(r'\\', "/", "newpath=\"" + newpath + "\""))
        print(re.sub(r'\\', "/", "allHours=" + str(allHours)))

        atpkProcess(srcFile, tifShpFile, newpath, allHours)
