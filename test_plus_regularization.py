import calendar
import multiprocessing
import os.path
import re, pathlib
import subprocess
import sys
import socket

import pyinterpolate
import pyproj
import shapely
from shapely.examples.geoms import polygon, Point

from jinlinfangScript.plotFunction import *

try:
    from osgeo import gdal
    from osgeo import ogr
except Exception:
    import gdal
    import ogr
import pandas, geopandas
from geopandas import points_from_xy
import numpy as np
from pyinterpolate import Blocks, VariogramCloud, TheoreticalVariogram, kriging, Deconvolution, PointSupport
import os

global workerNo, processNo, taskPattern
workerNo = 1
processNo = 1
taskPattern = "2\d{3}"

'''
reference Page：https://pyinterpolate.readthedocs.io/en/latest/usage/tutorials/Blocks%20to%20points%20Ordinary%20Kriging%20interpolation%20%28Intermediate%29.html
'''
os.environ['OPENBLAS_NUM_THREADS'] = '1'


def convertProj(p):
    print("转换{},")
    wgs84 = pyproj.CRS('EPSG:4326')
    utm = pyproj.CRS('EPSG:3395')
    project = pyproj.Transformer.from_crs(wgs84, utm, always_xy=True).transform
    # for e in polygon
    utm_p = shapely.ops.transform(project, p)
    return utm_p


def atpkOKProcessParallel(srcFile, tifShpFile, newpath, monthHours):
    print(tifShpFile)
    print(srcFile)
    newRoot = os.path.dirname(newpath)
    yearmonth = re.search(r'(\d+).HDF.tif.clip.tif$', srcFile).group(1)
    smoothedFile = re.sub(r'.TRMM1KmAfterQHshpMaskAndATPK.tif',
                          '.TRMM1KmAfterQHshpMaskAndATPK.smoothedPointMatrixOKatpk.csv',
                          newpath)

    cmd = "wc -l " + smoothedFile
    returnVal = 0
    if os.path.isfile(smoothedFile):
        try:
            returnVal = subprocess.check_output(cmd, shell=True).decode()
            returnVal = str(returnVal).strip().split(r" ")
            returnVal = int(returnVal[0])
        except Exception:
            print("its platform is not linux ")

    condition = (not os.path.isfile(smoothedFile)) or (not (returnVal != 0 and returnVal >= 1272601))

    if not condition: return "don't need to trans"

    modelFilePrefix = re.sub(r'.TRMM1KmAfterQHshpMaskAndATPK.tif', '.TRMM1KmAfterQHshpMaskAndATPK.model.',
                             newpath)
    smoothedShp = re.sub(r'.TRMM1KmAfterQHshpMaskAndATPK.tif',
                         '.TRMM1KmAfterQHshpMaskAndATPK.smoothedPointMatrixOKatpk.shp',
                         newpath)
    # 打开面tif数据
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
    tdf.geometry = tdf["geometry"]  # 1116*3的数据框 0.25 degree的点数据

    # 获取block dataframe Read areal data
    gdf = geopandas.GeoDataFrame.from_file(tifShpFile, encoding='gb18030')
    gdf.columns = ["TrmmValue", "geometry"]
    filterGdf = gdf[gdf["TrmmValue"] >= 0]
    filterGdf["TrmmValue"] = filterGdf["TrmmValue"] * monthHours
    filterGdf["geometry_lonlat"] = filterGdf["geometry"]
    filterGdf = filterGdf.to_crs(crs='EPSG:3395')
    filterGdf["polygonIndex"] = filterGdf.index

    blocks = pyinterpolate.Blocks()
    blocks.from_geodataframe(filterGdf, "TrmmValue", geometry_col='geometry', index_col=None)
    areal_centroids = blocks.data[['centroid_x', 'centroid_y', 'TrmmValue']].values
    '''
    
    from pyinterpolate import Blocks, read_txt, calc_point_to_point_distance, VariogramCloud
    from pyinterpolate.variogram.empirical.experimental_variogram import calculate_semivariance
    distances = calc_point_to_point_distance(areal_centroids[:,:-1])
    maximum_range = np.max(distances) / 2
    743357.5182533629
    '''


    # 使用1000+个面的中心点做point——support
    psdf = pandas.DataFrame(areal_centroids, columns=["X", "Y", "TrmmValue"])

    point_support = PointSupport()
    psGdf = geopandas.GeoDataFrame(data=psdf, geometry=geopandas.points_from_xy(psdf.X, psdf.Y), crs='EPSG:3395')

    point_support.from_geodataframes(point_support_dataframe=psGdf, blocks_dataframe=filterGdf,
                                     point_support_geometry_col="geometry", point_support_val_col="TrmmValue",
                                     blocks_geometry_col="geometry", blocks_index_col="polygonIndex",
                                     use_point_support_crs=True)

    # 正则化解卷积训练半方差
    reg_mod = Deconvolution(verbose=True)
    #分十段
    reg_mod.fit(agg_dataset=blocks, point_support_dataset=point_support, agg_step_size=75000, agg_max_range=750000,
                variogram_weighting_method="closest", model_types="basic")
    initialVarFigPath = modelFilePrefix + "initialVarFig.png"
    initialVarFigTitle = str(yearmonth)+"_"
    deviationFigPath = modelFilePrefix + "deviationFig.png"
    deviationFigTitle = str(yearmonth)+"_Deviation: regularized semivariogram vs theoretical model"
    weightFigPath = modelFilePrefix + "weightFig.png"
    weightFigTitle= str(yearmonth)+"_weight change over each iteration"
    finalVarFigPath = modelFilePrefix + "finalVarFig.png"
    finalVarFigTitle = str(yearmonth)+"_"

    plot_variograms(reg_mod, initialVarFigPath,initialVarFigTitle)
    reg_mod.transform(max_iters=5)
    plot_deviations(reg_mod, deviationFigPath,deviationFigTitle)
    plot_weights(reg_mod, weightFigPath,weightFigTitle)
    plot_variograms(reg_mod, finalVarFigPath,finalVarFigTitle)
    reg_mod.export_model(modelFilePrefix + "regularizedModel.json")

    # 制作空值的点 127W
    print("{}: 创建新的栅格".format(yearmonth))
    degree_unit = 0.00925925933
    lonS = trmm_x_origin - degree_unit
    lonE = trmm_x_origin + pixel_width * srcArr.shape[1] + degree_unit
    latS = trmm_y_origin + degree_unit
    latE = trmm_y_origin - (pixel_width * srcArr.shape[0]) - degree_unit

    newWidth = int((lonE - lonS) / degree_unit) + 1
    newHeight = int((latS - latE) / degree_unit) + 1
    lonArr = [e * degree_unit + lonS for e in range(newWidth + 1)]
    latArr = [e * degree_unit + latE for e in range(newHeight + 1)]
    latArr.reverse()

    # print("总共点数：{}".format(len(lonArr) * len(latArr)))
    print("{}: 总共点数：{}".format(yearmonth, len(lonArr) * len(latArr)))
    posindexLst = []

    posLst = []
    for lon in lonArr:
        for lat in latArr:
            posLst.append([lon, lat, Point(lon, lat)])

    print("{}: 得到poslst".format(yearmonth))
    posDf = pandas.DataFrame(posLst, columns=["lon", "lat", "geometry"])
    # print("开始制作geodataframe")
    print("{}: 开始制作geodataframe".format(yearmonth))
    gdf_pts = geopandas.GeoDataFrame(posDf, crs='EPSG:4326')
    gdf_pts["geometry_lonlat"] = gdf_pts["geometry"]

    print("{}: 开始转换X,Y坐标".format(yearmonth))
    gdf_pts = gdf_pts.to_crs(crs='EPSG:3395')

    print("{}: 转换X,Y坐标成功".format(yearmonth))
    gdf_pts['x'] = gdf_pts.geometry.x
    gdf_pts['y'] = gdf_pts.geometry.y
    print(gdf_pts.head())

    if condition:
        semi_model = reg_mod.final_theoretical_model

        print("{} 没有已经跑完的插值文件 或者文件不合格 现在开始进行插值".format(yearmonth))
        # Build a map of interpolated values
        preds_cols = []
        errs_cols = []

        print("{}: model {} 插值开始 ".format(yearmonth, 'p ' + semi_model.name[:3] + ' ' + str(
            len(semi_model.lags))))
        try:
            kriged = kriging(observations=areal_centroids,
                             theoretical_model=semi_model,
                             points=gdf_pts[['x', 'y']].values,
                             neighbors_range=None,
                             no_neighbors=8,
                             use_all_neighbors_in_range=True,
                             allow_approx_solutions=True,
                             number_of_workers=workerNo
                             )

            # Interpolate missing values and uncertainty
            pred_col_name = 'p ' + semi_model.name[:3] + ' ' + str(len(semi_model.lags))
            uncertainty_col_name = 'e ' + semi_model.name[:3] + ' ' + str(len(semi_model.lags))
            gdf_pts[pred_col_name] = kriged[:, 0]
            gdf_pts[uncertainty_col_name] = kriged[:, 1]
            preds_cols.append(pred_col_name)
            errs_cols.append(uncertainty_col_name)

            for pcol1 in preds_cols:
                print('Column:', pcol1)
                for pcol2 in preds_cols:
                    if pcol1 == pcol2:
                        pass
                    else:
                        mad = gdf_pts[pcol1] - gdf_pts[pcol2]
                        mad = np.abs(np.mean(mad))
                        print(f'Mean Absolute Difference with {pcol2} is {mad:.4f}')
        except RuntimeError:
            print("{} 's {} ,{} model interpolation run error ,maybe it's singularity matrix ".format(yearmonth,
                                                                                                      semi_model.name,
                                                                                                      semi_model.lags))
        print("{}: model {} 插值结束 ".format(yearmonth,
                                              'p ' + semi_model.name[:3] + ' ' + str(len(semi_model.lags))))
        print("{}: 插值结束 开始写入磁盘".format(yearmonth))
        gdf_pts.to_csv(smoothedFile, sep="\t")
        print("{}: 写入插值结果结束".format(yearmonth))
        return "success"


'''
    # 开始 建立模型
    
    print("{}: 开始Create analysis parameters".format(yearmonth))

    maximum_range = 700000
    number_of_lags = [4, 8, 16]
    step_sizes = [maximum_range / x for x in number_of_lags]
    variogram_clouds = []
    for step_size in step_sizes:
        print("开始处理step size： {}".format(step_size))
        vc = VariogramCloud(input_array=areal_centroids, step_size=step_size, max_range=maximum_range + step_size)
        variogram_clouds.append(vc)
    # Now remove outliers from each lag
    print("{} :Now remove outliers from each lag".format(yearmonth))

    _ = [vc.remove_outliers(method='iqr', iqr_lower_limit=3, iqr_upper_limit=1.5, inplace=True) for vc in
         variogram_clouds]

    print("开始Create a theoretical semivariogram model")

    theoretical_semivariograms = []
    count = 1
    for idx, vc in enumerate(variogram_clouds):
        print(f'Semivariance calculated for {vc.lags} lags.')
        print('')
        # Calculate experimental model
        exp_model = vc.calculate_experimental_variogram()


        # Assign experimental model and data to TheoreticalSemivariogram

        theo_semi = TheoreticalVariogram()
        theo_semi.autofit(experimental_variogram=exp_model)

        theoretical_semivariograms.append(theo_semi)
        theo_semi.plot()
        theo_semi.to_json(modelFilePrefix + str(count) + ".json")
        count = count + 1
        print('')
        print('Model parameters:')
        print('Model type:', theo_semi.name)
        print('Nugget:', theo_semi.nugget)
        print('Sill:', theo_semi.sill)
        print('Range:', theo_semi.rang)
        print('Model error:', theo_semi.rmse)
        print('')
        print('#####')

   

    if condition:
        print("{} 没有已经跑完的插值文件 或者文件不合格 现在开始进行插值".format(yearmonth))
        # Build a map of interpolated values
        preds_cols = []
        errs_cols = []

        for semi_model in theoretical_semivariograms:
            print("{}: model {} 插值开始 ".format(yearmonth, 'p ' + semi_model.name[:3] + ' ' + str(
                len(semi_model.lags))))
            try:
                kriged = kriging(observations=areal_centroids,
                                 theoretical_model=semi_model,
                                 points=gdf_pts[['x', 'y']].values,
                                 neighbors_range=None,
                                 no_neighbors=8,
                                 use_all_neighbors_in_range=True,
                                 allow_approx_solutions=True,
                                 number_of_workers=workerNo
                                 )

                # Interpolate missing values and uncertainty
                pred_col_name = 'p ' + semi_model.name[:3] + ' ' + str(len(semi_model.lags))
                uncertainty_col_name = 'e ' + semi_model.name[:3] + ' ' + str(len(semi_model.lags))
                gdf_pts[pred_col_name] = kriged[:, 0]
                gdf_pts[uncertainty_col_name] = kriged[:, 1]
                preds_cols.append(pred_col_name)
                errs_cols.append(uncertainty_col_name)

                for pcol1 in preds_cols:
                    print('Column:', pcol1)
                    for pcol2 in preds_cols:
                        if pcol1 == pcol2:
                            pass
                        else:
                            mad = gdf_pts[pcol1] - gdf_pts[pcol2]
                            mad = np.abs(np.mean(mad))
                            print(f'Mean Absolute Difference with {pcol2} is {mad:.4f}')
            except RuntimeError:
                print("{} 's {} ,{} model interpolation run error ,maybe it's singularity matrix ".format(yearmonth,
                                                                                                          semi_model.name,
                                                                                                          semi_model.lags))
            print("{}: model {} 插值结束 ".format(yearmonth,
                                                  'p ' + semi_model.name[:3] + ' ' + str(len(semi_model.lags))))
        print("{}: 插值结束 开始写入磁盘".format(yearmonth))
        gdf_pts.to_csv(smoothedFile, sep="\t")
        # gdf_pts.to_file(smoothedShp)
        print("{}: 写入插值结果结束".format(yearmonth))

        # if not os.path.isfile(smoothedShp):
        #     gdf_pts_fromFile = geopandas.GeoDataFrame().from_file()
        #
        #     gdf_pts.to_file(smoothedShp, sep="\t")

        return "success"

'''

if __name__ == '__main__':

    taskDic = {

        "localhost.localdomain": ("^200[\d]", 36, 1, "/mnt/sata3_storage/jinlinfang/wpThesis/"),
        "sc153": ("^201\d", 48, 4, "/data/linfang.jin/sc_okr/reorganizedData")
    }
    pcName = socket.gethostname()
    import argparse

    dataRoot = r"C:\Users\crawl\Desktop\workForWP\wangpeng-master-thesis\reorganizedRawData"
    parallel = False

    parser = argparse.ArgumentParser(description='命令行中传入一个数字')
    parser.add_argument('--rootDir', '-r', type=str, required=False, help='根目录')
    parser.add_argument("-parallel", action='store_true', required=False, help='是否并行')
    parser.add_argument("--workerNo", '-w', type=int, required=False, help='插值单任务线程数')
    parser.add_argument("--processNo", '-c', type=int, required=False, help='并行进程数')
    parser.add_argument("--taskPattern", '-t', type=str, required=False, help='处理任务年份')

    args = parser.parse_args()
    if args.rootDir:
        dataRoot = r'' + args.rootDir
    else:
        if pcName in taskDic.keys():
            dataRoot = taskDic[pcName][3]

    if args.workerNo:
        workerNo = args.workerNo

    else:
        if pcName in taskDic.keys():
            workerNo = taskDic[pcName][2]

    if args.processNo:
        processNo = args.processNo
        pass
    else:
        if pcName in taskDic.keys():
            processNo = taskDic[pcName][1]
    if args.taskPattern:
        taskPattern = args.taskPattern

    else:
        if pcName in taskDic.keys():
            taskPattern = taskDic[pcName][1]

    parallel = args.parallel
    outDir = os.path.join(dataRoot, "TRMMtifATPKplusRegularization")
    tifSrcDir = os.path.join(dataRoot, "TRMMtifMaskedByQH")
    tifShpDir = os.path.join(dataRoot, "TRMMtif2ShpPolygon")

    srcFileDic = dict()
    paramsInfoLst = []
    for d in os.listdir(tifSrcDir):
        if not re.search(r'^\d+$', d):
            continue
        for f in os.listdir(os.path.join(tifSrcDir, d)):
            if not re.search(r'HDF.tif.clip.tif$', f):
                continue
            srcFile = os.path.join(tifSrcDir, d, f)
            yearmonth = re.search(r'(\d+).HDF.tif.clip.tif', f).group(1)
            srcFileDic[int(yearmonth)] = (srcFile, d, f)

    sortedSrcFileDic = dict(sorted(srcFileDic.items()))

    print(dataRoot)
    print(taskPattern)
    print(processNo)
    print(workerNo)

    for yearmonth in sortedSrcFileDic.keys():
        if not re.search(r'' + taskPattern, str(yearmonth)):
            continue

        (srcFile, d, f) = sortedSrcFileDic[yearmonth]
        yearmonth = str(yearmonth)
        year = int(str(yearmonth)[0:4])
        month = int(str(yearmonth)[4:])

        newName = re.sub(r'.HDF.tif.clip.tif', '.TRMM1KmAfterQHshpMaskAndATPK.tif', f)
        newSubDir = os.path.join(outDir, d)
        newpath = os.path.join(newSubDir, newName)
        pathlib.Path(newSubDir).mkdir(parents=True, exist_ok=True)

        shpName = yearmonth + ".TRMMPolygonShp.shp"
        tifShpFile = os.path.join(tifShpDir, d, shpName)
        allHours = 24 * calendar.monthrange(year, month)[1]

        # print(re.sub(r'\\', "/", "srcFile=\"" + srcFile + "\""))
        # print(re.sub(r'\\', "/", "tifShpFile=\"" + tifShpFile + "\""))
        # print(re.sub(r'\\', "/", "newpath=\"" + newpath + "\""))
        # print(re.sub(r'\\', "/", "monthHours=" + str(allHours)))

        paramsInfoLst.append((srcFile, tifShpFile, newpath, allHours))
        if not parallel:
            atpkOKProcessParallel(srcFile, tifShpFile, newpath, allHours)

    # #并行化处理
    rLst = []
    if parallel:
        pool = multiprocessing.Pool(processNo)
        processLst = [pool.apply_async(atpkOKProcessParallel, args=param) for param in paramsInfoLst]
        for p in processLst:
            r = p.get()
            rLst.append(r)
    print(rLst)
