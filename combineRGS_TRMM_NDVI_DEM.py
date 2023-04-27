import calendar
import datetime
import json
import os.path
import re, pathlib
import shutil
from osgeo import gdal
from shapely.geometry import Point
import geopandas
import numpy

try:
    from osgeo import gdal
    from osgeo import ogr
except Exception:
    import gdal
    import ogr
import pandas


def parseTif(tif):
    infoDic = dict()
    dr = gdal.Open(tif)
    transform = dr.GetGeoTransform()
    infoDic["x_origin"] = transform[0]
    infoDic["y_origin"] = transform[3]
    infoDic["widthRes"] = transform[1]
    infoDic["heightRes"] = transform[5]

    srcArr = None
    try:
        srcArr = dr.ReadAsArray()
    except Exception:
        pass
    infoDic["arr"] = srcArr
    return infoDic


def collectRGSinfo(froot):
    print(froot)
    sumInfo = dict()
    yearDic = dict()
    for f in os.listdir(froot):
        m = re.search(r'qh_(\d+).txt', f)
        if not m:
            continue
        yearmonth = m.group(1)
        year = str(yearmonth)[0:4]
        month = str(yearmonth)[4:]
        fpath = os.path.join(froot, f)
        df = pandas.read_csv(fpath, sep="\t", index_col=0)
        df = df[["rgs"]]
        tempDic = dict()
        for index, row in df.iterrows():
            tempDic[index] = row["rgs"]

        if not year in sumInfo.keys():
            sumInfo[year] = {yearmonth: tempDic}
        else:
            sumInfo[year][yearmonth] = tempDic
        if not year in yearDic.keys():
            # {year:{id: rgssum}}
            yearDic[year] = tempDic
        else:
            for sid in tempDic.keys():
                yearDic[year][sid] = tempDic[sid] + yearDic[year][sid]

    for year in sumInfo.keys():
        sumInfo[year]["yearlySumAmount"] = yearDic[year]

    return sumInfo


def list_all_files(rootdir):
    import os
    _files = []
    # 列出文件夹下所有的目录与文件
    list = os.listdir(rootdir)
    for i in range(0, len(list)):
        # 构造路径
        path = os.path.join(rootdir, list[i])
        # 判断路径是否为文件目录或者文件
        # 如果是目录则继续递归
        if os.path.isdir(path):
            _files.extend(list_all_files(path))
        if os.path.isfile(path):
            _files.append(path)
    return _files


def mergeRgsFile(files, sumFile):
    cols = ["stationID", "stationLat", "stationLon", "stationEle", "YEAR", "MONTH", "rgs"]
    tpltDf = None
    for f in files:
        tempDf = pandas.read_csv(f, sep="\t", index_col=0, header=0)
        if tpltDf is None:
            tpltDf = tempDf
            tpltDf["stationID"] = tpltDf.index
            continue
        for index, row in tempDf.iterrows():
            if not index in tpltDf.index:
                print("")
                continue
            tpltDf.loc[index, "rgs"] = tpltDf.loc[index, "rgs"] + row["rgs"]

    tpltDf = tpltDf[cols]
    tpltDf.to_csv(sumFile, sep="\t")


def mergeTrmmFile(files, sumFile):
    for f in files:
        pass

    return


def mergeNdviFile(files, sumFile):
    yearTifLst = files
    sumTif = sumFile
    sumGtf = ()
    sumPrj = None
    arrLst = []
    dt = None

    for tif in yearTifLst:
        # 去每月的nvdi数据 进行求平均
        dr = gdal.Open(tif)
        b = dr.GetRasterBand(1)
        dt = b.DataType

        transform = dr.GetGeoTransform()
        sumGtf = transform

        trmm_x_origin = transform[0]
        trmm_y_origin = transform[3]
        pixel_width = transform[1]
        pixel_height = transform[5]
        prj = dr.GetProjection()
        sumPrj = prj

        srcArr = None
        try:
            srcArr = dr.ReadAsArray()
        except Exception:
            pass
        arrLst.append(srcArr)
    print("start to find mean")

    yearMeanArr = numpy.apply_along_axis(numpy.max, axis=0, arr=numpy.array(arrLst))
    print("finish find mean")
    xSize = yearMeanArr.shape[1]
    ySize = yearMeanArr.shape[0]

    # print(yearMeanArr)

    tiffDriver = gdal.GetDriverByName("Gtiff")
    out_ds = tiffDriver.Create(sumTif, xSize, ySize, 1, dt)
    out_ds.SetProjection(sumPrj)
    out_ds.SetGeoTransform(sumGtf)
    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(yearMeanArr)

    out_ds.FlushCache()
    out_ds.BuildOverviews("average", [2, 4, 8, 16, 32])
    del out_ds


def mergeCsvFile(files, sumFile):
    pass


def collectFiles(fileRoot=None, monthlyPattern=None, yearlyPattern=None, dataType=None):
    fileInfo = dict()
    allFiles = list_all_files(fileRoot)
    yearDic = dict()
    eg_ext = ".tif"
    for f in allFiles:
        fname = os.path.basename(f)

        mm = re.search(r"" + monthlyPattern, fname)
        if mm:
            yearmonth = mm.group(1)
            fileInfo[yearmonth] = f
            eg_ext = os.path.splitext(f)[1]
            year = yearmonth[:4]
            if not year in yearDic.keys():
                yearDic[year] = [f]
            else:
                yearDic[year].append(f)
            continue

        if yearlyPattern:
            ym = re.search(r"" + yearlyPattern, fname)
            if ym:
                fileInfo[ym.group(1)] = f

    if not yearlyPattern:
        for y in yearDic.keys():
            fLst = yearDic[y]
            yearFile = os.path.join(fileRoot, y + "_merged" + eg_ext)
            fileInfo[y] = yearFile
            if dataType in ["rgs"]:
                mergeRgsFile(fLst, yearFile)
                continue
            if dataType in ["tif"]:
                mergeTrmmFile(fLst, yearFile)
                continue
            if dataType in ["csv"]:
                mergeCsvFile(fLst, yearFile)
                continue
            if dataType in ["ndvi"]:
                mergeNdviFile(fLst, yearFile)
                continue
    return fileInfo


def obtainInfo(t, rgsBasicDf, rgsFiles, demRoot, ndviFiles, trmm25Files,
               trmm1Files):
    cols = ["stationId", "lon", "lat", "time", "rgs", "ndvi", "trmm25", "trmm1", "slope", "elevation", "aspect"]
    finalDf = pandas.DataFrame(columns=cols)
    demSlopeFile = os.path.join(demRoot, "qh_1km_slope.tif")
    demElevationFile = os.path.join(demRoot, "qh_1km_elevation.tif")
    demAspectFile = os.path.join(demRoot, "qh_1km_aspect.tif")

    slopeData = []
    elevationData = []
    aspectData = []

    trmm25Data = []
    trmm1Data = []
    ndviData = []

    ndviFile = ndviFiles[t]
    rgsFile = rgsFiles[t]
    trmm25File = trmm25Files[t]
    trmm1File = trmm1Files[t]


def parseDemFile(demDir):
    demSlopeFile = os.path.join(rootDEM1kmDataDir, "qh_1km_slope.tif")
    demElevationFile = os.path.join(rootDEM1kmDataDir, "qh_1km_elevation.tif")
    demAspectFile = os.path.join(rootDEM1kmDataDir, "qh_1km_aspect.tif")


    pass


if __name__ == '__main__':

    rootDir = r"C:\Users\crawl\Desktop\workForWP\wangpeng-master-thesis\reorganizedRawData"
    rootRGSdataDir = os.path.join(rootDir, "RGSData-afterExceptionProcessAndStationDeletion")
    eg_rgs_file = os.path.join(rootRGSdataDir, "qh_200001.txt")
    rootTRMM25kmdataDir = os.path.join(rootDir, "TRMMtifMaskedByQH")
    rootTRMM1kmdataDir = os.path.join(rootDir, "TRMMtifATPKplusRegularization")
    rootNDVI1kmDataDir = os.path.join(rootDir, "NDVI_1KM")
    rootDEM1kmDataDir = os.path.join(rootDir, "DEM_1KM_raw")
    demSlopeFile = os.path.join(rootDEM1kmDataDir, "qh_1km_slope.tif")
    demElevationFile = os.path.join(rootDEM1kmDataDir, "qh_1km_elevation.tif")
    demAspectFile = os.path.join(rootDEM1kmDataDir, "qh_1km_aspect.tif")

    '''
    stationid  stationlon  stationlat X Y geometry  
    '''
    rgsBasicData = pandas.read_csv(eg_rgs_file, sep="\t", index_col=0)
    glst = []
    for index, row in rgsBasicData.iterrows():
        glst.append(Point(row["stationLon"], row["stationLat"]))
    rgsBasicData["geometry"] = glst
    rgsBasicDataGdf = geopandas.GeoDataFrame(rgsBasicData, crs='EPSG:4326')
    rgsBasicDataGdf["lonlat"] = rgsBasicDataGdf["geometry"]
    rgsBasicDataGdf = rgsBasicDataGdf.to_crs(crs='EPSG:3395')
    rgsBasicDataGdf["stationID"] = rgsBasicDataGdf.index
    selectCols = ["stationID", 'lonlat', 'geometry']
    rgsBasicDataGdf = rgsBasicDataGdf[selectCols]
    rgsDataDic = collectRGSinfo(rootRGSdataDir)
    ndviFiles = collectFiles(fileRoot=rootNDVI1kmDataDir, monthlyPattern="(\d{6}).NVDI.1KM.tif",
                             yearlyPattern="^(\d{4})_merged.tif$", dataType="ndvi")
    rgsFiles = collectFiles(fileRoot=rootRGSdataDir, monthlyPattern="qh_(\d{6}).txt", dataType="rgs")
    # 2000.yealySumPrecipitation.HDF.tif.clip.tif
    trmm25Files = collectFiles(fileRoot=rootTRMM25kmdataDir,
                               monthlyPattern="(\d{6}).HDF.tif.clip.tif$",
                               yearlyPattern="^(\d{4}).yealySumPrecipitation.HDF.tif.clip.tif$",
                               dataType="tif")
    # 2000.yealySumPrecipitation.TRMM1KmAfterQHshpMaskAndATPK.smoothedPointMatrixOKatpk.csv
    trmm1Files = collectFiles(fileRoot=rootTRMM1kmdataDir,
                              monthlyPattern="(\d{4,6}).TRMM1KmAfterQHshpMaskAndATPK.smoothedPointMatrixOKatpk.csv$",
                              yearlyPattern="^(\d{4}).yealySumPrecipitation.TRMM1KmAfterQHshpMaskAndATPK.smoothedPointMatrixOKatpk.csv$",
                              dataType="csv")

    targetMonthLst = ["200001"]
    targetYearLst = ["2000"]
    demDic=parseDemFile(rootDEM1kmDataDir)

    for t in targetMonthLst:
        if not re.search(r'^\d{6}$', t):
            continue
        obtainInfo(t, rgsBasicDataGdf, rgsFiles, rootDEM1kmDataDir, ndviFiles, trmm25Files, trmm1Files)

    for t in targetYearLst:
        if not re.search(r'^\d{4}$', t):
            continue
        obtainInfo(t, rgsBasicDataGdf, rgsFiles, rootDEM1kmDataDir, ndviFiles, trmm25Files, trmm1Files)
