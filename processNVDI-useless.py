
import calendar
import os.path
import re, pathlib
import sys

import numpy
import pyinterpolate
try:
    from osgeo import gdal
    from osgeo import ogr
except Exception:
    import gdal
    import ogr
import pandas, geopandas
from geopandas import points_from_xy
from pyinterpolate import build_experimental_variogram, Blocks, PointSupport, Deconvolution, TheoreticalVariogram, smooth_blocks




dataRoot=r'C:\Users\crawl\Desktop\workForWP\wangpeng-master-thesis\reorganizedRawData'
if sys.platform =="linux":
    dataRoot="/mnt/sata3_storage/jinlinfang/wpThesis"

outDir=os.path.join(dataRoot,"NDVIresampleTo1KM")
tifSrcDir = os.path.join(dataRoot,"NDVI")
# tifShpDir = os.path.join(dataRoot,"TRMMtif2ShpPolygon")

mergeDic=dict()


for d in os.listdir(tifSrcDir):
    for f in os.listdir(os.path.join(tifSrcDir, d)):
        m= re.search(r'(\d+).(h2[56]).hdf$', f)
        if not m:
            continue
        srcFile = os.path.join(tifSrcDir,d,f)

        yearmonth = m.group(1)
        blockType=m.group(2)

        if yearmonth not in mergeDic.keys():
            mergeDic[yearmonth][blockType] = srcFile
            # outMergeNVDI = os.path.join()
            # mergeDic[yearmonth]["outMergeNVDI"] = os.path.join(outDir , outMergeNVDIName)
            # mergeDic[yearmonth]["outMergeNVDIto1KM"] = os.path.join(outDir , outMergeNVDIto1KMName)
            # mergeDic[yearmonth]["outMergeNVDIto0.25Degree"] = os.path.join(outDir , outMergeNVDIto25DegreeName)

        # mergeDic[blockType] =
        year = int(yearmonth[0:4])
        month = int(yearmonth[4:])
        newName = re.sub(r'.h2[56].hdf$', '.TRMM1KmAfterQHshpMaskAndATPK.tif', f)
        newSubDir = os.path.join(outDir, d)
        pathlib.Path(newSubDir).mkdir(parents=True, exist_ok=True)
        newpath = os.path.join(newSubDir, newName)
        srcFile = os.path.join(tifSrcDir, d, f)
        shpName = yearmonth+"TRMMPolygonShp.shp"
        # tifShpFile=os.path.join(tifShpDir,d,shpName)
        allHours = 24 * calendar.monthrange(year, month)[1]

        print(re.sub(r'\\', "/", "srcFile=\"" + srcFile + "\""))
        # print(re.sub(r'\\', "/", "tifShpFile=\"" + tifShpFile + "\""))
        print(re.sub(r'\\', "/", "newpath=\"" + newpath + "\""))
        # print(re.sub(r'\\', "/", "allHours=" + str(allHours) ))

        #start
        ds = gdal.Open(srcFile)
        dsNDVI = ds.GetSubDatasets()[0][0]
        rasterNDVI = gdal.Open(dsNDVI)
        nvdiArr = rasterNDVI.ReadAsArray()
        Metadata = ds.GetMetadata()




