import calendar
import datetime
import json
import os.path
import re, pathlib
import shutil

import pandas as pd
import shapefile
from osgeo import gdal, osr

try:
    from osgeo import gdal
    from osgeo import ogr
except Exception:
    import gdal
    import ogr

import numpy as np
import pandas, geopandas

from geopandas import points_from_xy

from pyinterpolate import TheoreticalVariogram
from pyinterpolate import Blocks, PointSupport
from pyinterpolate import smooth_blocks

outDir = r"C:\Users\crawl\Desktop\workForWP\wangpeng-master-thesis\reorganizedRawData\TRMMtif2ShpPolygonWithHours"
srcDir = r"C:\Users\crawl\Desktop\workForWP\wangpeng-master-thesis\reorganizedRawData\TRMMtifMaskedByQH"

from pyproj import CRS
from pyproj import Transformer

crs_WGS84 = CRS.from_epsg(4326)  # WGS84地理坐标系
crs_WebMercator = CRS.from_epsg(3395)  # Web墨卡托投影坐标系
import osgeo.osr as osr  # 使用gdal方法


def point_transform(source_ref, target_ref, x, y):
    # 创建目标空间参考
    spatialref_target = osr.SpatialReference()
    spatialref_target.ImportFromEPSG(target_ref)  # 转换坐标
    # 创建原始空间参考
    spatialref_source = osr.SpatialReference()
    spatialref_source.ImportFromEPSG(source_ref)
    # 原始坐标
    # 构建坐标转换对象，用以转换不同空间参考下的坐标
    trans = osr.CoordinateTransformation(spatialref_source, spatialref_target)
    coordinate_after_trans = trans.TransformPoint(x, y)
    return coordinate_after_trans  # 一个二维向量

infoDic=dict()

for d in os.listdir(srcDir):
    for f in os.listdir(os.path.join(srcDir, d)):
        if not re.search(r'HDF.tif.clip.tif$', f):
            continue
        yearmonth = re.search(r'(\d+).HDF.tif.clip.tif', f).group(1)
        year = int(yearmonth[0:4])
        month = int(yearmonth[4:])
        newName = re.sub(r'.HDF.tif.clip.tif', '.TRMMPolygonShp.shp', f)
        xyShpName = re.sub(r'.HDF.tif.clip.tif', '.TRMMPolygonShp.lonlat.shp', f)
        perfectShpName = re.sub(r'.HDF.tif.clip.tif', '.TRMMPolygonShp.perfect.shp', f)
        newSubDir = os.path.join(outDir, d)
        pathlib.Path(newSubDir).mkdir(parents=True, exist_ok=True)
        newpath = outshp = os.path.join(newSubDir, newName)
        perfectShp=os.path.join(newSubDir,perfectShpName)
        # finalShp = os.path.join(newSubDir, xyShpName)
        srcFile = raster = os.path.join(srcDir, d, f)
        print(re.sub(r'\\', "/", "srcFile=\"" + srcFile + "\""))
        print(re.sub(r'\\', "/", "newpath=\"" + newpath + "\""))
        allHours = 24 * calendar.monthrange(year, month)[1]



        inraster = gdal.Open(raster)  # 读取路径中的栅格数据
        inband = inraster.GetRasterBand(1)  # 这个波段就是最后想要转为矢量的波段，如果是单波段数据的话那就都是1
        prj = osr.SpatialReference()
        prj.ImportFromWkt(inraster.GetProjection())  # 读取栅格数据的投影信息，用来为后面生成的矢量做准备

        drv = ogr.GetDriverByName("ESRI Shapefile")
        if os.path.exists(outshp):  # 若文件已经存在，则删除它继续重新做一遍
            drv.DeleteDataSource(outshp)
        Polygon = drv.CreateDataSource(outshp)  # 创建一个目标文件
        Poly_layer = Polygon.CreateLayer(raster[:-4], srs=prj, geom_type=ogr.wkbMultiPolygon)  # 对shp文件创建一个图层，定义为多个面类
        newField = ogr.FieldDefn('value', ogr.OFTReal)  # 给目标shp文件添加一个字段，用来存储原始栅格的pixel value
        Poly_layer.CreateField(newField)

        gdal.FPolygonize(inband, None, Poly_layer, 0)  # 核心函数，执行的就是栅格转矢量操作
        Polygon.SyncToDisk()
        Polygon = None
        infoDic[outshp] = [allHours,perfectShp]

for shp in infoDic.keys():


    allHours = infoDic[shp][0]
    perfectShp = infoDic[shp][1]
    print(shp)
    print(allHours)
    print(perfectShp)

    gdf = geopandas.GeoDataFrame.from_file(shp, encoding='gb18030')
    gdf.columns = ["TrmmValue", "geometry"]
    filterGdf = gdf[gdf["TrmmValue"] >= 0]
    filterGdf["TrmmValue"] = filterGdf["TrmmValue"] * allHours
    filterGdf["geometry_lonlat"] = filterGdf["geometry"]
    filterGdf = filterGdf.to_crs(crs='EPSG:3395')
    filterGdf = filterGdf[["geometry","TrmmValue"]]
    filterGdf.to_file(perfectShp)


