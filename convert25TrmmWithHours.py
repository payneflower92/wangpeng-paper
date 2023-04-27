import calendar
import os.path
import re, pathlib
try:
    from osgeo import gdal
    from osgeo import ogr
except Exception:
    import gdal
    import ogr

import geopandas


outDir = r"C:\Users\crawl\Desktop\workForWP\wangpeng-master-thesis\reorganizedRawData\TRMMtifMaskedByQHWithHours"
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
        newName = re.sub(r'.HDF.tif.clip.tif', '.withHours.tif', f)
        newSubDir = os.path.join(outDir, d)
        pathlib.Path(newSubDir).mkdir(parents=True, exist_ok=True)
        newpath = outshp = os.path.join(newSubDir, newName)
        srcFile = raster = os.path.join(srcDir, d, f)
        print(re.sub(r'\\', "/", "srcFile=\"" + srcFile + "\""))
        print(re.sub(r'\\', "/", "newpath=\"" + newpath + "\""))
        allHours = 24 * calendar.monthrange(year, month)[1]

        inraster = gdal.Open(raster)  # 读取路径中的栅格数据
        inband = inraster.GetRasterBand(1)  # 这个波段就 是最后想要转为矢量的波段，如果是单波段数据的话那就都是1
        prj = osr.SpatialReference()
        prj.ImportFromWkt(inraster.GetProjection())  # 读取栅格数据的投影信息，用来为后面生成的矢量做准备

