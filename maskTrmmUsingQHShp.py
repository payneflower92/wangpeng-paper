import calendar
import datetime
import json
import os.path
import re, pathlib
import shutil,numpy
import numpy as np
from osgeo import gdal, gdal_array

try:
    from osgeo import gdal
    from osgeo import ogr
except Exception:
    import gdal
    import ogr
import pandas

try:
    import Image
    import ImageDraw
except:
    from PIL import Image, ImageDraw

# _*_ coding: utf-8 _*_


import operator
from osgeo import gdal, gdal_array, osr
import shapefile


def image2Array(i):
    """
    将一个Python图像库的数组转换为一个gdal_array图片
    """
    a = gdal_array.numpy.fromstring(i.tobytes(), 'b')
    a.shape = i.im.size[1], i.im.size[0]
    return a


def world2Pixel(geoMatrix, x, y):
    """
    使用GDAL库的geomatrix对象((gdal.GetGeoTransform()))计算地理坐标的像素位置
    """
    ulx = geoMatrix[0]
    uly = geoMatrix[3]
    xDist = geoMatrix[1]
    yDist = geoMatrix[5]
    rtnX = geoMatrix[2]
    rtnY = geoMatrix[4]
    pixel = int((x - ulx) / xDist)
    line = int((uly - y) / abs(yDist))
    return (pixel, line)


# 用于裁剪的栅格数据
raster = 'C:/Users/crawl/Desktop/workForWP/wangpeng-master-thesis/reorganizedRawData/TiffFromTRMMAfterqinghaiShp0.25Degree/2000/TRMM0.25AfterQHshpMask_200001.tif'
# 用于裁剪的多边形shp文件
shp = "C:/Users/crawl/Desktop/codeProject/wangpengMasterThesis/reOrganizedData/BOUNDARY_SHP/Qinghai.shp"

m = re.search(r'^TRMM0.25AfterQHshpMask_(\d{6}).tif', os.path.basename(raster))
yearmonth = None
if m:
    yearmonth = m.group(1)
year = int(yearmonth[0:4])
month = int(yearmonth[4:])
allHours = 24 * calendar.monthrange(year, month)[1]
# 裁剪后的栅格数据
output = os.path.join(os.path.dirname(raster), os.path.basename(raster) + ".gdalMosaicQh.dataframe.tif")
# outArray = []
# 将数据源作为gdal_array载入
srcArray = gdal_array.LoadFile(raster)
# 同时载入gdal库的图片从而获取geotransform
srcImage = gdal.Open(raster)
geoTrans = srcImage.GetGeoTransform()
# # 使用PyShp库打开shp文件
# r = shapefile.Reader("{}".format(shp))
# # 将图层扩展转换为图片像素坐标
# minX, minY, maxX, maxY = r.bbox
# ulX, ulY = world2Pixel(geoTrans, minX, maxY)
# lrX, lrY = world2Pixel(geoTrans, maxX, minY)
# # 计算新图片的尺寸
# pxWidth = int(lrX - ulX)+1
# pxHeight = int(lrY - ulY)+1
# clip = srcArray
# # 为图片创建一个新的geomatrix对象以便附加地理参照数据
# geoTrans = list(geoTrans)
# geoTrans[0] = minX
# geoTrans[3] = maxY
# # 在一个空白的8字节黑白掩膜图片上把点映射为像元绘制市县
# # 边界线
# pixels = []
# for p in r.shape(0).points:
#     pixels.append(world2Pixel(geoTrans, p[0], p[1]))
# rasterPoly = Image.new("L", (pxWidth, pxHeight), 1)
# # 使用PIL创建一个空白图片用于绘制多边形
# rasterize = ImageDraw.Draw(rasterPoly)
# rasterize.polygon(pixels, 0)
# # 使用PIL图片转换为Numpy掩膜数组
# mask = image2Array(rasterPoly)
# # 根据掩膜图层对图像进行裁剪
# clip = gdal_array.numpy.choose(mask, (clip, 0)).astype(gdal_array.numpy.float32)
#
# for x_index in range(mask.shape[0]):
#     for y_index in range(mask.shape[0]):
#         if mask[x_index, y_index] == 1:
#             continue
#         lon = 0
#         lat = 0
#         trmmValue = clip[x_index, y_index] * allHours
#         outArray.append([lon, lat, trmmValue])

# gdal_array.SaveArray(clip, "{}.tif".format(output),
#                      format="GTiff", prototype=raster)

from osgeo import gdal,gdalconst
shppath = r'D:\Africa\Africa_city.shp'
tifpath = r'D:\regionImg\VNL_2012Africa.tif'
outtif1 = r'D:\Africa\Africa_FID0.tif'
cutlineWhere = 'FID = 2485'
ds = gdal.Warp(
output,    #裁剪后图像保存的完整路径（包括文件名）
raster,    #待裁剪的影像完整路径（包括文件名）
format='GTiff', # 保存图像的格式
cutlineDSName=shp, # 矢量文件的完整路径
cropToCutline=True, # 保证裁剪后影像大小跟矢量文件的图框大小一致（设置为False时，结果图像大小会跟待裁剪影像大小一样，则会出现大量的空值区域）
# cutlineWhere=cutlineWhere #矢量文件筛选条件,
#dstNodata=0
)




