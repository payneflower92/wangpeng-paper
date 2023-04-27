import datetime
import json
import os.path
import re,pathlib
import shutil

import numpy
from osgeo import gdal


rootDir = r"C:\Users\crawl\Desktop\U盘\5_scale"
newRootDir=r'C:\Users\crawl\Desktop\workForWP\wangpeng-master-thesis\reorganizedRawData\NDVI_1KM'
mirrorDic = dict()

ndviFileDic = dict()
years = range(2000, 2020)
months = range(1, 13)
ndviDir = rootDir
for f in os.listdir(ndviDir):
    #QH_WGS84_MOD13A3.A2020336.NDVI.tif
    m = re.search(r'^QH_WGS84_MOD13A3.A(\d{7}).NDVI.tif$', f.strip())
    if m:
        tackleDate = m.group(1)
        year = int(tackleDate[0:4])
        dateNo = int(tackleDate[4:7])
        firstDay = datetime.datetime(year, 1, 1)
        wantedDay = firstDay + datetime.timedelta(dateNo - 1)
        wantedDayYear = year
        wantedDayBeforeMonth = int(datetime.datetime.strftime(wantedDay, "%m")) - 1
        if wantedDayBeforeMonth == 0:
            wantedDayYear = wantedDayYear - 1
            wantedDayBeforeMonth = 12
        keyStr = "%04d" % wantedDayYear + "%02d" % wantedDayBeforeMonth
        ndviFileDic[keyStr] = os.path.join(ndviDir, f)
mirrorDic = dict()
yearDic = dict()
for year in years:
    if year not in yearDic.keys():
        yearDic[year] = dict()
    for month in months:
        timeStr = '%04d' % year + '%02d' % month
        newOutDir = os.path.join(newRootDir, '%04d' % year)
        pathlib.Path(newOutDir).mkdir(parents=True, exist_ok=True)
        destPath = os.path.join(newOutDir, timeStr + ".NVDI.1KM.tif")
        if timeStr in ndviFileDic.keys():
            # shutil.copyfile(ndviFileDic[timeStr], destPath)
            mirrorDic[destPath] = ndviFileDic[timeStr]
            yearDic[year][month] = destPath
        else:
            print("不正常的情況：{}]" % [timeStr])

# 处理年份

for year in years:
    yearTifLst = yearDic[year].values()
    sumTif = os.path.join(os.path.join(newRootDir,str(year), "{}.meanNVDI1KM.tif".format(year) ))
    sumGtf = ()
    sumPrj = None
    arrLst=[]
    dt=None

    for tif in yearTifLst:
        # 去每月的nvdi数据 进行求平均
        dr = gdal.Open(tif)
        b=dr.GetRasterBand(1)
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

    yearMeanArr=numpy.apply_along_axis(numpy.mean,axis=0,arr=numpy.array(arrLst))
    print("finish find mean")
    xSize=yearMeanArr.shape[1]
    ySize=yearMeanArr.shape[0]

    # print(yearMeanArr)

    tiffDriver = gdal.GetDriverByName("Gtiff")
    out_ds = tiffDriver.Create(sumTif,xSize,ySize,1,dt )
    out_ds.SetProjection(sumPrj)
    out_ds.SetGeoTransform(sumGtf)
    out_band=out_ds.GetRasterBand(1)
    out_band.WriteArray(yearMeanArr)

    out_ds.FlushCache()
    out_ds.BuildOverviews("average",[2,4,8,16,32])
    del out_ds











