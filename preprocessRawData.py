# -*- coding:utf-8 -*-
# author:jinlinfang
# datetime:2023/2/18/0018 21:56
# software: PyCharm
import datetime
import json
import os.path
import re,pathlib
import shutil
from osgeo import gdal


def reorganizeTRMMdata(rawDataRoot, outDataRoot, dataType):
    pass


if __name__ == '__main__':
    years = range(2000, 2020)
    months = range(1, 13)
    dataType = {"TRMM": ["TRMM"], "NDVI": [], "DEM": [], "RGS": []}
    rawDataRoot = "D:\\Data\\TRMM_Qinghai\\"
    outDataRoot = "D:\\workForWP\\wangpeng-master-thesis\\reorganizedRawData"

    trmmFileDic=dict()
    trmmDir=os.path.join(rawDataRoot,"TRMM","trmm_hdf")
    for f in os.listdir(trmmDir):
        m=re.search(r'^3B43.(\d{6})\S+.HDF$', f.strip())
        if m:
            if m.group(1) in trmmFileDic.keys():
                print("原始TRMM有2份:" + m.group(1))
            trmmFileDic[m.group(1)] = os.path.join(trmmDir,f)
    mirrorDic=dict()
    for year in years:
        for month in months:
            timeStr='%04d' % year +'%02d' % month
            newOutDir=os.path.join(outDataRoot,"TRMM", '%04d' % year )
            pathlib.Path(newOutDir).mkdir(parents=True, exist_ok=True)
            destPath=os.path.join(newOutDir,timeStr+".HDF")
            if timeStr in trmmFileDic.keys():
                shutil.copyfile(trmmFileDic[timeStr], destPath)
                mirrorDic[destPath] = trmmFileDic[timeStr]
            else:
                print("不正常的情況："+ timeStr)
    with open(os.path.join(outDataRoot,"TRMM","文件来源.json"),'w')  as jf:
        json.dump(mirrorDic,jf)

    ndviFileDic = dict()
    ndviDir = os.path.join(rawDataRoot, "MOD13A3")
    for f in os.listdir(ndviDir):
        m = re.search(r'^MOD13A3.A(\d{7}).(h25|h26)v\S+.hdf$', f.strip())
        if m:
            tackleDate=m.group(1)
            sliceNo = m.group(2)
            year=int(tackleDate[0:4])
            dateNo=int(tackleDate[4:7])
            firstDay=datetime.datetime(year,1,1)
            wantedDay=firstDay+datetime.timedelta(dateNo-1)
            wantedDayYear=year
            wantedDayBeforeMonth=int(datetime.datetime.strftime(wantedDay,"%m")) -1
            if wantedDayBeforeMonth==0:
                wantedDayYear=wantedDayYear-1
                wantedDayBeforeMonth=12
            keyStr="%04d" % wantedDayYear + "%02d" %wantedDayBeforeMonth

            if (keyStr,sliceNo) in ndviFileDic.keys():
                print("原始ndvi有2份:{}, {}"  % [keyStr,sliceNo])
            ndviFileDic[(keyStr,sliceNo)] =os.path.join(ndviDir,f)

    mirrorDic = dict()
    for year in years:
        for month in months:
            for sliceType in ["h25", "h26"]:
                timeStr = '%04d' % year + '%02d' % month
                newOutDir = os.path.join(outDataRoot, "NDVI", '%04d' % year)
                pathlib.Path(newOutDir).mkdir(parents=True, exist_ok=True)
                destPath = os.path.join(newOutDir, timeStr +"." + sliceType + ".hdf")
                if (timeStr,sliceType) in ndviFileDic.keys():
                    shutil.copyfile(ndviFileDic[(timeStr,sliceType)], destPath)
                    mirrorDic[destPath] = ndviFileDic[(timeStr,sliceType)]
                else:
                    print("不正常的情況：{},{}" % [timeStr,sliceType])
    with open(os.path.join(outDataRoot, "NDVI", "文件来源.json"), 'w') as jf:
        json.dump(mirrorDic, jf)


    #拼接DEM数据  参考网址https://blog.csdn.net/qq_39085138/article/details/120457154

    import re
    import subprocess
    flist = [os.path.join("D:\\Data\\TRMM_Qinghai\\GDEM_V2", e) for e in os.listdir("D:\\Data\\TRMM_Qinghai\\GDEM_V2")
             if re.search(r'ASTGTM2_\S+_dem.tif$', e)]
    path_utility = "C:\\Users\\payneflower\\anaconda3\\envs\\pythonEnv37\\Scripts\\gdal_merge.py"
    pathout = "testmergedem.tif"
    filein = " ".join(flist)
    options = "COMPRESS=LZW"
    cmd = 'python {0} -co {3} -o {1} {2}'.format(path_utility, pathout, filein, options)

    #重采样DEM数据 90m->1km  https://github.com/xbr2017/PyGdal_resample/blob/master/code/resample.py

    in_ds = gdal.Open('nat_color.tif')
    out_rows = int(in_ds.RasterYSize / 50)
    out_columns = int(in_ds.RasterXSize / 50)
    num_bands = in_ds.RasterCount

    gtiff_driver = gdal.GetDriverByName('GTiff')
    out_ds = gtiff_driver.Create('nat_color_resampled.tif',
                                 out_columns, out_rows, num_bands)

    out_ds.SetProjection(in_ds.GetProjection())
    geotransform = list(in_ds.GetGeoTransform())
    geotransform[1] *= 50
    geotransform[5] *= 50
    out_ds.SetGeoTransform(geotransform)

    data = in_ds.ReadRaster(
        buf_xsize=out_columns, buf_ysize=out_rows)
    out_ds.WriteRaster(0, 0, out_columns, out_rows, data)
    out_ds.FlushCache()




    #处理rgs





















