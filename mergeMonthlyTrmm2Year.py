import calendar
import os.path
import re, pathlib

import numpy

try:
    from osgeo import gdal
    from osgeo import ogr,osr
except Exception:
    import gdal
    import ogr


import os

import fiona
import rasterio
import rasterio.mask

if __name__ == '__main__':


    dataRoot = r"C:\Users\crawl\Desktop\workForWP\wangpeng-master-thesis\reorganizedRawData"
    tifSrcDir = os.path.join(dataRoot, "TRMMtifMaskedByQH")
    outDir = os.path.join(dataRoot, "TRMMtifMaskedByQH/yearly")

    srcFileDic=dict()

    for d in os.listdir(tifSrcDir):
        if not re.search(r'^\d+$', d):
            continue
        for f in os.listdir(os.path.join(tifSrcDir, d)):
            if not re.search(r'HDF.tif.clip.tif$', f):
                continue
            srcFile = os.path.join(tifSrcDir, d, f)
            year = re.search(r'(\d{4})\d{2}.HDF.tif.clip.tif', f).group(1)
            if  int(year) not in srcFileDic.keys():
                srcFileDic[int(year)] = [srcFile]
            else:
                srcFileDic[int(year)].append(srcFile)

    sortedSrcFileDic = dict(sorted(srcFileDic.items()))

    outTiffFileLst = []
    for year in sortedSrcFileDic.keys():
        yearTiffLst=sortedSrcFileDic[year]
        sumArr=numpy.zeros((31,56), dtype = float, order = 'C')
        geoTransform = None
        projection = None

        for tiff in yearTiffLst:
            yearmonth = re.search(r'(\d+).HDF.tif.clip.tif', os.path.basename(tiff)).group(1)
            month = int(yearmonth[4:])
            allHours = 24 * calendar.monthrange(year, month)[1]
            raster=tiff

            inraster = gdal.Open(raster)  # 读取路径中的栅格数据
            inband = inraster.GetRasterBand(1)  # 这个波段就是最后想要转为矢量的波段，如果是单波段数据的话那就都是1

            rawArr = inband.ReadAsArray()# 读取栅格数据的投影信息，用来为后面生成的矢量做准备
            arr=numpy.where(rawArr>=-1,rawArr,0.000)
            sumArr=sumArr+arr*allHours
            if not geoTransform:
                geoTransform = inraster.GetGeoTransform()
            if not projection:
                projection = inraster.GetProjectionRef()
        # sumArr = sumArr/(len(yearTiffLst))
        pathlib.Path(outDir).mkdir(parents=True,exist_ok=True)
        outYearTiff = os.path.join(outDir, "{}.yealySumPrecipitation.HDF.tif.clip.tif".format(year))

        rows, cols = sumArr.shape

        #开始写tiff
        driver = gdal.GetDriverByName("GTiff")
        out_ds = driver.Create(outYearTiff, cols, rows, 1, gdal.GDT_Float32)
        out_ds.SetProjection(projection)
        out_ds.SetGeoTransform(geoTransform)
        band = out_ds.GetRasterBand(1)
        band.WriteArray(sumArr)
        band.FlushCache()
        band.ComputeStatistics(False)
        out_ds=None
        outTiffFileLst.append(outYearTiff)





        # dataset = driver.Create(outYearTiff, cols, rows, 1,    gdal.GDT_Float64)
        # dataset.SetGeoTransform(geoTransform)
        # dataset.SetProjection(projection)
        #
        #
        # band = dataset.GetRasterBand(1)  # GetRasterBand is not zero indexed
        # band.WriteArray(sumArr)  # Numpy is zero indexed
        # band.FlushCache()
        #
        # del dataset


    qhShp=r'C:\Users\crawl\Desktop\workForWP\wangpeng-master-thesis\reorganizedRawData\qinghaiShpBoundary\Qinghai.shp'

    with fiona.open(qhShp, "r") as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]
    for unmaskedTiff in outTiffFileLst:
        finalTiff = re.sub(r"yealySumPrecipitation.HDF.tif.clip.tif", "yealySumPrecipitation.qhShpMasked.tif", unmaskedTiff)
        with rasterio.open(unmaskedTiff) as src:
            out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
            out_meta = src.meta
            out_meta.update({"driver": "GTiff",
                             "height": out_image.shape[1],
                             "width": out_image.shape[2],
                             "transform": out_transform})

            with rasterio.open(finalTiff, "w", **out_meta) as dest:
                dest.write(out_image)



