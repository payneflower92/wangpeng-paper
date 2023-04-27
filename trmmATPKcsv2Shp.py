import multiprocessing
import re, pathlib

import sys, argparse
from datetime import datetime


import pandas


try:
    from osgeo import gdal
    from osgeo import ogr, osr
except Exception:
    import gdal
    import ogr

import pandas as pd
import os,numpy
import fiona
import rasterio
import rasterio.mask
global dataRoot
global workerNo, processNo, taskPattern
workerNo = 1
processNo = 1

def csv2shp(csv_path, shp_path, layerName):
    # 解决中文字符问题
    gdal.SetConfigOption("GDAL_FILENAME_IS_UTF8", "NO")
    gdal.SetConfigOption("SHAPE_ENCODING", "")

    # 设置空间参考,4326代表WGS84
    SpatialReference = osr.SpatialReference()
    SpatialReference.ImportFromEPSG(4326)

    # 新建DataSource,Layer
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.CreateDataSource(shp_path)
    layer = data_source.CreateLayer(layerName, SpatialReference, ogr.wkbPoint)

    # 读取csv文件
    csv_df = pd.read_csv(csv_path, sep="\t", index_col=0)
    # csv所有列名,即shp的字段名
    filed_names = list(csv_df)[6:]

    # layer添加上述字段
    for field_name in filed_names:
        print(str(csv_df[field_name].dtypes))
        if ("int" in str(csv_df[field_name].dtypes)):
            field = ogr.FieldDefn(field_name, ogr.OFTInteger)
            field.SetWidth(32)
        elif ("float" in str(csv_df[field_name].dtypes)):
            field = ogr.FieldDefn(field_name, ogr.OFTReal)
            field.SetWidth(32)
            field.SetPrecision(5)
        else:
            field = ogr.FieldDefn(field_name, ogr.OFTString)
            field.SetWidth(32)
        layer.CreateField(field)

    # 从layer中读取相应的feature类型，并创建feature
    featureDefn = layer.GetLayerDefn()
    feature = ogr.Feature(featureDefn)

    # 设定几何形状
    point = ogr.Geometry(ogr.wkbPoint)

    # 循环输入字段属性值
    for i in range(len(csv_df)):
        for j in range(len(filed_names)):
            if ("int" in str(csv_df[filed_names[j]].dtypes)):
                feature.SetField(filed_names[j], int(csv_df.iloc[i, j + 6]))
            elif ("float" in str(csv_df[filed_names[j]].dtypes)):
                feature.SetField(filed_names[j], float(csv_df.iloc[i, j + 6]))
            else:
                feature.SetField(filed_names[j], str(csv_df.iloc[i, j + 6]))

        # 设置几何信息部分
        # 利用经纬度创建点,X为经度,Y为纬度(我的数据第5列是经度,第6列是纬度)
        point.AddPoint(float(csv_df.iloc[i, 0]), float(csv_df.iloc[i, 1]))
        feature.SetGeometry(point)

        # 将feature写入layer
        layer.CreateFeature(feature)

    # 从内存中清除 ds，将数据写入磁盘中
    data_source.Destroy()


def shapeClip(
        trmmShpFile,
        qhShpFile,
        outShp):
    """
    矢量裁剪
    :param trmmShpFile: 要裁剪的矢量文件
    :param qhShpFile: 掩膜矢量文件
    :param outShp: 裁剪后的矢量文件保存目录
    :return:
    """
    ogr.RegisterAll()
    gdal.SetConfigOption("GDAL_FILENAME_IS_UTF8", "YES")
    # 载入要裁剪的矢量文件
    print("# 载入要裁剪的矢量文件")
    baseData = ogr.Open(trmmShpFile)
    baseLayer = baseData.GetLayer()
    spatial = baseLayer.GetSpatialRef()
    geomType = baseLayer.GetGeomType()
    baseLayerName = baseLayer.GetName()
    # 载入掩膜矢量文件
    print("载入掩膜矢量文件")
    maskData = ogr.Open(qhShpFile)
    maskLayer = maskData.GetLayer()
    maskLayerName = maskLayer.GetName()
    # 生成裁剪后的矢量文件
    outLayerName = maskLayerName + "_Clip_" + baseLayerName

    gdal.SetConfigOption("SHAPE_ENCODING", "GBK")
    driver = ogr.GetDriverByName("ESRI Shapefile")
    outData = driver.CreateDataSource(outShp)
    outLayer = outData.CreateLayer(outLayerName, spatial, geomType)
    print("开始裁剪矢量文件")
    baseLayer.Clip(maskLayer, outLayer,["PROMOTE_TO_MULTI=YES",])
    print("裁剪矢量文件结束")
    outData.Release()
    baseData.Release()
    maskData.Release()


def csv2Tiff(csvFile, tifFile,srcTif):
    ptsDf = pandas.read_csv(csvFile, sep="\t", index_col=0)
    # 打开面tif数据
    dr = gdal.Open(srcTif)
    transform = dr.GetGeoTransform()
    trmm_x_origin = transform[0]
    trmm_y_origin = transform[3]
    pixel_width = transform[1]
    pixel_height = transform[5]
    projection=dr.GetProjectionRef()

    srcArr = None
    try:
        srcArr = dr.ReadAsArray()
    except Exception:
        pass

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
    targetArr=numpy.zeros((2, len(latArr),len(lonArr)),dtype=numpy.float32)

    count=0
    for i in range(len(lonArr)):
        for j in range(len(latArr)):
            targetArr[0,j,i] = ptsDf.iat[count,6] if ptsDf.iat[count,6]>=0 else 0
            # if ptsDf.iat[count,6]<0:
            #     print(ptsDf.iat[count,6])
            targetArr[1, j, i] = ptsDf.iat[count, 7]
            count = count+1

    # # get GeoTiff driver by
    # image_type = 'GTiff'
    # driver = gdal.GetDriverByName(image_type)
    #
    # # passing the filename, x and y direction resolution, no. of bands, new raster.
    # new_raster = driver.Create(tifFile, len(lonArr), len(latArr), 2, gdal.GDT_Float32)
    #
    # # transforms between pixel raster space to projection coordinate space.
    # new_raster.SetGeoTransform((lonArr[0], degree_unit, 0, latArr[0], 0, -degree_unit))
    #
    # new_raster.SetProjection(dr.GetProjectionRef())
    # new_raster.FlushCache()
    #
    # # get required raster band.
    # band = new_raster.GetRasterBand(1)
    # # band = out_ds.GetRasterBand(1)
    # band.WriteArray(targetArr[0,:,:])
    # band.FlushCache()
    # # get required raster band.
    # band2 = new_raster.GetRasterBand(2)
    # # band = out_ds.GetRasterBand(1)
    # band2.WriteArray(targetArr[1,:,:])
    # band2.FlushCache()
    # band2.ComputeStatistics(False)
    # # new_raster.ComputeStatistics(False)
    # new_raster=None

    # 开始写tiff
    driver = gdal.GetDriverByName("GTiff")
    out_ds = driver.Create(tifFile, len(lonArr), len(latArr), 2, gdal.GDT_Float32)
    out_ds.SetProjection(projection)
    geoTransform=(lonArr[0], degree_unit, 0, latArr[0], 0, -degree_unit)
    out_ds.SetGeoTransform(geoTransform)
    band = out_ds.GetRasterBand(1)
    arr=targetArr[0,:,:]
    band.WriteArray(arr)
    band.FlushCache()
    band.ComputeStatistics(False)

    band = out_ds.GetRasterBand(2)
    arr=targetArr[1,:,:]
    band.WriteArray(arr)
    band.FlushCache()
    band.ComputeStatistics(False)


    out_ds = None

    qhShp=os.path.join(dataRoot, 'qinghaiShpBoundary/Qinghai.shp')

    with fiona.open(qhShp, "r") as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]

    finalTiff = re.sub(r"TRMM1KmAfterQHshpMaskAndATPK.smoothedPointMatrixOKatpk.raw.tif", "TRMM1KmAfterQHshpMaskAndATPK.smoothedPointMatrixOKatpk.qhcliped.tif", tifFile)
    with rasterio.open(tifFile) as src:
        out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
        out_meta = src.meta
        out_meta.update({"driver": "GTiff",
                         "height": out_image.shape[1],
                         "width": out_image.shape[2],
                         "transform": out_transform})

        with rasterio.open(finalTiff, "w", **out_meta) as dest:
            dest.write(out_image)

def processParallel(srcFile, unclipShpFile, layerName, unclipTiffFile,egTif):
    start_t = datetime.now()
    print("now begin to process {} ".format(layerName))
    if not os.path.isfile(unclipShpFile):
        csv2shp(srcFile, unclipShpFile, layerName)
    end1_t = datetime.now()
    elapsed_sec = (end1_t - start_t).total_seconds()
    print("转csv到shp 共消耗: " + "{:.2f}".format(elapsed_sec) + " 秒")

    if not os.path.isfile(unclipTiffFile):
        csv2Tiff(srcFile, unclipTiffFile, egTif)
        end2_t = datetime.now()
        elapsed_sec = (end2_t - start_t).total_seconds()
        print("转csv到tif 共消耗: " + "{:.2f}".format(elapsed_sec) + " 秒")
    return "scucess"




if __name__ == '__main__':
    dataRoot = r"C:\Users\crawl\Desktop\workForWP\wangpeng-master-thesis\reorganizedRawData"
    if sys.platform == "linux":
        dataRoot = "/home/cloudam/wpThesis"
    egTif= os.path.join(dataRoot,"TRMMtifMaskedByQH/2000/200001.HDF.tif.clip.tif")
    parallel = False

    parser = argparse.ArgumentParser(description='命令行中传入一个数字')
    parser.add_argument('--rootDir', '-r', type=str, required=False, help='根目录')
    parser.add_argument("--processNo", '-c', type=int,  required=False, help='并行进程数')
    parser.add_argument("-parallel", action='store_true', required=False, help='是否并行')
    parser.add_argument("--workerNo", '-w', type=int, required=False, help='插值单任务线程数')




    args = parser.parse_args()
    if args.rootDir:
        dataRoot = r'' + args.rootDir
    if args.processNo:
        processNo=args.processNo
    if args.workerNo:
        workerNo = args.workerNo


    csvDir = os.path.join(dataRoot, "TRMMtifATPKplusRegularization")
    outDir = os.path.join(dataRoot, "TRMMtifATPKplusRegularizationTransferToShp1")
    qhShpFile= os.path.join(dataRoot, "qinghaiShpBoundary","Qinghai.shp")

    srcFileDic = dict()
    paramsInfoLst = []
    for d in os.listdir(csvDir):
        if not re.search(r'^(yearly|\d+)$', d):
            continue
        for f in os.listdir(os.path.join(csvDir, d)):
            if not re.search(r'TRMM1KmAfterQHshpMaskAndATPK.smoothedPointMatrixOKatpk.csv$', f):
                continue
            srcFile = os.path.join(csvDir, d, f)
            yearmonth = re.search(r'(\S+).TRMM1KmAfterQHshpMaskAndATPK.smoothedPointMatrixOKatpk.csv', f).group(1)
            srcFileDic[yearmonth] = (srcFile, d, f)

    sortedSrcFileDic = dict(sorted(srcFileDic.items(),reverse=True))
    print(sortedSrcFileDic)

    print(dataRoot)

    for yearmonth in sortedSrcFileDic.keys():
        # print("{}：开始转换csv到shp".format(yearmonth))

        (srcFile, d, f) = sortedSrcFileDic[yearmonth]
        yearmonth = str(yearmonth)
        # year = int(str(yearmonth)[0:4])
        # month = int(str(yearmonth)[4:])

        unclipShpName = re.sub(r'TRMM1KmAfterQHshpMaskAndATPK.smoothedPointMatrixOKatpk.csv',
                         'TRMM1KmAfterQHshpMaskAndATPK.smoothedPointMatrixOKatpk.raw.shp', f)
        qhclipShpName = re.sub(r'TRMM1KmAfterQHshpMaskAndATPK.smoothedPointMatrixOKatpk.csv',
                               'TRMM1KmAfterQHshpMaskAndATPK.smoothedPointMatrixOKatpk.qhShpClip.shp', f)
        unclipTiffName = re.sub(r'TRMM1KmAfterQHshpMaskAndATPK.smoothedPointMatrixOKatpk.csv',
                         'TRMM1KmAfterQHshpMaskAndATPK.smoothedPointMatrixOKatpk.raw.tif', f)
        newSubDir = os.path.join(outDir, d)
        unclipShpFile = os.path.join(newSubDir, unclipShpName)
        unclipTiffFile = os.path.join(newSubDir, unclipTiffName)
        qhclipShp = os.path.join(newSubDir,qhclipShpName)
        pathlib.Path(newSubDir).mkdir(parents=True, exist_ok=True)

        paramsInfoLst.append((srcFile, unclipShpFile, "ATPK_" + yearmonth, unclipTiffFile,egTif ))
        processParallel(srcFile, unclipShpFile, "ATPK_" + yearmonth, unclipTiffFile,egTif)
        print()



    # #并行化处理
    # rLst = []
    #
    # pool = multiprocessing.Pool(processNo)
    # processLst = [pool.apply_async(processParallel, args=param) for param in paramsInfoLst]
    # for p in processLst:
    #     r = p.get()
    #     rLst.append(r)
    # print(rLst)





