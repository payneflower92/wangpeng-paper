# -*- coding:utf-8 -*-
# author:jinlinfang
# datetime:2023/2/18/0018 22:35
# software: PyCharm

import pathlib

from osgeo import gdal



def openTRMMdata(fpath):
    pass
    # D:\wangpengMasterThesis\下载的相关文件\TRMM_Qinghai\TRMM\trmm_hdf





if __name__ == '__main__':
    fpath = "D:/workForWP/wangpeng-master-thesis/reorganizedRawData/TRMM/2001/200101.HDF"
    hdfObj = gdal.Open(fpath)
    metadata1 = hdfObj.GetMetadata()

    # from arcgis.gis import GIS

