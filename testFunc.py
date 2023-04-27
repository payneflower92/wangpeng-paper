import re
import unittest

from combineRGS_TRMM_NDVI_DEM import *


class MyTestCase(unittest.TestCase):
    def testOpenDEM(self):
        # root=r'C:\Users\crawl\Desktop\workForWP\wangpeng-master-thesis\reorganizedRawData\RGSData-afterExceptionProcessAndStationDeletion'
        # info = collectRGSinfo(root)
        # print(len(info))
        slope = r'C:\Users\crawl\Desktop\workForWP\wangpeng-master-thesis\reorganizedRawData\DEM_1KM_raw\qh_1km_slope.tif'
        ele = r'C:\Users\crawl\Desktop\workForWP\wangpeng-master-thesis\reorganizedRawData\DEM_1KM_raw\qh_1km_elevation.tif'
        aspect = r'C:\Users\crawl\Desktop\workForWP\wangpeng-master-thesis\reorganizedRawData\DEM_1KM_raw\qh_1km_aspect.tif'
        parseDemData(slope,ele,aspect)


if __name__ == '__main__':
    unittest.main()
