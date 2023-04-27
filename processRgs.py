import datetime
import json
import os.path
import re,pathlib
import shutil
from osgeo import gdal

try:
   from osgeo import gdal
   from osgeo import ogr
except Exception:
   import gdal
   import ogr
import pandas

rootDir="D:\\wangpeng\\2000-2019\\qh_rgs\\"
outDir="D:\\workForWP\\wangpeng-master-thesis\\reorganizedRawData\\RGSdata"
pathlib.Path(outDir).mkdir(parents=True,exist_ok=True)
rgsFileLst=[os.path.join(rootDir, e) for e in os.listdir(rootDir) if re.search(r'qh_\d+.shp$', e)]

allStationDic=dict()
for rgsFile in  rgsFileLst:
    year = re.search(r'(\d{4})\d{2}\.shp$',os.path.basename(rgsFile)).group(1)
    outName=re.search(r'(\S+)\.shp',os.path.basename(rgsFile)).group(1) +".txt"

    gdal.SetConfigOption("GDAL_FILENAME_IS_UTF8","NO")
    gdal.SetConfigOption("SHAPE_ENCODING","")


    ogr.RegisterAll()
    ds = ogr.Open(rgsFile, 0)
    ds.GetLayerCount()
    ds.GetMetadata()
    ds.GetLayer(0)
    ly = ds.GetLayer(0)

    stationDic=dict()
    for f in ly:
        keys = f.keys()
        currentKey = ""
        for k in keys:
            v = f.GetField(k)
            if k=="stationNum":
                currentKey=v
                stationDic[v] = dict()
                if not currentKey in allStationDic.keys():
                    allStationDic[currentKey]={year:[]}
                else:
                    if not year in allStationDic[v].keys():
                        allStationDic[currentKey][year] = []
                continue
            else:
                stationDic[currentKey][k] = v
                if k=="rgs":
                    allStationDic[currentKey][year].append(v)
                elif k not in  allStationDic[currentKey].keys() and k not in ["YEAR","MONTH"]:
                    allStationDic[currentKey][k] = v


    pandas.DataFrame.from_dict(stationDic)
    df = pandas.DataFrame.from_dict(stationDic)
    # df.transpose()
    df = df.transpose()
    print( df.shape[0])
    # df.to_csv(os.path.join(outDir,outName),sep="\t")



# print(allStationDic)

#增加从处理后的文件中集合rgs数据的步骤

processMonthlyDataRoot = "D:\\workForWP\\wangpeng-master-thesis\\reorganizedRawData\\RGSData-afterExceptionProcess"
del2StationDir="D:\\workForWP\\wangpeng-master-thesis\\reorganizedRawData\\RGSData-afterExceptionProcessAndStationDeletion"
dataFileLst=[os.path.join(processMonthlyDataRoot, e) for e in os.listdir(processMonthlyDataRoot) if re.search(r'qh_\d+.txt$', e)]
allStationDic=dict()
i=1
for f in dataFileLst:
    m=re.search(r'qh_(\d{4})\d{2}.txt$',os.path.basename(f))
    year=m.group(1)
    del2stationFile=os.path.join(del2StationDir,os.path.basename(f))
    df=pandas.read_csv(f,sep="\t",index_col=0)
    df=df.loc[~df.index.isin([52854,52941])]
    df.to_csv(del2stationFile, sep="\t")

    for index , record in df.iterrows():
        if index not in allStationDic.keys():
            allStationDic[index] = {year:[]}
        elif year not in allStationDic[index].keys():
            allStationDic[index][year] = []

        for k in record.keys():
            if k == "rgs":
                allStationDic[index][year].append(record[k])
            elif k not in  allStationDic[index].keys() and k not in ["YEAR","MONTH"]:
                allStationDic[index][k] = record[k]
filteredYearStationDic=dict()
for sno in allStationDic.keys():
    if sno not in filteredYearStationDic:
        filteredYearStationDic[sno] = dict()
    for k in allStationDic[sno]:
        if re.search(r'\d{4}',k):
            if len(allStationDic[sno][k]) ==12:
                filteredYearStationDic[sno][k]=sum(allStationDic[sno][k])
df2=pandas.DataFrame.from_dict(filteredYearStationDic).transpose()
df2.to_csv(os.path.join(del2StationDir,"yearrgs.txt"),sep="\t")


