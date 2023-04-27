# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.




# Press the green button in the gutter to run the script.
import os


def processEveryMonthData(rootDir, outDir):

    convertTRMM(rootDir, outDir)
    convertNDVI(rootDir, outDir)
    convertDEM(rootDir, outDir)
    convertRGS(rootDir, outDir)

    applyAtpk(outDir)
    extractInfoForRGS(outDir)
    applyGWR(outDir)


def processEveryYearData(rootDir, outDir,monthLst):
    mergeYearTRMMdata()
    mergeYearNDVIdata()
    mergeRGSdata(rootDir,monthLst)


    pass


if __name__ == '__main__':
    years=range(2000,2020)
    months=range(1,13)
    rawDataRoot=""
    resultDataRoot=""

    for year in years:
        rootDirCurrentYear = os.path.join(rawDataRoot, str(year) )
        outDirCurrentYear = os.path.join(resultDataRoot,"年份水平上", str(year))
        print(rootDirCurrentYear)
        print(outDirCurrentYear)
        processEveryYearData(rootDirCurrentYear, outDirCurrentYear,months)


        for month in months:
            rootDirCurrentTime=os.path.join(rawDataRoot,str(year),str(month))
            outDirCurrentTime=os.path.join(resultDataRoot,"月份水平上",str(year),str(month))
            print(rootDirCurrentTime)
            print(outDirCurrentTime)
            processEveryMonthData(rootDirCurrentTime,outDirCurrentTime)