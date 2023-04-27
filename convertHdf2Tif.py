import arcpy, re
from arcpy import env
import os.path

sourceDir = "D:\\workForWP\\wangpeng-master-thesis\\reorganizedRawData\\TRMM\\2000"
destDir = "D:\\workForWP\\wangpeng-master-thesis\\reorganizedRawData\\TiffFromTRMM\\2000"
# wgs84Prj = "D:/workForWP/wangpeng-master-thesis/WGS 1984.prj"

# pathlib.Path(destDir).mkdir(parents=True, exist_ok=True)
# if not  os.path.isdir(destDir):
#     os.mkdir(destDir)
arcpy.CheckOutExtension("Spatial")

# wgs84PrjObj = arcpy.SpatialReference()
# wgs84PrjObj.create(wgs84Prj)

env.workspace = sourceDir
env.scratchWorkspace = sourceDir
hdfLst = arcpy.ListRasters("*", "HDF")
for hdf in hdfLst:
    targetName = re.search(r'(\S+)\.HDF$', os.path.basename(hdf)).group(1) + ".tif"
    tempData = arcpy.ExtractSubDataset_management(hdf, os.path.join(destDir, targetName), "1")
    arcpy.Rotate_management(os.path.join(destDir, targetName), os.path.join(destDir, targetName) + ".flip.tif",
                            angle=270)
    # arcpy.DefineProjection_management(in_dataset=os.path.join(destDir, targetName) + ".flip.tif", coor_system=wgs84Prj)
    # '''
    #
    # '''
    arcpy.ProjectRaster_management(in_dataset=os.path.join(destDir, targetName) + ".flip.tif",
                                   out_raster=os.path.join(destDir, targetName) + ".projected.tif",
                                   resampling_type="NEAREST",
                                   Registration_Point=(), vertical=False
                                   )

