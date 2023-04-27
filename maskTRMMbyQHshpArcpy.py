#author CUIT_RS191_TCY
import arcpy
import os
import glob

arcpy.CheckOutExtension('Spatial')
arcpy.gp.overwriteOutput = 1

#arcpy.env.workspace = 'G:/Ice_phenology/2002_LST/test/'
input_file = 'G:/Test/'
output_file = 'G:/Test/Clip/'
shp_file = 'G:/Test/Test.shp'
prefix = '.tif'
tail = os.path.splitext(shp_file)[0][36:]#这里可自行修改

if not os.path.exists(output_file):
    os.mkdir(output_file)
#rasters = arcpy.ListRasters('*',".tiff")
file_all = os.listdir(input_file)

for file_i in file_all:
    if file_i.endswith(prefix):
        file_name = input_file + file_i
        output_name = output_file + os.path.splitext(file_i)[0] + tail + '_Clip.tif'
        arcpy.gp.ExtractByMask_sa(file_name,shp_file,output_name)
        print('The clip of {} has completed!'.format(file_i))
print('The project has completed! Congratulation!!!')

for file_i in glob.glob(os.path.join(output_file,'*.xml')):
    os.remove(file_i)
for file_i in glob.glob(os.path.join(output_file,'*.tfw')):
    os.remove(file_i)
