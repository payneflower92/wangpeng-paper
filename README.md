# wangpeng-paper

本项目记录了毕业论文写过的一些代码，主要是对遥感数据处理及相关插值算法实现

主要包括两部分工作：
1、对原始文件处理，包括格式转换（hdf转tif、tif转shp、tif转csv）；
单位修正（mm/hr转mm/month、TRMM\NDVI缺省值修改、固定倍率）；图像裁剪、拼接、投影（tif、shp）；
栅格计算（如dem求取aspect和slope）
2、算法构建，包括Area-to-Point Kriging和GWRK等，ATPK主要参考了Pyinterpolate包、GWRK主要参考了MGWR包，
感谢两包作者的无私奉献；上述两包都能在github上找到
