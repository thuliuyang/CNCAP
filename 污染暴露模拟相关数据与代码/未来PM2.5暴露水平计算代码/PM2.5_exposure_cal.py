# -*- coding: utf-8 -*-
#
# Copyright © 2023 Department of Earth System Science Tsinghua University
#
########################################################################
#
# Version is 1.0
#
# 用于将空气质量模型模拟得到的未来PM2.5浓度与TAP浓度数据进行同化，并结合人口数据计算不同地区的人口暴露。
#
########################################################################


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os 
from netCDF4 import Dataset
import netCDF4 as nc
import math
from scipy import interpolate
from scipy.interpolate import Rbf,bisplrep,bisplev,griddata
import pylab as pl
import matplotlib as mpl


#-----------------------------------------------------------基准函数------------------------------------------------------------------------------
def pop_weighted_pm_region(PM_data,pop_data,region_mask):
    #用于计算不同区域的PM2.5暴露
    #PM_data：PM2.5浓度网格数据
    #pop_data：人口网格数据
    #region_mask：不同区域的掩膜数据（区域为为1，区域外为0）
    
    pop_weighted_pm = (PM_data * pop_data * region_mask).sum() / (pop_data * region_mask).sum()
    
    return pop_weighted_pm

def read_mask(mask_file):
    #用于读取不同地区的掩膜数据
    #mask_file：掩膜文件路径

    mask_data = Dataset(Prov_file)
    region_mask = (mask_data.variables['mask'][:,:]).data 
    
    return region_mask

def read_conc(conc_file):
    #用于读取浓度数据
    #conc_file：浓度文件路径
    nc_obj_tap = Dataset(conc_file)
    conc = nc_obj_tap.variables['PM2.5'][:,:].data.round(2) 


    return conc
    
def Interp_ass(conc_base,conc_sim,conc_tap,points_cmaq):
    #用于对模拟数据插值同化
    #conc_base_file：基准年模拟的文件
    #conc_sim_file： 未来年模拟的文件
    #conc_tap_file： TAP的浓度文件
    #points_cmaq:模拟的区域经纬度

    # 针对模拟结果进行计算插值，与人口数据匹配
    xnew = np.arange(70,136,0.1) #x
    ynew = np.arange(15,55,0.1)#y
    xx, yy = np.meshgrid(xnew, ynew) 
    conc_interp_base = griddata(points_cmaq, conc_base.reshape(-1), (xx, yy), method='nearest')
    conc_interp_cmaq = griddata(points_cmaq, conc_sim.reshape(-1), (xx, yy), method='nearest')
    #进行数据同化
    #计算比值
    conc_ratio = conc_interp_cmaq/conc_interp_base       
    #选出与tap对应的网格
    mask = conc_tap.copy()
    mask[mask>0] =1
    #只保留tap有数的中国地区
    conc_cmaq_scale_tap_out = conc_tap * conc_ratio * mask 
    
    return conc_cmaq_scale_tap_out    
            

def nc_output(out_dir,outfile,conc_data):
    #用于将同化后的数据输出成nc文件
    #out_dir：输出路径
    #outfile：输出文件
    #conc_data： 输出的浓度数据
        
    of_out = out_dir + outfile
    of_nc_grid = nc.Dataset(of_out, 'w', format='NETCDF3_CLASSIC') #'NETCDF3_CLASSIC'
    lat = of_nc_grid.createDimension('lat', int(400))
    lon = of_nc_grid.createDimension('lon', int(660))
    grid_nc = of_nc_grid.createVariable('PM25','f4',('lat','lon'))
    grid_nc.unit = 'μg/m3'

    #输出文件
    grid_data = conc_data.copy()
    grid_nc[:,:]= grid_data[:,:]  
    of_nc_grid.by = 'CNCAP Team'    

    return None 
           
        

#-----------------------------------------------------------输入输出数据信息------------------------------------------------------------------------------
#需要运行的情景名称，样例中以ref情景为例
scenario_list = ['ref']
#需要进行运行的年份
year_list = [2050]

#不同地区的掩膜网格路径
BTH_file = './Input data/china_mask/China_226_mask_0.1deg.nc'     #京津冀2+26城市地区，区域内为1，区域外为0
FW_file = './Input data/china_mask/China_FW_mask_0.1deg.nc'       #汾渭地区，区域内为1，区域外为0
YRD_file = './Input data/china_mask/China_YRD_mask_0.1deg.nc'     #长三角地区，区域内为1，区域外为0
Prov_file = './Input data/china_mask/Prov_Boundary.01deg.nc'      #省级地区掩膜，各省不同的代码
#读取不同地区的掩膜网格数据
BTH_mask= read_mask(BTH_file)
FW_mask = read_mask(FW_file)
YRD_mask = read_mask(YRD_file)
Prov_mask = read_mask(Prov_file)

#模拟的区域信息文件，读取模型输出结果的经纬度
cmaq=Dataset('./Input data/GRIDCRO2D_China_d01.nc')
lon_cmaq=cmaq.variables['LON'][0,0,:,:].data.reshape(-1)
lat_cmaq=cmaq.variables['LAT'][0,0,:,:].data.reshape(-1)
points_cmaq = np.array(list(zip(lon_cmaq,lat_cmaq)))

#输出文件路径

#输出PM2.5污染暴露浓度数据
PM_conc_filepath = './Output/PM25.csv'
out_dir = "./Output/"
        

#------------------------------------------循环计算不同情景与不同年份的PM2.5污染暴露-----------------------------------------------------------------

for scenario in scenario_list:
    for year in year_list:
    
        #读取人口网格数据，以SSP1情景下的2050年为例
        if year == 2050:
            pop_file = r'./Input data/Population/SSP1_China_population_2050_01deg.nc' 
            pop_data_sr = Dataset(pop_file)
            pop_data = (pop_data_sr.variables['pop'][:,:]).data 
            pop_data[pop_data<0] =0
        
        #读取基准年模拟结果数据
        conc_base = read_conc('./Input data/Simulation/2020/Yearlymean.d01.nc')
        conc_cmaq = read_conc('./Input data/Simulation/'+scenario+'/'+str(year)+'/Yearlymean.d01.nc')
        

        
        #读取TAP数据
        nc_obj_tap = Dataset('./Input data/TAP/Annualmean_base_TAP_PM25_2020_01deg.nc')
        conc_tap = nc_obj_tap.variables['PM25'][:,:].data.round(2)
        
        #模拟数据同化
        conc_cmaq_scale_tap_out = Interp_ass(conc_base,conc_cmaq,conc_tap,points_cmaq)

        #输出同化后的文件
        outfile = scenario+"_"+ str(year) + "_" + "TAP_scale_cmaq.01x01.nc"
        nc_output(out_dir,outfile,conc_cmaq_scale_tap_out)
        
        
        #计算不同区域的污染暴露
        #新建一个浓度表
        data_conc_out = pd.DataFrame(np.zeros([4,4]),columns = ['year','scenario','region','conc']) 
        data_conc_out.loc[:,'year'] = year
        data_conc_out.loc[:,'scenario'] = scenario
        data_conc_out.loc[0,'region'] = 'China'
        data_conc_out.loc[0,'conc'] = (conc_cmaq_scale_tap_out * pop_data).sum() / pop_data.sum() 
        data_conc_out.loc[1,'region'] = 'BTH'
        data_conc_out.loc[1,'conc'] = pop_weighted_pm_region(conc_cmaq_scale_tap_out,pop_data,BTH_mask)
        data_conc_out.loc[2,'region'] = 'FW'
        data_conc_out.loc[2,'conc'] = pop_weighted_pm_region(conc_cmaq_scale_tap_out,pop_data,FW_mask)               
        data_conc_out.loc[3,'region'] = 'YRD'
        data_conc_out.loc[3,'conc'] = pop_weighted_pm_region(conc_cmaq_scale_tap_out,pop_data,YRD_mask)
              
        #把浓度数据存放到一个表中
        if (scenario == scenario_list[0])&(year == year_list[0]):
            data_out_all = data_conc_out.copy()
        else:
            data_out_all = data_out_all.append(data_conc_out)

#输出文件
data_out_all.to_csv(PM_conc_filepath)


