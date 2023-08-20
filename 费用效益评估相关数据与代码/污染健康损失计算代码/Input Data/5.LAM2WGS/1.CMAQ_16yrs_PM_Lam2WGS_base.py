# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 21:03:35 2016

@author: YIXUAN
"""

#regrid cmaq-derived PM2.5 concentration into 0.1 degree
#for 6 simulations of national-10 sensitivity simulation 

#speeded up as compare to previous approach (nested loop)
#yzheng 2017/10/20

import numpy as np
from netCDF4 import Dataset
import netCDF4 as nc
from pyproj import Proj
import pyproj
import datetime
import os
import json

if_cal_ind = 0

casename = 'fix_mics2010'

idir_cmaq = '/Users/YIXUAN/Research/2017/1.2000-2015/2.Simulation_Related/3.AQ_Sim/3.Output/1.Sim_Result/annualmean/ncep_' + casename + '/'
odir_cmaq = '/Users/YIXUAN/Research/2017/1.2000-2015/2.Simulation_Related/3.AQ_Sim/3.Output/2.Sim_Result_WGS/annualmean/ncep_' + casename + '/'

idir_gridind = '/Users/YIXUAN/Research/phd_dissertation/C4/3.Output/3.Attri_Conc2Measures/1.senSim_WGS/'

species = ['PM25','NO3','SO4','NH4','BC','OC','OTHERPMFINE']

years = np.arange(2005,2006)
#ctrl_year = 2015
#ctrl_syr = '%4d' % ctrl_year


if (os.path.isdir(odir_cmaq) == False):
    os.makedirs(odir_cmaq)

#使用反演的2015年PM2.5浓度作为mask
if_mask = '/Users/YIXUAN/Research/2016/2.BME_MAP/3.output/6.XT_Best/2.2Stage/fit8/Optimal PM2.5 AnnualAvg 2013.nc' 
omask = Dataset(if_mask)['PM25'][:]
mask_ind = np.where(omask.data == -9999)


wgs84=Proj("+init=EPSG:4326")
lcc = Proj(proj='lcc', lat_1=25, lat_2=40,lon_0=110,lat_0=34, ellps='WGS84')
i_CMAQxorig = -3096000 #m
i_CMAQyorig = -2286000
i_CMAQsize = 36000
i_CMAQcols = 172
i_CMAQrows = 127
o_CMAQxorig = 70 #deg
o_CMAQyorig = 15 
o_CMAQsize = 0.1
o_CMAQncols = 660
o_CMAQnrows = 400



### make high resolution x y matrix of CMAQ simulation output



t_CMAQsize = 1000

scaler = i_CMAQsize/t_CMAQsize
t_CAMQcols = scaler*i_CMAQcols
t_CAMQrows = scaler*i_CMAQrows

if if_cal_ind:
    
    t_x_CMAQ = np.linspace(i_CMAQxorig+t_CMAQsize/2,i_CMAQxorig + i_CMAQsize*i_CMAQcols - t_CMAQsize/2,t_CAMQcols)
    t_y_CMAQ = np.linspace(i_CMAQyorig+t_CMAQsize/2,i_CMAQyorig + i_CMAQsize*i_CMAQrows - t_CMAQsize/2,t_CAMQrows)
    
    x_tmp,y_tmp = np.meshgrid(t_x_CMAQ,t_y_CMAQ)
    t_lon_CMAQ,t_lat_CMAQ = pyproj.transform(lcc,wgs84,x_tmp,y_tmp)
    
    def cal_lonind(val):
        return float(int((val-o_CMAQxorig)/o_CMAQsize))
    def cal_latind(val):
        return float(int((val-o_CMAQyorig)/o_CMAQsize))
    
    v_cal_lonind = np.vectorize(cal_lonind)  # or use a different name if you want to keep the original f
    v_cal_latind = np.vectorize(cal_latind)
    
    t_ind_lon = v_cal_lonind(t_lon_CMAQ)  
    t_ind_lat = v_cal_latind(t_lat_CMAQ)
    
    t_ind_lon[t_ind_lon >= o_CMAQncols ] = np.nan
    t_ind_lon[t_ind_lon < 0 ] = np.nan
    
    t_ind_lat[t_ind_lat >= o_CMAQncols ] = np.nan
    t_ind_lat[t_ind_lat < 0 ] = np.nan
             
    #t_ind = t_ind_lat * t_CAMQcols + t_ind_lon
    
    #将索引转换为一维数组
    t_ind_lon_1 = np.reshape(t_ind_lon,t_CAMQrows * t_CAMQcols)
    t_ind_lat_1 = np.reshape(t_ind_lat,t_CAMQrows * t_CAMQcols)
    t_ind_1 = t_ind_lat_1 * o_CMAQncols + t_ind_lon_1 #!!!!!!!!!!o_CMAQncols
    
    #将0.1度目标网格的索引也转换为1维的
    o_ind = np.linspace(0,o_CMAQncols*o_CMAQnrows-1,o_CMAQncols*o_CMAQnrows)
    
    #定义输出字典
    dict_ijind = {}
    for i in o_ind:
        dict_ijind[i] = []
    
    #计算索引;遍历目标矩阵的下标使用np.where函数来查找对应CMAQ数据网格索引的过程太慢了，因此遍历CMAQ数据网格索引，直接赋值给目标队列
    for i in range(0,len(t_ind_1)):
        print str(i) + '/'+str(len(t_ind_1))
        ky = t_ind_1[int(i)]
        if dict_ijind.has_key(ky):
            dict_ijind[ky].append(i)
    
    #去除没有有效索引值的key
    for i in o_ind:
        if len(dict_ijind[i]) == 0:
            dict_ijind.pop(i)
        
    #output dict to json, avoid repeat indexing in the future
    with open(odir_cmaq + '1km_2_01deg_index.json', 'w') as fp:
        json.dump(dict_ijind, fp)
    
    
#read json from output file    
else:
    with open(idir_gridind + '1km_2_01deg_index.json') as data_file:    
        dict_ijind = json.load(data_file)


#@@@@@for control scenarios

def cal_oinds(inds):
    return np.mod(inds,o_CMAQncols),(inds/o_CMAQncols)
def cal_tinds(inds):
    return np.mod(inds,t_CAMQcols),(inds/t_CAMQcols)

for year in years:
    
    syr = str(year)
    if_cmaq = idir_cmaq+'Annualmean.'+syr+'.nc' 
    if os.path.isfile(if_cmaq):
        i_CMAQ = Dataset(if_cmaq)
        o_arr = {} #store result for each species
        for spe in species:
            print year,spe
            i_PM = i_CMAQ[spe][:][:]
            o_PM = np.zeros((o_CMAQnrows,o_CMAQncols))
            o_PM[:] = np.nan
            #o_count = np.zeros((o_CMAQnrows,o_CMAQncols))
            tmp_PM = np.repeat(np.repeat(i_PM,scaler,axis = 0),scaler,axis=1)
            
            #tmp_PM[np.isnan(t_ind_lon)] = np.nan
            #tmp_PM[np.isnan(t_ind_lat)] = np.nan
                  
            for o_ind in dict_ijind.keys():
                inds = dict_ijind[o_ind]
                o_ind = int(float(o_ind))

                inds_x,inds_y = cal_tinds(np.array(inds))
                o_lon,o_lat = cal_oinds(np.array(o_ind))
                #print o_lon,o_lat
                o_PM[o_lat,o_lon] = np.nanmean(tmp_PM[inds_y,inds_x])
            
            #o_PM[mask_ind] = np.nan
            o_arr[spe] = o_PM.copy()
            
        i_CMAQ.close()
            
        of_PM = odir_cmaq+'Annualmean.'+syr+'.01deg.nc'        
        if os.path.isfile(of_PM):
            os.remove(of_PM)
        of_nc = nc.Dataset(of_PM, 'w', format='NETCDF4')
        
        ncols = of_nc.createDimension('ncols', o_CMAQncols)
        nrows = of_nc.createDimension('nrows', o_CMAQnrows)
        
        for spe in species:
            o_spe = of_nc.createVariable(spe,'f4',('nrows','ncols'))
            o_spe.unit = 'ug/m3'
            o_spe[:] = o_arr[spe]
    
        #write global attribute
        of_nc.by = 'Yixuan Zheng'
        of_nc.creattime = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        of_nc.coordinatesys = 'WGS84'
        of_nc.xllcorner = o_CMAQncols
        of_nc.yllcorner = o_CMAQnrows
        of_nc.pixelsize = o_CMAQsize
        of_nc.close()

        

            
            
    
            
