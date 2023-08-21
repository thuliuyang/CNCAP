# -*- coding: utf-8 -*-
#
# Copyright © 2023 Department of Earth System Science Tsinghua University
#
########################################################################
#
# Version is 1.0
#
# 用于计算未来PM2.5污染暴露相关过早死亡人数。
#   --输入数据包含未来人口分布，人口结构，PM2.5浓度，基准死亡率
#   --输出包括全国、分省及重点区域过早死亡人数
#
########################################################################



import pandas as pd
import numpy as np
import datetime
import yaml
import os
import scipy as sp
from scipy import stats
import time
from netCDF4 import Dataset
import netCDF4 as nc
import bisect

np.seterr(all='ignore') #break when warning occur



#-----------------------------------------------基础函数------------------------------------------------------------------

def mean_confidence_interval(data, confidence=0.95):
    # 计算数据的均值和标准误差
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), stats.sem(a)
    # 计算置信区间的上下界
    h = se * sp.stats.t.ppf((1+confidence)/2., n-1)
    return m, m-h, m+h


def causes_allage_mort(idir_IERPar,tmrel_data,cause,cause_file,conc1):
    """
    计算与年龄段无关疾病的过早死亡

    Parameters
    ----------
    idir_IERPar : 参数路径
    tmrel_data : 暴露响应曲线数据
    cause : 疾病类型
    cause_file : 疾病文件
    conc1 : 污染暴露浓度

    Returns
    -------
    af_pm_mean : 归因系数的均值
    af_pm_lb : 归因系数的下界
    af_pm_hb : 归因系数的上界

    """
    pm_1d = np.reshape(conc1,(-1,1))
    pm_data = pd.DataFrame(pm_1d)
    
    #读取相对风险曲线数据
    rr_curve = pd.read_csv(idir_IERPar+cause_file+'.csv',index_col=0)
    rr_curve.index = rr_curve.index.map(lambda x: float("%.2g" % x))

    tmrel_rr = tmrel_data.copy()
    for ii in range(1000):
        row_tmrel = bisect.bisect_left(rr_curve.index, tmrel_data.iloc[ii,][0]) 
        tmrel_rr.iloc[ii,0] = rr_curve.iloc[row_tmrel,ii]
        
    colnames = rr_curve.columns[0:1000]
    
    #处理PM2.5数据，计算相对风险
    pm_data.index = pm_data[0].map(lambda x: float("%.2g" % x))

    pm_data['sort'] = range(264000)
    pm_data_rr = pd.merge(pm_data,rr_curve,how = "left",left_index=True,right_index=True)
    pm_data_rr.sort_values(by = 'sort',inplace=True)
    
    pm_data_rr.loc[:,colnames] = 1-1/(pm_data_rr.loc[:,colnames] / tmrel_rr.T.values)
    pm_data_rr = pm_data_rr.loc[:,colnames]

    pm_data_rr[pm_data_rr<0]=0

    pm_data_rr['mean'] = np.mean(pm_data_rr.loc[:,colnames],1)
    pm_data_rr['lb'] = np.percentile(pm_data_rr.loc[:,colnames],5,1)
    pm_data_rr['hb'] = np.percentile(pm_data_rr.loc[:,colnames],95,1)
    
    #计算归因系数
    af_pm_mean = np.reshape(pm_data_rr['mean'].values,(400,660))
    af_pm_lb  = np.reshape(pm_data_rr['lb'].values,(400,660))
    af_pm_hb  = np.reshape(pm_data_rr['hb'].values,(400,660))
    
    return af_pm_mean,af_pm_lb,af_pm_hb 

def causes_agespecific_mort(idir_IERPar,tmrel_data,cause,cause_file,conc1,age):
    """
    计算与年龄段相关疾病的过早死亡

    Parameters
    ----------
    idir_IERPar : 参数路径
    tmrel_data : 暴露响应曲线数据
    cause : 疾病类型
    cause_file : 疾病文件
    conc1 : 污染暴露浓度
    age : 年龄段

    Returns
    -------
    af_pm_mean : 归因系数的均值
    af_pm_lb : 归因系数的下界
    af_pm_hb : 归因系数的上界

    """
    pm_1d = np.reshape(conc1,(-1,1))
    pm_data = pd.DataFrame(pm_1d)
    
    #读取相对风险曲线数据
    rr_curve = pd.read_csv(idir_IERPar+cause_file+str(age)+'.csv',index_col=0)
    rr_curve.index = rr_curve.index.map(lambda x: float("%.2g" % x))

    tmrel_rr = tmrel_data.copy()
    for ii in range(1000):
        row_tmrel = bisect.bisect_left(rr_curve.index, tmrel_data.iloc[ii,][0]) 
        tmrel_rr.iloc[ii,0] = rr_curve.iloc[row_tmrel,ii]
        
    colnames = rr_curve.columns[0:1000]
    
    #处理PM2.5数据，计算相对风险
    pm_data.index = pm_data[0].map(lambda x: float("%.2g" % x))

    pm_data['sort'] = range(264000)
    pm_data_rr = pd.merge(pm_data,rr_curve,how = "left",left_index=True,right_index=True)
    pm_data_rr.sort_values(by = 'sort',inplace=True)
    
    pm_data_rr.loc[:,colnames] = 1-1/(pm_data_rr.loc[:,colnames] / tmrel_rr.T.values)
    pm_data_rr = pm_data_rr.loc[:,colnames]

    pm_data_rr[pm_data_rr<0]=0

    pm_data_rr['mean'] = np.mean(pm_data_rr.loc[:,colnames],1)
    pm_data_rr['lb'] = np.percentile(pm_data_rr.loc[:,colnames],5,1)
    pm_data_rr['hb'] = np.percentile(pm_data_rr.loc[:,colnames],95,1)
    
    #计算归因系数
    af_pm_mean = np.reshape(pm_data_rr['mean'].values,(400,660))
    af_pm_lb  = np.reshape(pm_data_rr['lb'].values,(400,660))
    af_pm_hb  = np.reshape(pm_data_rr['hb'].values,(400,660))
    
    return af_pm_mean,af_pm_lb,af_pm_hb       



def disarrange(a, axis=-1):
    """
    在给定的轴上原地打乱a的顺序。
    参数：
        a: 要打乱顺序的数组
        axis: 要打乱顺序的轴，默认为最后一个轴
    返回：
        无返回值，但会原地修改数组`a`的顺序
    """
    b = a.swapaxes(axis, -1)
    # 将指定轴和最后一个轴进行交换，得到交换后的数组`b`
    # `b`是`a`的一个视图，所以`a`也会被原地打乱。
    
    shp = b.shape[:-1]
    # 获取除了最后一个轴之外的所有轴的形状
    
    for ndx in np.ndindex(shp):
        # 遍历除了最后一个轴之外的所有轴的索引
        # 通过np.ndindex生成一个迭代器
    
        np.random.shuffle(b[ndx])
        # 在`b`的当前索引位置上打乱顺序
    
    return

#-----------------------------------------------输入输出数据设置------------------------------------------------------------------
#基础参数路径
idir_bmr = './Input Data/1.Parameter/'
idir_IERPar = './Input Data/1.Parameter/IER_2019/MRBRT/mrbrt/'
#人口数据路径
idir_dmg = './Input Data/3.population/'

tmrel_draws = './Input Data/1.Parameter/IER_2019/MRBRT/tmrel_draws.csv'
tmrel_data = pd.read_csv(tmrel_draws)
#浓度数据路径
idir_conc = './Input Data/2.PM2.5_con/'
#浓度数据路径
odir = './Output/'

#区域网格数据
#yaml-based provincial code and name mapping file
if_pname = './Input Data/5.LAM2WGS/prov_code_name_map.yml'
with open(if_pname,'rb') as f:
    i_pname =yaml.safe_load(f)['province']
#i_pname = yaml.load(open(if_pname))
itbl_pname = pd.DataFrame.from_records(data = i_pname)
itbl_pname = itbl_pname.set_index(itbl_pname['regionId'])
itbl_pname = itbl_pname[itbl_pname.index<70] #exclude Taiwan
itbl_pname.set_index('regionNameE',inplace=True)
region_maskname = './Input Data/4.china_mask/'
itbl_jjj = Dataset(region_maskname+'China_226_mask_0.1deg.nc')
mask_jjj = itbl_jjj['mask'][:][:]
itbl_yrd = Dataset(region_maskname+'China_YRD_mask_0.1deg.nc')
mask_yrd = itbl_yrd['mask'][:][:]
itbl_fw = Dataset(region_maskname+'China_FW_mask_0.1deg.nc')
mask_fw = itbl_fw['mask'][:][:]
        

#要计算的情景
#PM2.5浓度情景
scenarios = ['reference']
#所使用的人口情景
ssp_sce = ['SSP1']

product = 'Optimal PM2.5'
#计时
stc = time.process_time()


if (os.path.isdir(odir)==False):
    os.makedirs(odir)


#人口网格数据与省级掩膜，netcdf-based pop and china-mask
idir_pop = './Input Data/3.population/' #POP
if_mask = './Input Data/4.china_mask/Prov_Boundary.01deg.nc' #mask

#设置要预测的疾病相关信息
causes_all = ['COPD','LC','LRI','T2_DM','IHD', 'STROKE']

causes_allage = ['COPD','LC','LRI','T2_DM'] #
causes_allage_filename = {'COPD':'resp_copd','LC':'neo_lung','LRI':'lri','T2_DM':'t2_dm'}
causes_agespecific = ['IHD', 'STROKE']
causes_agespecific_name = {'IHD':'cvd_ihd_', 'STROKE':'cvd_stroke_'}

#准备性别-年龄相关参数
ages = list(np.arange(25,85,5))#list(np.arange(0,85,5))# [0] + list(np.arange(25,85,5))
genders = ['female','male']
tcolumns = []
agecolumns = []
for ai in range(len(ages)):
    age_str = '%02d'  % ages[ai] 
    agecolumns.append(age_str)
for gender in genders:
        for ai in range(len(ages)):
            age_str = '%02d'  % ages[ai] 
            tcolumns.append(gender+'.'+age_str) 

tcolumns_LRI = []
agecolumns_LRI = ['00']
ages_LRI = [0]
for gender in genders:
        for ai in range(len(ages_LRI)):
            age_str = '%02d'  % ages_LRI[ai] 
            tcolumns_LRI.append(gender+'.'+age_str) 
agecolumns_all =  agecolumns


#运行年份设置

years = np.array([2050])

#核算输出范围
percencase = ['Mean','LB','HB']

#未来基准死亡率预测方法
future = ['WHO'] 

#-----------------------------------循环采用不同的基准死亡率预测方法预测未来PM2.5污染相关过早死亡------------------------------------------------------------------
for proj_mor in future:
   ssp_i = 0
   for scenario in scenarios:
       print(scenario)
       

       if (os.path.isdir(odir)==False):
           os.makedirs(odir)
           
       #读取人口数据
       is_pop = {}
       for year in years:
           if_pop = idir_pop + ssp_sce[ssp_i] +'_China_population_'+str(year)+'_01deg.nc'  
           id_pop = Dataset(if_pop)
           i_pop = id_pop['pop'][:][:]
           is_pop[year] = i_pop.copy()
       
       #读取网格数据
       id_mask = Dataset(if_mask)
       i_mask = id_mask['mask'][:][:]
       i_mask[i_mask==71] = 0 #exclude Taiwan
       id_prov = np.unique(i_mask[:])
       id_prov = id_prov[id_prov > 0]

       ttag = '%4d' % year  
       #=============读取基准死亡率=================
       itbl_bmrs = {}    
       for year in years:
           syr = str(year)
           if_bmr = idir_bmr + 'Cause-specific_Deaths_Rate_China_GBD2019.xls'
           itbl_bmrs[year] = pd.read_excel(if_bmr,str(year)+'_'+proj_mor,index_col = 0) 
       
       #==========读取人口年龄结构============

       if_popstru = idir_dmg + 'POPstru_Ratio_China_SSP1.csv' 
       print(if_pop)
       print(if_popstru)

       itbl_popstru  = pd.read_csv(if_popstru,index_col=0)
       
       #==========读取输出数据的标准模板==============
       
       if_pwpm = './Input Data/' + 'output_demo.csv'
       itbl_pwpm = pd.read_csv(if_pwpm,index_col=0)
       itbl_pwpm['regionId'] = itbl_pname['regionId']
       
       #init tables for output statistics
       ##national
       cn_stat = {}
       for cause in (['TOTAL']+causes_all):
           cn_stat[cause] = pd.DataFrame(columns=years,index=percencase)
       
       cn_prov = {}
       for cause in (['TOTAL'] + causes_all):
           tcols = ['regionName']
           for year in years:
               for pcase in percencase:
                   tcols.append(cause+'_'+str(year)+'_'+pcase)
                   
               cn_prov_tmp = pd.DataFrame(index = itbl_pwpm.index,columns = tcols)
           
               cn_prov_tmp['regionName'] =itbl_pwpm.index
               cn_prov_tmp.set_index(['regionName'],inplace = True)
               cn_prov[cause] = cn_prov_tmp.copy()    
               
       cn_region = {}
       for cause in (['TOTAL'] + causes_all):
           tcols = ['regionName']
           for year in years:
               for pcase in percencase:
                   tcols.append(cause+'_'+str(year)+'_'+pcase)
               cn_region_tmp = pd.DataFrame(index = ['JJJ','YRD','FW'],columns = tcols)
               cn_region_tmp['regionName'] = ['JJJ','YRD','FW']
               cn_region_tmp.set_index(['regionName'],inplace = True)
               cn_region[cause] = cn_region_tmp.copy()  

       
       #循环进行污染相关死亡核算
       for year in years:
           
           syr = str(year)
           if_pm = idir_conc+ scenario + '_' + syr +'_exposure_01deg.nc'  
           id_pm = Dataset(if_pm)
           i_pm = id_pm['PM25'][:][:]  
           i_pm[np.isnan(i_pm)] = 0
           i_pm[i_mask == 0 ] = 0

           
           itbl_bmr = itbl_bmrs[year]
           i_pop = is_pop[year]
           i_pop.data[i_pop.data<0] = 0
           i_pop.data[i_mask == 0 ] = 0
           prov_num = int((np.shape(itbl_pwpm))[0])
           
       
       #******************************生成数据存储模板***********************************************
           stc = time.process_time()
           tot_array_mean = np.zeros([prov_num,4+len(agecolumns_all)*2]) # 32 provinces, 38 endpoints (COPD+LC+LRI+DM+(IHD+STROKE)*17),
           tot_array_lb = np.zeros([prov_num,4+len(agecolumns_all)*2]) # 32 provinces, 39 endpoints (COPD+LC+LRI+(IHD+STROKE+NCD_LRI)*12),
           tot_array_hb = np.zeros([prov_num,4+len(agecolumns_all)*2]) # 32 provinces, 39 endpoints (COPD+LC+LRI+(IHD+STROKE+NCD_LRI)*12),
               
           region_death_mean = np.zeros([3,4+len(agecolumns_all)*2])     # BTH,YRD,FW
           region_death_lb = np.zeros([3,4+len(agecolumns_all)*2])     # BTH,YRD,FW
           region_death_hb = np.zeros([3,4+len(agecolumns_all)*2])     # BTH,YRD,FW
               
           cause_count = 0
       
       #======计算所有疾病相关健康损失========
       #------不分年龄段的健康损失-------       
           for cause in  causes_allage:
               print('processing ' + syr+ ' ' + cause)
               af_sst = datetime.datetime.now()
               
               #创建空表存放数据
               cause_mort_prov_mean = np.zeros([prov_num]) #to store result for each endpoint
               cause_mort_prov_lb = np.zeros([prov_num])        
               cause_mort_prov_hb = np.zeros([prov_num])    
               
               AF_mean,AF_lb,AF_hb =  causes_allage_mort(idir_IERPar,tmrel_data,cause,causes_allage_filename[cause],i_pm)
               dr_cause = (itbl_bmr.loc[cause,tcolumns]*itbl_popstru.loc[year,tcolumns]).sum()
               
               #核算全国过早死亡
               cause_mort_grid_nation_mean = AF_mean *i_pop.data * dr_cause
               cause_mort_grid_nation_lb = AF_lb *i_pop.data * dr_cause
               cause_mort_grid_nation_hb = AF_hb *i_pop.data * dr_cause
               
               #核算各省过早死亡
               for pid,prov in enumerate(itbl_pwpm.index[0:-1]):
                   ind_arr = np.where((i_mask==itbl_pwpm.loc[prov,'regionId'])& (i_pm>0))
                   
                   cause_mort_prov_mean[pid] = cause_mort_grid_nation_mean[ind_arr].sum()	
                   cause_mort_prov_lb[pid] = cause_mort_grid_nation_lb[ind_arr].sum()            
                   cause_mort_prov_hb[pid] = cause_mort_grid_nation_hb[ind_arr].sum()
                   
                   cause_mort_prov_mean[(prov_num-1)] = cause_mort_prov_mean[0:-1].sum() #sum to national total
                   cause_mort_prov_lb[(prov_num-1)] = cause_mort_prov_lb[0:-1].sum() #sum to national total            
                   cause_mort_prov_hb[(prov_num-1)] = cause_mort_prov_hb[0:-1].sum() #sum to national total            
                   
               ept_count = cause_count
               print(ept_count)
               
               #区域过早死亡核算
               region_death_mean[0,ept_count] = (cause_mort_grid_nation_mean * mask_jjj).sum()
               region_death_mean[1,ept_count] = (cause_mort_grid_nation_mean * mask_yrd).sum()
               region_death_mean[2,ept_count] = (cause_mort_grid_nation_mean * mask_fw).sum()
        
               region_death_lb[0,ept_count] = (cause_mort_grid_nation_lb * mask_jjj).sum()
               region_death_lb[1,ept_count] = (cause_mort_grid_nation_lb * mask_yrd).sum()
               region_death_lb[2,ept_count] = (cause_mort_grid_nation_lb * mask_fw).sum()
       
               region_death_hb[0,ept_count] = (cause_mort_grid_nation_hb * mask_jjj).sum()
               region_death_hb[1,ept_count] = (cause_mort_grid_nation_hb * mask_yrd).sum()
               region_death_hb[2,ept_count] = (cause_mort_grid_nation_hb * mask_fw).sum()
               
               tot_array_mean[:,ept_count] = cause_mort_prov_mean[:]
               tot_array_lb[:,ept_count] = cause_mort_prov_lb[:]        
               tot_array_hb[:,ept_count] = cause_mort_prov_hb[:]        
               af_set = datetime.datetime.now()
               #print (af_set-af_sst).seconds
               cause_count = cause_count+1
       
   
       
       #------分年龄段的疾病健康损失计算-------    

           cause_count = 0
           #分年龄段计算
           for cause in  causes_agespecific:
               print('processing ' + syr+ ' ' + cause)
               cause_mort_prov_age_mean = np.zeros([prov_num,len(agecolumns)])
               cause_mort_prov_age_lb = np.zeros([prov_num,len(agecolumns)])        
               cause_mort_prov_age_hb = np.zeros([prov_num,len(agecolumns)])
               
               cause_mort_grid_age_mean = np.zeros([400,660,len(agecolumns)])
               cause_mort_grid_age_lb = np.zeros([400,660,len(agecolumns)])        
               cause_mort_grid_age_hb = np.zeros([400,660,len(agecolumns)])        
               for ai in range(0,len(agecolumns)):
                   age=agecolumns[ai]
                   AF_mean,AF_lb,AF_hb = causes_agespecific_mort(idir_IERPar,tmrel_data,cause,causes_agespecific_name[cause],i_pm,age)
                   for gender in genders:
                       unit_bmr = itbl_bmr.at[cause,gender+'.'+age]
                       unit_dmg = itbl_popstru.loc[year,gender+'.'+age]
                       cause_mort_grid_age_mean[:,:,ai] = cause_mort_grid_age_mean[:,:,ai]+AF_mean*i_pop.data*unit_bmr*unit_dmg
                       cause_mort_grid_age_lb[:,:,ai] = cause_mort_grid_age_lb[:,:,ai]+AF_lb*i_pop.data*unit_bmr*unit_dmg
                       cause_mort_grid_age_hb[:,:,ai] = cause_mort_grid_age_hb[:,:,ai]+AF_hb*i_pop.data*unit_bmr*unit_dmg
                       
               #分省核算
               for pid,prov in enumerate(itbl_pwpm.index[0:-1]):
                   ind_arr = np.where((i_mask==itbl_pwpm.loc[prov,'regionId'])& (i_pm>0))
                   for ai in range(0,len(agecolumns)):
                       cause_mort_prov_age_mean[pid,ai] = np.sum(cause_mort_grid_age_mean[:,:,ai][ind_arr])
                       cause_mort_prov_age_lb[pid,ai] = np.sum(cause_mort_grid_age_lb[:,:,ai][ind_arr])
                       cause_mort_prov_age_hb[pid,ai] = np.sum(cause_mort_grid_age_hb[:,:,ai][ind_arr])
         
               age_count = 0
               for ai in range(0,len(agecolumns)):
                   cause_mort_prov_age_mean[(prov_num-1),ai] = cause_mort_prov_age_mean[0:-1,ai].sum()#sum to national total
                   cause_mort_prov_age_lb[(prov_num-1),ai] = cause_mort_prov_age_lb[0:-1,ai].sum()#sum to national total
                   cause_mort_prov_age_hb[(prov_num-1),ai] = cause_mort_prov_age_hb[0:-1,ai].sum()#sum to national total
               
                   ept_count = cause_count*len(agecolumns)+age_count+4
                   
                   #全国核算
                   tot_array_mean[:,ept_count] = cause_mort_prov_age_mean[:,ai]
                   tot_array_lb[:,ept_count] = cause_mort_prov_age_lb[:,ai]
                   tot_array_hb[:,ept_count] = cause_mort_prov_age_hb[:,ai]
                   
                   #分区域核算
                   region_death_mean[0,ept_count] = (cause_mort_grid_age_mean[:,:,ai] * mask_jjj).sum()
                   region_death_mean[1,ept_count] = (cause_mort_grid_age_mean[:,:,ai] * mask_yrd).sum()
                   region_death_mean[2,ept_count] = (cause_mort_grid_age_mean[:,:,ai] * mask_fw).sum()
       
                   region_death_lb[0,ept_count] = (cause_mort_grid_age_lb[:,:,ai] * mask_jjj).sum()
                   region_death_lb[1,ept_count] = (cause_mort_grid_age_lb[:,:,ai] * mask_yrd).sum()
                   region_death_lb[2,ept_count] = (cause_mort_grid_age_lb[:,:,ai] * mask_fw).sum()
   #    
                   region_death_hb[0,ept_count] = (cause_mort_grid_age_hb[:,:,ai] * mask_jjj).sum()
                   region_death_hb[1,ept_count] = (cause_mort_grid_age_hb[:,:,ai] * mask_yrd).sum()
                   region_death_hb[2,ept_count] = (cause_mort_grid_age_hb[:,:,ai] * mask_fw).sum()
                   
                   print('ept_count: ',ept_count)
                   age_count =  age_count + 1
               cause_count = cause_count+1          


           #====================汇总数据=======================
           #summary for endpoints
           causes_ind_all = np.array([(0,1,2,3,4,21),(1,2,3,4,21,38)])
           
           for re in range(0,3):
               reg = ['JJJ','YRD','FW'][re]
               for cid in range(0,len(causes_all)):
                   cause = causes_all[cid]
                   tar_cols = [cause + '_'+syr+'_'+i for i in  percencase]
                   
                   region_epts_mean = region_death_mean[re,causes_ind_all[0,cid]:causes_ind_all[1,cid]]
                   region_epts_lb = region_death_lb[re,causes_ind_all[0,cid]:causes_ind_all[1,cid]]
                   region_epts_hb = region_death_hb[re,causes_ind_all[0,cid]:causes_ind_all[1,cid]]
                   
                   cn_region[cause].loc[reg,tar_cols[0]] = np.round(region_epts_mean.sum(axis=0))
                   cn_region[cause].loc[reg,tar_cols[1]] = np.round(region_epts_lb.sum(axis=0))
                   cn_region[cause].loc[reg,tar_cols[2]] = np.round(region_epts_hb.sum(axis=0))
               
               
           for prov in itbl_pwpm.index:
               provind = np.where(itbl_pwpm.index == prov)
               prov_arr_mean = np.squeeze(tot_array_mean[provind,:])
               prov_arr_lb = np.squeeze(tot_array_lb[provind,:])
               prov_arr_hb = np.squeeze(tot_array_hb[provind,:])

               prov_sim_mean = prov_arr_mean.copy() 
               prov_sim_lb = prov_arr_lb.copy() 
               prov_sim_hb = prov_arr_hb.copy() 
               
               prov_sim_shuffle_mean = prov_sim_mean.copy()
               prov_sim_shuffle_lb = prov_sim_lb.copy()
               prov_sim_shuffle_hb = prov_sim_hb.copy()
       
               #disarrange(prov_sim_shuffle,axis=0)
               prov_sim_sum_mean = prov_sim_shuffle_mean[0:38].sum()
               prov_sim_sum_lb = prov_sim_shuffle_lb[0:38].sum()
               prov_sim_sum_hb = prov_sim_shuffle_hb[0:38].sum()  
               
                        
               tar_cols = ['TOTAL_'+syr+'_'+i for i in  percencase]
               
               cn_prov['TOTAL'].loc[prov,tar_cols[0]] = np.round(prov_sim_sum_mean)         
               cn_prov['TOTAL'].loc[prov,tar_cols[1]] = np.round(prov_sim_sum_lb) 
               cn_prov['TOTAL'].loc[prov,tar_cols[2]] = np.round(prov_sim_sum_hb) 
       
               
               for cid in range(0,len(causes_all)):
                   cause = causes_all[cid]
                   
                   tar_cols = [cause + '_'+syr+'_'+i for i in  percencase]
                   
                   prov_epts_mean = prov_sim_shuffle_mean[causes_ind_all[0,cid]:causes_ind_all[1,cid]].copy()
                   prov_epts_lb = prov_sim_shuffle_lb[causes_ind_all[0,cid]:causes_ind_all[1,cid]].copy()
                   prov_epts_hb = prov_sim_shuffle_hb[causes_ind_all[0,cid]:causes_ind_all[1,cid]].copy()
                   #disarrange(prov_epts,axis=1)
                   prov_epts_sum_mean = prov_epts_mean.sum(axis=0)
                   prov_epts_sum_lb = prov_epts_lb.sum(axis=0)
                   prov_epts_sum_hb = prov_epts_hb.sum(axis=0)
                   cn_prov[cause].loc[prov,tar_cols[0]] = np.round(prov_epts_sum_mean)     
                   cn_prov[cause].loc[prov,tar_cols[1]] = np.round(prov_epts_sum_lb)
                   cn_prov[cause].loc[prov,tar_cols[2]] = np.round(prov_epts_sum_hb)


       #输出省级、区域结果
       of_mort_cn = odir +  scenario +'_Mort_Uncertainty_'+ttag+'_prov_agespecific_'+proj_mor+'.xlsx'
       of_mort_cn_re = odir +  scenario +'_Mort_Uncertainty_'+ttag+'_region_agespecific_'+proj_mor+'.xlsx'

       writer = pd.ExcelWriter(of_mort_cn) 
       writer_reg = pd.ExcelWriter(of_mort_cn_re) 
       for cause in (['TOTAL']+causes_all):
           cn_prov[cause].to_excel(writer,cause)
           cn_region[cause].to_excel(writer_reg,cause)
       writer.save() 
       writer_reg.save()
       ssp_i = ssp_i +1
       
       #输出结果
       tcols = []
       of_mort_cn = odir + scenario +'_Mort_Uncertainty_'+ttag+'_China_agespecific_'+proj_mor+'.xlsx'
       
       writer = pd.ExcelWriter(of_mort_cn) 
       
       for pcase in percencase:
           otbl_causespecific = pd.DataFrame(index=causes_all,columns=years)
           
           for cause in causes_all+['TOTAL']:
               for year in years:
                   otbl_causespecific.loc[cause,year] = cn_prov[cause].loc['Total',cause+'_'+str(year)+'_'+pcase]
           
           otbl_causespecific.to_excel(writer,pcase)
       writer.save()


       

       
