# -*- coding: utf-8 -*-
"""
Created on Sun Jun 12 19:34:33 2016

@author: YIXUAN
"""

#calculate PM2.5 attributable mortality based on GEMM NCD-LRI model
#generate parameter table of GEMM at the beginning to ensure each province applied same parameters. 

#yzheng,2018/10/05

#wrl,2020/04/01

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
import matplotlib.pyplot as plt

np.seterr(all='ignore') #break when warning occur

#calculate distribution of AF reduction during 2013-2017
#conc1 and conc is the arrays storing pm2.5 for all counties in 2013 and 2017

def cal_theta(theta,se_theta):
    return np.random.normal(loc=theta,scale=se_theta,size=1000)

def causes_agespecific_mort(itbl_all, theta_sample, conc1, age, threshold = 2.4):
    
    #GEMM(z)=exp{θlog(z/α+1)/(1+exp{-(z-μ)/ν})}, where z=max(0, PM2.5-2.4μg/m3)
    
    
    #cause = 'COPD'        
    itbl_causes = itbl_all[itbl_all['age'] == age]
    
    #theta = float(itbl_causes['theta'])  #calculation based on nparray will be faster
    #se_theta = float(itbl_causes['se_theta'])
    alpha = float(itbl_causes['alpha'])
    mu = float(itbl_causes['mu'])
    nu = float(itbl_causes['nu'])
    
    #generate 1000 samples of theta
    #theta_sample = np.random.normal(loc=theta,scale=se_theta,size=1000)
    
    def f(x):
        #print x
        x=np.max([0,x-threshold])
        return np.exp(theta_sample*np.log(x/alpha+1)/(1+np.exp(-(x-mu)/nu)))
    
#    def f_se_hb(x):
#        #print x
#        x=np.max([0,x-threshold])
#        return np.exp((theta+1.96*se_theta)*np.log(x/alpha+1)/(1+np.exp(-(x-mu)/nu)))
#    def f_se_lb(x):
#        #print x
#        x=np.max([0,x-threshold])
#        return np.exp((theta-1.96*se_theta)*np.log(x/alpha+1)/(1+np.exp(-(x-mu)/nu)))
    #st = datetime.datetime.now()
    
    conc_arr_rsp1 = np.reshape(conc1,[1,np.shape(conc1)[0]])
    #conc_arr_rsp2 = np.reshape(conc2,[1,np.shape(conc2)[0]])
    
    rr = np.apply_along_axis(f, 0, conc_arr_rsp1) #ip_pm)
    #rr_hb = np.apply_along_axis(f_se_hb, 0, conc_arr_rsp1) #ip_pm)
    #rr_lb = np.apply_along_axis(f_se_lb, 0, conc_arr_rsp1) #ip_pm)    
    
    af = (rr-1)/rr
    #af_hb = (rr_hb-1)/rr_hb      
    #af_lb = (rr_lb-1)/rr_lb      
         
    return af #,af_hb,af_lb #mean,lb,hb #,mean_,lb_,hb_ 

def CIstat(x):
    return np.mean(x), np.percentile(x, 2.5),np.percentile(x, 97.5) 
    #USE INTERPOLATION WILL ALTER THE PERCENTILE RESULT AT SOME EXTENT
    #return np.mean(x), np.percentile(x, 2.5,interpolation='nearest'),np.percentile(x, 97.5,interpolation='nearest') 
def disarrange(a, axis=-1):
    """
    Shuffle `a` in-place along the given axis.

    Apply numpy.random.shuffle to the given axis of `a`.
    Each one-dimensional slice is shuffled independently.
    """
    b = a.swapaxes(axis, -1)
    # Shuffle `b` in-place along the last axis.  `b` is a view of `a`,
    # so `a` is shuffled in place, too.
    shp = b.shape[:-1]
    for ndx in np.ndindex(shp):
        np.random.shuffle(b[ndx])
    return


styr = 2017
edyr = 2017
years = np.arange(styr,edyr+1)
dmg_year = 2017
gbd_year = 2017
sgbd_year = str(gbd_year)

#idir_GEMMPar = '/Users/YIXUAN/Research/Models/5.GEMM/3.GEMM_Parameters/'
idir_bmr = '/Users/ruiliwu/Desktop/Research/2018/GEMM/'

idir_dmg = idir_bmr

idir_conc = '/Users/ruiliwu/Desktop/Research/2018/GEMM/CONC/'
odir = '/Users/ruiliwu/Desktop/Research/2018/GEMM/Mort_NCD-LRI/'

#netcdf-based pop and china-mask
idir_pop = '/Users/ruiliwu/Desktop/Research/2018/GEMM/GPW_POP/' #POP
if_mask = '/Users/ruiliwu/Desktop/Research/2018/GEMM/GPW_POP/Prov_Boundary.01deg.noTW.nc' #mask

if_pop = idir_pop + 'GPW_POP_01deg_2017_predicted_china.nc'
id_pop = Dataset(if_pop)
i_pop = id_pop['pop'][:][:]
    
id_mask = Dataset(if_mask)
i_mask = id_mask['mask'][:][:]
i_mask[i_mask==71] = 0 #exclude Taiwan

id_prov = np.unique(i_mask[:])
id_prov = id_prov[id_prov > 0]

ttag = '%4d' % styr + '-' + '%4d' % edyr 
causes_agespecific = ['NCD-LRI']
causes_all = causes_agespecific
genders = ['Female','Male']
percencase = ['Mean','LB','HB']
ages = list(np.arange(25,85,5))# [0] + list(np.arange(25,85,5))

tcolumns = []
agecolumns = []
for ai in range(len(ages)):
    age_str = '%02d'  % ages[ai] 
    agecolumns.append(age_str)
for gender in genders:
        for ai in range(len(ages)):
            age_str = '%02d'  % ages[ai] 
            tcolumns.append(gender+'.'+age_str) 

agecolumns_all = agecolumns


if (os.path.isdir(odir)==False):
    os.makedirs(odir)

#yaml-based provincial code and name mapping file
if_pname = '/Users/ruiliwu/Desktop/Research/2018/GEMM//prov_code_name_map.yml'
i_pname = yaml.load(open(if_pname))['province']
itbl_pname = pd.DataFrame.from_records(data = i_pname)
itbl_pname = itbl_pname.set_index(itbl_pname['regionId'])
itbl_pname = itbl_pname[itbl_pname.index<70] #exclude Taiwan

itbl_pname.set_index('regionNameE',inplace=True)

#=============read bmr tables=================

#if_bmr = idir_bmr + 'Cause-specific_Deaths_Rate_China_GBD' + sgbd_year + '_1990-' + sgbd_year + '_byYears.xls'
if_bmr = idir_bmr + 'NCD-LRI_Deaths_Rate_China_GBD' + sgbd_year + '_1990-' + sgbd_year + '_byYears.xls'
itbl_bmr = pd.read_excel(if_bmr,str(dmg_year),index_col = 0) 

#==========read national demographic table=========
#----------
if_pop = idir_dmg + 'POPstru_China_GBD' + sgbd_year + '_1990-' + sgbd_year + '.csv' #need total population for scaling with the 2010 pop
if_popstru = idir_dmg + 'POPstru_Ratio_China_GBD' + sgbd_year + '_1990-' + sgbd_year + '.csv' #need total population for scaling with the 2010 pop

itbl_pop  = pd.read_csv(if_pop,index_col=0)
itbl_popstru  = pd.read_csv(if_popstru,index_col=0)


#==========read PWPM2.5 concentrations as a template to generate output statistics============

#total population in 2010 is also included

if_pwpm = '/Users/ruiliwu/Desktop/Research/2018/GEMM/Annual_PopWeighted_PMExposure_baseNnoCtrl_2002-2012.csv'
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


#==========read gemm parameters =============
if_IERPar = odir + '/GEMM_CRCurve_parameters_1000samples_NCD-LRI.xlsx' #25-80 for IHD and STROKE, 
#if_IERPar = '/Users/ruiliwu/Desktop/Research/2018/GEMM/GEMM_CRCurve_parameters.xlsx'
itbls_all = {}
itbls_theta = {}
for cause in causes_all:
    
    itbl_all = pd.read_excel(if_IERPar,cause)
    itbls_all[cause] = itbl_all.copy()
    itbls_all[cause+'_theta'] = pd.read_excel(if_IERPar,cause+'_theta')
#    itbls_all[cause] = pd.read_excel(if_IERPar,cause)
    
    
#==============================================main start====================================

#casename='SR_ground_mar_conc_pSO4'

#casenames=['inmap_reduce50','inmap_reduce70','inmap_reduce90','cmaq_reduce70','cmaq_reduce90']
casenames = ['']
#main start
for year in years:
    year = year
    syr = str(year)

    if_pm = idir_conc+ 'wgs_mean.2017.01deg.'+casename+'.nc'
    id_pm = Dataset(if_pm)
    i_pm = id_pm['pSO4'][:][:]  
    i_pm[i_mask == 0 ] = 0
       
#*************************************************************************************************

    #tot_array for one year
    tot_array = np.zeros([264000,len(agecolumns_all)*len(causes_agespecific)]) #1000 simulations, 32 provinces, 12 endpoints (NCD-LRI)*12),
    cause_count = 0

#======calculate Mort for causes_agespecific========
#------Age-specific BMR-------     
    #calculate RR for causes_agespecific
    for cause in  causes_agespecific:
        itbl_all = itbls_all[cause]
        itbl_all_theta = itbls_all[cause + '_theta']
        print('processing ' + syr+ ' ' + cause)
   
# store mort array -wrl
# calculate the national by grid 
    
        ind_arr = np.where(i_mask>=0)
        conc_arr = i_pm[ind_arr]
        pop_arr = i_pop[ind_arr]  

        cause_mort_grid_age3 = np.zeros([1000,264000,len(agecolumns_all)*len(causes_agespecific)])

        age_count = 0           
        for ai in range(0,len(agecolumns)):
            age=agecolumns[ai]
            AF = causes_agespecific_mort(itbl_all,itbl_all_theta[float(age)].values,conc_arr,float(age))
            cause_mort_grid_age = np.zeros([1000,len(pop_arr)])
            cause_mort_grid_age2 = np.zeros([len(pop_arr)])
            for gender in genders:
            #save time by 1/4, yz, 2017/05/07
                unit_bmr = itbl_bmr.at[cause,gender+'.'+age]
                unit_dmg = itbl_popstru.loc[dmg_year,gender+'.'+age]
                for sim in range(0,1000):
                    cause_mort_grid_age[sim,:] = cause_mort_grid_age[sim,:]+AF[sim,:]*pop_arr*unit_bmr*unit_dmg

            cause_mort_grid_age3[:,:,ai] =cause_mort_grid_age
        
            ept_count = cause_count*len(agecolumns)+age_count
            tot_array[:,ept_count] = cause_mort_grid_age3[:,:,ai].sum(axis=0)
#            print 'ept_count: ',ept_count
            age_count =  age_count + 1
        cause_count = cause_count+1 
        
        prov_arr = np.squeeze(tot_array[:,:])
        prov_sim = prov_arr.copy()
        
        prov_sim_shuffle = prov_sim.copy()
        disarrange(prov_sim_shuffle,axis=0)
        prov_sim_sum = prov_sim_shuffle.sum(axis=1)
                
        tar_cols = ['TOTAL_'+syr+'_'+i for i in  percencase]
        
        cn_prov['TOTAL'].loc[prov,tar_cols] = np.round(np.apply_along_axis(CIstat, 0, prov_sim_sum))
        
        causes_ind_all = np.array([(0),(12)])
        for cid in range(0,len(causes_all)):
            cause = causes_all[cid]
            
            tar_cols = [cause + '_'+syr+'_'+i for i in  percencase]
            
            prov_epts = prov_sim_shuffle[:,causes_ind_all[0]:causes_ind_all[1]].copy()
            disarrange(prov_epts,axis=1)
            prov_epts_sum = prov_epts.sum(axis=1)
            cn_prov[cause].loc[prov,tar_cols] = np.round(np.apply_along_axis(CIstat, 0, prov_epts_sum))      
            
# output data
        
        of_PM = odir+'Grid_Mort_Array_'+ casename+syr +str(dmg_year)+'_NCD-LRI.nc'        
        if os.path.isfile(of_PM):
            os.remove(of_PM)
        of_nc = nc.Dataset(of_PM, 'w', format='NETCDF4')
    
        ncols = of_nc.createDimension('ncols', 660)
        nrows = of_nc.createDimension('nrows', 400)    

        MORT = of_nc.createVariable('MORT','f4',('nrows','ncols'))
        MORT.unit = 'people'
        MORT[:] = np.reshape(cn_prov,[400,660])

# write global attribute
        
        of_nc.by = 'Ruili Wu'
        of_nc.creattime = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        of_nc.causes = '(NCD-LRI) *12 '
        of_nc.ages = agecolumns
        of_nc.close()
        print('Good girl!')