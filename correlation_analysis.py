# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 10:09:09 2021

@author: scalvinib
"""
#Functions for correlation analysis between CT parameters and folding rate


import numpy as np
import scipy.stats
import scipy as sy
import matplotlib.pyplot as plt
from xlwt import Workbook
import seaborn as sns
import pandas as pd
#=============================================================================

def line(x, a , b):
    y= a*x+b
    return y

def correlate (var1, var2):
    var1=np.asarray(var1.astype(np.float16))
    var2=np.asarray(var2.astype(np.float16))
    corr=scipy.stats.pearsonr(var1,var2)
    return corr
 

def lin_fit(x,y):    
    guess1= 0.0
    guess2=0.0
    guess=[guess1, guess2]

    errfunc2 = lambda p, Series_two, y: ( line(x, *p) - y)**2
    optim, success= sy.optimize.leastsq(errfunc2, guess[:], args=( x, y))
    fit=line(x, optim[0] , optim[1])
    return fit

def fractions_fit(data_two, data_multi):
    color_graph= 'b'
    color_graph_2= 'orange'
    #sns.set_style('darkgrid')
    series_corr=correlate(data_two["Series"],data_two['Folding rate'])
    series_pvalue_two="%.3f" % round(series_corr[1], 3)
    series_corr_two="%.2f" % round(series_corr[0], 2)
    parallel_corr=correlate(data_two["Parallel"],data_two['Folding rate'])
    parallel_pvalue_two="%.3f" % round(parallel_corr[1], 3)
    parallel_corr_two="%.2f" % round(parallel_corr[0], 2)
    cross_corr=correlate(data_two["Cross"],data_two['Folding rate'])
    cross_pvalue_two="%.3f" % round(cross_corr[1], 3)  
    cross_corr_two="%.2f" % round(cross_corr[0], 2)    

    series_corr=correlate(data_multi["Series"],data_multi['Folding rate'])
    series_pvalue_multi="%.3f" % round(series_corr[1], 3)
    series_corr_multi="%.2f" % round(series_corr[0], 2)
    parallel_corr=correlate(data_multi["Parallel"],data_multi['Folding rate'])
    parallel_pvalue_multi="%.3f" % round(parallel_corr[1], 3)
    parallel_corr_multi="%.2f" % round(parallel_corr[0], 2)
    cross_corr=correlate(data_multi["Cross"],data_multi['Folding rate'])
    cross_pvalue_multi="%.3f" % round(cross_corr[1], 3)
    cross_corr_multi="%.2f" % round(cross_corr[0], 2)

    fig, ax = plt.subplots(nrows=1, ncols=3)
    fig.set_figheight(4)
    fig.set_figwidth(15)
    ax[0].set_title('Series', fontsize=20)
    ax[0].scatter(data_two["Series"],data_two['Folding rate'],label= 'corr_two= {}, p= {}'.format(series_corr_two, series_pvalue_two), color= color_graph)
    ax[0].scatter(data_multi["Series"],data_multi['Folding rate'],label= 'corr_multi= {}, p={}'.format(series_corr_multi, series_pvalue_multi), color=color_graph_2)
    fit_Series_multi=lin_fit(data_multi["Series"],data_multi['Folding rate'])
    fit_Series_two=lin_fit(data_two["Series"],data_two['Folding rate'])
    ax[0].legend()
    ax[0].plot(data_multi["Series"],fit_Series_multi, color=color_graph_2)
    ax[0].plot(data_two["Series"],fit_Series_two,color= color_graph)
    ax[0].set_xlabel('Series (%)')
    ax[0].set_ylabel('Folding rate (ln kf)')

    
    ax[1].set_title('Parallel', fontsize=20)
    ax[1].scatter(data_two["Parallel"],data_two['Folding rate'], label= 'corr_two= {}, p= {}'.format(parallel_corr_two, parallel_pvalue_two), color= color_graph)
    ax[1].scatter(data_multi["Parallel"],data_multi['Folding rate'], label= 'corr_multi= {}, p= {}'.format(parallel_corr_multi, parallel_pvalue_multi), color= color_graph_2)
    fit_Parallel_multi=lin_fit(data_multi["Parallel"],data_multi['Folding rate'])
    fit_Parallel_two=lin_fit(data_two["Parallel"],data_two['Folding rate'])
    ax[1].legend()
    ax[1].plot(data_two["Parallel"],fit_Parallel_two, color= color_graph)
    ax[1].plot(data_multi["Parallel"],fit_Parallel_multi,color= color_graph_2)
    ax[1].set_xlabel('Parallel (%)')
  
    ax[2].set_title('Cross',fontsize=20)
    ax[2].scatter(data_two["Cross"],data_two['Folding rate'], label= 'corr_two= {}, p= {}'.format(cross_corr_two, cross_pvalue_two), color= color_graph)
    ax[2].scatter(data_multi["Cross"],data_multi['Folding rate'], label= 'corr_multi= {}, p= {}'.format(cross_corr_multi, cross_pvalue_multi), color= color_graph_2)
    fit_Cross_multi=lin_fit(data_multi["Cross"],data_multi['Folding rate'])
    fit_Cross_two=lin_fit(data_two["Cross"],data_two['Folding rate'])
    ax[2].legend()  
    ax[2].plot(data_two["Cross"],fit_Cross_two,color= color_graph)
    ax[2].plot(data_multi["Cross"],fit_Cross_multi,color= color_graph_2)
    ax[2].set_xlabel('Cross (%)')
    
    correlation_vec_two=[series_corr_two, series_pvalue_two, parallel_corr_two, parallel_pvalue_two, cross_corr_two, cross_pvalue_two]  
    correlation_vec_multi=[series_corr_multi, series_pvalue_multi, parallel_corr_multi,parallel_pvalue_multi, cross_corr_multi, cross_pvalue_multi] 
    return correlation_vec_two,correlation_vec_multi
        
def remove_pdbs(database, pdb_strings):  
    n_pdbs=len(pdb_strings)
    for t in range(n_pdbs):
        
        database=database[database['PDB Id']!=pdb_strings[t]]
                                          
    N_contacts = "N contacts" in database
    if N_contacts:
        
        database_zeros=database[database['N contacts']==0]
        if (database_zeros.empty == False):
            print('ERROR: NUMBER OF CONTACTS=0')
            print(database_zeros) 
            print('number of proteins with zero contacts:')
            print(len(database_zeros))
        database= database[database['N contacts']!=0]     
    else:
        print('No contacts')
    return database   


def save_results(corr_two,corr_multi,filename):
    wb = Workbook()   
    sheet1 = wb.add_sheet('Sheet 1') 
    sheet1.write(0,0,'series corr')
    sheet1.write(0,1,'series pvalue')
    sheet1.write(0,2,'parallel corr')
    sheet1.write(0,3, 'parallel pvalue')
    sheet1.write(0,4,'cross corr')
    sheet1.write(0,5,'cross pvalue')
    sheet1.write(0,6,'folder')

    for t in range(len(corr_two)+1):
        if(t == len(corr_two)):
            sheet1.write(1,t,'two')
            sheet1.write(2,t,'multi')
        else:
            sheet1.write(1,t,corr_two[t])
            sheet1.write(2,t,corr_multi[t])
    wb.save(filename)
    
def normalize(x):
        return (x - np.min(x)) / np.ptp(x)   
    
    
def set_layout(SMALL_SIZE=15,MEDIUM_SIZE=18,BIGGER_SIZE=21 ):

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)
    
def corr_matrix(higher, ave, lower):
    higher_matrix=higher.to_numpy()
    ave_matrix=ave.to_numpy()
    lower_matrix=lower.to_numpy()
    matrix=[lower_matrix,ave_matrix,higher_matrix]
    #print(matrix[1])
    
    dim_row=lower_matrix.shape[0]
    dim_col=lower_matrix.shape[1]-1
    corr_matrix=np.zeros((dim_row*3,dim_col))
    corr_coeff=np.zeros((dim_col, int(dim_col/2)))
    pvalue=np.zeros((dim_col, int(dim_col/2)))
    
    for t in range(3):
        corr_matrix[t*dim_row:dim_row*(t+1), 0:dim_col]=np.copy(matrix[t][0:dim_row, 0:dim_col])
    for t in range(3):    
        corr_coeff[:,t]=corr_matrix[:,t*2]
        pvalue[:,t]=corr_matrix[:,(t*2)+1]
  
    pvalue=1-pvalue
    pvalue[pvalue>=0.95]=1
    pvalue[pvalue<0.95]=0.0
    corr_corrected=corr_coeff*pvalue
    return corr_corrected


def correlation_map_database(lower_lim, upper_lim, data_two, data_multi, corr_quantity):
    data_multi_lowerCO=data_multi[data_multi['Contact Order']<lower_lim]
    data_two_lowerCO=data_two[data_two['Contact Order']<lower_lim]

    data_multi_ave=data_multi[(data_multi['Contact Order']>lower_lim) &(data_multi['Contact Order']<upper_lim)]
    data_two_ave=data_two[(data_two['Contact Order']>lower_lim) &(data_two['Contact Order']<upper_lim)]

    data_multi_higherCO=data_multi[data_multi['Contact Order']>upper_lim]
    data_two_higherCO=data_two[data_two['Contact Order']>upper_lim]
    
    list_db=[data_two_lowerCO,data_multi_lowerCO,data_two_ave,data_multi_ave,data_two_higherCO,data_multi_higherCO]
    
    matrix=np.zeros((len(list_db),1))
    pvalue=np.zeros((len(list_db),1))
    
    for t in range(len(list_db)):        
        CO_corr=correlate(list_db[t][corr_quantity],list_db[t]['ln kf'])
        pvalue[t]="%.2f" % round(CO_corr[1], 2)
        matrix[t]="%.2f" % round(CO_corr[0], 2)
        
    pvalue=1-pvalue
    pvalue[pvalue>=0.95]=1
    pvalue[pvalue<0.95]=0.0
    corr_corrected=matrix*pvalue
    
    return corr_corrected


def barplot_fractions(data, upper_lim, lower_lim):
    data_lowerCO=data[data['Contact Order']<lower_lim]
    data_ave=data[(data['Contact Order']>lower_lim) &(data['Contact Order']<upper_lim)]
    data_higherCO=data[data['Contact Order']>upper_lim]
    
    CO_df=[data_lowerCO,data_ave,data_higherCO]
    new_df=[]
    for t in range(len(CO_df)):
        series=np.array(CO_df[t]['Series'])
        parallel=np.array(CO_df[t]['Parallel'])
        cross=np.array(CO_df[t]['Cross'])

        tot_len=len(series)+len(parallel)+len(cross)
        nu_top_vec=np.zeros(tot_len)
        nu_top_vec[0:len(series)]=np.copy(series)
        nu_top_vec[len(series):len(series)+len(parallel)]=np.copy(parallel)
        nu_top_vec[len(series)+len(parallel):tot_len]=np.copy(cross)
        vec_fraction=['']*tot_len
        vec_fraction[0:len(series)]=np.copy(['Series']*len(series))
        vec_fraction[len(series):len(series)+len(parallel)]=np.copy((['Parallel']*len(parallel)))
        vec_fraction[len(series)+len(parallel):tot_len]=np.copy(['Cross']*len(cross))

        data_new= {'Top_num':nu_top_vec, 'Top_fract':vec_fraction}
        df= pd.DataFrame(data=data_new)
        new_df.append(df)
        
        #Make plot
    ig, ax = plt.subplots(1, 3, figsize=(14,3))
    sns.set(style="whitegrid", color_codes=True)
    pal = sns.color_palette("Set2", 3)
    labels=["Lower CO","Average CO","Upper CO"]
    for t in range(len(new_df)):
        sns.barplot(x="Top_fract", y="Top_num", hue="Top_fract", data=new_df[t], dodge = False,
                        palette=pal, ax=ax[t] )
        ax[t].set_ylim(0.0,0.9)
        ax[t].set_title(labels[t], fontsize=20)
        ax[t].legend(fontsize=14)
        
    return ig
    