import pandas as pd
import numpy as np

import os
from os.path import join, exists

from scipy.stats import spearmanr
import math

from sklearn.preprocessing import StandardScaler
import argparse

eCON = '/data/bswift-1/Pessoa_Lab/eCON'
SCR_betapath = join(eCON,'SCR_new/dataset/results_ShockUncensored/{subj}/{group}lable/splitted_regs/noProx/buttonPress/shock_IM/{subj}_bucket_LSS.1D')
roi_betapath = join(eCON,'dataset/results_ShockUncensored/{subj}/{group}lable/splitted_regs/ROI_final/noProx/buttonPress/shock_IM/{subj}_betas_3dLSS.1D') # final (new) rois
yoked = pd.read_excel(join(eCON,'SCR_new/scripts/CON_yoked_table.xlsx'))
yoked = yoked.query('use == 1')

# load ROI information
rois = pd.read_csv(join(eCON,'ROI_masks/EmoproxII_ROIs_final/readme'),sep='\t',index_col='Index')
rois = rois['ROI'].to_dict()

def get_zcorr_df(df,rois):
    rba_df = pd.DataFrame()
    for roi in rois:
        for i, row in yoked.iterrows():
            con_roi = df[(df['Pair']=='P{:02d}'.format(i)) & (df['Group']=='control')][roi].values
            uncon_roi = df[(df['Pair']=='P{:02d}'.format(i)) & (df['Group']=='uncontrol')][roi].values

            con_SCR = df[(df['Pair']=='P{:02d}'.format(i)) & (df['Group']=='control')]['SCR'].values
            uncon_SCR = df[(df['Pair']=='P{:02d}'.format(i)) & (df['Group']=='uncontrol')]['SCR'].values

            con_rp,_ = spearmanr(con_roi,con_SCR)
            uncon_rp,_ = spearmanr(uncon_roi,uncon_SCR)

            con_zp = math.atanh(con_rp)
            uncon_zp = math.atanh(uncon_rp)

            # Y is the difference (uncon-con) of z-transformed correlation coefficient
            tmp_df = pd.DataFrame(['P{:02}'.format(i), roi,
                                   uncon_zp-con_zp,uncon_zp,con_zp],
                                  index=['Pair','ROI','Y','uncontrol','control']).T
            rba_df = pd.concat([rba_df,tmp_df],axis=0,ignore_index=True)

    return rba_df

def make_df():
    # load trail-by-trial SCR and ROI stressor responses
    df = pd.DataFrame()
    for i, row in yoked.iterrows():
        for kind in ['control','uncontrol']:
            scr_beta = np.loadtxt(SCR_betapath.format(subj=row[kind],group=kind))
            roi_beta = np.loadtxt(roi_betapath.format(subj=row[kind],group=kind)).T
            all_betas = np.concatenate((scr_beta[:,None],roi_beta),axis=1)

            tmp_df = pd.DataFrame(all_betas,columns=['SCR']+list(rois.values()))
            tmp_df['SubjID'] = row[kind]
            tmp_df['Group'] = kind
            tmp_df['Pair'] = 'P{:02d}'.format(i)
            df = pd.concat([df,tmp_df],axis =0)


    # Create a contrast (uncon - con) df named paired_df
    paired_df = pd.DataFrame()
    for i,eff in enumerate(list(rois.values())+['SCR'],start=1):
        paired_df['{} Diff'.format(eff)] = df[eff][df['Group']=='uncontrol'] - df[eff][df['Group']=='control']

    #paired_df.rename(columns={'rBNST':'rBNST_diff','SCR':'SCR_diff'},inplace=True)
    df.reset_index(inplace=True)

    rba_df = get_zcorr_df(df,list(rois.values()))
    rba_df['Y'] = rba_df['Y'].astype(float)
    rba_df['uncontrol'] = rba_df['uncontrol'].astype(float)
    rba_df['control'] = rba_df['control'].astype(float)
    
    # load covariates
    bpd = pd.read_csv('DATA/button_presses.txt',sep='\t')
    STAI = pd.read_csv('DATA/STAI_scores.txt',sep='\t')

    # combine covariates
    cov_df = pd.merge(bpd,STAI,how='left')
    cov_df.set_index('Pair',inplace=True)

    # Create contrast df by subtracting control from uncontrol
    control = cov_df[cov_df['Group'] == 'control'][["buttPress","TRAIT","STATE"]]
    uncontrol = cov_df[cov_df['Group'] == 'uncontrol'][["buttPress","TRAIT","STATE"]]

    cov_df = uncontrol.subtract(control)
    cov_df.rename(columns={"TRAIT":"TRAITdiff",
                              "STATE":"STATEdiff",
                              "buttPress":"BPdiff"},
                     inplace=True)
    cov_df['TRAITmean'] = uncontrol.add(control)['TRAIT']/2 
    cov_df['STATEmean'] = uncontrol.add(control)['STATE']/2

    cov_df = cov_df[cov_df.index.isin(rba_df.Pair)]

    # Standardize covariates 
    ss = StandardScaler()
    cov_df_std = pd.DataFrame(ss.fit_transform(cov_df),columns=cov_df.columns,index=cov_df.index)
    cov_df_std.reset_index(inplace=True)
    rba_df = rba_df.merge(cov_df_std)
    
    return rba_df[['Pair','ROI','Y',
                   'uncontrol','control',
                   'TRAITmean','TRAITdiff',
                   'STATEmean','STATEdiff',
                   'BPdiff']]
        
if __name__ == '__main__':
    '''
    This script extracts every participants trial stressor responses from 24 ROIs, standardizes
    covariates and returns a dataframe ready for BML analysis.
    '''
    
    parser = argparse.ArgumentParser(description='Extracts every participants stressor responses from 24 ROIs, standardizes covariates and returns a dataframe ready for BML analysis.')
    
    parser.add_argument('-o','--output',type=str,
                        help='output path to save the output df')
    parser.add_argument('--overwrite',type=int,default=0,
                        help='overwrite the existing file? (0 or 1)')
    
    args = parser.parse_args()
    
    if exists(args.output+'/uncon_v_con_ROI_SCR_zscorr.txt') and not args.overwrite:
        raise OSError(args.output+"/uncon_v_con_ROI_SCR_zscorr.txt: File exists. Use '-overwrite 1' to overwrite")
    else:
        os.makedirs(args.output,exist_ok=True)
        df = make_df()
        df.to_csv(args.output+'/uncon_v_con_ROI_SCR_zscorr.txt',sep='\t',
                  index=False,float_format='%.5f')
    