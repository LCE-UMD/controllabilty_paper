from sklearn.preprocessing import StandardScaler
import pandas as pd
import numpy as np
import os
from os.path import join, isdir, dirname, basename, exists
from os import mkdir
import glob as glob
import argparse

# Project folder path
eCON = '/data/bswift-1/Pessoa_Lab/eCON'

# Path beta estimate text files for every subject
bucket_path = join(eCON,'dataset/results_ShockUncensored/{subj}/{group}/splitted_regs/emoproxII_ROIs_final_new/noProx/buttonPress/{subj}_bucket_REML_clean.1D')

# load ROI information
rois = pd.read_csv(join(eCON,'ROI_masks/EmoproxII_ROIs_final/readme_new.txt'),
                   sep='\t',index_col='Index')
rois = rois['ROI'].to_dict()
print('Total number of ROIs: ',len(rois.keys()))

# load yoked participant info
yoked = pd.read_excel(join(eCON,'onsetdir/CON_yoked_table.xlsx'))
yoked = yoked.query('use == 1')

def combine_beta(beta,t):
    weights = (t/beta)**2
    combined_beta = np.sum(weights*beta,axis=1)/np.sum(weights,axis=1)
    combined_var = 1/np.sum(weights,axis=1)
    return combined_beta, np.sqrt(combined_var)


# Standardize the four STAI scores
def standardize(df):
    ss = StandardScaler()
    df_norm = pd.DataFrame(ss.fit_transform(df),
                           columns=df.columns,
                           index=df.index)
    return df_norm


def make_df(args):
    main_df = pd.DataFrame()
    for i, row in yoked.iterrows():
        pair = i
        for kind in ['control','uncontrol']:
            subj = row[kind]
            nruns = np.sum(row.loc['run0':'run5'].astype(bool))
            group = kind+'lable'
            data = np.loadtxt(bucket_path.format(subj=subj,group=group))
            shock_est = data[:,25:][:,:nruns*2]
            beta = shock_est[:,::2]
            t = shock_est[:,1::2]
            combined_beta, combined_sd = combine_beta(beta=beta,t=t)

            beta_df = pd.DataFrame(combined_beta,index=list(rois.values())).T
            beta_df['Pair'] = 'P{:02d}'.format(i)
            beta_df['Subject'] = subj
            beta_df['Group'] = kind
            beta_df = pd.melt(beta_df,id_vars=['Pair','Subject','Group'],
                              var_name='ROI',value_name='beta')

            sd_df = pd.DataFrame(combined_sd,index=list(rois.values())).T
            sd_df['Pair'] = 'P{:02d}'.format(i)
            sd_df['Subject'] = subj
            sd_df['Group'] = kind
            sd_df = pd.melt(sd_df,id_vars=['Pair','Subject','Group'],
                            var_name='ROI',value_name='sd')

            tmp_df = pd.merge(beta_df,sd_df)

            main_df = pd.concat([main_df,tmp_df],axis=0,ignore_index=True)

    main_df['cond'] = main_df['Group'].apply(lambda kind: 0.5 if kind=='uncontrol' else -0.5)

    # load button_presses.txt 
    button_press = pd.read_csv(args.output+'/button_presses.txt',sep='\t',
                               index_col=['Pair','Subject','Group'])
    # load STAI_scores.txt
    STAI = pd.read_csv(args.output+'/STAI_scores.txt',sep='\t',
                       index_col=['Subject','Group'])

    behavior_df = pd.merge(button_press.reset_index(),STAI.reset_index())

    # Merge ROI responses and behavioral dfs into single df
    df = pd.merge(main_df,behavior_df)
    df.set_index(['Pair','ROI'],inplace=True)

    control = df[df['Group']=='control'][['beta','buttPress','TRAIT','STATE']]
    uncontrol = df[df['Group']=='uncontrol'][['beta','buttPress','TRAIT','STATE']]

    final_df = uncontrol.subtract(control)
    final_df.rename(columns={'TRAIT':'TRAITdiff','STATE':'STATEdiff'},inplace=True)
    final_df['TRAITmean'] = uncontrol.add(control)['TRAIT']/2 
    final_df['STATEmean'] = uncontrol.add(control)['STATE']/2

    final_df = final_df[['beta']].join(standardize(final_df[['TRAITmean','TRAITdiff',
                                                             'STATEmean','STATEdiff',
                                                             'buttPress']]))
    return final_df

if __name__ == '__main__':
    '''
    This script extracts every participants stressor responses from 24 ROIs, standardizes
    covariates and returns a dataframe ready for BML analysis.
    '''
    
    parser = argparse.ArgumentParser(description='Extracts every participants stressor responses from 24 ROIs, standardizes covariates and returns a dataframe ready for BML analysis.')
    
    parser.add_argument('-o','--output',type=str,
                        help='output path to save the output df')
    parser.add_argument('--overwrite',type=int,default=0,
                        help='overwrite the existing file? (0 or 1)')
    
    args = parser.parse_args()
    
    if exists(args.output+'/uncon_v_con_stressor.txt') and not args.overwrite:
        raise OSError(args.output+"/uncon_v_con_stressor.txt: File exists. Use '-overwrite 1' to overwrite")
    else:
        os.makedirs(args.output,exist_ok=True)
        df = make_df(args)
        df.to_csv(args.output+'/uncon_v_con_stressor.txt',sep='\t',
                  index=True,float_format='%.6f')
    