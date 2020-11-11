from os.path import join, exists, basename
import os
import pandas as pd
import numpy as np
import nilearn as nil
from nilearn.masking import apply_mask, unmask
import argparse

# Define path to the shock beta map of every subject
eCON = '/data/bswift-1/Pessoa_Lab/eCON'
voxelwise_path = join(eCON,'dataset/results_ShockUncensored/{subj}/{group}lable/splitted_regs/shock_analysis/noProx/buttonPress/{subj}_shock_beta.nii.gz')

# load the csv file that contains the list of yoked-participaths
yoked = pd.read_excel(join(eCON,'onsetdir/CON_yoked_table.xlsx'))
yoked = yoked.query('use == 1')


def getSubMask(indx,mask):
    '''
    Takes in the mask and ROI index, and outputs a binary (0s and 1s) 
    mask with only that ROI.
    '''
    mask_idx = np.where(mask.get_data() == indx)
    sub_mask = np.zeros_like(mask.get_data())
    sub_mask[mask_idx] = 1
    sub_mask_img = nil.image.new_img_like(mask,sub_mask)
    return sub_mask_img

def main(args):
    mask = nil.image.load_img(args.mask)
    # Total number of ROIs
    print('Total number of ROIs: ',len(np.unique(mask.get_data())[1:]))
    
    df = pd.DataFrame(columns=['Subj','ROI','VOX','Y'])
    for indx in np.unique(mask.get_data())[1:]:
        submask = getSubMask(indx,mask)
        print('roi{:02d}'.format(int(indx)))
        print('Total number of voxels: ',(submask.get_data() == 1).sum())
        #plot_glass_brain(submask)
        for ii,row in yoked.iterrows():
            uncon = nil.image.load_img(voxelwise_path.format(subj=row['uncontrol'],
                                                             group='uncontrol'))
            con = nil.image.load_img(voxelwise_path.format(subj=row['control'],
                                                             group='control'))

            uncon_beta = apply_mask(uncon,submask)
            con_beta = apply_mask(con,submask)
            Y = uncon_beta - con_beta
            tmp_df = pd.DataFrame(Y,columns=['Y'],
                                  index=['{name}_{k:03d}'.format(name='roi{:02d}'.format(int(indx)),k=k) for k in range(Y.shape[0])])
            tmp_df.index.name = 'VOX'
            tmp_df.reset_index(inplace=True)
            tmp_df['ROI'] = 'roi{:02d}'.format(int(indx))
            tmp_df['Subj'] = 'P{:02d}'.format(ii)
            df = pd.concat([df,tmp_df],ignore_index=True)
        
    # Load the covariate file
    cov_list= ["Pair","TRAITmean","TRAITdiff","STATEmean","STATEdiff","BPdiff"]
    cov = pd.read_csv('data/ROIwise/uncon_v_con_stressor.txt',
                      sep='\t')[cov_list].drop_duplicates()
    cov.reset_index(drop=True,inplace=True)
    
    # merge the covariate columns to the dataframe.
    df = df.merge(cov,left_on='Subj',right_on='Pair')
    
    return df
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description = "Extract's voxelwise stressor response for the given ROI mask")
    
    parser.add_argument('-m','--mask',default='data/mask/right_insula_10ROIs.nii.gz',
                      type=str,help='path/to/ROI-mask')
    parser.add_argument('-o','--output',type=str,
                      default='data/voxelwise',help='path/to/output/folder')
    parser.add_argument('--overwrite',help='overwirte: 0 (default) or 1')
    
    args = parser.parse_args()
    
    maskname = basename(args.mask)
    output_path = join(args.output,'uncon_v_con_{}.txt').format(maskname.split('.')[0])
    
    if exists(output_path) and not args.overwrite:
        raise OSError(output_path+": File exists. Use '-overwrite 1' to overwrite")
    else:
        os.makedirs(args.output,exist_ok=True)
        df = main(args)
        df.to_csv(output_path,sep='\t',index=False,float_format='%.6f')