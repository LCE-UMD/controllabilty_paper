from os.path import dirname, basename, join, exists
import pandas as pd
import numpy as np
import argparse
import nilearn as nil
from nilearn.masking import apply_mask, unmask

def combinePos(ROI,VOX):
    '''
    Combines the two levels (ROI and VOX) of random effect posteriors.
    
    Parameters
    ----------
    ROI: ROI level random effect posteriors (pd.DataFrame)
    VOX: VOX (voxel) level random effect posteriors (pd.DataFrame)
    
    Returns
    -------
    posteriors: combined (ROI + VOX) random effects posteriors (pd.DataFrame)
    '''
    
    posteriors = pd.DataFrame()
    for roi in ROI.columns:
        roi_pos = ROI[roi].values
        vox_pos = VOX[[col for col in VOX.columns if roi in col]].values
        
        tmp_df = pd.DataFrame(roi_pos.reshape(-1,1) + vox_pos,
                              columns=[col for col in VOX.columns if roi in col])
        posteriors = pd.concat([posteriors,tmp_df],axis=1)
    return posteriors

def add_mean(posteriors,fixef,effect):
    '''
    Combines fixed and random effect posteriors
    
    Parameters
    ----------
    posteriors: combined random effect posteriors (pd.DataFrame)
    fixef: fixed effect posteriors (pd.DataFrame)
    effect: Intercept, TRAITmean, TRAITdiff, STATEmean, STATEdiff, or BPdiff (string)
    
    Return
    ------
    pd.DataFrame: (number of iterations x VOX)
    '''
    return posteriors.add(fixef[effect],axis=0)

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

def main(posterior_dir, args):
    # Laod fixed effect posteriors
    fixef = pd.read_csv(join(posterior_dir,'POP.txt'),sep='\t')
    # Load random effect posteriors
    ROI_Intercept = pd.read_csv(join(posterior_dir,'ROI_Intercept.txt'),sep='\t')
    VOX_Intercept = pd.read_csv(join(posterior_dir,'VOX_Intercept.txt'),sep='\t')
    
    posteriors = combinePos(ROI_Intercept,VOX_Intercept)
    posteriors = add_mean(posteriors,fixef,'Intercept')
    
    # Load the insula mask
    mask = nil.image.load_img(args.mask)
    roi_indx = np.unique(mask.get_data())[1:]
    print('ROI indexes: ',roi_indx)
    
    # P+ of every voxel is computed and rendered on a brain template for visualization
    probMask = nil.image.new_img_like(mask,np.zeros_like(mask.get_data()))
    for ii, n in enumerate(roi_indx):
        n_vox = (mask.get_data() == n).sum()
        roi_mask = getSubMask(n,mask)
        print('ROI index:',n)
        print('Num of voxels in data table:',n_vox)
        print('Number of vox in the roi mask:',(roi_mask.get_data() == 1).sum())
        p_plus = []
        for jj, m in enumerate(range(n_vox)):
            vox = posteriors['roi{:02d}_{:03d}'.format(int(n),m)]
            p_plus.append(np.mean(vox > 0))

        rendered = unmask(p_plus,roi_mask)
        probMask = nil.image.math_img("img1 + img2",img1=probMask,img2=rendered)
    
    return probMask

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'combine the voxelwise fixed and random effect posteriors and compute the posterior probabilities (`P+`) and render the onto the MNI brain:')
    
    parser.add_argument('-m','--mask',default='data/masks/right_insula_10ROIs.nii.gz',
                      type=str,help='path/to/ROI-mask')
    parser.add_argument('-o','--output',type=str,
                      default='results/voxelwise',help='path/to/output/folder')
    parser.add_argument('--overwrite',help='overwirte: 0 (default) or 1')
    
    args = parser.parse_args()
    
    # Mask name without extension
    maskname = basename(args.mask).split('.')[0]
    posterior_dir = join(args.output,'uncon_v_con_'+maskname)
    
    if exists(join(posterior_dir,maskname+'_Pmap.nii.gz')) and not args.overwrite:
        raise OSError(join(posterior_dir,maskname+'_Pmap.nii.gz')+": File exists. Use '-overwrite 1' to overwrite")
    else:
        probMask = main(posterior_dir, args)
        probMask.to_filename(join(posterior_dir,maskname+'_Pmap.nii.gz'))