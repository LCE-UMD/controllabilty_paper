# Controllability over stressor decreases responses in key threat-related brain areas.

## Requirements
To install Python requirements:
```
pip install -r requirements.txt 
```
R requirements:
> R version 3.6.0 (or higher)  
> Libriaries: brms (version 2.14.0), tidyverse (version 1.3.0)  

## Bayesian multilevel analysis at the level of region of interest
---

#### __Preprocessing__  
Get total number button-presses and STAI scores for all participants.

To extract total number of button presses:
```
$ python 00-make_buttonPress_df.py  
```
output: `data/button_presses.txt`  

To get STAI (state and trait) scores:
```
$ python 01-make_STAI_df.py  
```
output: `data/STAI_scores.txt`  

To extract stressor response from 24 ROIs along with standardized
covariates:
```
$ python 02-extract_stressor_resp.py  
```
output: `data/uncon_v_con_stressor.txt`  

The output contains following standardized covariates: 
- TRAITmean: average trait scores of yoked participants.  
- TRAITdiff: trait difference of the yoked participants (uncontrollable - contollable).  
- STATEmean: average state scores of the yoked participants.  
- STATEdiff: state difference of the yoked participants.  

#### __BML analysis__  

To run BML analysis
```
$ Rscript 03-uncon_v_con_stressor.r
```
output: (R workspace) `uncon_v_con_stressor/results.RData`  
The BML output is saved as results/uncon_v_con_stressor/results.RDtata.

To extract posteriors:
```
$ Rscript 04-uncon_v_con_stressor_posteriors.r
```
output: (posterior) `results/uncon_v_con_stressor/*.txt`  

## Brain-skin conductance correlation
---
#### __Preprocessing__  
To extract brain-skin conductance correlations for the 24 ROIs 
```
$ python 05-make_ROI_SCR_df.py  
```
output: `data/uncon_v_con_ROI_SCR_zscorr.txt`  

#### __BML analysis__
To run BML analysis on the left and right BST-skin conductance correlations:
```
$ Rscript 06-uncon_v_con_ROI_SCR_zscorr.r
```

output:  
    (R workspace) `results/uncon_v_con_ROI_SCR/results.RData`  
    (posterior) `results/uncon_v_con_ROI_SCR/uncon_v_con_lBST_SCR_corr.txt`  
    (posterior) `results/uncon_v_con_ROI_SCR/uncon_v_con_rBST_SCR_corr.txt`  

## Bayesian multilevel analysis of insula voxels
---
#### __Preprocessing__
**Right Insula**: To extract voxelwise stressor responses for right insula voxels:
```
$ python 07-insulaVoxelwiseStressorResp.py \
    --mask data/masks/right_insula_10ROIs.nii.gz \
    --output data/voxelwise  
```
output: `data/voxelwise/uncon_v_con_right_insula_10ROIs.txt`.  

`right_insula_10ROIs.nii.gz` is a right insula mask with 10 sub-divisions. 
                            The sub-divisions were created using k-means clustering,
                            and limiting the total number of clusters to roughly 10.
                            
**Left Insual**: To extract voxelwise stressor responses for insula insula voxels:
```
$ python 07-insulaVoxelwiseStressorResp.py \
    --mask data/masks/left_insula_11ROIs.nii.gz \
    --output data/voxelwise  
```
output: `data/voxelwise/uncon_v_con_left_insula_11ROIs.txt`  

`left_insula_11ROIs.nii.gz` is a left insula mask with 11 sub-divisions.  

#### __BML analysis__

To run BML analysis on __right insula__:
```
$ Rscript 08-uncon_v_con_insula_voxelwise_stressor.r \
    data/voxelwise/uncon_v_con_right_insula_10ROIs.txt \
    results/voxelwise  
```
outputs in `results/voxelwise/uncon_v_con_right_insula_10ROIs/`: 
- R workspace: `results.RData`  
- Fixed effect posterior: `POP.txt`
- Random effect posterior: `ROI_*.txt` & `VOX_*.txt`

To run BML analysis on __left insula__:
```
$ Rscript 08-uncon_v_con_insula_voxelwise_stressor.r \
    data/voxelwise/uncon_v_con_left_insula_11ROIs.txt \
    results/voxelwise  
```
outputs in `results/voxelwise/uncon_v_con_left_insula_11ROIs/`: 
- R workspace: `results.RData`  
- Fixed effect posterior: `POP.txt`
- Random effect posterior: `ROI_*.txt` & `VOX_*.txt`