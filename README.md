# Controllability over stressor decreases responses in key threat-related brain areas.

## Requirements
To install Python requirements:
```
pip install -r requirements.txt 
```
R requirements:
> R version 3.6.0 (or higher)  
> Libriaries: brms (version 2.14.0), tidyverse (version 1.3.0)  

## Data
```
    |-- data
        |-- behavioral
            |-- button_presses.txt               <- Total number of button presses.  
            |-- STAT_score.txt                   <- State and Trait scores.
        |-- CON_yoked_table.xlsx                 <- Participant IDs and usable runs.
        |-- CON_yoked_table_SCR.xlsx
        |-- stressor_response_estimates
            |-- Fig4.txt                         <- stressor response estimates
                                                    for 24 ROIs (paper: Figure 4).
            |-- Fig6-left_PI.1D                  <- stressor response estimates 
                                                    in left posterior insula.
            |-- Fig6-right_PI.1D                 <- stressor response estimates
                                                    in right posterior insula.
            |-- Fig7-PCC.txt                     <- stressor response estimates
                                                    in posterior cingulate cortex.
            |-- Fig8-MFG                         <- stressor response estimates
                                                    in middle frontal gyrus.
        |-- masks
            |-- emoproxII_ROIs_final.nii.gz      <- 24 ROI mask.
            |-- emoproxII_ROIs_final_info.txt    <- ROI information.
            |-- left_insula_11ROIs.nii.gz        <- left insula mask with 11 sub-ROIs.
            |-- right_insula_10ROIs.nii.gz       <- right insula mask with 10 sub-ROIs.
        |-- ROIwise
            |-- uncon_v_con_stressor.txt         <- Brain (24 ROIs) stressor response.
            |-- uncon_v_con_ROI_SCR_zscorr.txt   <- Brain-skin conductance response (SCR)
                                                    correlation.
        |-- subjects
            |-- skin_conductance
                |-- CON???_bucket_LSS.1D         <- participant's trial-by-trail
                                                    SCR response to stresssor.
            |-- stressor_canonical
                |-- CON???_bucket_REML.1D        <- participant's ROI response to stressor.
            |-- stressor_canonical_voxelwise
                |-- CON???_shock_beta.nii.gz     <- participant's whole-brain response 
                                                    to stressor.
            |-- stressor_trial_by_trail          
                |-- CON???_betas_3dLSS.1D        <- participant's trial-by-trail
                                                    ROI response to stresssor.
```  

## Bayesian multilevel analysis at the level of region of interest
---

### __Preprocessing__  

To extract stressor response from 24 ROIs along with standardized
covariates:
```
$ python 02-extract_stressor_resp.py  
```
Output: `data/ROIwise/uncon_v_con_stressor.txt`  

The output contains following standardized covariates: 
- TRAITmean: average trait scores of yoked participants.  
- TRAITdiff: trait difference of the yoked participants (uncontrollable - contollable).  
- STATEmean: average state scores of the yoked participants.  
- STATEdiff: state difference of the yoked participants.  
- BPdiff: difference in the total number of button-presses in the yoked participants (uncontrollable - contollable).  

### __BML analysis__  

To run BML analysis
```
$ Rscript 03-uncon_v_con_stressor.r
```
Output: (R workspace) `results/ROIwise/uncon_v_con_stressor/results.RData`  

To extract posteriors from the saved R workspace:
```
$ Rscript 04-uncon_v_con_stressor_posteriors.r
```
Outputs in `results/ROIwise/uncon_v_con_stressor/`:
- ROI posteriors for __uncontrollability > controllability__ contrast (Paper: `Figure 3`): `Intercept_post.txt`  

## Brain-skin conductance correlation
---
### __Preprocessing__  
To extract brain-skin conductance correlations for the 24 ROIs 
```
$ python 05-make_ROI_SCR_df.py  
```
Output: `data/ROIwise/uncon_v_con_ROI_SCR_zscorr.txt`  

### __BML analysis__
To run BML analysis on the left and right BST-skin conductance correlations:
```
$ Rscript 06-uncon_v_con_ROI_SCR_zscorr.r
```

Outputs in `results/ROIwise/uncon_v_con_ROI_SCR/`:  
- (R workspace): `results.RData`  
- Posterior for Left BST & skin conductance correlation (Paper: `Figure 5`): `uncon_v_con_lBST_SCR_corr.txt`  
- Posterior for Right BST & skin conductance correlation (Paper: `Figure 5`): `uncon_v_con_rBST_SCR_corr.txt`  

## Bayesian multilevel analysis of insula voxels
---
### __Preprocessing__
**Right Insula**: To extract voxelwise stressor responses for right insula voxels:
```
$ python 07-insulaVoxelwiseStressorResp.py \
    --mask data/masks/right_insula_10ROIs.nii.gz \
    --output data/voxelwise  
```
Output: `data/voxelwise/uncon_v_con_right_insula_10ROIs.txt`.  

`right_insula_10ROIs.nii.gz` is a right insula mask with 10 sub-divisions. 
                            The sub-divisions were created using k-means clustering,
                            and limiting the total number of clusters to roughly 10.
                            
**Left Insual**: To extract voxelwise stressor responses for insula insula voxels:
```
$ python 07-insulaVoxelwiseStressorResp.py \
    --mask data/masks/left_insula_11ROIs.nii.gz \
    --output data/voxelwise  
```
Output: `data/voxelwise/uncon_v_con_left_insula_11ROIs.txt`  

`left_insula_11ROIs.nii.gz` is a left insula mask with 11 sub-divisions.  

### __BML analysis__

To run BML analysis on __right insula__:
```
$ Rscript 08-uncon_v_con_insula_voxelwise_stressor.r \
    data/voxelwise/uncon_v_con_right_insula_10ROIs.txt \
    results/voxelwise  
```
Outputs in `results/voxelwise/uncon_v_con_right_insula_10ROIs/`: 
- R workspace: `results.RData`  
- Fixed effect posterior: `POP.txt`
- Random effect posterior: `ROI_*.txt` & `VOX_*.txt`

To run BML analysis on __left insula__:
```
$ Rscript 08-uncon_v_con_insula_voxelwise_stressor.r \
    data/voxelwise/uncon_v_con_left_insula_11ROIs.txt \
    results/voxelwise  
```
Outputs in `results/voxelwise/uncon_v_con_left_insula_11ROIs/`:
- R workspace: `results.RData`  
- Fixed effect posterior: `POP.txt`
- Random effect posterior: `ROI_*.txt` & `VOX_*.txt`

To combine the voxelwise fixed and random effect posteriors and compute the posterior probabilities (`P+`) and render onto the MNI brain for visualization:  
For __right insula__,  
```
$ python 09-uncon_v_con_insula_voxelwise_stressor_posteriors.py \
    -m data/masks/right_insula_10ROIs.nii.gz \
    -o results/voxelwise  
```

For __left insula__,
```
$ python 09-uncon_v_con_insula_voxelwise_stressor_posteriors.py \
    -m data/masks/left_insula_11ROIs.nii.gz \
    -o results/voxelwise
```

Outputs:
- Right insula P+ map (Paper: `Figure 6`): `results/voxelwise/uncon_v_con_right_10ROIs/right_insula_10ROIs_Pmap.nii.gz`  
- Left insula P+ map (Paper: `Figure 6`): `results/voxelwise/uncon_v_con_left_11ROIs/left_insula_11ROIs_Pmap.nii.gz`  

