# Controllability over stressor decreases responses in key threat-related brain areas.

## Bayesian multilevel analysis at the level of region of interest
---

#### __Preprocessing__  
First get the behavioral and STAI scores that with form covarites in the model.

To extract total number of button presses:
```
$ python 00-make_buttonPress_df.py --output data/ --overwrite 1  
output: data/button_presses.txt
```

To gets STAI scores:
```
$ python 01-make_STAI_df.py --output data/ --overwrite 1
output: data/STAI_scores.txt
```

To extract stressor response from 24 ROIs along with standardized
covariates:
```
$ python 02-extractStressorResp.py --output data/ --overwrite 1
output: data/uncon_v_con_stressor.txt
```

#### __BML analysis__  
Requirements to run BML: 
> R version 3.6.0 (or higher)  
> brms (version 2.14.0)
> tidyverse (version 1.3.0)

To run BML analysis
```
$ Rscript 03-uncon_v_con_stressor.r
output: (R workspace) uncon_v_con_stressor/results.RData
```
The BML output is saved as results/uncon_v_con_stressor/results.RDtata.

To extract posteriors:
```
$ Rscript 04-uncon_v_con_stressor_posteriors.r
output: (posterior) results/uncon_v_con_stressor/*.txt
```

## Brain-skin conductance correlation
---
#### __Preprocessing__  
To extract brain-skin conductance correlations for the 24 ROIs 
```
$ python 05-make_ROI_SCR_df.py -o data/ 
output: data/uncon_v_con_ROI_SCR_zscorr.txt
```

#### __BML analysis__
To run BML analysis on the left and right BST-skin conductance correlations:
```
$ Rscript 06-uncon_v_con_ROI_SCR_zscorr.r
output: (R workspace) results/uncon_v_con_ROI_SCR/results.RData  
        (posterior) results/uncon_v_con_ROI_SCR/uncon_v_con_lBST_SCR_corr.txt  
        (posterior) results/uncon_v_con_ROI_SCR/uncon_v_con_rBST_SCR_corr.txt  
```
