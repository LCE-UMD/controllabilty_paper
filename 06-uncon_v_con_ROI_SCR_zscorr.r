library(brms)
library(tidyverse) # needed for data manipulation.
library(parallel)

# Create a result directory
outdir <- "results/ROIwise/uncon_v_con_ROI_SCR"
if (!dir.exists(outdir)){dir.create(outdir,recursive = TRUE)}

df <- read.table('data/ROIwise/uncon_v_con_ROI_SCR_zscorr.txt',header = TRUE,sep='\t')

################################ left BST ##################################
# Filter out all ROIs except left BST
dataTable <- filter(df,ROI=='L BST')
head(dataTable)

# Model
mod = '1 + TRAITmean + TRAITdiff + STATEmean + STATEdiff + BPdiff'
modelForm = paste('Y ~',mod)
priorRBA <- get_prior(formula = modelForm,data=dataTable,family = 'student')
priorRBA$prior[1] <- "student_t(3,0,10)"
priorRBA$prior[7] <- "student_t(3,0,10)"
priorRBA$prior[8] <- "gamma(3.325,0.1)"
priorRBA$prior[9] <- "student_t(3,0,10)"
priorRBA

iterations <- 10000
chains <- 4
SCALE <- 1
ns <- iterations*chains/2

EOI = 'Intercept,TRAITmean,TRAITdiff,STATEmean,STATEdiff,BPdiff'

# number of ROIs
nR <- nlevels(dataTable$ROI)

# number of sigfigs to show on the table
nfigs <- 4

lBST <- brm(modelForm,
          data=dataTable,
          chains = chains,
          family = 'student',
          prior = priorRBA,
          inits=0, iter=iterations, 
          control = list(adapt_delta = 0.99, max_treedepth = 15))

############################## right BST ####################################
dataTable <- filter(df,ROI=='R BST')
head(dataTable)

rBST <- brm(modelForm,
          data=dataTable,
          chains = chains,
          family = 'student',
          prior = priorRBA,
          inits=0, iter=iterations, 
          control = list(adapt_delta = 0.99, max_treedepth = 15))

#############################################################################
getPosteriors <- function(model){
    data <- fixef(model,summary = FALSE)
}

write.csv(getPosteriors(lBST),paste0(outdir,"/uncon_v_con_lBST_SCR_corr.txt"))
write.csv(getPosteriors(rBST),paste0(outdir,"/uncon_v_con_rBST_SCR_corr.txt"))

# Save the results as a RData file
save.image(file=paste0(outdir,"/results.RData"))
