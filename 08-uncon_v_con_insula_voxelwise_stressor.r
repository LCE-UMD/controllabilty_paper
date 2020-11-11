## Insula voxelwise BML
library(brms)
library(dplyr)
library(parallel)

args <- commandArgs(TRUE)
datapath <- args[1]
outdir <- paste0(args[2],'/',tools::file_path_sans_ext(basename(datapath)))

dataset <- read.table(datapath,header=TRUE,sep='\t')

# Convert VOX, ROI and Subj columns into factors
dataset$VOX <- factor(dataset$VOX)
dataset$ROI <- factor(dataset$ROI)
dataset$Subj <- factor(dataset$Subj)

# Print number of levels in each grouping variable
print(paste("Number of Subjects:-",nlevels(dataset$Subj)))
print(paste("Number of ROIs:-",nlevels(dataset$ROI)))
print(paste("Number of Voxels:-",nlevels(dataset$VOX)))
# Print number of cores allocated
print(paste0('Number of cores available: ', detectCores(all.tests = FALSE, logical = TRUE)))

# Set nuber of cores to use. 
# To run each chain on a single core, set number of core to 4
print(getOption("mc.cores"))
options(mc.cores = parallel::detectCores())
print(paste0('Number of cores allocated',getOption("mc.cores")))

# Number of iterations for the MCMC sampling
iterations <- 20000
# Number of chains to run
chains <- 4
SCALE <- 1
ns <- iterations*chains/2
# number of sigfigs to show on the table
nfigs <- 4

mod = '1 + TRAITmean + TRAITdiff + STATEmean + STATEdiff + buttPress'
modelForm = paste('Y ~', mod,'+ (1 | gr(Subj, dist= "student")) + (',mod,'| gr(ROI, dist="student")) + (',mod,'|gr(VOX, dist="student"))')
print('Model Formula:')
print(modelForm)

# get dafualt priors for data and model
priorRBA <- get_prior(formula = modelForm,data=dataset,family = 'student')
# You can assign prior or your choice to any of the parameter in the table below. 
# For example. If you want to assign a student_t(3,0,10) prior to all parameters of class b, 
# the following line does that for you. Parameters in class b are the population effects (cond, STATE and TRAIT)
priorRBA$prior[1] <- "student_t(3, 0, 10)"
priorRBA$prior[7] <- "lkj(2)"
priorRBA$prior[10:12] <- "gamma(3.325,0.1)"
priorRBA$prior[13] <- "student_t(3, 0, 10)"
priorRBA$prior[14] <- "gamma(3.325,0.1)"
priorRBA$prior[15] <- "student_t(3, 0, 10)"
priorRBA$prior[32] <- "student_t(3, 0, 10)"

print("")
# Print the table with priors
print(priorRBA)

# Create a result directory
if (!dir.exists(outdir)){dir.create(outdir,recursive = TRUE)}

# Generate the Stan code for our own reference
stan_code <- make_stancode(modelForm,data=dataset,chains = chains,family = 'student',prior = priorRBA)
cat(stan_code,file = paste0(outdir,"/stancode.stan"),sep='\n')

# Following run the BML model
fm <- brm(modelForm,
          data=dataset,
          chains = chains,
          family = 'student',
          prior = priorRBA,
          inits=0, iter=iterations, 
          control = list(adapt_delta = 0.99, max_treedepth = 15))

# Extract posteriors for fixed (aa) and random effects (bb)
aa <- fixef(fm, summary = FALSE)/SCALE # Population-Level Estimates
bb <- lapply(ranef(fm, summary = FALSE), `/`, SCALE) # Extract Group-Level (or random-effect)

write.table(as.data.frame(aa),file=paste(outdir,'/POP.txt'),sep='\t',
            col.names = TRUE, row.names = FALSE)

for (col in colnames(aa)){
  roi_eff <- data.frame(bb[['ROI']][,,col])
  write.table(roi_eff,paste0(outdir,'/ROI_',col,'.txt'),sep='\t',
              col.names = TRUE, row.names = FALSE)
  
  vox_eff <- data.frame(bb[['VOX']][,,col])
  write.table(vox_eff,paste0(outdir,'/VOX_',col,'.txt'),sep='\t',
            col.names = TRUE, row.names = FALSE)
}

# Save the results as a RData file
save.image(file=paste0(outdir,'/results.RData'))
