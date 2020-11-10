library(brms)
library(tidyverse) # needed for data manipulation.

df <- read.table('DATA/uncon_v_con_stressor.txt',header = TRUE,sep = "\t")
head(df)

iterations <- 10000
chains <- 4
SCALE <- 1
ns <- iterations*chains/2

df$Pair<- factor(df$Pair)
df$ROI <- factor(df$ROI)
# number of ROIs
print(paste0('Number of ROI: ',nlevels(df$ROI)))
print(paste0('Number of pairs: ',nlevels(df$Pair)))
print(paste0('Number of cores available: ', detectCores(all.tests = FALSE, logical = TRUE)))
# number of sigfigs to show on the table
nfigs <- 4

head(df)

# Set nuber of cores to use. 
# To run each chain on a single core, set number of core to 4
print(getOption("mc.cores"))
options(mc.cores = parallel::detectCores())
print(getOption("mc.cores"))

mod = ' 1 + TRAITmean + TRAITdiff + STATEmean + STATEdiff + buttPress'
modelForm = paste('beta ~',mod,'+ (1 | gr(Pair, dist= "student")) + (',mod,'| gr(ROI, dist="student"))')
priorRBA <- get_prior(formula = modelForm,data=df,family = 'student')
priorRBA$prior[1] <- "student_t(3,0,10)"
priorRBA$prior[7] <- "lkj(2)"
priorRBA$prior[9:10] <- "gamma(3.325,0.1)"
priorRBA$prior[12] <- "gamma(3.325,0.1)"
print(priorRBA)

# Create a result directory
outdir <- "results/uncon_v_con_stressor"
if (!dir.exists(outdir)){dir.create(outdir)}

# Generate the Stan code for our own reference
stan_code <- make_stancode(modelForm,
                           data=df,
                           chains = chains,
                           family = 'student',
                           prior = priorRBA)
cat(stan_code,file = paste0(outdir,"/stancode.stan"),sep='\n')

# Following run the BML model
fm <- brm(modelForm,
          data=df,
          chains = chains,
          family = 'student',
          prior = priorRBA,
          inits=0, iter=iterations, 
          control = list(adapt_delta = 0.99, max_treedepth = 15))

# Save the results as a RData file
save.image(file=paset0(outdir,"/results.RData"))
