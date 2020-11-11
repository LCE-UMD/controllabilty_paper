require(brms)

outdir <- 'results/ROIwise/uncon_v_con_stressor/'

# Load the BML output image
load(paste0(outdir,'results.RData'))

EOIq <- unlist(lapply(strsplit(mod,'\\+')[[1]],trimws))
if(!('Intercept' %in% EOIq)) EOIq <- c('Intercept', EOIq)
EOIq <- EOIq[!grepl('1', EOIq)]

print('Fixed effect model terms:')
print(EOIq)

# Extract posteriors for fixed (aa) and random effects (bb)
aa <- fixef(fm, summary = FALSE)/SCALE # Population-Level Estimates
bb <- lapply(ranef(fm, summary = FALSE), `/`, SCALE) # Extract Group-Level (or random-effect)


# Sum the fixed and random effect (only ROI) posterior. Following function does this.
# The function adds the "the poaterior of the global intercept (or slope) with the 
# intercept (or slope) posterior of each roi. bb[['ROI']][,,tm] is a matrix with 
# columns containing posteriors for every ROI
# aa[,tm] is a single column posterior for the global intercept/slope.
# intercept or slope is defined by tm, where tm is the index of the model term.
psROI <- function(aa, bb, tm) {
  ps <- apply(bb[['ROI']][,,tm], 2, '+', aa[,tm])
  return(ps)
}

# generates ridge plots for effects (Intercept, cond, TRAIT and STATE), and saves them.
for (ii in 1:length(EOIq)){
  posteriors <- psROI(aa, bb, EOIq[ii])
  write.table(posteriors,paste0(outdir,EOIq[ii],'_post.txt'),sep='\t',
              col.names = TRUE, row.names = FALSE)
}