# IMPORT LIBRARIES:

library(tidyverse)
library(ggridges)
library(ggthemes)
library(cowplot)
library(here)

rois <- c("L aMCC", "L basolateral Amygdala",	"L central/medial amygdala",	"L dorsal anterior Insula",	"L ventral anterior Insula",	"L BST",	"L mid Hippocampus",	"L mid posterior Insula",	"L PAG",	"L Thalamus",	"PCC/Precuneus",	"PCC",	"R aMCC",	"R basolateral Amygdala",	"R central/medial Amygdala",	"R dorsal anterior Insula",	"R ventral anterior Insula",	"R BST",	"R mid posterior Insula",	"R PAG",	"R posterior Hippocampus",	"R Thalamus",	"anterior vmPFC",	"posterior vmPFC")

##### FIGURE 3: POSTERIOR PLOT #####
# this plot is larger than others so did not use function.

print("Starting figure 3...")

# import & transform to long form
df <- read.table(here("results/ROIwise/uncon_v_con_stressor","Intercept_post.txt"), sep = "", header = T)
df$X <- NULL
colnames(df) <- rois
iterations <- length(df[,1])
df.long <- df %>% gather(ROI)
colnames(df.long) <- c("ROI", "Y")
df.long <- df.long %>%
  mutate(index = rep(1:length(rois), each = iterations)) %>%
  group_by(ROI) %>%
  mutate(mean = mean(Y)) %>%
  mutate(p = ((sum(Y>0)/iterations))) %>%
  mutate(p.plot = p) %>%
  mutate(p.plot = replace(p.plot, p.plot > 0.15 & p.plot < 0.85, NA))

# summary dataframe for P+ labels
P <- df.long %>%
  select(ROI, index, mean, p) %>%
  unique() %>%
  arrange(p)

# Because this plot is larger it has different settings than 90mm plots
intercept <- ggplot(df.long, aes(x = Y, y = as.numeric(reorder(index,p)), group = ROI, fill = p.plot)) +
  coord_cartesian(xlim = c(-0.14, 0.23)) +
  geom_density_ridges(quantile_lines = TRUE,
                      quantiles = 2,
                      scale = 2.3,
                      rel_min_height = .01,
                      color = "#404040",
                      size = .75) +
  geom_vline(xintercept = 0, alpha = .85, color = "black", size = .8) +
  scale_y_continuous(breaks = 1:length(P$ROI),
                     expand = c(0,0.1),
                     labels = P$ROI,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(P$ROI),
                                         labels = format(round(P$p, 3),nsmall = 2))) +
  scale_x_continuous(breaks = c(-0.1, 0, 0.1, 0.2), labels = c("-0.10", "0", "0.10", "0.20")) +
  scale_fill_gradientn(limits = c(0,1),
                       colors = c("blue","cyan",
                                  "gray","yellow","red"),
                       values = c(0,0.15,
                                  0.150000001, 0.85,
                                  0.850000001, 1.0),
                       breaks = c(0, 0.15, 0.85, 1)) +
  scale_color_gradientn(limits = c(0,1),
                        colors = c("blue","cyan",
                                   "gray","yellow","red"),
                        values = c(0,0.15,
                                   0.150000001, 0.85,
                                   0.850000001, 1.0),
                        breaks = c(0, 0.15, 0.85, 1)) +
  guides(fill = guide_colorbar(barwidth = 1.25,
                               barheight = 8,
                               nbin = 50,
                               frame.colour = "black",
                               frame.linewidth = 1.25,
                               ticks.colour = "black")) +
  theme_stata() +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.grid.major.y = element_line(color = "lightgray", size = .75),
    plot.title = element_text(size = 12.5, margin = unit(c(0,0.1,.25,02),"cm"), face = "plain", hjust = 0.5),
    legend.title = element_text(size = 11, hjust = 0),
    legend.text = element_text(size = 8, angle = 0),
    legend.position = c(.85,.31),
    legend.background = element_rect(),
    legend.box.background = element_rect(colour = "black", size = .5),
    text = element_text(family = "Helvetica"),
    axis.title = element_text(size = 16),
    axis.line = element_line(size = .25),
    axis.text.y = element_text(size= 9, color = "black", margin = unit(c(0,-0.05,0,0.05),"cm"), angle = 0, vjust = 0),
    axis.text.y.right = element_text(size = 9, color = "black",margin = unit(c(0,0,0,-0.05),"cm"), angle = 0),
    axis.text.x = element_text(size = 9, color = "black", margin = unit(c(0.04,0,0,0),"cm")),
    axis.ticks.x = element_line(size = .5),
    axis.ticks.length=unit(.2, "cm"),
    axis.ticks.y = element_blank()) +
  labs(
    x = NULL,
    y = NULL,
    title = "Uncontrollable vs. Controllable stressor",
    fill = "P+")

plot <- ggdraw(intercept) +
  draw_text("P+", x = .94, y = .92, size= 11) +
  draw_text("controllable > uncontrollable", x = .31, y = 0.015, size = 11) +
  draw_text("uncontrollable > controllable", x = .74, y = 0.015, size = 11)

ggsave(here('plots','Fig3_intercept_posteriors.pdf'), plot = plot, dpi = 600, height = 120, width = 180, units = "mm")

print("Done with Figure 3, onto the next...")



# Fig 5. BST


lBST <- read.table(here('results/ROIwise/uncon_v_con_ROI_SCR','uncon_v_con_lBST_SCR_corr.txt'), header = TRUE, sep = ',')
lBST$X <- NULL
lBST$ROI <- rep("lBST", length(lBST[,1]))

rBST <- read.table(here('results/ROIwise/uncon_v_con_ROI_SCR','uncon_v_con_rBST_SCR_corr.txt'), header = TRUE, sep = ',')
rBST$X <- NULL
rBST$ROI <- rep("rBST",length(rBST[,1]))

BST <- rbind(lBST,rBST)


BST <- BST %>%
  select("Intercept","ROI") %>%
  mutate(Measure = rep("Controllability", length(BST[,1]))) %>%
  mutate(index = rep(1:2, each = iterations)) %>%
  group_by(ROI) %>%
  mutate(p.val = rep((sum(Intercept > 0)/iterations)))

labels <- BST %>%
  select("ROI","p.val","index") %>%
  unique()


BST.control <- BST %>%
  ggplot(aes(x = Intercept, y = as.numeric(index), group = ROI, fill = p.val)) +
  coord_cartesian(xlim = c(-0.13, 0.34), clip = 'off') +
  geom_density_ridges(quantile_lines = TRUE,
                      quantiles = 2,
                      scale = 2.7,
                      rel_min_height = .01,
                      color = "#404040",
                      size = .35) +
  geom_vline(xintercept = 0, alpha = .85, color = "black", size = .45) +
  scale_x_continuous(breaks = c(-0.1, 0, 0.1, 0.2, 0.3),
                     labels = c("-0.10", "0", "0.10", "0.20", "0.30")) +
  scale_y_continuous(breaks = 1:length(labels$index),
                     expand = c(0,0),
                     labels = labels$ROI,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(labels$index),
                                         labels = format(round(labels$p.val, 3),nsmall = 2))) +
  scale_fill_gradientn(limits = c(0,1),
                       colors = c("blue","cyan",
                                  "gray","yellow","red"),
                       values = c(0,0.15,
                                  0.150000001, 0.85,
                                  0.850000001, 1.0),
                       breaks = c(0, 0.15, 0.85, 1)) +
  guides(fill = guide_colorbar(barwidth = .5,
                               barheight = 3,
                               nbin = 50,
                               frame.colour = "black",
                               frame.linewidth = .75,
                               ticks.colour = "black")) +
  theme_stata() +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.grid.major.y = element_line(color = "lightgray", size = .25),
    plot.title = element_text(size = 7.5, margin = unit(c(0,0.1,.1,02),"cm"), face = "plain", hjust = 0.5),
    legend.title = element_text(size = 5, hjust = 0),
    legend.text = element_text(size = 4, angle = 0),
    legend.position = c(.87,.72),
    legend.background = element_rect(),
    legend.box.background = element_rect(colour = "black", size = .05),
    text = element_text(family = "Helvetica"),
    axis.title = element_text(size = 16),
    axis.line = element_line(size = .15),
    axis.text.y = element_text(size= 5, color = "black", margin = unit(c(0,-0.05,0,0.05),"cm"), angle = 0, vjust = 0),
    axis.text.y.right = element_text(size = 5, color = "black",margin = unit(c(0,0,0,-0.05),"cm"), angle = 0),
    axis.text.x = element_text(size = 5, color = "black", margin = unit(c(0.04,0,0,0),"cm")),
    axis.ticks.x = element_line(size = .25),
    axis.ticks.length=unit(.12, "cm"),
    axis.ticks.y = element_blank()) +
  labs(
    x = NULL,
    y = NULL,
    title = "Controllability",
    fill = "P+")



plot <- ggdraw(BST.control) +
  draw_text("P+", x = .93, y = .9, size= 6) +
  draw_text("controllable > uncontrollable", x = .31, y = 0.015, size = 6) +
  draw_text("uncontrollable > controllable", x = .78, y = 0.015, size = 6)

ggsave(here('plots','Fig5_BST_SCR.pdf'), plot = plot, dpi = 600, height = 45, width = 70, unit = "mm")



##### SUPPLEMENTARY FIGURE 1: Covariate plots #####



# create a function to generate plots

func <- function(file,xMin,xMax,breaks,title){
  df <- utils::read.table(here('results/ROIwise/uncon_v_con_stressor',file), fill = TRUE, sep = "", header = T)
  print(paste("Dataset",file, "imported!"))
  rois <- c("L aMCC", "L basolateral Amygdala",	"L central/medial amygdala",	"L dorsal anterior Insula",	"L ventral anterior Insula",	"L BST",	"L mid Hippocampus",	"L mid posterior Insula",	"L PAG",	"L Thalamus",	"PCC/Precuneus",	"PCC",	"R aMCC",	"R basolateral Amygdala",	"R central/medial Amygdala",	"R dorsal anterior Insula",	"R ventral anterior Insula",	"R BST",	"R mid posterior Insula",	"R PAG",	"R posterior Hippocampus",	"R Thalamus",	"anterior vmPFC",	"posterior vmPFC")
  df$X <- NULL
  colnames(df) <- rois
  iterations <- length(df[,1])
  df.long <- df %>% gather(ROI)
  colnames(df.long) <- c("ROI", "Y")
  df.long <- df.long %>%
    mutate(index = rep(1:length(rois), each = iterations)) %>%
    group_by(ROI) %>%
    mutate(mean = mean(Y)) %>%
    mutate(p = ((sum(Y>0)/iterations))) %>%
    mutate(p.plot = p) %>%
    mutate(p.plot = replace(p.plot, p.plot > 0.15 & p.plot < 0.85, NA))


  print("Creating summary dataset...")
  P <- df.long %>%
    dplyr::select(ROI, index, mean, p) %>%
    unique() %>%
    dplyr::arrange(p)


  print(paste("Plotting",title,"..."))
  plot <- ggplot(df.long, aes(x = Y, y = as.numeric(reorder(index,p)),
                              group = ROI,
                              fill = df.long$p.plot)) +
    coord_cartesian(xlim = c(xMin, xMax)) +
    geom_density_ridges(quantile_lines = TRUE,
                        quantiles = 2,
                        scale = 2.5,
                        rel_min_height = .01,
                        color = "#404040",
                        size = .2) +
    geom_vline(xintercept = 0, alpha = .85, color = "black", size = .3) +
    scale_y_continuous(breaks = 1:length(P$ROI),
                       expand = c(0,0.1),
                       labels = P$ROI,
                       sec.axis = sec_axis(~.,
                                           breaks = 1:length(P$ROI),
                                           labels = format(round(P$p, 3),nsmall = 2))) +
    scale_x_continuous(breaks = breaks,
                       labels = format((breaks), nsmall = 2)) +
    scale_fill_gradientn(limits = c(0,1),
                         colors = c("blue","cyan",
                                    "#505154","yellow","red"),
                         values = c(0,0.15,
                                    0.150000001, 0.85,
                                    0.850000001, 1.0),
                         breaks = c(0, 0.15, 0.85, 1)) +
    labs(x = NULL, y = NULL, title = title) +
    theme_stata() +
    theme(
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.grid.major.y = element_line(color = "lightgray", size = .25),
      plot.title = element_text(size = 8, margin = unit(c(0,0.1,.1,02),"cm"), face = "plain", hjust = 0.5),
      legend.title = element_text(size = 6, hjust = 0),
      legend.text = element_text(size = 5, angle = 0),
      legend.background = element_rect(),
      legend.box.background = element_rect(colour = "black", size = .25),
      legend.position = 'none',
      text = element_text(family = "Helvetica"),
      axis.title = element_text(size = 16),
      axis.line = element_line(size = .25),
      axis.text.y = element_text(size= 5.25, color = "black", margin = unit(c(0,-0.05,0,0.05),"cm"), angle = 0, vjust = 0),
      axis.text.y.right = element_text(size = 5.25, color = "black",margin = unit(c(0,0,0,-0.05),"cm"), angle = 0),
      axis.text.x = element_text(size = 5.25, color = "black", margin = unit(c(0.04,0,0,0),"cm")),
      axis.ticks.x = element_line(size = .25),
      axis.ticks.length=unit(.12, "cm"),
      axis.ticks.y = element_blank())

  plot <- ggdraw(plot) +
    draw_text("P+", x = .94, y = .9, size= 6.5) +
    draw_text("controllable > uncontrollable", x = .31, y = 0.015, size = 6.5) +
    draw_text("uncontrollable > controllable", x = .78, y = 0.015, size = 6.5)



  output.name <- paste('Supplementary_',title,"_posteriors.pdf")

  ggsave(here('plots',output.name), plot = plot, dpi = 600, height = 60, width = 90, units = "mm")
  print('done!')

}


### Panel a - State anxiety: Mean

file <- 'STATEmean_post.txt'
xMin <- -0.09
xMax <- 0.13
breaks <-  c(-0.1,-0.05, 0, 0.05, 0.1)
title <- "State anxiety: mean"

func(file, xMin, xMax, breaks, title)


### Panel b - State anxiety: difference

data <- 'STATEdiff_post.txt'
xMin <- -0.16
xMax <- 0.15
breaks <-  c(-0.1, 0, 0.1)
title <- "State anxiety: difference"

func(data, xMin, xMax, breaks, title)


### Panel c - Trait anxiety: mean

data <- 'TRAITmean_post.txt'
xMin <- -0.09
xMax <- 0.07
breaks <- c(-0.05, 0, 0.05)
title <- "Trait anxiety: mean"

func(data, xMin, xMax, breaks, title)


### Panel d - Trait anxiety: difference

data <- 'TRAITdiff_post.txt'
xMin <- -0.09
xMax <- 0.07
breaks <- c(-0.05, 0, 0.05)
title <- "Trait anxiety: mean"

func(data, xMin, xMax, breaks, title)


### Panel e - Button presses


data <- 'BPdiff_post.txt'
xMin <- -0.20
xMax <- 0.18
breaks <- c(-0.20, -0.10, 0, 0.10, 0.20)
title <- "Button Presses"

func(data,xMin,xMax,breaks,title)

print("All posterior plots generated!")
