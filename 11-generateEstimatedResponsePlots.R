library(tidyverse)
library(here)
# plyr is also required but not front loaded because of conflicts with here package


df <- read.table(here('data/stressor_response_estimates','Fig4.txt'), sep = ',', header = TRUE)
x.axis.breaks <- seq(from = 0, to = 13.75, by = 3)
x.axis.labs <- c("0.00",  "3.00",  "6.00", " 9.00", "12.00")
x.axis.expand <- c(0.0,0.0)
y.axis.expand <- c(0.04, 0)


# Run these sections to be able to manually manipulate facet axes! It is very dumb that this still is an issue, but this works. This solution is from: https://fishandwhistle.net/post/2018/modifying-facet-scales-in-ggplot2/
# Modifying specific scales
# The problem of modifying specific scales has long been a requested feature in ggplot2, but it’s a bit tricky to implement and hard to think of exactly the syntax should go. The solution I use below is far from perfect but gets the job done, and probably will for a while. It hinges on modifying the init_scales() method of the underlying Facet class. From theGitHub source:
#Essentially, when the plot is built, init_scales() copies whatever scale you’ve specified for the x and y aesthetics (or the default one if you haven’t) over all the panels. As long as you know what panel you’re looking to modify, you could create a specification that overrides an x or y scale at a particular panel. There’s probably better ways to go about this, but one might be this:


scale_override <- function(which, scale) {
  if (!is.numeric(which) || (length(which) != 1) || (which%%1 != 0)) {
    stop("which must be an integer of length 1")
  }

  if (is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
    stop("scale must be an x or y position scale")
  }

  structure(list(which = which, scale = scale), class = "scale_override")
}

# Next, we need to implement a version of init_scales() that looks for scale
# overrides and replaces the default scale with the overridden one.
# It<e2><80><99>s a bit verbose, but because init_scales() isn<e2><80><99>t
# really well-documented it<e2><80><99>s hard to predict what gets called with
# what. I<e2><80><99>m overriding FacetWrap here, but it could easily be
# FacetGrid.

CustomFacetWrap <- ggproto("CustomFacetWrap", FacetWrap, init_scales = function(self,
  layout, x_scale = NULL, y_scale = NULL, params) {
  # make the initial x, y scales list
  scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale,
    params)

  if (is.null(params$scale_overrides))
    return(scales)

  max_scale_x <- length(scales$x)
  max_scale_y <- length(scales$y)

  # ... do some modification of the scales$x and scales$y here based on
  # params$scale_overrides
  for (scale_override in params$scale_overrides) {
    which <- scale_override$which
    scale <- scale_override$scale

    if ("x" %in% scale$aesthetics) {
      if (!is.null(scales$x)) {
        if (which < 0 || which > max_scale_x)
          stop("Invalid index of x scale: ", which)
        scales$x[[which]] <- scale$clone()
      }
    } else if ("y" %in% scale$aesthetics) {
      if (!is.null(scales$y)) {
        if (which < 0 || which > max_scale_y)
          stop("Invalid index of y scale: ", which)
        scales$y[[which]] <- scale$clone()
      }
    } else {
      stop("Invalid scale")
    }
  }

  # return scales
  scales
})


# Lastly, we need a constructor function. Unfortunately, facet_wrap()
# doesn<e2><80><99>t let you specify a Facet class, so we have to hack the result
# of facet_wrap() so we can use the syntax that we are used to with the function.

facet_wrap_custom <- function(..., scale_overrides = NULL) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)

  # sanitize scale overrides
  if (inherits(scale_overrides, "scale_override")) {
    scale_overrides <- list(scale_overrides)
  } else if (!is.list(scale_overrides) || !all(vapply(scale_overrides, inherits,
    "scale_override", FUN.VALUE = logical(1)))) {
    stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  }

  facet_super$params$scale_overrides <- scale_overrides

  ggproto(NULL, CustomFacetWrap, shrink = facet_super$shrink, params = facet_super$params)
}


print("Custom facet wrap scale override loaded in! Now plot y-axes can be specified on an individual basis!")
print("Thanks https://fishandwhistle.net/post/2018/modifying-facet-scales-in-ggplot2/ !!")

print("-----------------------------------------------------------------------------------------------------")
print("Starting figure 4 panels...")


#### Figure 4: Estimated stressor responses


# Panel a

# create average for each time point in each ROI for both groups (control/uncontrol)
df <- df %>% group_by(ROI, Time, Group) %>%
  dplyr::summarize(Mean = mean(Response, na.rm=TRUE))

# 8 ROIs separated for panel A (those with highest effect for uncontrol > control)
P.ROIs <- df %>% group_by(ROI) %>%
  filter(ROI == "L lBST" |
           ROI == "R BST" |
           ROI == "L Anterior dorsal insula" |
           ROI == "L Thalamus (emoproxI shock intersection)" |
           ROI == "R Amygdala (BLBM)" |
           ROI == "R Anterior dorsal insula" |
           ROI == "L Amygdala (CeME)" |
           ROI == "R Thalamus (emoproxI shock intersection)")


P.ROIs$ROI <- factor(P.ROIs$ROI)

# change ordering to correspond to P+ values
P.ROIs$ROI <- factor(P.ROIs$ROI, levels = c("L lBST",
                                            "R BST",
                                            "L Anterior dorsal insula",
                                            "L Thalamus (emoproxI shock intersection)",
                                            "R Amygdala (BLBM)",
                                            "R Anterior dorsal insula",
                                            "L Amygdala (CeME)",
                                            "R Thalamus (emoproxI shock intersection)"))

# change ROI names to something plotting appropriate
P.ROIs$ROI <- plyr::revalue(P.ROIs$ROI, c("L lBST" = "L BST",
                                    "R BST" = "R BST",
                                    "L Anterior dorsal insula" = "L dorsal anterior Insula",
                                    "L Thalamus (emoproxI shock intersection)"  = "L Thalamus",
                                    "R Amygdala (BLBM)" = "R basolateral Amygdala",
                                    "R Anterior dorsal insula" = "R dorsal anterior Insula",
                                    "L Amygdala (CeME)" = "L central/medial Amygdala",
                                    "R Thalamus (emoproxI shock intersection)" = "R Thalamus"))

# Manually insert P+ values as characters
P <- c("P+ = 0.998",
       "P+ = 0.994",
       "P+ = 0.990",
       "P+ = 0.939",
       "P+ = 0.925",
       "P+ = 0.904",
       "P+ = 0.901",
       "P+ = 0.883")

P <- rep(P, each = length(unique(P.ROIs$Time))*2)
P.ROIs <- P.ROIs %>% arrange(ROI) %>% add_column(P)


# Panel b: positive-going ROIs

# filter positive-going ROIs that are not in Panel A
Pos <- df %>% group_by(ROI) %>% filter(max(Mean) > 0.25 & ROI ==  "R Amygdala (CeME)" |
                                         ROI == "L aMCC" |
                                         ROI == "R aMCC" |
                                         ROI == "L Anterior ventral insula (anterior pole)" |
                                         ROI == "R Anterior ventral insula (anterior pole)" |
                                         ROI == "L Mid posterior insula" |
                                         ROI == "R Mid posterior insula" |
                                         ROI == "L PAG" |
                                         ROI == "R PAG" |
                                         ROI == "L Amygdala (BLBM)")

Pos$ROI <- factor(Pos$ROI)
Pos$ROI <- factor(Pos$ROI, levels = c("L aMCC",
                                      "R aMCC",
                                      "L Anterior ventral insula (anterior pole)",
                                      "R Anterior ventral insula (anterior pole)",
                                      "L Mid posterior insula",
                                      "R Mid posterior insula",
                                      "L Amygdala (BLBM)",
                                      "R Amygdala (CeME)",
                                      "L PAG",
                                      "R PAG"))

Pos$ROI <- plyr::revalue(Pos$ROI, c("L aMCC" = "L aMCC",
                              "R aMCC" = "R aMCC",
                              "L Anterior ventral insula (anterior pole)" = "L ventral anterior Insula",
                              "R Anterior ventral insula (anterior pole)" = "R ventral anterior Insula",
                              "L Mid posterior insula" = "L mid/posterior Insula",
                              "R Mid posterior insula" = "R mid/posterior Insula",
                              "L Amygdala (BLBM)" = "L basolateral Amygdala",
                              "R Amygdala (CeME)" = "R central/medial Amygdala",
                              "L PAG" = "L PAG",
                              "R PAG" = "R PAG"))

# Panel c: Negative going ROIs

# filter negative-going ROIs, but do not include any that are in Panel a
Neg <- df %>% group_by(ROI) %>%
  filter(max(Mean) < 0.25 & ROI != "R Amygdala (BLBM)" & ROI != "L Amygdala (BLBM)")

Neg$ROI <- factor(Neg$ROI)

# specify order
Neg$ROI <- factor(Neg$ROI,levels = c("vmPFC (from out lab)",
                                     "vmPFC (new 5mm sphere from Hartley paper)",
                                     "PCC (1)",
                                     "PCC (2)",
                                     "L Mid Hippocampus",
                                     "R Posterior Hippocampus"))

# rename to something plot appropriate
Neg$ROI <- plyr::revalue(Neg$ROI, c("vmPFC (from out lab)" = "anterior vmPFC",
                              "vmPFC (new 5mm sphere from Hartley paper)" = "posterior vmPFC",
                              "PCC (1)" = "PCC/Precuneus",
                              "PCC (2)" = "PCC",
                              "L Mid Hippocampus" = "L mid Hippocampus",
                              "R Posterior Hippocampus" = "R posterior Hippocampus"))


print("datasets ready, plotting commences")

# theme for Figure 4 plots
theme <- theme_bw() +
  theme(
    text = element_text(family = "Helvetica"),
    strip.background = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text.x = element_text(size = 6, color = "black", hjust = 0.5),
    aspect.ratio = .9,
    legend.position = 'none',
    axis.ticks = element_line(size = .4),
    axis.ticks.length=unit( .05, "cm"),
    axis.title = element_text(size = 6.5),
    axis.text.x = element_text(size = 4.5,margin = unit(c(0.04,0,0,0), "cm")),
    axis.text.y = element_text(size = 4.5, color = "black", margin = unit(c(0,0.05,0,0), "cm")),
    panel.border = element_rect(colour = "black", fill=NA, size=.5),
    plot.margin=grid::unit(c(0,.5,0,0), "cm"),
    panel.spacing.x=unit(0.2, "lines"),
    panel.spacing.y = unit(0.1,"lines"))



# panel a plotting
P.plot <- ggplot(P.ROIs, aes(x = Time, y = Mean, color = Group, label = P)) +
  geom_label(x = 9.2, y = .86, color = "black", size = 1.5, alpha = .65) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "dimgray", size = .5) +
  geom_point(size = .3) +
  geom_line(size = .7) +
  facet_wrap(~ROI, ncol = 4, scale = "free", labeller = label_wrap_gen(width = 20, multi_line = TRUE)) +
  scale_y_continuous(limits = c(-0.3, 1),
                     breaks = seq(from = -0.3, to = 1, by = .3),
                     labels = c("-0.30", "0.00", "0.30", "0.60","0.90"),
                     expand = y.axis.expand) +
  scale_x_continuous(breaks = x.axis.breaks,
                     labels = c("0.00", "3.00", "6.00", "9.00", "12.00"),
                     expand = x.axis.expand) +
  scale_colour_manual(values = c("#0043B7", "#F5652D"), labels = c("Controllable", "Uncontrollable")) +
  labs(x ="Time (seconds)", y = "% signal change", color = NULL, title =NULL) +
  theme


ggsave(here("plots","Fig4-PanelA_8ROIs.pdf"),dpi = 600, height = 57, width = 115, unit = "mm")

print("Panel 4a done.")

# panel b plotting

pos.plot <- ggplot(Pos, aes(x = Time, y = Mean, color = Group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "dimgray", size = .5) +
  geom_point(size = .3) +
  geom_line(size = .7) +
  facet_wrap(~ROI, scale = "free", ncol = 5, labeller = label_wrap_gen(width = 17, multi_line = TRUE)) +
  scale_x_continuous(breaks = x.axis.breaks,
                     labels = x.axis.labs,
                     expand = x.axis.expand) +
  scale_y_continuous(limits = c(-0.4, 1.3),
                     breaks = c( -0.4, 0, 0.4, 0.8, 1.2),
                     labels = c("-0.40", "0.00", "0.40", "0.80", "1.20"),
                     expand = y.axis.expand) +
  scale_colour_manual(values = c("#0043B7", "#F5652D"), labels = c("controllable", "uncontrollable")) +
  labs(x = "Time (seconds)", y = "% signal change", title = NULL, color = NULL) +
  theme

ggsave(here('plots','Fig4-PanelB_posROIs.pdf'), plot = pos.plot ,dpi = 600, height = 62, width = 150, unit = "mm")
print("Panel 4b done.")

# panel c plotting
neg.plot <- ggplot(Neg, aes(x = Time, y = Mean, color = Group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "dimgray", size = .5) +
  geom_point(size = .3) +
  geom_line(size = .7) +
  facet_wrap_custom(~ROI, scales = "free", ncol = 6, labeller = label_wrap_gen(width = 12, multi_line = TRUE),
                    scale_overrides = list(
                      scale_override(1, scale_y_continuous(limits = c(-0.6, 0.1),
                                                           breaks = c(-0.6, -0.4, -0.2, 0),
                                                           labels = c("-0.60", "-0.40", "-0.20", "0.00"))),
                      scale_override(2, scale_y_continuous(limits = c(-0.6, 0.1),
                                                           breaks = c(-0.6, -0.4, -0.2, 0),
                                                           labels = c("-0.60", "-0.40", "-0.20", "0.00"))),
                      scale_override(3, scale_y_continuous(limits = c(-0.6, 0.1),
                                                           breaks = c(-0.6, -0.4, -0.2, 0),
                                                           labels = c("-0.60", "-0.40", "-0.20", "0.00"))),
                      scale_override(4, scale_y_continuous(limits = c(-1, 0.1),
                                                           breaks= c(-1,-0.6,-0.2,0.2),
                                                           labels = c("-1.00", "-0.60", "-0.20", "0.20"))),
                      scale_override(5, scale_y_continuous(limits = c(-0.6, 0.1),
                                                           breaks = c(-0.6, -0.4, -0.2, 0),
                                                           labels = c("-0.60", "-0.40", "-0.20", "0.00"))),
                      scale_override(6, scale_y_continuous(limits = c(-0.6, 0.1),
                                                           breaks = c(-0.6, -0.4, -0.2, 0),
                                                           labels = c("-0.60", "-0.40", "-0.20", "0.00")))
                    )) +
  scale_x_continuous(breaks = x.axis.breaks,
                     labels = c("0.00", "3.00", "6.00", "9.00", "12.00"),
                     expand = x.axis.expand) +
  scale_y_continuous(breaks = c(-0.60, -.40, -.20, 0, .20),
                     labels = c("-0.60","-0.40","-0.20","0.00","0.20"),
                     expand = y.axis.expand) +
  scale_colour_manual(values = c("#0043B7", "#F5652D"), labels = c("Controllable", "Uncontrollable")) +
  theme +
  labs(x = "Time (seconds)", y = "% signal change", title = NULL, legend = NULL) +

  ggsave(here("plots","Fig4-PanelC_Neg_ROIs.pdf"),dpi = 600, height=38, width = 158, unit = "mm")

print("Panel 4c done.")

# theme for small (45mm) plots

theme <- theme_bw() +
  theme(
    text = element_text(family = "Helvetica"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(size = .25),
    aspect.ratio = .85,
    legend.position = 'none',
    plot.title = element_text(size = 6, hjust = 0.5, margin = unit(c(0,0,0.06,0),"cm")),
    axis.text.x = element_text(size = 5, margin = unit(c(0.07,0,0,0),"cm")),
    axis.title.x = element_text(size = 6),
    axis.text.y = element_text(size = 5, color = "black", margin = unit(c(0,0.05,0,0), "cm")),
    axis.title.y = element_text(size = 6),
    axis.ticks = element_line(size = .3),
    axis.ticks.length=unit(.05, "cm"),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.border = element_rect(colour = "black", fill=NA, size=.3))


print("Onto Figure 6!")
# Figure 6 Insula plots

# import data
l.pINS <- read.csv(here('data/stressor_response_estimates','Fig6-left_PI.1D'), sep = '', header = FALSE,
                   col.names = c("Controllable", "Uncontrollable"))


r.pINS <- read.csv(here('data/stressor_response_estimates','Fig6-right_PI.1D'), sep = '', header = FALSE,
                   col.names = c("Controllable", "Uncontrollable"))



# panel a - left posterior Insula ROI

l.pINS <- l.pINS %>%
  gather(Group, Mean) %>%
  mutate(Time = rep(l.pINS$Time <- seq(from = 0, to = 13.75, by = 1.25),2)) %>%
  ggplot(aes(x = Time, y = Mean, group = Group, color = Group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "dimgray", size = .75) +
  geom_point(size = .5, shape = 20) +
  geom_line(size = .75) +
  coord_cartesian(ylim = c(-0.2,.8)) +
  scale_x_continuous(breaks = x.axis.breaks,
                     labels = x.axis.labs,
                     expand = x.axis.expand) +
  scale_y_continuous(breaks = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8),
                     labels = c("-0.20", "0.00", "0.20", "0.40", "0.60","0.80"),
                     expand = y.axis.expand) +
  scale_colour_manual(values = c("#0043B7", "#F5652D"), labels = c("controllable", "uncontrollable")) +
  labs(x = 'Time (seconds)', y = '% signal change', color = NULL, title = "L posterior Insula") +
  theme

ggsave(here("plots","Fig6-PanelA_l-INS.pdf"),height = 45, width = 45, unit = "mm", dpi = 600)

print("Panel 6a is done.")

# panel b - right posterior Insula ROI


r.pINS <- r.pINS %>%
  gather(Group, Mean) %>%
  mutate(Time = rep(Time <- seq(from = 0, to = 13.75, by = 1.25),2)) %>%
  ggplot(aes(x = Time, y = Mean, group = Group, color = Group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "dimgray", size = .75) +
  geom_point(size = .5, shape = 20) +
  geom_line(size = .75) +
  coord_cartesian(ylim = c(-0.2,.8)) +
  scale_x_continuous(breaks = x.axis.breaks,
                     labels = x.axis.labs,
                     expand = x.axis.expand) +
  scale_y_continuous(breaks = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8),
                     labels = c("-0.20", "0.00", "0.20", "0.40", "0.60", "0.80"),
                     expand = y.axis.expand) +
  scale_colour_manual(values = c("#0043B7", "#F5652D")) +
  labs(x = 'Time (seconds)', y = '% signal change', color = NULL, title = "R posterior Insula") +
  theme_bw() +
  theme_bw() +
  theme

ggsave(here("plots","Fig6-PanelB_r-INS.pdf"),height = 45, width = 45, unit = "mm", dpi = 600)

print("Panel 6b is done.")
print("Onto figure 7!")

# Figure 7 plot

df.shock <- read.table(here('data/stressor_response_estimates','Fig7-PCC.txt'), fill = TRUE, sep=",", header = TRUE)
df.shock <- df.shock %>% gather(ROI, Response, 3:23)
df.shock %>% filter(ROI == "R.Posterior.cingulate..cortex") %>%  # Filter out r PCC & plot
  mutate(ROI = recode(ROI, R.Posterior.cingulate..cortex = "r PCC")) %>%
  ggplot(aes(x = Time, y = Response, color = Group)) +
  coord_cartesian(ylim = c(-0.6,.2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "dimgray", size = .75) +
  geom_point(size = .5, shape = 20) +
  geom_line(size = .75)  +
  scale_x_continuous(breaks = x.axis.breaks,
                     labels = x.axis.labs,
                     expand = x.axis.expand) +
  scale_y_continuous(breaks = c(-0.60, -.40, -.20, 0, .20),
                     labels = c("-0.60","-0.40","-0.20","0.00","0.20"),
                     expand = y.axis.expand) +
  scale_colour_manual(values = c("#0043B7", "#F5652D"), labels = c("Controllable", "Uncontrollable")) +
  labs(x = "Time (seconds)", y = "% signal change", color = NULL, title = "R PCC") +
  theme

ggsave(here("plots",'Fig7-PCC_plot.pdf'), height = 45, width = 45, unit = "mm", dpi = 600)

print("Figure 7 done, onto 8!")

# Figure 8: Exploratory analysis plot

explore <- read.table(here('data/stressor_response_estimates','Fig8-MFG.txt'), sep = ',', header = TRUE)

explore %>% filter(ROI == "SupMidGyrus") %>%
  mutate(ROI = recode(ROI, SupMidGyrus = "medial Frontal Gyrus")) %>%
  ggplot(aes(x = Time, y = Response, color = Group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "dimgray", size = .75) +
  geom_point(size = .5, shape = 20) +
  geom_line(size = .75)  +
  coord_cartesian(ylim = c(-0.6,.2)) +
  scale_x_continuous(breaks = x.axis.breaks,
                     labels = x.axis.labs,
                     expand = x.axis.expand) +
  scale_y_continuous(breaks = c(-0.60, -.40, -.20, 0, .20),
                     labels = c("-0.60","-0.40","-0.20","0.00","0.20"),
                     expand = y.axis.expand) +
  scale_colour_manual(values = c("#0043B7", "#F5652D")) +
  labs(x = 'Time (seconds)', y = '% signal change', color = NULL, title = "medial Frontal Gyrus") +
  theme


ggsave(here("plots","Fig8-MFG.pdf"), height = 45, width = 45, unit = "mm", dpi = 600)

print("Figure 8 done.")


print("All estimated response plots have been generated!")
