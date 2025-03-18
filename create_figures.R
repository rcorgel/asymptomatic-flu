################################################################################
# File Name: create_figures                                                    #
#                                                                              #
# Purpose:   Create figures from analyses.                                     #
# Project:   Asymptomatic Influenza Associations                               #
################################################################################

####################
# 1. SET-UP SCRIPT #
####################

# Start with a clear environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(lubridate)
library(assertr)
library(mvnfast)
library(corrplot)
library(MCMCpack)
library(car)
library(Hmisc)
library(cowplot)
library(reshape2)
library(moments)
library(sf)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Documents/asymptomatic-flu-proj/')

#####################
# 2. CREATE FIGURES #
#####################

# Load data for figures
regress_dat <- read.csv('out/regress_dat.csv')
data_frame_all <- read.csv('tmp/data_frame_all.csv')
combo_list <- read.csv('tmp/combo_list.csv')
traces_dist <- read.csv('tmp/traces_dist.csv')
data_frame_bio <- read.csv('tmp/data_frame_bio.csv')
data_frame_adm <- read.csv('tmp/data_frame_adm.csv')
bar_plot <- read.csv('tmp/bar_plot.csv')
waic_forward <- read.csv('tmp/waic_forward.csv')
data_frame_for <- read.csv('tmp/data_frame_for.csv')

# Figure 1 #

# Load maps
usa_albers_state <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_state.rds'))   # convert to sf
usa_albers_county <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_county.rds')) # convert to sf

# Join data
usa_albers_county <- left_join(usa_albers_county, regress_dat, by = c('GEOID' ='county_fips'))

# Create map
map <- ggplot() +
  geom_sf(data = usa_albers_county, aes(fill = percent_asymp_flu_scale/100, group = GEOID), 
          color= 'black', linewidth = 0.15) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), 
          fill = '#FFFFFF00', color= 'black', linewidth = 0.4) +
  scale_fill_gradient('Percent\n', low = "white", high = "#347dab") + 
  ggtitle('Spatial Patterns') + 
  theme_void() + theme(legend.position = 'right',
                       legend.text = element_text(size = 16),
                       legend.title = element_text(size = 20),
                       axis.title = element_text(size=20),
                       legend.margin=margin(),
                       strip.background = element_blank(),
                       legend.spacing.y = unit(0.25, 'cm'),
                       legend.key.height = unit(1.5, 'cm'),
                       strip.text = element_text(size = 16),
                       plot.title = element_text(size=20, hjust = 0.5))

# Create density
density <- ggplot(data = regress_dat) +
  geom_density(aes(percent_asymp_flu_scale /100), color = '#347dab',
               fill = '#347dab', alpha = 0.65) + 
  ggtitle('County Distribution') + 
  geom_vline(xintercept=mean(regress_dat$percent_asymp_flu_scale /100), 
             lty=2, linewidth = 1) + 
  ylab("Density") + xlab("Percent") + xlim(0, 1) + 
  theme_minimal() +  # use a white background
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 16),
        axis.text = element_text(size=16),
        axis.text.y = element_blank(),
        axis.title = element_text(size=20),
        legend.position = "bottom",
        legend.box="vertical",
        legend.margin=margin(),
        strip.background = element_blank(),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.key.size = unit(1, 'cm'),
        strip.text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0.5))

# Combine figure and save
figure_1 <- cowplot::plot_grid(map, density,
                               nrow = 2, rel_widths = c(1, 1),
                               labels = c('', ''),
                               label_size = 26, hjust = 0) 

ggsave('./figs/figure_1_p.jpg', plot = figure_1, height = 7, width = 6)

# Figure 2 #

bio_plot <- ggplot(data=data_frame_bio[data_frame_bio$Variables != 'Intercept', ], 
                   aes(x=Variables, y= `50%`, ymin= `2.5%`, ymax= `97.5%`)) +
  geom_pointrange(aes(), color = '#5a9374', size = 1, linewidth = 1) + 
  ggtitle('Immune/Biology Model') + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Median (2.5%-97.5%)") +
  ylim(-1.6, 1.6) + 
  theme_minimal() +  # use a white background
  theme(legend.position = 'bottom') +
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 20),
        axis.text = element_text(size=16),
        axis.title = element_text(size=20),
        legend.position = "bottom",
        legend.box="vertical",
        legend.margin=margin(),
        strip.background = element_blank(),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.key.size = unit(1, 'cm'),
        strip.text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0.5))


adm_plot <- ggplot(data=data_frame_adm[data_frame_adm$Variables != 'Intercept', ], 
                   aes(x=Variables, y= `50%`, ymin= `2.5%`, ymax= `97.5%`)) +
  geom_pointrange(aes(), color = '#9086ba', size = 1, linewidth = 1) + 
  ggtitle('Burden/SES Model') + 
  ylim(-1.6, 1.6) + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Median (2.5%-97.5%)") +
  theme_minimal() +  # use a white background
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 20),
        axis.text = element_text(size=16),
        axis.title = element_text(size=20),
        legend.position = "bottom",
        legend.box="vertical",
        legend.margin=margin(),
        strip.background = element_blank(),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.key.size = unit(1, 'cm'),
        strip.text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0.5))

bar_plot <- ggplot(data=bar_data, aes(x=Model, y=WAIC, fill = Model)) + 
  geom_bar(stat="identity", position="identity", color = 'black', 
           width = 0.5, alpha = 1) + 
  theme_classic() +
  scale_fill_manual('', values = c('#5a9374', '#9086ba'),
                    breaks = c('Immune/Biology', 'Burden/SES')) + 
  ylab('WAIC') + xlab("") + ggtitle('Model Comparison') +
  coord_cartesian(ylim = c(24000, 25000)) +
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 20),
        axis.text = element_text(size=16),
        axis.title = element_text(size=20),
        legend.position = 'none',
        legend.box="vertical",
        legend.margin=margin(),
        strip.background = element_blank(),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.key.size = unit(1, 'cm'),
        strip.text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0.5))

# Combine figure and save
figure_2 <- cowplot::plot_grid(adm_plot, bio_plot, ggplot + theme_void(), bar_plot,
                               nrow = 1, rel_widths = c(1, 1, 0.1, 1),
                               labels = c('', '', ''),
                               label_size = 26, hjust = 0) 

ggsave('./figs/figure_2.jpg', plot = figure_2, height = 5, width = 20)

figure_2 <- cowplot::plot_grid(adm_plot, bio_plot,
                               nrow = 2, rel_widths = c(1, 1),
                               labels = c('', ''),
                               label_size = 26, hjust = 0) 

ggsave('./figs/figure_2_p.jpg', plot = figure_2, height = 7, width = 6)

# Figure 3 #

forward_line <- ggplot() + geom_line(aes(y = waic_forward$all_waic, x = 1:15), color = '#c36272', 
                                     linewidth = 2) + 
  xlab('Step Number') + ggtitle('WAIC by Forward Selection Step') +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14)) +
  ylab('WAIC') + theme_minimal() +  # use a white background
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 20),
        axis.text = element_text(size=16),
        axis.title = element_text(size=20),
        legend.position = "bottom",
        legend.box="vertical",
        legend.margin=margin(),
        strip.background = element_blank(),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.key.size = unit(1, 'cm'),
        strip.text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0.5))

forward_plot <- ggplot(data=data_frame_for[data_frame_for$Variables != 'Intercept', ], 
                       aes(x=Variables, y= `50%`, ymin= `2.5%`, ymax= `97.5%`)) +
  geom_pointrange(aes(), color = '#c36272', size = 1, linewidth = 1) + 
  ggtitle('Forward Model') + 
  ylim(-1.6, 1.6) + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Median (2.5%-97.5%)") +
  theme_minimal() +  # use a white background
  theme(legend.position = 'bottom') +
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 20),
        axis.text = element_text(size=16),
        axis.title = element_text(size=20),
        legend.position = "bottom",
        legend.box="vertical",
        legend.margin=margin(),
        strip.background = element_blank(),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.key.size = unit(1, 'cm'),
        strip.text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0.5))

# Combine figure and save
figure_3 <- cowplot::plot_grid(forward_line, forward_plot,
                               nrow = 1, rel_widths = c(1, 1),
                               labels = c('', ''),
                               label_size = 26, hjust = 0) 

ggsave('./figs/figure_3.jpg', plot = figure_3, height = 5, width = 20)

figure_3 <- cowplot::plot_grid(forward_line, forward_plot,
                               nrow = 2, rel_widths = c(1, 1),
                               labels = c('', ''),
                               label_size = 26, hjust = 0) 

ggsave('./figs/figure_3_p.jpg', plot = figure_3, height = 7, width = 6)

# Figure 4 #

density <- ggplot(data = combo_list) +
  geom_density(aes(WAIC), color = '#C7A939',
               fill = '#C7A939', alpha = 0.65) + 
  ggtitle('WAIC Distribution') +
  ylab("Density") + xlab("WAIC") + xlim(24450, 25800) +
  theme_minimal() +  # use a white background
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 16),
        axis.text = element_text(size=16),
        axis.text.y = element_blank(),
        axis.title = element_text(size=20),
        legend.position = "bottom",
        legend.box="vertical",
        legend.margin=margin(),
        strip.background = element_blank(),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.key.size = unit(1, 'cm'),
        strip.text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0.5))

all_plot <- ggplot(data=data_frame_all[data_frame_all$Variables != 'Intercept', ], 
                   aes(x=Variables, y= `50%`, ymin= `2.5%`, ymax= `97.5%`)) +
  geom_pointrange(aes(), color = '#C7A939', size = 1, linewidth = 1) + 
  ggtitle('All Combos Model') + 
  ylim(-1.6, 1.6) + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Median (2.5%-97.5%)") +
  theme_minimal() +  # use a white background
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 20),
        axis.text = element_text(size=16),
        axis.title = element_text(size=20),
        legend.position = "bottom",
        legend.box="vertical",
        legend.margin=margin(),
        strip.background = element_blank(),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.key.size = unit(1, 'cm'),
        strip.text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0.5))

# Combine figure and save
figure_4 <- cowplot::plot_grid(density, all_plot,
                               nrow = 2, rel_widths = c(1, 1),
                               labels = c('', ''),
                               label_size = 26, hjust = 0) 

ggsave('./figs/figure_4_p.jpg', plot = figure_4, height = 7, width = 6)

# Figure A #

traces_plot <- ggplot(data = regress_dat) +
  geom_density(data = traces_dist, 
               aes(value /100, group = variable), color = 'gray', 
               linewidth = 0.05, alpha = 0.01) + 
  ggtitle('Predicted County Distributions') + 
  geom_vline(xintercept=mean(regress_dat$percent_asymp_flu_scale /100), 
             lty=2, linewidth = 1) + 
  ylab("Density") + xlab("Percent") + xlim(0, 1) + 
  theme_minimal() +  # use a white background
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 16),
        axis.text = element_text(size=16),
        axis.text.y = element_blank(),
        axis.title = element_text(size=20),
        legend.position = "bottom",
        legend.box="vertical",
        legend.margin=margin(),
        strip.background = element_blank(),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.key.size = unit(1, 'cm'),
        strip.text = element_text(size = 16),
        plot.title = element_text(size=20, hjust = 0.5))

# Save
ggsave('./figs/figure_a.jpg', plot = traces_plot, height = 4, width = 8)

################################################################################
################################################################################


