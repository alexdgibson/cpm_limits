# figures that will be included in the final article
# included in the protocol
# created March 2025

# load libraries
library(ggplot2)
library(cowplot)
library(dplyr)


### Plots for the discrimination ###
# create the forest plots for the cumulative meta-analyses
# cardiology 
cardio_forest <- ggplot()+
  # geom_point(aes(x = , y = ), shape = 15)+ # commented out for the protocol
  # geom_linerange(aes(y = ,
  #                    xmin = ,
  #                    xmax = ))+
  # geom_linerange(data = ,aes(y = 0,
  #                                          xmin = ,
  #                                          xmax = ))+
  theme_classic()+
  labs(y = "Study")+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  geom_vline(aes(xintercept = 0.5), linetype = "longdash")+
  scale_x_continuous(limits = c(0.4,1))+
  scale_y_continuous(limits = c(2,100), labels = c("100", "75", "50", "25", "2"))


# obstetrics and gynaecology 
obgyn_forest <- ggplot()+
  # geom_point(aes(x = , y = ), shape = 15)+ # commented out for the protocol
  # geom_linerange(aes(y = ,
  #                    xmin = ,
  #                    xmax = ))+
  # geom_linerange(data = ,aes(y = 0,
  #                                          xmin = ,
  #                                          xmax = ))+
  theme_classic()+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  geom_vline(aes(xintercept = 0.5), linetype = "longdash")+
  scale_x_continuous(limits = c(0.4,1))+
  scale_y_continuous(limits = c(2,100), labels = c("100", "75", "50", "25", "2"))


# gasteroenterology 
gastro_forest <- ggplot()+
  # geom_point(aes(x = , y = ), shape = 15)+ # commented out for the protocol
  # geom_linerange(aes(y = ,
  #                    xmin = ,
  #                    xmax = ))+
  # geom_linerange(data = ,aes(y = 0,
  #                                          xmin = ,
  #                                          xmax = ))+
  theme_classic()+
  labs(x = "AUC Summary Estimate",
       y = "Study")+
  theme(text = element_text(size = 12),
        panel.background = element_rect(colour = "black", fill = NA, linewidth = 1))+
  geom_vline(aes(xintercept = 0.5), linetype = "longdash")+
  scale_x_continuous(limits = c(0.4,1))+
  scale_y_continuous(limits = c(2,100), labels = c("100", "75", "50", "25", "2"))


# neurology 
neuro_forest <- ggplot()+
  # geom_point(aes(x = , y = ), shape = 15)+ # commented out for the protocol
  # geom_linerange(aes(y = ,
  #                    xmin = ,
  #                    xmax = ))+
  # geom_linerange(data = ,aes(y = 0,
  #                                          xmin = ,
  #                                          xmax = ))+
  theme_classic()+
  labs(x = "AUC Summary Estimate")+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  geom_vline(aes(xintercept = 0.5), linetype = "longdash")+
  scale_x_continuous(limits = c(0.4,1))+
  scale_y_continuous(limits = c(2,100), labels = c("100", "75", "50", "25", "2"))


# join each plot together with cowplot
plot_grid(cardio_forest, obgyn_forest, gastro_forest, neuro_forest,
          labels = c("A","B", "C", "D"),
          label_size = 10)

#save the plot
ggsave(file = "03_figures/protocol_exmaple_discr_forest_layout.jpg",
       width = 8,
       height = 6,
       dpi = 300)







# create the plots of the change in meta-analysis
# cardiology
cardio_forest_change <- ggplot()+
# geom_line(linewidth = 1)+ # commented out for the protocol
  geom_hline(aes(yintercept = 0.0), linetype = "longdash")+
theme_classic()+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  labs(y = "Change in AUC Summary Estimate")+
  scale_y_continuous(limits = c(-0.1, 0.1))+
  scale_x_continuous(limits = c(0,100))
  

# obstetrics and gynocology
obgyn_forest_change <- ggplot()+
  # geom_line(linewidth = 1)+ # commented out for the protocol
  geom_hline(aes(yintercept = 0.0), linetype = "longdash")+
  theme_classic()+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  scale_y_continuous(limits = c(-0.1, 0.1))+
  scale_x_continuous(limits = c(0,100))


# gasteroenterfolgy
gastro_forest_change <- ggplot()+
  # geom_line(linewidth = 1)+ # commented out for the protocol
  geom_hline(aes(yintercept = 0.0), linetype = "longdash")+
  theme_classic()+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  labs(x = "Meta-Analysis Step",
       y = "Change in AUC Summary Estimate")+
  scale_y_continuous(limits = c(-0.1, 0.1))+
  scale_x_continuous(limits = c(0,100))


# neurology
neuro_forest_change <- ggplot()+
  # geom_line(linewidth = 1)+ # commented out for the protocol
  geom_hline(aes(yintercept = 0.0), linetype = "longdash")+
  theme_classic()+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  labs(x = "Meta-Analysis Step")+
  scale_y_continuous(limits = c(-0.1, 0.1))+
  scale_x_continuous(limits = c(0,100))


#join the plots together with cowplot
plot_grid(cardio_forest_change, obgyn_forest_change, gastro_forest_change, neuro_forest_change,
          labels = c("A","B", "C", "D"),
          label_size = 10)


# save the plot
ggsave(file = "03_figures/protocol_exmaple_discr_forest_change_layout.jpg",
       width = 8,
       height = 6,
       dpi = 300)




# line plot of the change in the summary estimate with simulated study outcomes and sample sizes
# cardiology
cardio_sim_line <- ggplot()+
#geom_line(aes(x = cstat, y = Mean, colour = as.factor(sample)))+
theme_classic()+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
labs(y = "AUC Summary Estiamte",
     colour = "Sample Size")+
  scale_y_continuous(limits = c(0.5,1))+
  scale_x_continuous(limits = c(0.5,1))


# obstetrics and gynecology
obgyn_sim_line <- ggplot()+
  #geom_line(aes(x = cstat, y = Mean, colour = as.factor(sample)))+
  theme_classic()+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  labs(colour = "Sample Size")+
  scale_y_continuous(limits = c(0.5,1))+
  scale_x_continuous(limits = c(0.5,1))


# gastroenterology
gastro_sim_line <- ggplot()+
  #geom_line(aes(x = cstat, y = Mean, colour = as.factor(sample)))+
  theme_classic()+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  labs(x = "Simulated Study AUC-Value",
       y = "AUC Summary Estiamte",
       colour = "Sample Size")+
  scale_y_continuous(limits = c(0.5,1))+
  scale_x_continuous(limits = c(0.5,1))


# neurology
neuro_sim_line <- ggplot()+
  #geom_line(aes(x = cstat, y = Mean, colour = as.factor(sample)))+
  theme_classic()+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  labs(x = "Simulated Study AUC-Value",
       colour = "Sample Size")+
  scale_y_continuous(limits = c(0.5,1))+
  scale_x_continuous(limits = c(0.5,1))


# join them together in a plot grid
plot_grid(cardio_sim_line, obgyn_sim_line, gastro_sim_line, neuro_sim_line,
          labels = c("A","B", "C", "D"),
          label_size = 10)


# save the plot
ggsave(file = "03_figures/protocol_exmaple_discr_sim_change_layout.jpg",
       width = 8,
       height = 6,
       dpi = 300)








### Plots for the calibration ###

# create the forest plots for the cumulative meta-analyses
# cardiology 
cardio_forest_cal <- ggplot()+
  # geom_point(aes(x = , y = ), shape = 15)+ # commented out for the protocol
  # geom_linerange(aes(y = ,
  #                    xmin = ,
  #                    xmax = ))+
  # geom_linerange(data = ,aes(y = 0,
  #                                          xmin = ,
  #                                          xmax = ))+
  theme_classic()+
  labs(y = "Study")+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  geom_vline(aes(xintercept = 0.5), linetype = "longdash")+
  scale_x_continuous(limits = c(0,3))+
  scale_y_continuous(limits = c(2,100), labels = c("100", "75", "50", "25", "2"))


# obstetrics and gynaecology 
obgyn_forest_cal <- ggplot()+
  # geom_point(aes(x = , y = ), shape = 15)+ # commented out for the protocol
  # geom_linerange(aes(y = ,
  #                    xmin = ,
  #                    xmax = ))+
  # geom_linerange(data = ,aes(y = 0,
  #                                          xmin = ,
  #                                          xmax = ))+
  theme_classic()+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  geom_vline(aes(xintercept = 0.5), linetype = "longdash")+
  scale_x_continuous(limits = c(0,3))+
  scale_y_continuous(limits = c(2,100), labels = c("100", "75", "50", "25", "2"))


# gasteroenterology 
gastro_forest_cal <- ggplot()+
  # geom_point(aes(x = , y = ), shape = 15)+ # commented out for the protocol
  # geom_linerange(aes(y = ,
  #                    xmin = ,
  #                    xmax = ))+
  # geom_linerange(data = ,aes(y = 0,
  #                                          xmin = ,
  #                                          xmax = ))+
  theme_classic()+
  labs(x = "O:E Summary Estimate",
       y = "Study")+
  theme(text = element_text(size = 12),
        panel.background = element_rect(colour = "black", fill = NA, linewidth = 1))+
  geom_vline(aes(xintercept = 0.5), linetype = "longdash")+
  scale_x_continuous(limits = c(0,3))+
  scale_y_continuous(limits = c(2,100), labels = c("100", "75", "50", "25", "2"))


# neurology 
neuro_forest_cal <- ggplot()+
  # geom_point(aes(x = , y = ), shape = 15)+ # commented out for the protocol
  # geom_linerange(aes(y = ,
  #                    xmin = ,
  #                    xmax = ))+
  # geom_linerange(data = ,aes(y = 0,
  #                                          xmin = ,
  #                                          xmax = ))+
  theme_classic()+
  labs(x = "O:E Summary Estimate")+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  geom_vline(aes(xintercept = 0.5), linetype = "longdash")+
  scale_x_continuous(limits = c(0,3))+
  scale_y_continuous(limits = c(2,100), labels = c("100", "75", "50", "25", "2"))


# join each plot together with cowplot
plot_grid(cardio_forest_cal, obgyn_forest_cal, gastro_forest_cal, neuro_forest_cal,
          labels = c("A","B", "C", "D"),
          label_size = 10)

# save the plot
ggsave(file = "03_figures/protocol_exmaple_cal_forest_layout.jpg",
       width = 8,
       height = 6,
       dpi = 300)




# create the plots of the change in meta-analysis
# cardiology
cardio_forest_change_cal <- ggplot()+
  # geom_line(linewidth = 1)+ # commented out for the protocol
  geom_hline(aes(yintercept = 0.0), linetype = "longdash")+
  theme_classic()+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  labs(y = "Change in O:E Summary Estimate")+
  scale_y_continuous(limits = c(-0.2, 0.2))+
  scale_x_continuous(limits = c(0,100))


# obstetrics and gynocology
obgyn_forest_change_cal <- ggplot()+
  # geom_line(linewidth = 1)+ # commented out for the protocol
  geom_hline(aes(yintercept = 0.0), linetype = "longdash")+
  theme_classic()+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  scale_y_continuous(limits = c(-0.1, 0.1))+
  scale_x_continuous(limits = c(0,100))


# gastroenterology
gastro_forest_change_cal <- ggplot()+
  # geom_line(linewidth = 1)+ # commented out for the protocol
  geom_hline(aes(yintercept = 0.0), linetype = "longdash")+
  theme_classic()+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  labs(x = "Meta-Analysis Step",
       y = "Change in O:E Summary Estimate")+
  scale_y_continuous(limits = c(-0.1, 0.1))+
  scale_x_continuous(limits = c(0,100))


# neurology
neuro_forest_change_cal <- ggplot()+
  # geom_line(linewidth = 1)+ # commented out for the protocol
  geom_hline(aes(yintercept = 0.0), linetype = "longdash")+
  theme_classic()+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  labs(x = "Meta-Analysis Step")+
  scale_y_continuous(limits = c(-0.1, 0.1))+
  scale_x_continuous(limits = c(0,100))


#join the plots together with cowplot
plot_grid(cardio_forest_change_cal, obgyn_forest_change_cal, gastro_forest_change_cal, neuro_forest_change_cal,
          labels = c("A","B", "C", "D"),
          label_size = 10)


# save the plot
ggsave(file = "03_figures/protocol_exmaple_cal_forest_change_layout.jpg",
       width = 8,
       height = 6,
       dpi = 300)






# line plot of the change in the summary estimate with simulated study outcomes and sample sizes
# cardiology
cardio_sim_line_cal <- ggplot()+
  #geom_line(aes(x = cstat, y = Mean, colour = as.factor(sample)))+
  theme_classic()+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  labs(y = "O:E Summary Estiamte",
       colour = "Sample Size")+
  scale_y_continuous(limits = c(0.5,1))+
  scale_x_continuous(limits = c(0.5,1))


# obstetrics and gynecology
obgyn_sim_line_cal <- ggplot()+
  #geom_line(aes(x = cstat, y = Mean, colour = as.factor(sample)))+
  theme_classic()+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  labs(colour = "Sample Size")+
  scale_y_continuous(limits = c(0.5,1))+
  scale_x_continuous(limits = c(0.5,1))


# gastroenterology
gastro_sim_line_cal <- ggplot()+
  #geom_line(aes(x = cstat, y = Mean, colour = as.factor(sample)))+
  theme_classic()+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  labs(x = "Simulated Study O:E Ratio",
       y = "O:E Summary Estiamte",
       colour = "Sample Size")+
  scale_y_continuous(limits = c(0.5,1))+
  scale_x_continuous(limits = c(0.5,1))


# neurology
neuro_sim_line_cal <- ggplot()+
  #geom_line(aes(x = cstat, y = Mean, colour = as.factor(sample)))+
  theme_classic()+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  labs(x = "Simulated Study O:E Ratio",
       colour = "Sample Size")+
  scale_y_continuous(limits = c(0.5,1))+
  scale_x_continuous(limits = c(0.5,1))


# join them together in a plot grid
plot_grid(cardio_sim_line_cal, obgyn_sim_line_cal, gastro_sim_line_cal, neuro_sim_line_cal,
          labels = c("A","B", "C", "D"),
          label_size = 10)


# save the plot
ggsave(file = "03_figures/protocol_exmaple_cal_sim_change_layout.jpg",
       width = 8,
       height = 6,
       dpi = 300)






# combine the plots together that would show both the discrimination and calibration by specialty
plot_grid(cardio_forest,cardio_forest_change,cardio_sim_line,
          cardio_forest_cal, cardio_forest_change_cal, cardio_sim_line_cal,
          labels = c("A", "B", "C", "D", "E", "F"),
          label_size = 10)

ggsave(file = "03_figures/protocol_exmaple_figure_layout.jpg",
       width = 12,
       height = 8,
       dpi = 300)









### create the plot for the cumulative research activity over time 

ggplot()+
  #geom_line(aes(x = cstat, y = Mean, colour = as.factor(sample)))+
  theme_classic()+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  labs(y = "Cummulative Articles",
       x = "Date")+
  scale_y_continuous(limits = c(0,150))+
  scale_x_continuous(limits = c(2000,2025))

# save the plot
ggsave(file = "03_figures/protocol_exmaple_cummulative_activity.jpg",
       width = 4,
       height = 3,
       dpi = 300)



# create a plot for the number of predictors vs the performance of the AUC

#cardiology
cardio_predictor <- ggplot()+
  #geom_line(aes(x = cstat, y = Mean, colour = as.factor(sample)))+
  theme_classic()+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  labs(y = "AUC Value")+
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(limits = c(1,10))

# obstetrics and gynecology
obgyn_predictor <- ggplot()+
  #geom_line(aes(x = cstat, y = Mean, colour = as.factor(sample)))+
  theme_classic()+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(limits = c(1,10))

# gasteroenterology
gastro_predictor <- ggplot()+
  #geom_line(aes(x = cstat, y = Mean, colour = as.factor(sample)))+
  theme_classic()+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  labs(y = "AUC Value",
       x = "Number of Predictors")+
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(limits = c(1,10))

# neurology
neuro_predictor <- ggplot()+
  #geom_line(aes(x = cstat, y = Mean, colour = as.factor(sample)))+
  theme_classic()+
  theme(text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  labs(x = "Number of Predictors")+
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(limits = c(1,10))



# join the predictor plots together
plot_grid(cardio_predictor, obgyn_predictor, gastro_predictor, neuro_predictor,
          labels = c("A", "B", "C", "D"),
          label_size = 10)


# save the paneled plot together 
ggsave(file = "03_figures/protocol_exmaple_predictors.jpg",
       width = 8,
       height = 6,
       dpi = 300)
