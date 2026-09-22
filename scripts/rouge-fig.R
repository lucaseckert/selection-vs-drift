#### ROUGE FIGURE ####

## constants
source(here::here("scripts/global.R"))

#### DEMOGRAPHY ####

Ni <- 5000  # initial size
Nb <- 50    # bottleneck size
Ng <- 1     # generations at bottleneck
Nt <- 4     # total generations

## recovery
N_recovery_gens <- Nt - Ng
growth_factor <- (Ni / Nb)^(1 / N_recovery_gens)
recovery_N <- round(Nb * (growth_factor ^ seq_len(N_recovery_gens)))
recovery_N[N_recovery_gens] <- Ni

## pop trajectory
pop_size <- c(
  Ni,           # starting point
  rep(Nb, Ng),  # bottleneck phase
  recovery_N    # recovery phase
)

df <- data.frame(
  N = c(rep(5000,7),pop_size, rep(5000,4)),
  year = 2007:2022
)

pop_plot<-ggplot(df, aes(x = year, y = N)) +
  annotate("rect",
           xmin = 2014.5, xmax = 2015.5,
           ymin = -Inf, ymax = Inf,
           alpha = 0.5, fill = "grey50") +
  annotate("rect",
           xmin = 2015.5, xmax = 2017.5,
           ymin = -Inf, ymax = Inf,
           alpha = 0.5, fill = "grey75") +
  geom_line(linewidth = 0.5) +
  ylim(0, 5500) +
  labs(title="Population Size")+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size=16, hjust = 0.5))

#### PHENOTYPE ####

## data
data<-read.csv("tidy-data/rouge_morpho_data.csv")

## pic
img_plated <- readPNG("figures/pics/plated.png")
grob_plated <- rasterGrob(img_plated, interpolate = TRUE)
img_unplated <- readPNG("figures/pics/unplated.png")
grob_unplated <- rasterGrob(img_unplated, interpolate = TRUE)

## plate year
plates<-data %>% mutate(plated=if_else(lp_no>0,T,F)) %>% 
  filter(site=="lake") %>% 
  mutate(year=if_else(year==2012,2013,year)) %>% 
  group_by(year) %>% 
  summarize(plated_freq = mean(plated),
            unplated_freq = 1 - plated_freq,
            n = n()) %>%
  ungroup()

## plotting
plot<-ggplot(plates, aes(x=year, y=unplated_freq)) +
  annotation_custom(
    grob_plated,
    xmin = 2007.5, xmax = 2013,
    ymin = 0.225, ymax=0.325)+
  annotation_custom(
    grob_unplated,
    xmin = 2007.5, xmax = 2013,
    ymin = 0.325, ymax=0.425)+
  annotate("rect",
           xmin = 2014.5, xmax = 2015.5,
           ymin = -Inf, ymax = Inf,
           alpha = 0.5, fill = "grey50") +
  annotate("rect",
           xmin = 2015.5, xmax = 2017.5,
           ymin = -Inf, ymax = Inf,
           alpha = 0.5, fill = "grey75") +
  # annotate("segment",
  #          x = 2013.5, xend = 2014.25,
  #          y = 0.15, yend = 0.15,
  #          arrow = arrow(length = unit(0.25, "cm")),
  #          linewidth = 0.8) +
  annotate("text",
           x = 2015, y = 0.55,
           label = "Drought",
           size=6,
           hjust = 1,
           angle = 90) +
  # annotate("segment",
  #          x = 2018.5, xend = 2017.75,
  #          y = 0.25, yend = 0.25,
  #          arrow = arrow(length = unit(0.25, "cm")),
  #          linewidth = 0.8) +
  annotate("text",
           x = 2016.5, y = 0.55,
           label = "Recovery",
           size=6,
           hjust = 1,
           angle = 90) +
  geom_line(linetype="dashed")+
  geom_point(size=4)+
  labs(x="Year", y="Unplated Phenotype Frequency")+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=16),
        axis.text = element_text(size=14))+
  inset_element(pop_plot, left = 0.7, bottom = 0.45, right = 0.975, top = 0.95)

## saving
save_fig(
  plot = plot,
  filename = "figures/timeline.tiff",
  width = 0.6,
  height = 0.4)
