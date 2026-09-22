#### SIMULTIONS ####

## constants
source(here::here("scripts/global.R"))

#### FIXED PARAMETERS ####
N_initial  <- 5000
N_final    <- 5000
N_years    <- 3
P_phenotype_initial  <- 0.49
P_phenotype_observed <- 0.87
N_replicates <- 10000

#### VARIABLE PARAMETERS ####
N_bottleneck_vals <- seq(5, 50, 5)
s_vals            <- seq(0, 0.75, 0.05)
dominance_vals    <- c("dominant", "recessive")
selection_scenarios <- c(FALSE, TRUE)  # FALSE = bottleneck only, TRUE = continuous

#### FUNCTION: phenotype -> allele frequency ####
get_p_from_pheno <- function(P_pheno, dominance_type) {
  if (dominance_type == "dominant") {
    1 - sqrt(1 - P_pheno)
  } else {
    sqrt(P_pheno)
  }
}

#### SIMULATION FUNCTION ####
simulate_bottleneck <- function(p0, Nb, Ni, Nt, s, dom,
                                selection_in_recovery = FALSE,
                                n_reps = 10000) {
  
  p <- rep(p0, n_reps)
  
  # fitness
  if (dom == "dominant") {
    wAA <- 1; wAa <- 1; waa <- 1 - s
  } else {
    wAA <- 1; wAa <- 1 - s; waa <- 1 - s
  }
  
  #### BOTTLENECK (1 generation) ####
  q     <- 1 - p
  w_bar <- (p^2 * wAA) + (2 * p * q * wAa) + (q^2 * waa)
  p     <- (p^2 * wAA + p * q * wAa) / w_bar
  
  p <- rbinom(n_reps, size = 2 * Nb, prob = p) / (2 * Nb)
  
  #### RECOVERY ####
  N_recovery_gens <- Nt - 1
  
  if (N_recovery_gens > 0) {
    growth_factor <- (Ni / Nb)^(1 / N_recovery_gens)
    recovery_N <- round(Nb * (growth_factor ^ seq_len(N_recovery_gens)))
    recovery_N[N_recovery_gens] <- Ni
    
    for (N_curr in recovery_N) {
      
      if (selection_in_recovery) {
        q     <- 1 - p
        w_bar <- (p^2 * wAA) + (2 * p * q * wAa) + (q^2 * waa)
        p     <- (p^2 * wAA + p * q * wAa) / w_bar
      }
      
      p <- rbinom(n_reps, size = 2 * N_curr, prob = p) / (2 * N_curr)
    }
  }
  
  return(p)
}

#### PARAM GRID ####
params <- expand.grid(
  N_bottleneck = N_bottleneck_vals,
  s = s_vals,
  dominance = dominance_vals,
  selection_in_recovery = selection_scenarios,
  stringsAsFactors = FALSE
)

#### RUN SIMULATIONS ####
set.seed(42)

results <- params %>%
  mutate(
    p0 = map_dbl(dominance, ~ get_p_from_pheno(P_phenotype_initial, .x)),
    
    final_p = pmap(
      list(p0, N_bottleneck, s, dominance, selection_in_recovery),
      function(p, nb, sel, dom, sel_rec)
        simulate_bottleneck(
          p0 = p,
          Nb = nb,
          Ni = N_final,
          Nt = N_years,
          s = sel,
          dom = dom,
          selection_in_recovery = sel_rec,
          n_reps = N_replicates
        )
    ),
    
    final_pheno = map2(final_p, dominance, function(p_vec, dom) {
      if (dom == "dominant") {
        1 - (1 - p_vec)^2
      } else {
        p_vec^2
      }
    }),
    
    p_success = map_dbl(final_pheno, ~ mean(.x >= P_phenotype_observed)),
    mean_final_pheno = map_dbl(final_pheno, mean)
  ) %>%
  select(N_bottleneck, s, dominance, selection_in_recovery, p_success, mean_final_pheno)

#### PLOT ####

# labels for the left plot
lab_success <- tibble(
  dominance = c("dominant", "recessive"),
  label = c("A", "B"),
  x = 45,
  y = 0.7
)

# labels for the right plot
lab_pheno <- tibble(
  dominance = c("dominant", "recessive"),
  label = c("C", "D"),
  x = 45,
  y = 0.7
)

success<-results %>% filter(selection_in_recovery == FALSE) %>%
  mutate(label="Success Probability") %>% 
  ggplot(aes(x=N_bottleneck, y=s, z=p_success))+
  geom_contour_fill(bins=20) +
  scale_fill_viridis_c(option = "magma",
                       guide = guide_colorbar(
                         barheight = unit(0.5, "cm"),
                         barwidth = unit(5, "cm")))+
  labs(x = expression("Bottleneck Size (" * N[e] * ")"),
       y = expression("Selection Coefficient (" * italic(s) * ")"),
       fill = "Success Probability",
       title = "Success Probability") +
  ne_s_layer()+
  facet_grid(rows=vars(dominance))+
  geom_text(
    data = lab_success,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1,
    fontface = "bold", size = 6, color="white")+
  theme(panel.grid = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        plot.title = element_text(size=14, hjust=0.5),
        legend.text = element_text(size=10),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=14, face="plain"),
        legend.ticks = element_line(color="black"),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0))


pheno<-results %>% filter(selection_in_recovery == FALSE) %>%
  ggplot(aes(x=N_bottleneck, y=s, z=mean_final_pheno-0.49))+
  geom_contour_fill(bins=20) +
  scale_fill_viridis_c(option = "mako",
                       guide = guide_colorbar(
                         barheight = unit(0.5, "cm"),
                         barwidth = unit(5, "cm")))+
  labs(x = expression("Bottleneck Size (" * N[e] * ")"),
       y = expression("Selection Coefficient (" * italic(s) * ")"),
       fill = expression(Delta~"Phenotype Frequency"),
       title = expression(Delta~"Phenotype Frequency")) +
  ne_s_layer()+
  facet_grid(rows=vars(dominance))+
  geom_text(
    data = lab_pheno,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1,
    fontface = "bold", size = 6, color="white")+
  theme(panel.grid = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        plot.title = element_text(size=14, hjust=0.5),
        legend.text = element_text(size=10),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=14, face="plain"),
        legend.ticks = element_line(color="black"),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0))

plot<-plot_grid(success, pheno, ncol=2)

save_fig(
  plot = plot,
  filename = "figures/new-sims.tiff",
  width = 0.6,
  height = 0.5
)

#### VARYING INITIAL FREQUENCY ####

## initial frequency
P_phenotype_initial_vals <- c(0.39, 0.49, 0.59)

## grid
params2 <- expand.grid(
  P_phenotype_initial = P_phenotype_initial_vals,
  N_bottleneck = N_bottleneck_vals,
  s = s_vals,
  dominance = dominance_vals,
  selection_in_recovery = selection_scenarios,
  stringsAsFactors = FALSE
)

## simulate
results2 <- params2 %>%
  mutate(
    ## convert initial phenotype frequency to allele frequency
    p0 = map2_dbl(
      P_phenotype_initial,
      dominance,
      ~ get_p_from_pheno(.x, .y)
    ),
    ## simulate bottleneck
    final_p = pmap(
      list(p0, N_bottleneck, s, dominance, selection_in_recovery),
      function(p, nb, sel, dom, sel_rec)
        simulate_bottleneck(
          p0 = p,
          Nb = nb,
          Ni = N_final,
          Nt = N_years,
          s = sel,
          dom = dom,
          selection_in_recovery = sel_rec,
          n_reps = N_replicates
        )
    ),
    ## convert final allele frequency back to phenotype frequency
    final_pheno = map2(
      final_p,
      dominance,
      function(p_vec, dom) {
        if (dom == "dominant") {
          1 - (1 - p_vec)^2
        } else {
          p_vec^2
        }
      }
    ),
    ## probability of reaching observed phenotype frequency
    p_success = map_dbl(
      final_pheno,
      ~ mean(.x >= P_phenotype_observed)
    ),
    ## mean final phenotype frequency
    mean_final_pheno = map_dbl(
      final_pheno,
      mean
    )
  ) %>%
  select(
    P_phenotype_initial,
    N_bottleneck,
    s,
    dominance,
    selection_in_recovery,
    p_success,
    mean_final_pheno
  )

## plot
success2<-results2 %>% filter(selection_in_recovery == FALSE) %>%
  mutate(label="Success Probability") %>% 
  ggplot(aes(x=N_bottleneck, y=s, z=p_success))+
  geom_contour_fill(bins=20) +
  scale_fill_viridis_c(option = "magma",
                       guide = guide_colorbar(
                         barheight = unit(0.5, "cm"),
                         barwidth = unit(5, "cm")))+
  labs(x = expression("Bottleneck Size (" * N[e] * ")"),
       y = expression("Selection Coefficient (" * italic(s) * ")"),
       fill = "Success Probability",
       title = "Success Probability") +
  ne_s_layer()+
  facet_grid(cols=vars(dominance), rows=vars(P_phenotype_initial))+
  theme(panel.grid = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        plot.title = element_text(size=14, hjust=0.5),
        legend.text = element_text(size=10),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=14, face="plain"),
        legend.ticks = element_line(color="black"),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0))

## saving
save_fig(
  plot = success2,
  filename = "figures/p0-sims.tiff",
  width = 0.6,
  height = 0.8)

results2 %>% filter(selection_in_recovery == FALSE) %>% view()
results2 %>% filter(selection_in_recovery == FALSE, P_phenotype_initial==0.39) %>% view()
