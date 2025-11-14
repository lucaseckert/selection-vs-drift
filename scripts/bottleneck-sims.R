#### Bottleneck Simulations ####

## packages
library(tidyverse)
library(purrr)

#### FIXED PARAMETERS ####
N_initial <- 5000             # final census size
N_years <- 7                  # total generations (bottleneck + recovery)
P_phenotype_initial <- 0.49   # initial phenotype frequency
P_phenotype_observed <- 0.87  # observed final phenotype frequency
N_replicates <- 10000         # number of independent simulations

#### VARIABLE PARAMETERS ####
N_bottleneck_vals <- c(5, 10, 15, 20, 30, 50) # size of population at bottleneck
N_generations_vals <- c(1, 2, 3)              # duration of bottleneck in generations
s_vals <- c(0, 0.05, 0.1, 0.2, 0.5)           # selection coefficient of phenotype
dominance_vals <- c("dominant", "recessive")  # dominance of observed phenotype

#### FUNCTION: phenotype to allele frequency ####
get_p_from_pheno <- function(P_pheno, dominance_type) {
  if (dominance_type == "dominant") {
    p <- 1 - sqrt(1 - P_pheno)
  } else {
    p <- 1 - sqrt(P_pheno)
  }
  return(p)
}

#### FUNCTION: simulate drift + selection + exponential recovery ####
simulate_bottleneck <- function(p_allele_start, 
                                N_bottleneck, 
                                N_initial, 
                                N_bottleneck_gens, 
                                N_total_gens,
                                s,
                                dominance) {
  
  p <- p_allele_start
  
  ### BOTTLENECK GENERATIONS (SELECTION + DRIFT)
  for (g in 1:N_bottleneck_gens) {
    
    ### 1. SELECTION BEFORE DRIFT
    if (s > 0) {
      if (dominance == "dominant") {
        # fitness: AA=1+s, Aa=1+s, aa=1
        p <- (p * (1 + s)) / (1 + s * p)
      } else {
        # recessive phenotype favored: aa=1+s, Aa=1, AA=1
        # allele a has freq q = 1-p
        q <- 1 - p
        w_bar <- 1 + s * q^2     # mean fitness
        q_prime <- q * (1 + s * q) / w_bar
        p <- 1 - q_prime
      }
    }
    
    ### 2. DRIFT
    p <- rbinom(1, size = 2 * N_bottleneck, prob = p) / (2 * N_bottleneck)
  }
  
  ### RECOVERY GENERATIONS (DRIFT ONLY)
  N_recovery_gens <- N_total_gens - N_bottleneck_gens
  if (N_recovery_gens > 0) {
    growth_factor <- (N_initial / N_bottleneck)^(1 / N_recovery_gens)
    N_values <- round(N_bottleneck * (growth_factor ^ seq_len(N_recovery_gens)))
    
    for (N_gen in N_values) {
      p <- rbinom(1, size = 2 * N_gen, prob = p) / (2 * N_gen)
    }
  }
  
  return(p)
}

#### PARAMETER GRID ####
params <- expand.grid(
  N_bottleneck = N_bottleneck_vals,
  N_generations = N_generations_vals,
  dominance = dominance_vals,
  s = s_vals,
  stringsAsFactors = FALSE
)

#### RUN SIMULATIONS ####
set.seed(42)
results <- params %>%
  mutate(
    p0 = map_dbl(dominance, ~ get_p_from_pheno(P_phenotype_initial, .x)),
    p_target = map_dbl(dominance, ~ get_p_from_pheno(P_phenotype_observed, .x)),
    
    final_p = pmap(
      list(p0, N_bottleneck, N_generations, dominance, s),
      function(p0_i, Nb, Ng, dom, sel) {
        replicate(
          N_replicates,
          simulate_bottleneck(
            p_allele_start = p0_i,
            N_bottleneck = Nb,
            N_initial = N_initial,
            N_bottleneck_gens = Ng,
            N_total_gens = N_years,
            s = sel,
            dominance = dom
          )
        )
      }
    ),
    
    ### convert to phenotype
    final_pheno = map2(final_p, dominance, function(p_vec, dom) {
      q_vec <- 1 - p_vec
      if (dom == "dominant") {
        1 - (q_vec)^2 
      } else {
        (q_vec)^2 
      }
    }),
    
    n_success = map_int(final_pheno, ~ sum(.x >= P_phenotype_observed)),
    p_success = n_success / N_replicates,
    mean_final = map_dbl(final_pheno, mean),
    sd_final = map_dbl(final_pheno, sd)
  )

#### SUMMARY TABLE ####
summary_table <- results %>%
  mutate(
    q0 = 1 - p0,
    q_target = 1 - p_target,
    allele_increased = if_else((p_target - p0) >= (q_target - q0), "p", "q"),
    delta_frequency = if_else(allele_increased == "p",
                              p_target - p0,
                              q_target - q0)
  ) %>%
  select(N_bottleneck, N_generations, s, dominance,
         p0, q0, p_target, q_target, allele_increased, delta_frequency,
         p_success, mean_final, sd_final)

## view results
view(summary_table)

## plotting
ggplot(summary_table, aes(x=N_bottleneck, y=p_success, color=factor(N_generations)))+
  geom_line(linewidth = 1)+
  scale_color_viridis_d(option = "magma", begin = 0.15, end = 0.75, direction=-1) +
  labs(x="Bottleneck Size", y="Probability of Success", color="Bottleneck Duration")+
  facet_grid(rows=vars(dominance), cols=vars(s), labeller = labeller(s = function(x) paste0("s = ", x)))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        strip.text = element_text(size=14))
