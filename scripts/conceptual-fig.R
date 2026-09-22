#### CONCEPTUAL FIGURE ####

## constants
source(here::here("scripts/global.R"))

## parameters
n_reps<-1000         # reps
n_gen <- 5           # generations
N <- 100             # population size
h2 <- 0.3            # heritability
s <- -0.1            # selection differential 
initial_mean <- 10   # starting phenotype value
pheno_sd <- 1        # phenotypic sd

## results
final_means <- numeric(n_reps)

for (i in 1:n_reps) {
  current_mean <- initial_mean
  
  for (g in 1:n_gen) {
    ## selection
    selection_response <- h2 * s
    
    ## drift
    drift_variance <- pheno_sd^2 / N
    drift_effect <- rnorm(1, mean = 0, sd = sqrt(drift_variance))
    
    ## new mean
    current_mean <- current_mean + selection_response + drift_effect
  }
  final_means[i] <- current_mean
}

## plotting
hist(final_means, 
     breaks = 30, 
     col = "skyblue", 
     border = "white",
     main = "Distribution of Mean Phenotypes after 50 Generations",
     xlab = "Mean Phenotypic Value",
     ylab = "Frequency of Replicates")

## data
df <- data.frame(phenotype = final_means)
mu <- mean(df$phenotype)
sigma <- sd(df$phenotype)

## plot
plot<-ggplot(df, aes(x = phenotype)) +
  geom_area(stat = "function", 
            fun = dnorm, 
            args = list(mean = mu, sd = sigma), 
            fill = "grey80", 
            alpha = 0.5) +
  #geom_hline(yintercept = 0.75, color = "black", linewidth = 0.75) +
  geom_vline(xintercept = initial_mean, color = "black", linewidth = 0.75) +
  geom_vline(xintercept = 9.55, linetype="dashed", color = "black", linewidth = 0.75) +
  geom_vline(xintercept = 10.14, linetype="dashed", color = "black", linewidth = 0.75) +
  labs(x="Mean Phenotypic Value", y="Probability Density")+
  theme_classic()+
  theme(panel.grid=element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=16, family="sans"))

## saving
save_fig(
  plot = plot,
  filename = here::here("figures/conceptual-fig.tiff"),
  width = 0.5,
  height = 0.3
)
