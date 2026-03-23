rm(list=ls())
library(tidyverse)
library(sandwich)
library(lmtest)

set.seed(6)

# --- Parameters ---
n      = 100
J      = 50
N      = n * J
prob_H = 0.3
NREPS  = 1000
M      = sort(rep(1:J, n))

GRID = expand.grid(s = c(0.5, 0.9), USE_SYNTHETIC = c(TRUE, FALSE))

# --- Data generating process ---
sample_df = function(h0, s) {
  H2      = h0[M]
  D_tilde = rbinom(N, size = 1, prob = s)
  D       = D_tilde * H2
  df      = data.frame(D_tilde = D_tilde, D = D, cluster = M, H = H2)
  df      = df %>% mutate(Y = 1 + 0*D + 0.0*H + rnorm(N))
  return(df)
}

# --- Test statistic ---
tstat = function(df1, use_synthetic) {
  if (use_synthetic) {
    model = lm(Y ~ D_tilde + H, df1, weights = 1 - D)
    b     = coef(model)["H"]
    se    = sqrt(diag(vcovCL(model, cluster = ~cluster, type = "HC0")))
    return(b / se["H"])
  } else {
    model = lm(Y ~ D + H, df1, weights = 1 - D)
    return(coef(model)["H"])
  }
}

# --- Run simulation for one (s, USE_SYNTHETIC) combination ---
run_sim = function(s, use_synthetic) {
  rejections = logical(NREPS)
  for (irep in seq_len(NREPS)) {
    h0    = rbinom(J, size = 1, prob = prob_H)
    df    = sample_df(h0, s)
    tobs  = tstat(df, use_synthetic)
    tvals = replicate(100, {
      h_perm = sample(h0)
      tstat(df %>% mutate(H = h_perm[cluster]), use_synthetic)
    })
    rejections[irep] = mean(abs(tvals) > abs(tobs)) <= 0.05
    if(irep %% 10 ==0) {
      cat(sprintf("s=%.1f  synthetic=%-5s  irep=%3d  %%rej=%.1f\n",
                  s, use_synthetic, irep, 100 * mean(rejections[seq_len(irep)])))
    }
  }
  mean(rejections)
}

# --- 2x2 table ---
results = GRID %>%
  rowwise() %>%
  mutate(rejection_rate = run_sim(s, USE_SYNTHETIC)) %>%
  ungroup()

table_2x2 = results %>%
  mutate(USE_SYNTHETIC = ifelse(USE_SYNTHETIC, "Synthetic", "Standard")) %>%
  pivot_wider(names_from = s, values_from = rejection_rate,
              names_prefix = "s=") %>%
  rename(Method = USE_SYNTHETIC)

print(table_2x2)