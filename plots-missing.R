
library(data.table)
library(ggplot2)
library(patchwork)

coverage_pfi_mean = readRDS(sprintf("%s/coverage_pfi_mean.Rds", res_dir))
coverage_pdp_mean = readRDS(sprintf("%s/coverage_pdp_mean.Rds", res_dir))

# PFI
pfi100 <- ggplot(coverage_pfi_mean[n == 100, ],
       aes_string(x = "nrefits", y = "coverage", color = "impute")) + 
  facet_grid(problem ~ algorithm) + 
  geom_line() + 
  scale_y_continuous(sprintf("Confidence Interval %s", "Coverage"), limits = c(0, 1)) +
  scale_x_continuous("Number of Model Refits") +
  scale_color_discrete("Sampling Strategy") + 
  geom_hline(yintercept = .95, linetype = "dashed") + 
  ggtitle("PFI (n = 100)")

pfi1000 <- ggplot(coverage_pfi_mean[n == 1000, ],
       aes_string(x = "nrefits", y = "coverage", color = "impute")) + 
  facet_grid(problem ~ algorithm) + 
  geom_line() + 
  scale_y_continuous(sprintf("Confidence Interval %s", "Coverage"), limits = c(0, 1)) +
  scale_x_continuous("Number of Model Refits") +
  scale_color_discrete("Sampling Strategy") + 
  geom_hline(yintercept = .95, linetype = "dashed") + 
  ggtitle("PFI (n = 1000)")

# PDP
pdp100 <- ggplot(coverage_pdp_mean[n == 100, ],
                 aes_string(x = "nrefits", y = "coverage", color = "impute")) + 
  facet_grid(problem ~ algorithm) + 
  geom_line() + 
  scale_y_continuous(sprintf("Confidence Interval %s", "Coverage"), limits = c(0, 1)) +
  scale_x_continuous("Number of Model Refits") +
  scale_color_discrete("Sampling Strategy") + 
  geom_hline(yintercept = .95, linetype = "dashed") + 
  ggtitle("PDP (n = 100)")

pdp1000 <- ggplot(coverage_pdp_mean[n == 1000, ],
                  aes_string(x = "nrefits", y = "coverage", color = "impute")) + 
  facet_grid(problem ~ algorithm) + 
  geom_line() + 
  scale_y_continuous(sprintf("Confidence Interval %s", "Coverage"), limits = c(0, 1)) +
  scale_x_continuous("Number of Model Refits") +
  scale_color_discrete("Sampling Strategy") + 
  geom_hline(yintercept = .95, linetype = "dashed") + 
  ggtitle("PDP (n = 1000)")

# Plot together
pfi100 / pfi1000
ggsave("pfi_coverage_missings.pdf", width = 10, height = 20)
