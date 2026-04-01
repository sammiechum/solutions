# Title: "02_create_visualizations.R"
# Author: "Sammie Chum"
# Date: "2026-03-31"

# Load libraries & data --------------------------------------------------------
# install.packages("ggplot2")

library(pharmaverseadam)
library(ggplot2)
library(dplyr)

adae <- pharmaverseadam::adae 

# Plot 1: AE severity distribution by treatment --------------------------------
ggplot(adae, aes(x = ACTARM, fill = AESEV)) +
  geom_bar() +
  scale_fill_manual(
    values = c(
      "MILD"     = "#E24B4A",
      "MODERATE" = "#1D9E75",
      "SEVERE"   = "#378ADD"
    ),
    labels = c("Mild", "Moderate", "Severe"),
    name = "Severity/Intensity"
  ) +
  labs(
    title = "AE severity distribution by treatment",
    x = "Treatment Arm",
    y = "Count AEs"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right"
  )

# Save to png file
ggsave("data/ae_severity.png")


# Plot 2: Top 10 most frequent AEs  --------------------------------------------
# Preprocess data: filter to top 10, calculate percentage and CIs
top10_ae <- adae %>%
  filter(AETERM %in% (adae %>%
                        count(AETERM) %>%
                        slice_max(n, n = 10) %>%
                        pull(AETERM))) %>%
  group_by(AETERM) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(
    n_total = sum(n),
    percentage = n / n_total * 100
  ) %>%
  rowwise() %>%
  mutate(
    ci_lower = binom.test(n, n_total, conf.level = 0.95)$conf.int[1] * 100,
    ci_upper = binom.test(n, n_total, conf.level = 0.95)$conf.int[2] * 100
  ) %>%
  ungroup()

# Create plot
ggplot(top10_ae, aes(x = percentage, y = reorder(AETERM, percentage))) +
  geom_point(size = 3, colour = "black") +
  geom_errorbarh(
    aes(xmin = ci_lower, xmax = ci_upper),
    width = 0.3,
    colour = "black",
    linewidth = 0.8
  ) +
  labs(
    title = "Top 10 Most Frequent Adverse Events",
    x = "Percentage of Patients (%)",
    y = "Adverse Event",
    subtitle = "Error bars represent 95% Clopper-Pearson confidence intervals"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 10),
    plot.subtitle = element_text(hjust = 0, face = "italic")
  )

# Save plot to png file
ggsave(
  filename = "data/top10_ae.png",
  plot = last_plot(), width = 9, height = 6, dpi = 300, units = "in"
)
