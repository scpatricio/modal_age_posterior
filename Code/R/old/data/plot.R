setwd("~/Desktop/PhD/codes/mode_dist/data/")
mode_summary <- read.csv("~/Desktop/PhD/codes/mode_dist/data/mode_summary_new.txt", sep="")

library(ggplot2)

mode_summary %>% 
  filter(Sex != "Total") %>%
  ggplot(aes(x = Year, color = Sex)) +
  # Confidence interval shaded area with custom color and no border
  geom_ribbon(aes(ymin = quantile_2_5, ymax = quantile_97_5, fill = Sex), alpha = 0.4, color = NA) +
  # Line for the mean with differentiated linetype
  geom_line(aes(y = mean), size = .5) +
  # Line for the median with differentiated linetype
  # geom_point(aes(y = mode, linetype = "Mode"), size = .5) +
  # Facet by country
  facet_wrap(~ Country) +
  # Customize the legend for fill and linetype
  scale_fill_manual(values = c("Male" = "#2f82b4", "Female" = "#e34a1c"), name = "Sex") +
  scale_color_manual(values = c("Male" = "#2f82b4", "Female" = "#e34a1c"), name = "Sex") +
  scale_linetype_manual(values = c("Mean" = "solid", "Mode" = "dotted"), name = "Location") +
  # Labels and theme adjustments
  labs(x = "Year", y = "Mode", title = "Modal age at death") +
  theme_minimal() +
  theme(legend.position = "bottom")

mode_summary <- read.csv("~/Desktop/PhD/codes/mode_dist/data/mode_full.txt", sep="")

library(dplyr)

# Ensure prob is continuous and normalized
mode_summary %>%
  filter(Sex != "Total", prob > 0) %>%
  group_by(Sex, Contry, Year) %>% 
  mutate(prob = (prob - min(prob)) / (max(prob) - min(prob))) %>% 
  ungroup() %>%
  ggplot(aes(x = Year, y = mode, fill = Sex, alpha = prob)) +
  # Heatmap tiles
  geom_tile() +
  # Continuous alpha scale with explicit breaks
  scale_alpha_continuous(
    range = c(0, 1), 
    name = "Probability", 
    breaks = seq(0, 1, 0.2), 
    labels = scales::percent_format(accuracy = 1)
  ) +
  # Custom color for Sex
  scale_fill_manual(values = c("Male" = "#2f82b4", "Female" = "#e34a1c"), name = "Sex") +
  # Facet by country
  facet_wrap(~ Contry) +
  # Labels and theme adjustments
  labs(x = "Year", y = "Age", title = "Probability distribution of the modal age at death") +
  theme_minimal() +
  theme(legend.position = "bottom")



mode_summary %>%
  filter(Sex != "Total", prob > 0) %>%
  group_by(Sex, Contry, Year) %>% 
  mutate(prob = (prob - min(prob)) / (max(prob) - min(prob))) %>%
  ungroup() %>%
  ggplot(aes(x = Year, y = mode, fill = Sex, alpha = prob)) +
  geom_tile() +
  scale_alpha_continuous(
    # range = c(0, 1), 
    name = "Probability", 
    breaks = seq(0, 1, 0.1), 
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_fill_manual(
    values = c("Male" = "#2f82b4", "Female" = "#e34a1c"), 
    name = "Sex"
  ) +
  facet_wrap(~ Contry) +
  labs(
    x = "Year", 
    y = "Age", 
    title = "Probability distribution of the modal age at death"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  guides(alpha = guide_colourbar())  # Add this line

library(ggplot2)
library(dplyr)
library(scales)

custom_trans <- trans_new(
  name = "custom",
  transform = function(x) x^(1/2),  # Apply square root scaling
  inverse = function(x) x^2        # Inverse of square root scaling
)

# Aggregate mean values for the line
mean_values <- mode_summary %>%
  filter(Sex != "Total", prob > 0) %>%
  group_by(Sex, Contry, Year) %>%
  summarise(mean = sum(prob * mode), .groups = "drop")

# Prepare the data for the heatmap
df_plot <- mode_summary %>%
  filter(Sex != "Total", prob > 0.001) %>%
  group_by(Sex, Contry, Year) %>%
  mutate(
    prob = (prob - min(prob)) / (max(prob) - min(prob))
  ) %>%
  ungroup()

ggplot() +
  # Heatmap layer
  geom_tile(
    data = df_plot,
    aes(x = Year, y = mode, fill = Sex, alpha = prob)
  ) +
  # Line layer
  # geom_line(
  #   data = mean_values,
  #   aes(x = Year, y = mean, color = Sex),
  #   size = 1.2
  # ) +
  # Alpha scale for the heatmap
  scale_alpha_continuous(
    # range = c(0, 1),
    name = "Probability",
    labels = scales::percent_format(accuracy = 1),
    trans = custom_trans,
    guide = guide_colourbar()
  ) +
  # Fill scale for the heatmap
  scale_fill_manual(
    values = c("Male" = "#2f82b4", "Female" = "#e34a1c"),
    name = "Sex",
    guide = guide_legend(
      override.aes = list(alpha = 1)  # Show full opacity in the legend
    )
  ) +
  # Adjust y-axis to show intervals of 10
  scale_y_continuous(
    breaks = seq(0, 100, by = 5),
    expand = c(0, 0)  # Optional: remove extra padding on the axis
  ) +
  # Color scale for the line
  scale_color_manual(
    values = c("Male" = "#2f82b4", "Female" = "#e34a1c"), 
    name = "Sex"
  ) +
  # Facet by country
  facet_wrap(~ Contry, ncol = 2) +
  # Labels and titles
  labs(
    x = "Year",
    y = "Age",
    title = "Probability distribution of the modal age at death"
  ) +
  # Theme adjustments
  theme_minimal() +
  theme(legend.position = "bottom") +
  # Custom legends
  guides(
    alpha = guide_colourbar(order = 1),
    fill = guide_legend(order = 2),
    color = guide_legend(order = 2)  # Add a legend for the line
  )
 
