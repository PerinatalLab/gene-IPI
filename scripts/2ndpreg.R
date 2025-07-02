# second pregnancy model

library(data.table)
library(dplyr)
library(lubridate)
library(rlang)
library(tidyverse)
library(ggplot2)

setwd("~/scratch/agnes/gene-IPI")


filtered_data <- fread("results/phenotype/filtered_pregnancies_multiparous.csv")


required_cols <- c("M_ID_1724", "DeliveryDate", "ConceptionDate")

if (!all(required_cols %in% colnames(filtered_data))) {
  missing <- setdiff(required_cols, colnames(filtered_data))
  abort(paste("Hiányzó oszlop(ok):", paste(missing, collapse = ", ")))
}


first_two <- filtered_data %>%
  arrange(M_ID_1724, DeliveryDate) %>%
  group_by(M_ID_1724) %>%
  mutate(BirthOrder = row_number()) %>%
  filter(BirthOrder %in% c(1, 2)) %>%
  ungroup()

# IPI calculation
ipi_data <- first_two %>%
  group_by(M_ID_1724) %>%
  arrange(DeliveryDate) %>%
  mutate(
    PreviousDelivery = lag(DeliveryDate),
    IPI = as.numeric(difftime(ConceptionDate, PreviousDelivery, units = "days")) / 30.44
  ) %>%
  ungroup() %>%
  filter(BirthOrder == 2)

summary(ipi_data$IPI)

fwrite(ipi_data, "results/phenotype/ipi_first2births_multiparous.csv")


# IPI histogram
p_hist <- ggplot(ipi_data, aes(x = IPI)) +
  geom_histogram(binwidth = 2, color = "black", fill = "lightblue") +
  labs(title = "Distribution of interpregnancy interval (IPI)",
       x = "IPI (months)", y = "Count") +
  theme_minimal()
print(p_hist)

# density plot of IPI
p_density <- ggplot(ipi_data, aes(x = IPI)) +
  geom_density(fill = "skyblue", alpha = 0.6) +
  labs(title = "Density of interpregnancy interval (IPI)",
       x = "IPI (months)", y = "Density") +
  theme_minimal()
print(p_density)

# boxplot for IPI
p_box <- ggplot(ipi_data, aes(y = IPI)) +
  geom_boxplot(fill = "lightgreen", outlier.color = "red") +
  labs(title = "Boxplot of interpregnancy interval (IPI)",
       y = "IPI (months)") +
  theme_minimal()
print(p_box)

# QQ plot to check normality
qq_plot <- ggplot(ipi_data, aes(sample = IPI)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "QQ plot of IPI (normality check)") +
  theme_minimal()
print(qq_plot)





# full stat.
shapiro.test(ipi_data$IPI)
# not normally distributed

ipi_data <- ipi_data %>%
  mutate(
    log_IPI = log(IPI),
    sqrt_IPI = sqrt(IPI)
  )

# histogram: log(IPI)
p_log_hist <- ggplot(ipi_data, aes(x = log_IPI)) +
  geom_histogram(binwidth = 0.2, fill = "orchid", color = "black") +
  labs(title = "Histogram of log(IPI)", x = "log(IPI)", y = "Count") +
  theme_minimal()
print(p_log_hist)

# QQ plot: log(IPI)
p_log_qq <- ggplot(ipi_data, aes(sample = log_IPI)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "QQ plot of log(IPI)") +
  theme_minimal()
print(p_log_qq)

# histogram: sqrt(IPI)
p_sqrt_hist <- ggplot(ipi_data, aes(x = sqrt_IPI)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue", color = "black") +
  labs(title = "Histogram of sqrt(IPI)", x = "sqrt(IPI)", y = "Count") +
  theme_minimal()
print(p_sqrt_hist)

# QQ plot sqrt(IPI)
p_sqrt_qq <- ggplot(ipi_data, aes(sample = sqrt_IPI)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "QQ plot of sqrt(IPI)") +
  theme_minimal()
print(p_sqrt_qq)

print(shapiro.test(ipi_data$log_IPI))
print(shapiro.test(ipi_data$sqrt_IPI))

summary_stats <- summarise(ipi_data,
     sample_size = n(),
     min = min(IPI),
     q1 = quantile(IPI, 0.25),
     median = median(IPI),
     mean = mean(IPI),
     q3 = quantile(IPI, 0.75),
     max = max(IPI),
     sd = sd(IPI)
)

print(summary_stats)
