setwd("~/Sam")
library(readxl)
small2 <- read_excel("small2.xlsx")
data <- read_excel("PBRC GR Final 2020-04-21.xlsx")
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(lmerTest) 
library(emmeans)
library(gridExtra)
library(grid)

# Reshape Data set
data_long <- data %>% filter(!is.na(GRID) & GRID != "") %>%
  dplyr::select(GRID, Sex, Age, DP1, DP2, TCA, TCS, TCR, HDLA, HDLS, HDLR) %>%
  pivot_longer(cols = c(TCA, TCS, HDLA, HDLS), names_to = "lipid_treatment_raw", values_to = "Delta_Response") %>%
  mutate(lipid_type = if_else(str_detect(lipid_treatment_raw, "^TC"), "TCD", "HDLD"),
         treatment = as.factor(if_else(lipid_treatment_raw %in% c("TCA", "HDLA"), "AAD", "DASH")),
         # Classify chronological period by verifying source columns
         period = case_when(lipid_treatment_raw %in% c("TCA", "HDLA") & (DP1 == "AAD" | DP1 == "A") ~ "Period 1",
                            lipid_treatment_raw %in% c("TCS", "HDLS") & (DP1 == "S2D" | DP1 == "S") ~ "Period 1",
                            TRUE ~ "Period 2"),
    period = as.factor(period),
    trt_seq = as.factor(if_else(DP1 == "AAD" | DP1 == "A", "AAD-DASH", "DASH-AAD")),
    Baseline_Raw = if_else(lipid_type == "TCD", TCR, HDLR)) %>% arrange(GRID, lipid_type, period)

# EDA
## Demographic Information
print(table(data$Sex))
print(summary(data$Age))

p1 <- ggplot(data) + 
  geom_histogram(aes(x = Age), color = "black", fill = "pink", binwidth = 10) + 
  ylab("Number of Participants") + ggtitle("Age Range of Participants") + theme_minimal()

p2 <- ggplot(filter(data_long, lipid_type == "HDLD"), aes(x = Baseline_Raw)) + 
  geom_histogram(bins = 20, fill = "burlywood", color = "white") + 
  labs(title = "Distribution of Baseline HDL", x = "Baseline HDL", y = "Count") + theme_minimal()

p3 <- ggplot(filter(data_long, lipid_type == "TCD"), aes(x = Baseline_Raw)) + 
  geom_histogram(bins = 20, fill = "burlywood", color = "white") + 
  labs(title = "Distribution of Baseline TC", x = "Baseline TC", y = "Count") + theme_minimal()

print(p1); print(p2); print(p3)

## Raw Baseline Summary Statistics
print(summary(filter(data_long, lipid_type == "HDLD")$Baseline_Raw))
print(summary(filter(data_long, lipid_type == "TCD")$Baseline_Raw))

## Crossover trendlines using your original metrics
p4 <- data_long %>% filter(lipid_type == "HDLD") %>% group_by(treatment, period) %>% 
  summarise(mean_HDLD = mean(Delta_Response, na.rm = TRUE), 
            se_HDLD = sd(Delta_Response, na.rm = TRUE) / sqrt(n()), .groups = 'drop') %>% 
  ggplot(aes(x = period, y = mean_HDLD, color = treatment, group = treatment)) + 
  geom_point(size = 3) + geom_line(linewidth = 1) + 
  geom_errorbar(aes(ymin = mean_HDLD - se_HDLD, ymax = mean_HDLD + se_HDLD), width = 0.1) + 
  labs(title = "Mean Change in HDL (HDLD) by Treatment and Period", x = "Period", y = "Mean HDLD") + theme_minimal()

p5 <- data_long %>% filter(lipid_type == "TCD") %>% group_by(treatment, period) %>% 
  summarise(mean_TCD = mean(Delta_Response, na.rm = TRUE), se_TCD = sd(Delta_Response, na.rm = TRUE) / sqrt(n()), 
            .groups = 'drop') %>% 
  ggplot(aes(x = period, y = mean_TCD, color = treatment, group = treatment)) + 
  geom_point(size = 3) + geom_line(linewidth = 1) + 
  geom_errorbar(aes(ymin = mean_TCD - se_TCD, ymax = mean_TCD + se_TCD), width = 0.1) + 
  labs(title = "Mean Change in TC (TCD) by Treatment and Period", x = "Period", y = "Mean TCD") + theme_minimal()

p6 <- ggplot(filter(data_long, lipid_type == "HDLD"), aes(x = trt_seq, y = Delta_Response)) + 
  geom_boxplot(fill = "darkolivegreen3") + labs(title = "HDLD by Treatment Sequence", x = "Sequence", y = "HDLD") + 
  theme_minimal()
p7 <- ggplot(filter(data_long, lipid_type == "TCD"), aes(x = trt_seq, y = Delta_Response)) + 
  geom_boxplot(fill = "darkolivegreen3") + labs(title = "TCD by Treatment Sequence", x = "Sequence", y = "TCD") + 
  theme_minimal()

print(p4); print(p5); print(p6); print(p7)

# Create subsets and place to store Results
subsets <- list(TCD = "TC", HDLD = "HDL")
table2 <- data.frame()
table3 <- data.frame()
table4 <- data.frame()

# Results
for (resp in names(subsets)) {
  label_var <- subsets[[resp]]
  analysis_subset <- data_long %>% filter(lipid_type == resp)
  # Centering the baseline variable
  analysis_subset$Baseline_Cov <- as.numeric(scale(analysis_subset$Baseline_Raw, scale = FALSE))
  # Program 5
  f5 <- as.formula("Delta_Response ~ Baseline_Cov + trt_seq + period + treatment + (1 | GRID)")
  model5_sas <- lmer(f5, data = analysis_subset, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
  # Outputs
  m5_ano     <- anova(model5_sas, type = "3", ddf = "Satterthwaite")
  ls_means_5 <- emmeans(model5_sas, ~ treatment, level = 0.90, lmer.df = "Satterthwaite")
  m5_diff    <- pairs(ls_means_5, reverse = TRUE)
  # Create Type III Fixed Effects Table
  df_anova <- as.data.frame(m5_ano) %>%
    rownames_to_column(var = "Effect") %>%
    mutate(Response = label_var, `Numerator DF` = as.numeric(NumDF), `Denominator DF` = round(as.numeric(DenDF), 1),
           `F` = round(`F value`, 4), `P-value` = if_else(`Pr(>F)` < 0.0001, "<0.0001", as.character(round(`Pr(>F)`, 4))),
           Effect = case_when(Effect == "trt_seq" ~ "sequence", Effect == "period" ~ "period", 
                              Effect == "treatment"~ "treatment", Effect == "Baseline_Cov" ~ paste0(label_var, "R"),
                              TRUE ~ as.character(Effect))) %>%
    dplyr::select(Effect, Response, `Numerator DF`, `Denominator DF`, `F`, `P-value`)
  table4 <- rbind(table4, df_anova)
  # Create Adjusted Least Squares Means Table
  df_emm <- as.data.frame(ls_means_5) %>% 
    mutate(Response  = label_var, `LS Mean` = round(emmean, 4), `SE` = round(SE, 4),`DF` = round(df, 1), 
           `90% CI`  = paste0("[", round(lower.CL, 4), ", ", round(upper.CL, 4), "]")) %>%
    dplyr::select(Response, Treatment = treatment, `LS Mean`, SE, DF, `90% CI`)
  table2 <- rbind(table2, df_emm)
  # Create Differences in Least Squares Means Table
  df_pairs_ci <- as.data.frame(confint(m5_diff))
  df_pairs_p  <- as.data.frame(m5_diff)
  df_pairs <- df_pairs_ci %>% inner_join(df_pairs_p, by = c("contrast", "SE", "df", "estimate")) %>%
    mutate(Response = label_var, `Estimate` = round(estimate, 4), `SE` = round(SE, 4), `DF` = round(df, 1),
           `90% CI` = paste0("[", round(lower.CL, 4), ", ", round(upper.CL, 4), "]"),
           `T-value` = round(t.ratio, 4), `P-value` = if_else(p.value < 0.0001, "<0.0001", as.character(round(p.value, 4))),
           Comparison = contrast) %>%
    dplyr::select(Response, Comparison, Estimate, SE, DF, `90% CI`, `T-value`, `P-value`)
  table3 <- rbind(table3, df_pairs)}

# Export PDFs
export_pdf_table <- function(data_frame, file_name) {
  pdf(file_name, width = 8.5, height = 3.5)
  table_theme <- ttheme_default(core = list(bg_params = list(fill = c("white", "#F2F2F2")), fg_params = list(fontsize = 10)),
                                colhead = list(bg_params = list(fill = "darkolivegreen3"), 
                                               fg_params = list(col = "white", fontface = "bold", fontsize = 11)))
  table_grob <- tableGrob(data_frame, rows = NULL, theme = table_theme)
  grid.arrange(table_grob, ncol = 1)
  dev.off()}

export_pdf_table(table2, "table2.pdf")
export_pdf_table(table3, "table3.pdf")
export_pdf_table(table4, "table4.pdf")
