################################################################################
#Packages
################################################################################
library(car) #If you have never used R before, you will have to do something called install.packages(). For this package use this code "install.packages("car")"
library(ggplot2)
library(emmeans)
library(lme4)
library(multcomp)


################################################################################
#Loading in the data
################################################################################

#this sets the working directory
setwd() #This will need to be set to correctly look at the data
#this reads the data from the csv
datum <- read.csv("OSD.csv")
head(datum) #This shows you your data to double check it is the right dataset
datum$PopulationF = as.factor(datum$Population) #as.factor() makes sure R reads this as a categorical variable when graphing or analyzing
datum$SampleF = as.factor(datum$Sample)

################################################################################
# Visualize - Bar plots
################################################################################
ggplot(datum, aes(PopulationF, X8OHdG.10.dG , color = SampleF)) +
  geom_boxplot() +
  theme_classic()

add_section_break <- function(num_dashes = 100) {
  cat(rep("\n", 3))
  cat(rep("-", num_dashes), "\n", sep = "")  # Dynamically generate dashes
  cat(rep("\n", 3))
}

# Call the function whenever needed
add_section_break()


################################################################################
# Checking the Importance of Random Effects
################################################################################
# This section compares two models to evaluate the significance of random effects:
# - Model 1 (GLMM): Includes a random effect (1 | ID) for individual variability.
# - Model 2 (GLM): Excludes the random effect, using only fixed effects.
#
# The comparison involves:
# - Log-Likelihoods: Evaluates the fit of both models using their log-likelihoods.
# - Chi-Square Test: Uses a nested model comparison to test the significance 
#   of the random effect by calculating the likelihood ratio and associated 
#   p-value (Pr(>Chisq)).
#
# Interpretation:
# - If the p-value is significant (e.g., < 0.05), the random effect significantly 
#   improves the model fit, and the GLMM is preferred over the simpler GLM.
# - If the random effect does not significantly improve the fit, the simpler GLM 
#   may be sufficient.
################################################################################

# GLMM with random effect
results_glmm <- glmer(X8OHdG.10.dG ~ PopulationF * SampleF + 
                        (1 | ID), 
                      data = datum, 
                      family = Gamma(link = "log"))

# GLM without random effect
results_glm <- glm(X8OHdG.10.dG ~ PopulationF * SampleF, 
                   data = datum, 
                   family = Gamma(link = "log"))

# Compare models using likelihood ratio test
anova_results <- anova(results_glmm, results_glm, test = "Chisq")

# Print results
print(anova_results)

# Interpretation helper with color-coded output
if (anova_results$`Pr(>Chisq)`[2] < 0.05) {
  cat("\033[32mThe random effect significantly improves model fit (p < 0.05). The GLMM is preferred.\033[0m\n")  # Green text
} else {
  cat("\033[31mThe random effect does not significantly improve model fit (p >= 0.05). The GLM may suffice.\033[0m\n")  # Red text
}


################################################################################
#Stats LMER
################################################################################
# 1. Start with LMER assuming normality of residuals
# 2. Find the best model using backwards step-wise
# 3. Use the best model to assess for normality
################################################################################

# First model
results1.0 = lmer(X8OHdG.10.dG ~ PopulationF * SampleF + 
                 (1|ID), datum)
summary(results1.0)
Anova(results1.0)

add_section_break()

################################################################################
# Choosing between a GLMM and LMM: Normality Analysis of Response Variable
################################################################################
# This section assesses whether the response variable is normally distributed. 
# The results guide the choice between an LMM or a GLMM for analysis.
#
# 1. Histogram of residuals: Check if residuals are roughly symmetric and bell-shaped.
# 2. Q-Q plot of residuals: Assess if residuals align with the reference line.
# 3. Shapiro-Wilk test: Statistical test for normality of residuals. 
# 4. Residuals vs Fitted plot: Check for homoscedasticity (constant variance).
# Conclusion:
# - If the response variable does not meet the assumption of normality.
# - Based on this analysis, an LMM is not suitable, and a GLMM with a Gamma 
#   distribution (or another appropriate distribution) is preferred.
################################################################################

# Extract residuals from the model
residuals <- residuals(results1.0)

# 1. Histogram of residuals
hist(residuals, 
     main = "Histogram of Residuals", 
     xlab = "Residuals", 
     col = "lightblue", 
     breaks = 20)

# 2. Q-Q plot
qqnorm(residuals)
qqline(residuals, col = "red") 

# 3. Shapiro-Wilk test for normality
shapiro_results <- shapiro.test(residuals)

# Print Shapiro-Wilk test results
cat("\033[34mShapiro-Wilk Test Results:\033[0m\n")
print(shapiro_results)

# 4. Residuals vs Fitted plot
fitted_values <- fitted(results1.0)
plot(fitted_values, residuals, 
     main = "Residuals vs Fitted", 
     xlab = "Fitted Values", 
     ylab = "Residuals", 
     pch = 20)
abline(h = 0, col = "red", lwd = 2) 

# Interpretation Helper
cat("\033[33mInterpretation Helper:\033[0m\n")
cat("- Shapiro-Wilk p-value:\n")
if (shapiro_results$p.value < 0.05) {
  cat("\033[31m  Residuals deviate from normality (p < 0.05).\033[0m\n")
} else {
  cat("\033[32m  Residuals are consistent with normality (p >= 0.05).\033[0m\n")
}

add_section_break()


################################################################################
# Check Normality of Random Effects
################################################################################

# Extract random effects
random_effects <- ranef(results1.0, condVar = TRUE)

# Inspect the structure of the random effects
cat("\033[34mStructure of Random Effects:\033[0m\n")
str(random_effects)

# Extract random effects for the grouping factor 'ID'
re_values <- as.numeric(random_effects$ID[, "(Intercept)"])

# 1. Histogram of random effects
hist(re_values, 
     main = "Histogram of Random Effects (ID)", 
     xlab = "Random Effect Values", 
     col = "lightgreen", 
     breaks = 20)

# 2. Q-Q plot of random effects
qqnorm(re_values, main = "Q-Q Plot of Random Effects (ID)")
qqline(re_values, col = "blue")

# 3. Shapiro-Wilk test for random effects normality
shapiro_re_results <- shapiro.test(re_values)

# Print Shapiro-Wilk test results for random effects
cat("\033[34mShapiro-Wilk Test Results for Random Effects:\033[0m\n")
print(shapiro_re_results)

# Interpretation Helper for Random Effects
cat("\033[33mInterpretation Helper:\033[0m\n")
cat("- Shapiro-Wilk p-value for Random Effects:\n")

if (shapiro_re_results$p.value < 0.05) {
  cat("\033[31m  Random effects deviate from normality (p < 0.05).\033[0m\n")
} else {
  cat("\033[32m  Random effects are consistent with normality (p >= 0.05).\033[0m\n")
}


add_section_break()

#############################################################################
#Stats GLMER
#############################################################################
#- This analysis runs a GLMER and also  specifies a gamma 
#distribution if the residuals are not normally distributed
#############################################################################

# results1.1 (make sure to adjust this to whatever the correct results is)
results1.1 <- glmer(X8OHdG.10.dG ~ PopulationF:SampleF + PopulationF + SampleF 
                  + (1 | ID), 
                  data = datum, 
                  family = Gamma(link = "log"))
summary(results1.1)
Anova(results1.1)

add_section_break()

#############################################################################
# Multiple Comparisons for Interaction Term
#############################################################################
# This code performs post hoc analysis to evaluate pairwise differences 
# between the significant interaction term in the model. It uses estimated 
# marginal means (emmeans) and Tukey adjustment for multiple comparisons.

# How to interpret:
# 1. Blue bars: Represent 95% confidence intervals (CIs) for group means. 
#    Non-overlapping bars suggest potential differences but are not definitive.
# 2. Red arrows: Represent pairwise comparisons. A difference is statistically 
#    significant if the arrow does not cross zero.

#############################################################################

# Load necessary libraries
library(emmeans)
library(multcomp)
library(multcompView)

# Compute marginal means for the interaction term in your model
emmeans_interaction <- emmeans(results1.1, ~ PopulationF * SampleF)

# Perform pairwise comparisons with Tukey adjustment
pairwise_comparisons <- pairs(emmeans_interaction, adjust = "tukey")

# Display pairwise comparisons
print(pairwise_comparisons)

# Add compact letter display (cld) to marginal means
emmeans_interaction_cld <- cld(object = emmeans_interaction,
                               adjust = "Tukey",
                               Letters = letters,
                               alpha = 0.05)

# Display the emmeans with letters
print(emmeans_interaction_cld)

# Visualization with grouping letters
library(ggplot2)

#Kruskal-Wallis
emm_data$Group <- interaction(emm_data$PopulationF, emm_data$SampleF, sep = "_")

kruskal.test(emmean ~ Group, data = emm_data)

library(rstatix)

dunn_test_result <- dunn_test(emmean ~ Group, data = emm_data, p.adjust.method = "BH")
print(dunn_test_result)

ggplot(emm_data, aes(x = Group, y = emmean, fill = PopulationF)) +
  geom_bar(stat = "identity", position = position_dodge()) +  
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(width = 0.9)) +
  geom_text(aes(label = .group), position = position_dodge(width = 0.9), hjust = -0.3, vjust = -0.5) +
  ggtitle("Estimated Marginal Means with Interaction (PopulationF * SampleF)") +
  xlab("Group (PopulationF * SampleF)") +
  ylab("Estimated Marginal Means") +
  theme_minimal()


add_section_break()

