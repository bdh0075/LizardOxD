################################################################################
#Packages
################################################################################
library(car) #If you have never used R before, you will have to do something called install.packages(). For this package use this code "install.packages("car")"
library(ggplot2)
library(emmeans)
library(lme4)
library(multcomp)
library(multcompView)

################################################################################
# Loading in the data for warm-up temp
################################################################################

datum <- read.csv("OSTSExp2warmup.csv")
#This shows you your data to double check it is the right dataset
head(datum)
# Ensure categorical variables are properly set
datum$PopulationF = as.factor(datum$Population) #as.factor() makes sure R reads this as a categorical variable when graphing or analyzing
datum$ChallTrtF = as.factor(datum$ChallTrt)
# View the first few rows to verify the correct dataset
head(datum) 

add_section_break()

################################################################################
# Visualize - Scatter Plot
################################################################################

ggplot(datum, aes(x = ChallTemp, y = ChallDamage, color = PopulationF)) +
  geom_point(size = 3, alpha = 0.7) + # Scatter plot points
  theme_classic() +                  # Clean theme
  labs(
    x = "Challenge Temperature",     # Custom x-axis label
    y = "Challenge Damage",          # Custom y-axis label
    color = "Population"             # Legend title
  ) +
  theme(
    legend.position = "right",       # Place legend on the right
    text = element_text(size = 14)   # Increase text size
  )

add_section_break()

ggplot(datum, aes(x = Rate, y = ChallDamage, color = PopulationF)) +
  geom_point(size = 3, alpha = 0.7) + # Scatter plot points
  theme_classic() +                  # Clean theme
  labs(
    x = "Rate",     # Custom x-axis label
    y = "Challenge Damage",          # Custom y-axis label
    color = "Population"             # Legend title
  ) +
  theme(
    legend.position = "right",       # Place legend on the right
    text = element_text(size = 14)   # Increase text size
  )

add_section_break()


################################################################################
#Stats LM
################################################################################
# 1. Start with LM assuming normality of residuals
# 2. Find the best model using backwards step-wise
# 3. Use the best model to assess for normality
################################################################################

# Ensure sum contrasts for Type III ANOVA (important for categorical variables)
options(contrasts = c("contr.sum", "contr.poly"))

# Full model with all interactions
results_full <- lm(ChallDamage ~ PopulationF * Rate, data = datum)
summary(results_full)
Anova(results_full, type = "III")  # Type III tests for categorical variables

# Stepwise model selection (Backward elimination)
results_stepwise <- step(results_full, direction = "backward", test = "F")
summary(results_stepwise)
Anova(results_stepwise)

# Extract interaction F-statistic and DF
anova_results <- Anova(results_full, type = "III")
interaction_df <- anova_results["PopulationF:Rate", "Df"]
interaction_f_value <- anova_results["PopulationF:Rate", "F value"]

interaction_f_value  # F-value
interaction_df  # Numerator
summary(results_full)$df[2]  # Denominator


################################################################################
# Choosing between a GLMM and LM: Normality Analysis of Response Variable
################################################################################
# This section assesses whether the response variable is normally distributed. 
# The results guide the choice between an LM or a GLMM for analysis.
#
# 1. Histogram of residuals: Check if residuals are roughly symmetric and bell-shaped.
#    Deviations (skewness or outliers) suggest non-normality.
# 
# 2. Q-Q plot of residuals: Assess if residuals align with the reference line. 
#    Deviations, especially in tails, indicate non-normality.
# 
# 3. Shapiro-Wilk test: Statistical test for normality of residuals. 
#    A significant p-value (< 0.05) indicates residuals deviate from normality.
# 
# 4. Residuals vs Fitted plot: Check for homoscedasticity (constant variance). 
#    Random scatter is good; patterns (e.g., fanning) suggest issues like heteroscedasticity.
# 
# Conclusion:
# - If the response variable does not meet the assumption of normality.
# - Based on this analysis, an LM is not suitable, and a 
#   GLMM with a Gamma distribution (or another appropriate 
#   distribution) is preferred.
################################################################################

# Extract residuals from the model
residuals <- residuals(results_stepwise)

# 1. Histogram of residuals to assess the overall distribution
hist(residuals, 
     main = "Histogram of Residuals", 
     xlab = "Residuals", 
     col = "lightblue", 
     breaks = 20)

# 2. Q-Q plot to visually check if residuals follow a normal distribution
qqnorm(residuals)
qqline(residuals, col = "red")  # Add a reference line for normality

# 3. Shapiro-Wilk test for normality of residuals
shapiro_results <- shapiro.test(residuals)  # Store the test results

# Print Shapiro-Wilk test results
cat("\033[34mShapiro-Wilk Test Results:\033[0m\n")
print(shapiro_results)

# 4. Residuals vs Fitted plot to check for homoscedasticity
fitted_values <- fitted(results_stepwise)
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
# Create a scatter plot of ChallDamage versus Rate, with points colored 
################################################################################
# by PopulationF. This visualizes the raw data, allowing for observation 
# of differences between populations.
# Add individual data points to show the distribution of ChallDamage 
# at different Rates
# Overlay a smoothed curve using a Gamma GLM with a log link to model the 
# relationship between Rate and ChallDamage
# This highlights the trend while accounting for the distribution of the 
# response variable.
################################################################################  

ggplot(datum, aes(x = Rate, y = ChallDamage, color = PopulationF)) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = Gamma(link = "log"))) +
  theme_classic()


add_section_break()


