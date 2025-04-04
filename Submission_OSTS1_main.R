################################################################################
#Packages
################################################################################
library(car) #If you have never used R before, you will have to do something called install.packages(). For this package use this code "install.packages("car")"
library(ggplot2)
library(emmeans)
library(lme4)
library(multcomp)
library(multcompView)
library(performance)

################################################################################
#Loading in the data
################################################################################

#this sets the working directory
setwd() #This will need to be set in order to look at the data

#this reads the data from the csv
datum <- read.csv("OSTS1.0_8OHdG.csv")
head(datum) 

# Ensure categorical variables are properly set
datum$PopulationF <- factor(datum$Population)
datum$StatusF <- factor(datum$Status)
datum$TempF <- factor(datum$Temp)

# View the first few rows to verify the correct dataset
head(datum)

################################################################################
# Visualize - Bar plots
################################################################################
ggplot(filtered_data, aes(x = PopulationF, y = MR, color = StatusF, fill = TempF)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  theme_classic() +
  labs(x = "Population", y = "Metabolic Rate", 
       color = "Status", fill = "Temperature") +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

add_section_break <- function(num_dashes = 100) {
  cat(rep("\n", 3))
  cat(rep("-", num_dashes), "\n", sep = "")  # Dynamically generate dashes
  cat(rep("\n", 3))
}

# Call the function whenever needed
add_section_break()


################################################################################
# Stats GLM
################################################################################
# 1. Start with GLM assuming normality of residuals
# 2. Find the best model using backward stepwise selection
# 3. Use the best model to assess for normality
################################################################################

# Filter data so only 'Post' values
filtered_data <- datum[datum$StatusF == "Post", ]

# View the filtered data
print(filtered_data)

# Step 1: Fit the full GLM model
full_model <- glm(MR ~ PopulationF * TempF * Mass, data = filtered_data, family = Gamma(link = "log"))

# Print summary and ANOVA for the full model
summary(full_model)
Anova(full_model, type = 3, test = "F")

# Stepwise Model Simplification (Backward Selection)
current_model <- full_model

repeat {
  # Test dropping each term in turn
  drop_results <- drop1(current_model, test = "Chisq")
  
  # Find the term with the highest p-value (non-significant term)
  highest_p <- max(drop_results$`Pr(>Chi)`[-1], na.rm = TRUE)
  term_to_drop <- rownames(drop_results)[which.max(drop_results$`Pr(>Chi)`[-1]) + 1]
  
  # Check if the highest p-value is > 0.05
  if (highest_p > 0.05) {
    cat("\nDropping term:", term_to_drop, "(p =", highest_p, ")\n")
    new_formula <- reformulate(setdiff(all.vars(formula(current_model)), term_to_drop), response = "MR")
    current_model <- update(current_model, new_formula)
  } else {
    cat("\nAll terms significant. Stopping backward selection.\n")
    break
  }
}

# Step 2: Extract the final formula
final_formula <- formula(current_model)  # Save the best model formula

# Step 3: Fit and evaluate the final model
best_model <- glm(final_formula, data = filtered_data, family = Gamma(link = "log"))
summary(best_model)
Anova(best_model, type = 3)


################################################################################
# Visualize - Bar plots
################################################################################
ggplot(datum, aes(x = PopulationF, y = OxD, color = StatusF, fill = TempF)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  theme_classic() +
  labs(x = "Population", y = "Oxidative Damage", 
       color = "Status", fill = "Temperature") +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

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
results_glmm <- glmer(OxD ~ PopulationF * StatusF * TempF + (1 | ID), 
                      data = datum, 
                      family = Gamma(link = "log"),
                      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

# Check for warnings or issues during model fitting
if (isSingular(results_glmm, tol = 1e-5)) {
  warning("GLMM model is singular. Consider simplifying the model.")
}

# GLM without random effect
results_glm <- glm(OxD ~ PopulationF * StatusF * TempF, 
                   data = datum, 
                   family = Gamma(link = "log"))

# Compare models using nested model comparison (Chi-Square Test)
anova_results <- anova(results_glmm, results_glm, test = "Chisq")

# Print results
print(anova_results)

# Interpretation helper with color-coded output
if (anova_results$`Pr(>Chisq)`[2] < 0.05) {
  cat("\033[32mThe random effect significantly improves model fit (p < 0.05). The GLMM is preferred.\033[0m\n")  # Green text
} else {
  cat("\033[31mThe random effect does not significantly improve model fit (p >= 0.05). The GLM may suffice.\033[0m\n")  # Red text
}

add_section_break()


################################################################################
# Stats LMER
################################################################################
# 1. Start with LMER assuming normality of residuals
# 2. Find the best model using backward stepwise selection
# 3. Use the best model to assess for normality
################################################################################

# Step 1: Fit the full model
full_model <- lmer(OxD ~ PopulationF * StatusF * TempF + Mass + (1 | ID), data = datum)

# Stepwise backward selection
current_model <- full_model

repeat {
  # Perform a drop test
  drop_results <- drop1(current_model, test = "Chisq")
  
  # Identify the term with the **highest p-value**
  highest_p <- max(drop_results$`Pr(>Chi)`[-1], na.rm = TRUE)
  term_to_drop <- rownames(drop_results)[which.max(drop_results$`Pr(>Chi)`[-1]) + 1]
  
  # Stop if all remaining terms are significant
  if (highest_p <= 0.05) {
    cat("\nAll remaining terms are significant. Stopping backward selection.\n")
    break
  }
  
  # Special Rule: **Always drop the four-way and three-way interactions first if non-significant**
  if ("PopulationF:StatusF:TempF:Mass" %in% rownames(drop_results) && drop_results["PopulationF:StatusF:TempF:Mass", "Pr(>Chi)"] > 0.05) {
    term_to_drop <- "PopulationF:StatusF:TempF:Mass"
  } else if ("PopulationF:StatusF:TempF" %in% rownames(drop_results) && drop_results["PopulationF:StatusF:TempF", "Pr(>Chi)"] > 0.05) {
    term_to_drop <- "PopulationF:StatusF:TempF"
  }
  
  # Print term being removed
  cat("\nDropping term:", term_to_drop, "(p =", highest_p, ")\n")
  
  # Update the model without the removed term
  new_formula <- reformulate(setdiff(all.vars(formula(current_model)), term_to_drop), response = "OxD")
  current_model <- update(current_model, new_formula)
}

# Step 2: Extract the best model formula
final_formula <- formula(current_model)
print(paste("Best model formula:", deparse(final_formula)))

# Step 3: Fit and evaluate the final model
best_model <- lmer(final_formula, data = datum)
summary(best_model)
Anova(best_model, type = 3)




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
residuals <- residuals(best_model)

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
fitted_values <- fitted(best_model)
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
random_effects <- ranef(best_model, condVar = TRUE)

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
# Stats GLMER (Final Model)
#############################################################################
# This section runs the final GLMER model after model selection
# The fixed effects include only significant terms from backward selection.
#############################################################################

# Use the best (simplified) model from stepwise selection
summary(best_model)

# Type III ANOVA to assess significance of fixed effects
Anova(best_model, type = 3)

add_section_break()


#############################################################################
# Multiple Comparisons for Interaction Term - NEW
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

# Compute marginal means for the interaction term in your model
emmeans_interaction <- emmeans(best_model, ~ PopulationF * StatusF * TempF)

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

# Convert emmeans with letters to a data frame for ggplot
emm_data <- as.data.frame(emmeans_interaction_cld)

#############################################################################
# Visualization: Create a ggplot with shifted grouping letters
# this still does NOT work and needs debugged
#############################################################################

ggplot(emm_data, aes(x = PopulationF, y = emmean, fill = StatusF)) +
  geom_bar(stat = "identity", position = position_dodge()) +  # Bar chart
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, position = position_dodge(width = 0.9)) +  # Error bars
  geom_text(aes(label = .group), 
            position = position_dodge(width = 0.9), hjust = -0.3, vjust = -0.5) +  # Shifted letters
  ggtitle("Estimated Marginal Means with Grouping Letters") +
  xlab("PopulationF") +
  ylab("Estimated Marginal Means") +
  theme_minimal()

add_section_break()



R.version.string

RStudio.Version()$version