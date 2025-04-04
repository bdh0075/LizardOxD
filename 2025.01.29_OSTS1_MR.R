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
setwd("~/Google Drive/My Drive/Research/Manuscripts/Projects - In progress/OSTS lizards/2025.01 - JEB Submission/Data Analysis/Mark's better analysis")

#this reads the data from the csv
datum <- read.csv("OSTS1.0_8OHdG.csv")
head(datum) 

# Ensure categorical variables are properly set
datum$PopulationF <- factor(datum$Population)
datum$StatusF <- factor(datum$Status)
datum$TempF <- factor(datum$Temp)

# View the first few rows to verify the correct dataset
head(datum)

# Filter data so only 'Post' values
filtered_data <- datum[datum$StatusF == "Post", ]

# View the filtered data
print(filtered_data)


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




