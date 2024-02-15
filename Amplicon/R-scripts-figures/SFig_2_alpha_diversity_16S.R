# List of packages
packages <- c("phyloseq", "reshape2", "ggplot2")

# Function to check and install missing packages
check_and_install <- function(pkg){
  if (!require(pkg, character.only = TRUE)){
    install.packages(pkg, dependencies = TRUE)
  }
}

# Apply the function to each package
sapply(packages, check_and_install)

# available measurements: "Observed" "Chao1" "ACE" "Shannon" "Simpson" "InvSimpson" "Fisher"
########### Fisher
# Calculate alpha diversity (Fisher) and store it in physeq.SU object
# This step calculates the Fisher alpha diversity for each sample in the physeq.SU object.
alpha.object <- cbind(
  x = sample_data(physeq.SU),
  y = estimate_richness(physeq.SU, measures = 'Fisher')
)

# Data preparation (formatting)
# This step selects the necessary columns from the alpha diversity object and reshapes the data for further analysis.
selected_columns <- alpha.object[, c("x.Soil.Location", "x.Group", "x.Type", "Fisher")]
selected_columns2 <- melt(selected_columns)
names(selected_columns2) <- c("Soil_Location","Group" , "Type", "variable", "Fisher")

# Define the comparisons
# This step defines the pairs of groups that will be compared in the subsequent statistical tests.
my_comparisons <- list(
  c("SU.CL.BO", "CO.CL.BO"),
  c("SU.CL.YO", "CO.CL.YO"),
  c("SU.CY.BU", "CO.CY.BU"),
  c("SU.CY.YO", "CO.CY.YO"),
  c("SU.SC.HE", "CO.SC.HE"),
  c("SU.SC.SH", "CO.SC.SH"),
  c("SU.SL.AN", "CO.SL.AN"),
  c("SU.SL.BE", "CO.SL.BE"),
  c("SU.SL.SH", "CO.SL.SH")
)

# Initialise an empty data frame to store the results
# This step creates an empty data frame where the results of the statistical tests will be stored.
results <- data.frame()

# Perform t-tests for each pair of groups
# This loop performs a Wilcoxon test for each pair of groups defined in my_comparisons and stores the results in the results data frame.
for (i in seq_along(my_comparisons)) {
  group1_data <- selected_columns2$Fisher[selected_columns2$Group == my_comparisons[[i]][1]]
  group2_data <- selected_columns2$Fisher[selected_columns2$Group == my_comparisons[[i]][2]]
  
  wilcox_test_result <- wilcox.test(group1_data, group2_data)
  
  results <- rbind(results, data.frame(
    group1 = my_comparisons[[i]][1],
    group2 = my_comparisons[[i]][2],
    p.value = wilcox_test_result$p.value
  ))
}

# Adjust the p-values for multiple comparisons using the Benjamini-Hochberg procedure
# This step adjusts the p-values obtained from the Wilcoxon tests to account for multiple comparisons.
results$p.adjusted <- p.adjust(results$p.value, method = "BH")

# Add significance levels based on the adjusted p-values
# This step adds a column to the results data frame indicating the significance level of each test based on the adjusted p-values.
results$p.signif <- symnum(results$p.adjusted, corr = FALSE, na = FALSE,
                           cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                           symbols = c("***", "**", "*", ".", " "))

# print result
print(results)

# Defind colour
my_palette <- c("darkgoldenrod", "limegreen")

# Create the boxplot
p <- ggplot(data=selected_columns2, aes(x=Soil_Location, y= Fisher, fill = Type)) + 
  geom_boxplot(size = 0.5, 
               width = 0.825, 
               color = "grey20", 
               position = position_dodge(0.9)
  ) +
  scale_fill_manual(values = my_palette) +
  labs(x = element_blank(),                     
       y = "Alpha Diversity (Fisher)"          
  ) + 
  theme_classic() + 
  theme(
    text = element_text(size = 18, colour = "black"), 
    axis.ticks = element_line(colour = "black", size = 1.1),
    axis.line = element_line(colour = 'black', size = 1.1),
    axis.text.x = element_text(colour = "black", angle = 45, hjust = 1, 
                               size = 13, face = "bold"),
    axis.text.y = element_text(angle = 0, hjust = 0, colour = "black", 
                               size = 13, face = "bold"),
    axis.title.y = element_text(color = "black", size = 15, face = "bold"),
    legend.position = "right") +
  scale_y_continuous(breaks = seq(100, 900, by = 100), limits = c(100, 900))  

# Add the results of the comparisons to the plot
for (i in seq_len(nrow(results))) {
  
  # Set the y position for the label
  y_position <- 900
  
  # Add the label to the plot
  p <- p + annotate("text", x = i, 
                    y=y_position , 
                    label=round(results$p.adjusted[i], 4), 
                    size= 3.75, face = "bold")
}

# Create a plot for alpha diversity
pdf(file = "Fig5B_alpha_Fisher.pdf", width = 8,height = 5)

# Print the plot
print(p)

# Close the PDF device and save the plot to a file
dev.off()


########## Simpson
# Calculate alpha diversity (Simpson) and store it in physeq.SU object
# This step calculates the Simpson alpha diversity for each sample in the physeq.SU object.
alpha.object <- cbind(
  x = sample_data(physeq.SU),
  y = estimate_richness(physeq.SU, measures = 'Simpson')
)

# Data preparation (formatting)
# This step selects the necessary columns from the alpha diversity object and reshapes the data for further analysis.
selected_columns <- alpha.object[, c("x.Soil.Location", "x.Group", "x.Type", "Simpson")]
selected_columns2 <- melt(selected_columns)
names(selected_columns2) <- c("Soil_Location","Group" , "Type", "variable", "Simpson")

# Define the comparisons
# This step defines the pairs of groups that will be compared in the subsequent statistical tests.
my_comparisons <- list(
  c("SU.CL.BO", "CO.CL.BO"),
  c("SU.CL.YO", "CO.CL.YO"),
  c("SU.CY.BU", "CO.CY.BU"),
  c("SU.CY.YO", "CO.CY.YO"),
  c("SU.SC.HE", "CO.SC.HE"),
  c("SU.SC.SH", "CO.SC.SH"),
  c("SU.SL.AN", "CO.SL.AN"),
  c("SU.SL.BE", "CO.SL.BE"),
  c("SU.SL.SH", "CO.SL.SH")
)

# Initialise an empty data frame to store the results
# This step creates an empty data frame where the results of the statistical tests will be stored.
results <- data.frame()

# Perform t-tests for each pair of groups
# This loop performs a Wilcoxon test for each pair of groups defined in my_comparisons and stores the results in the results data frame.
for (i in seq_along(my_comparisons)) {
  group1_data <- selected_columns2$Simpson[selected_columns2$Group == my_comparisons[[i]][1]]
  group2_data <- selected_columns2$Simpson[selected_columns2$Group == my_comparisons[[i]][2]]
  
  wilcox_test_result <- wilcox.test(group1_data, group2_data)
  
  results <- rbind(results, data.frame(
    group1 = my_comparisons[[i]][1],
    group2 = my_comparisons[[i]][2],
    p.value = wilcox_test_result$p.value
  ))
}

# Adjust the p-values for multiple comparisons using the Benjamini-Hochberg procedure
# This step adjusts the p-values obtained from the Wilcoxon tests to account for multiple comparisons.
results$p.adjusted <- p.adjust(results$p.value, method = "BH")

# Add significance levels based on the adjusted p-values
# This step adds a column to the results data frame indicating the significance level of each test based on the adjusted p-values.
results$p.signif <- symnum(results$p.adjusted, corr = FALSE, na = FALSE,
                           cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                           symbols = c("***", "**", "*", ".", " "))

# print result
print(results)

# Defind colour
my_palette <- c("darkgoldenrod", "limegreen")

# Create the boxplot
p <- ggplot(data=selected_columns2, aes(x=Soil_Location, y= Simpson, fill = Type)) + 
  geom_boxplot(size = 0.5, 
               width = 0.825, 
               color = "grey20", 
               position = position_dodge(0.9)
  ) +
  scale_fill_manual(values = my_palette) +
  labs(x = element_blank(),                     
       y = "Alpha Diversity (Simpson)"          
  ) + 
  theme_classic() + 
  theme(
    text = element_text(size = 18, colour = "black"), 
    axis.ticks = element_line(colour = "black", size = 1.1),
    axis.line = element_line(colour = 'black', size = 1.1),
    axis.text.x = element_text(colour = "black", angle = 45, hjust = 1, 
                               size = 13, face = "bold"),
    axis.text.y = element_text(angle = 0, hjust = 0, colour = "black", 
                               size = 13, face = "bold"),
    axis.title.y = element_text(color = "black", size = 15, face = "bold"),
    legend.position = "right") +
  scale_y_continuous(breaks = seq(0.99, 1, by = 0.0025), limits = c(0.99, 1))  

# Add the results of the comparisons to the plot
for (i in seq_len(nrow(results))) {
  
  # Set the y position for the label
  y_position <- 1
  
  # Add the label to the plot
  p <- p + annotate("text", x = i, 
                    y=y_position , 
                    label=round(results$p.adjusted[i], 4), 
                    size= 3.75, face = "bold")
}

# Create a plot for alpha diversity
pdf(file = "Fig5B_alpha_Simpson.pdf", width = 8,height = 5)

# Print the plot
print(p)

# Close the PDF device and save the plot to a file
dev.off()


########## InvSimpson
# Calculate alpha diversity (InvSimpson) and store it in physeq.SU object
# This step calculates the Inverse Simpson alpha diversity for each sample in the physeq.SU object.
alpha.object <- cbind(
  x = sample_data(physeq.SU),
  y = estimate_richness(physeq.SU, measures = 'InvSimpson')
)

# Data preparation (formatting)
# This step selects the necessary columns from the alpha diversity object and reshapes the data for further analysis.
selected_columns <- alpha.object[, c("x.Soil.Location", "x.Group", "x.Type", "InvSimpson")]
selected_columns2 <- melt(selected_columns)
names(selected_columns2) <- c("Soil_Location","Group" , "Type", "variable", "InvSimpson")

# Define the comparisons
# This step defines the pairs of groups that will be compared in the subsequent statistical tests.
my_comparisons <- list(
  c("SU.CL.BO", "CO.CL.BO"),
  c("SU.CL.YO", "CO.CL.YO"),
  c("SU.CY.BU", "CO.CY.BU"),
  c("SU.CY.YO", "CO.CY.YO"),
  c("SU.SC.HE", "CO.SC.HE"),
  c("SU.SC.SH", "CO.SC.SH"),
  c("SU.SL.AN", "CO.SL.AN"),
  c("SU.SL.BE", "CO.SL.BE"),
  c("SU.SL.SH", "CO.SL.SH")
)

# Initialise an empty data frame to store the results
# This step creates an empty data frame where the results of the statistical tests will be stored.
results <- data.frame()

# Perform t-tests for each pair of groups
# This loop performs a Wilcoxon test for each pair of groups defined in my_comparisons and stores the results in the results data frame.
for (i in seq_along(my_comparisons)) {
  group1_data <- selected_columns2$InvSimpson[selected_columns2$Group == my_comparisons[[i]][1]]
  group2_data <- selected_columns2$InvSimpson[selected_columns2$Group == my_comparisons[[i]][2]]
  
  wilcox_test_result <- wilcox.test(group1_data, group2_data)
  
  results <- rbind(results, data.frame(
    group1 = my_comparisons[[i]][1],
    group2 = my_comparisons[[i]][2],
    p.value = wilcox_test_result$p.value
  ))
}

# Adjust the p-values for multiple comparisons using the Benjamini-Hochberg procedure
# This step adjusts the p-values obtained from the Wilcoxon tests to account for multiple comparisons.
results$p.adjusted <- p.adjust(results$p.value, method = "BH")

# Add significance levels based on the adjusted p-values
# This step adds a column to the results data frame indicating the significance level of each test based on the adjusted p-values.
results$p.signif <- symnum(results$p.adjusted, corr = FALSE, na = FALSE,
                           cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                           symbols = c("***", "**", "*", ".", " "))

# print result
print(results)

# Defind colour
my_palette <- c("darkgoldenrod", "limegreen")

# Create the boxplot
p <- ggplot(data=selected_columns2, aes(x=Soil_Location, y= InvSimpson, fill = Type)) + 
  geom_boxplot(size = 0.5, 
               width = 0.825, 
               color = "grey20", 
               position = position_dodge(0.9)
  ) +
  scale_fill_manual(values = my_palette) +
  labs(x = element_blank(),                     
       y = "Alpha Diversity (InvSimpson)"          
  ) + 
  theme_classic() + 
  theme(
    text = element_text(size = 18, colour = "black"), 
    axis.ticks = element_line(colour = "black", size = 1.1),
    axis.line = element_line(colour = 'black', size = 1.1),
    axis.text.x = element_text(colour = "black", angle = 45, hjust = 1, 
                               size = 13, face = "bold"),
    axis.text.y = element_text(angle = 0, hjust = 0, colour = "black", 
                               size = 13, face = "bold"),
    axis.title.y = element_text(color = "black", size = 15, face = "bold"),
    legend.position = "right") +
  scale_y_continuous(breaks = seq(100, 900, by = 100), limits = c(100, 900))  

# Add the results of the comparisons to the plot
for (i in seq_len(nrow(results))) {
  
  # Set the y position for the label
  y_position <- 900
  
  # Add the label to the plot
  p <- p + annotate("text", x = i, 
                    y=y_position , 
                    label=round(results$p.adjusted[i], 4), 
                    size= 3.75, face = "bold")
}

# Create a plot for alpha diversity
pdf(file = "Fig5B_alpha_InvSimpson.pdf", width = 8,height = 5)

# Print the plot
print(p)

# Close the PDF device and save the plot to a file
dev.off()
