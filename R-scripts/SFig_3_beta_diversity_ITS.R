# List of packages
packages <- c("phyloseq", "ggplot2")

# Function to check and install missing packages
check_and_install <- function(pkg){
  if (!require(pkg, character.only = TRUE)){
    install.packages(pkg, dependencies = TRUE)
  }
}

# Apply the function to each package
sapply(packages, check_and_install)

# available measurements [c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")]
##########  Simpson
# Calculate alpha diversity (Simpson) and store it in physeq.a object
alpha.object <- cbind(
  x = sample_data(physeq.a),
  y = estimate_richness(physeq.a, measures = 'Simpson')
)

#### pdf
# Open a PDF device to save the plot
pdf(file = "Alpha_ITS_Simpson.pdf", width = 8,height = 5)

# Create a plot for alpha diversity
# This step creates a scatter plot of the Simpson alpha diversity for each group in the physeq.a object.
ggplot(data = alpha.object, aes(x = x.Group, y = Simpson, color = x.Group)) + 
  scale_color_brewer(palette = "Paired") +
  theme_classic() + 
  labs(
    x = element_blank(),                     # No x-axis label
    y = "Alpha Diversity (Simpson)"          # y-axis label
  ) + 
  geom_point(alpha = 1, position = "jitter", size = 4) + 
  geom_boxplot(alpha = 0, colour = "black", size = 0.8)+ 
  theme(
    text = element_text(size = 18, colour = "black"), 
    axis.ticks = element_line(colour = "black", size = 1.1),
    axis.line = element_line(colour = 'black', size = 1.1),
    axis.text.x = element_text(colour = "black", angle = 45, 
                               hjust = 1, size = 13, face = "bold"),
    axis.text.y = element_text(angle = 0, hjust = 0, 
                               colour = "black", size = 13, face = "bold"),
    axis.title.y = element_text(color = "black", size = 15, face = "bold"),
    legend.position = "none"                 # Hide legend
  ) +  
  scale_y_continuous(breaks = seq(0.95, 0.99, by = 0.01), limits = c(0.95, 0.99))  

# Close the PDF device and save the plot to a file
dev.off()

########## InvSimpson
# Calculate alpha diversity (InvSimpson) and store it in physeq.a object
alpha.object <- cbind(
  x = sample_data(physeq.a),
  y = estimate_richness(physeq.a, measures = 'InvSimpson')
)

#### pdf
# Open a PDF device to save the plot
pdf(file = "Alpha_ITS_InvSimpson.pdf", width = 8,height = 5)

# Create a plot for alpha diversity
# This step creates a scatter plot of the Inverse Simpson alpha diversity for each group in the physeq.a object.
ggplot(data=alpha.object, aes(x=x.Group, y=InvSimpson, color=x.Group)) + 
  scale_color_brewer(palette="Paired") +
  theme_classic() + 
  labs(
    x=element_blank(),                     # No x-axis label
    y="Alpha Diversity (InvSimpson)"       # y-axis label
  ) + 
  geom_point(alpha=1, position="jitter", size=4) + 
  geom_boxplot(alpha=0, colour="black", size=0.8)+ 
  theme(
    text=element_text(size=18, colour="black"), 
    axis.ticks=element_line(colour="black", size=1.1),
    axis.line=element_line(colour='black', size=1.1),
    axis.text.x=element_text(colour="black", angle=45,
                             hjust=1, size=13, face="bold"),
    axis.text.y=element_text(angle=0, hjust=0,
                             colour="black", size=13, face="bold"),
    axis.title.y=element_text(color="black", size=15, face="bold"),
    legend.position="none"                 # Hide legend
  ) +
  scale_y_continuous(breaks=seq(20,90,by=10), limits=c(20,90))  

# Close the PDF device and save the plot to a file
dev.off()


##########  Fisher
# Calculate alpha diversity (Fisher) and store it in physeq.a object
# This step calculates the Fisher alpha diversity for each sample in the physeq.a object.
alpha.object <- cbind(
  x = sample_data(physeq.a),
  y = estimate_richness(physeq.a, measures = 'Fisher')
)

#### pdf
# Open a PDF device to save the plot
pdf(file = "Alpha_ITS_Fisher.pdf", width = 8,height = 5)

# Create a plot for alpha diversity
# This step creates a scatter plot of the Fisher alpha diversity for each group in the physeq.a object.
ggplot(data=alpha.object, aes(x=x.Group, y=Fisher, color=x.Group)) + 
  scale_color_brewer(palette="Paired") +
  theme_classic() + 
  labs(
    x=element_blank(),                     # No x-axis label
    y="Alpha Diversity (Fisher)"           # y-axis label
  ) + 
  geom_point(alpha=1, position="jitter", size=4) + 
  geom_boxplot(alpha=0, colour="black", size=0.8)+ 
  theme(
    text=element_text(size=18, colour="black"), 
    axis.ticks=element_line(colour="black", size=1.1),
    axis.line=element_line(colour='black', size=1.1),
    axis.text.x=element_text(colour="black", angle=45,
                             hjust=1, size=13, face="bold"),
    axis.text.y=element_text(angle=0, hjust=0,
                             colour="black", size=13, face="bold"),
    axis.title.y=element_text(color="black", size=15, face="bold"),
    legend.position="none"                 # Hide legend
  ) +
  scale_y_continuous(breaks=seq(20,70,by=10), limits=c(20,70))  

# Close the PDF device and save the plot to a file
dev.off()

