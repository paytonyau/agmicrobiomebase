
df = read.csv("16s_length_distribution.tsv", header = T, sep = ",")
df$SequenceLength <- apply(df, 1, function(x) nchar(x[['Sequence']]))
# Sort the data frame by SequenceLength
df <- df[order(df$SequenceLength), ]
# Calculate the frequency of each sequence length
df2 <- table(df$SequenceLength)
# Create a data frame for plotting
plot_df <- data.frame(SequenceLength = as.numeric(names(df2)), Frequency = as.numeric(df2))
# Calculate the accumulated percentage
plot_df$AccumulatedPercentage <- cumsum(plot_df$Frequency) / sum(plot_df$Frequency) * 100
# Filter the data frame
plot_df2 <- plot_df[plot_df$AccumulatedPercentage >= 0.5 & plot_df$AccumulatedPercentage <= 99.5, ]

# Load the ggplot2 package
library(ggplot2)

# Create a line graph
p = ggplot(plot_df2, aes(x = SequenceLength, y = AccumulatedPercentage)) +
  labs(x = "Sequence Length", y = "Accumulated Percentage (%)", 
       title = "Distribution of Accumulated Percentage of Sequence Length Frequency from 0.5% to 99.5%") +
  theme_classic() + 
  geom_line(aes(color = SequenceLength), size = 1.5) +
  theme(
    text = element_text(size = 16, colour = "black"), 
    axis.ticks = element_line(colour = "black", size = 1.1),
    axis.line = element_line(colour = 'black', size = 1.1),
    axis.text.x = element_text(colour = "black", angle = 0, hjust = 0.5, size = 13, face = "bold"),
    axis.text.y = element_text(colour = "black", angle = 0, hjust = 0.5, size = 13, face = "bold"),
    axis.title.y = element_text(color = "black", size = 20, face = "bold"), 
    axis.title.x = element_text(color = "black", size = 20, face = "bold"))

# Add vertical lines and text labels at the 2%, 50%, and 98% points
p + geom_vline(aes(xintercept = plot_df$SequenceLength[which.min(abs(plot_df$AccumulatedPercentage - 2))]), 
               linetype="dashed", size = 1.5, color = "wheat4") +
  geom_text(aes(x = plot_df$SequenceLength[which.min(abs(plot_df$AccumulatedPercentage - 2))], y = 2, 
                label = paste("2% (", plot_df$SequenceLength[which.min(abs(plot_df$AccumulatedPercentage - 2))], ")", sep = "")), vjust = 1, color = "wheat4") +
  geom_vline(aes(xintercept = plot_df$SequenceLength[which.min(abs(plot_df$AccumulatedPercentage - 50))]), linetype="dashed", size = 1.5, color = "forestgreen") +
  geom_text(aes(x = plot_df$SequenceLength[which.min(abs(plot_df$AccumulatedPercentage - 50))], y = 50, 
                label = paste("50% (", plot_df$SequenceLength[which.min(abs(plot_df$AccumulatedPercentage - 50))], ")", sep = "")), vjust = -1, color = "forestgreen") +
  geom_vline(aes(xintercept = plot_df$SequenceLength[which.min(abs(plot_df$AccumulatedPercentage - 98))]), linetype="dashed", size = 1.5, color = "coral4") +
  geom_text(aes(x = plot_df$SequenceLength[which.min(abs(plot_df$AccumulatedPercentage - 98))], y = 98, 
                label = paste("98% (", plot_df$SequenceLength[which.min(abs(plot_df$AccumulatedPercentage - 98))], ")", sep = "")), vjust = -1, color = "coral4") +
  scale_x_continuous(breaks = seq(min(plot_df$SequenceLength), max(plot_df$SequenceLength), by = 5))
