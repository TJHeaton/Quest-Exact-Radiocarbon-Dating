# Plot histogram of block sizes for tree-rings in IntCal20

NoAv_Extract <- read.csv("Datasets/TreeRings_IntCal20_NotMerged.csv", header = TRUE)
Individual_Measurement_Spans <- NoAv_Extract$calspan

hist(Individual_Measurement_Spans, freq = TRUE, ylab = "",
     col = 'skyblue3',
     breaks = seq(0.5, 50.5, by = 1), 
     xlab = "Number of annual rings", main = "Frequency",  las = 1, ylim = c(0, 6000))