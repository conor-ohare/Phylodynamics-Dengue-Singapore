
## install.packages("devtools")
## devtools::install_github("aezarebski/timtamslam")
library(timtamslamR)

## Set up the sequence data

setwd("/Users/conor/NUS/dengue/timtam/new_timtam")

x <- read_fasta("sample-sequences.fasta")

png("tutorial-2-sequence-dates.png", width=600, height=450)
plot_dates(x)
dev.off()

y <- rename_dates_to_times_a(x)

png("tutorial-2-sequence-times.png", width=600, height=450)
plot_times(y)
dev.off()

write_fasta(y, "sample-sequences-new-times.fasta")

## Set up the occurrence data

p <- get_present(x, y)

z <- read.csv("aggregated-occurrences-updated.csv")
z <- rename_time_series(p, z)
write.csv(z, "aggregated-occurrences-new-times.csv",
          row.names=FALSE)

sink("disaster-strings.txt")
print("Here are the disaster sizes:\n")
paste(z$size, sep = "", collapse = " ")
print("Here are the backward-times of the disasters:\n")
paste(z$bwd_times, sep = "", collapse = " ")
sink()
