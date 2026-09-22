
# Library to interest
library(ggplot2) # for nice figures

# Define directory to work
save_dir = "/shared/projects/2538_eb3i_n1_2025/atelier_scrnaseq/cours_intro_rmd/"

# Set as working directory
setwd(save_dir)
getwd() # check

# Create a dataframe
my_data = data.frame(A = c(1,1,2,4,4,4,6),
                     B = c(2,2,1,3,4,5,6))
my_data

# Explore dataframe
dim(my_data)

summary(my_data)

# Make a histogram
hist(my_data$B)

# Make a beautiful histogram
ggplot(my_data, aes(x = B)) +
  geom_histogram()

# Prepare output path to save file
filename_to_save = paste0(save_dir, "my_data.csv")

# Save data
write.csv(my_data, file = filename_to_save)
