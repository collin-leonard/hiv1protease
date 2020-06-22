##########################
## Initialize workspace ##
##########################

# There are only two attributes: a 8 letter string (octamer) and a label that tells whether this
# string represents a site in a peptide (or protein) where the HIV-1 protease cleaves it (+1 if yes, -1 if no). 
# The 8 letter string can also be viewed as 8 independent attributes. It is common to map the 8 letters
# to an orthogonal binary representation; a matrix with 160 elements.

libraries <- c("dplyr","tidyverse", "lubridate", "Matrix", 
               "dslabs", "caret", "reshape2", "randomForest", 
               "data.table", "gam", "ggplot2")

xfun::pkg_attach(libraries)
rm(libraries)

setwd("~/Projects/HIV")

##################################
## Load and manipulate raw data ##
##################################

load_dat <- function(name_vector){
  # Initialize empty data frame
  dat <- data.frame(octamer = character(), flag = factor())
  
  # For loop takes names of each file, extracts the data, and fills dat with it 
  for (name in name_vector) {
    dat_name <- paste("./data/",name, sep = "")
    temp_dat <- read.delim2(dat_name, header = FALSE,sep = ",", col.names = c("octamer","flag"))
    dat <- rbind(dat, temp_dat)
  }
  
 # Returns only the unique values of dat, (leaves behind octamers that produce conflicting results in different studies)  
 return(unique(dat, by = octamer))
}

expand_octamers <- function(data_mat){
  ## Create matrix of individual octamers
  my_frame <- as.data.frame(data_mat$octamer, stringsAsFactors = FALSE)
  split_frame <- as.data.frame(tstrsplit(my_frame[,1], '', type.convert = TRUE), stringsAsfactors = FALSE)
  name_list<- names(table(split_frame[,1]))
  
  ## Create blank matrix of octamer letters
  blank_mat <- matrix(0,nrow(data_mat)+1, length(name_list))
  blank_mat[1,] <- name_list
  colnames(blank_mat) <- name_list
  
  ## Fill in blank matrix with data 
  for (index in 1:nrow(split_frame)) {
    for (col in 1:ncol(split_frame)) {
      temp_col <- which(blank_mat[1,] == split_frame[index,col]) #
      y <- as.numeric(blank_mat[index+1,temp_col])
      blank_mat[index+1,temp_col] <- y + 1
    }
  }
  
  ## Add flags
  new_dat <- as.data.frame(blank_mat[-1,], stringsAsFactors=FALSE)
  new_dat <- cbind(data_mat$flag, new_dat)
  row_vec <- as.character(data_mat$octamer)
  rownames(new_dat) <- row_vec
  colnames(new_dat)[1] <- "flag"
  
  return(mutate_all(new_dat, function(x) as.numeric(as.character(x))))
}

data_names <- c("schillingData.txt", "impensData.txt", "746Data.txt", "1625Data.txt")

dat <- load_dat(data_names)

# Remove octamers that show up in different data sets, with different flag results
rep_ind <- which(duplicated(dat$octamer) == TRUE)
dat <- dat[-rep_ind,]

# Binarize flagging data (1 = cleavage, 0 = no cleavage)
dat$flag <- (dat$flag+1)/2

# Expand data to create octamer matrix  
exp_dat <- expand_octamers(dat)

# Clean up workspace
rm(data_names, rep_ind)

#########################
## Understand raw data ##
#########################

# Find mean

naive_mean <- mean(exp_dat$flag)

# With a mean of .171, the majority of cases exhibit no cleavage

# Calculate the correlations between each element
cor_dat <- cor(as.matrix(exp_dat))

# my_image function produces a heatmap using image(t(x))
my_image <- function(x, title, ... ){
  colors = rev(RColorBrewer::brewer.pal(9, "RdBu"))
  cols <- 1:ncol(x)
  rows <- 1:nrow(x)
  image(cols, rows, t(x[rev(rows),,drop=FALSE]), xaxt = "n", yaxt = "n",
        xlab="", ylab="",  col = colors, ...)
  abline(h=rows + 0.5, v = cols + 0.5)
  axis(side = 2, labels = FALSE, tick = FALSE, las = 2)
  title(main = title)
}

png(file = "./figs/heat_fig.png", width = 600, height = 400)
my_image(cor_dat, "Correlation between aminos and flag")
dev.off()

# cor_dat = the relationship between aminos and the flag result
cor_dat <- cor_dat[2:21,1]
cor_dat <- as.data.frame(cor_dat)

# The top 5 letters cor with flag (G with .187)
top_cors <- top_n(abs(cor_dat),5)

cor_dat <- cbind(rownames(cor_dat),cor_dat)
colnames(cor_dat) <- c("amino","cor")  
cor_dat <- cor_dat %>% mutate(amino = reorder(amino,cor))

# Plot cor relationships
p1 <- ggplot(cor_dat, aes(x = amino, cor)) +
  geom_point(color = "#00AFBB", size = 3) +
  geom_abline(slope = 0, intercept = 0, color = "grey") +
  ggtitle("Correlation between each Amino and the Resulting Cleavage") +
  scale_y_continuous(name = "Aminos") +
  scale_x_discrete(name = "Correlation") +
  theme_dark()

# Save file

png(file = "./figs/cor_fig.png", width = 600, height = 400)
p1
dev.off()

# Clean up workspace
rm(p1)

# Build data frame of sums of amino occurences

dat_1 <- exp_dat %>% filter(flag == 1)
dat_0 <- exp_dat %>% filter(flag == 0)

dat_1_sums <- as.data.frame(colSums(dat_1)[2:ncol(dat_1)])
dat_1_sums <- cbind(rownames(dat_1_sums), dat_1_sums,1)
colnames(dat_1_sums) <- c("amino","height","flag")

dat_0_sums <- as.data.frame(colSums(dat_0)[2:ncol(dat_0)])
dat_0_sums <- cbind(rownames(dat_0_sums), dat_0_sums,0)
colnames(dat_0_sums) <- c("amino","height","flag")

dat_sums <- full_join(dat_0_sums,dat_1_sums)
dat_sums$flag <- as.factor(dat_sums$flag) 

# Plot sums of occurences for the flag values
p1 <-ggplot(dat_sums, aes(amino, height))+ 
  geom_bar(stat = "identity", aes(fill = flag), position = "dodge2") +
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  xlab("Amino") + ylab("Number of Occurences") +
  ggtitle("Occurences of Aminos") +
  theme_dark()

png(file = "./figs/sums_fig.png", width = 600, height = 400)
p1
dev.off()

# Compare the ratios of occurences between

left_sums <- left_join(dat_0_sums[-3],dat_1_sums[-3], by="amino")
colnames(left_sums)[2:3] <- c("flag_0","flag_1")

left_sums <- left_sums %>% mutate(ratio = flag_1/flag_0) %>% mutate(resid = ratio - nrow(dat_1)/nrow(dat_0)) %>% mutate(amino = reorder(amino,resid))

# Plot Occurence residuals

p1 <- ggplot(left_sums, aes(x = amino, y =resid)) +
  geom_point(color = "#00AFBB", size = 3) +
  ggtitle("Difference between Expected and Observed Occurences in Aminos") +
  geom_abline(slope = 0, intercept = 0, color = "grey") +
  scale_y_continuous(name = "Residuals") +
  scale_x_discrete(name = "Aminos") +
  theme_dark()

png(file = "./figs/resids_fig.png", width = 600, height = 400)
p1
dev.off()

left_sums <- left_sums %>% mutate(resid = abs(resid))
top_resids <- top_n(left_sums, 5, wt = resid)

# Clean up workspace
rm(dat_0,dat_0_sums,dat_1,dat_1_sums,p1)

######################################
## Developing and Testing Algorithms##
######################################

# Create training and test sets from data set
exp_dat$flag <- factor(as.character(exp_dat$flag)) 

set.seed(719)
test_index <- createDataPartition(exp_dat$flag, times = 1, p = 0.2, list = FALSE)
test_set <- exp_dat[test_index,]
train_set <- exp_dat[-test_index,]

# Train 6 different models, and compare to determine the best for this application

models <- c("glm", "lda", "knn", "gamLoess", "multinom", "rf")

# Generate six models
fits <- lapply(models, function(model){ 
  print(model)
  train(flag ~ ., method = model, data = train_set)
}) 

names(fits) <- models

# Use fits to generate predictions
pred <- sapply(fits, function(object) 
  predict(object, newdata = test_set))

pred <- as.data.frame(pred)

# Accuracy of the models
acc <- colMeans(pred == test_set$flag)

# Calculate the F2 values for each model
F_vals <- cbind(models,vector(mode = "numeric", length = length(models)))
colnames(F_vals)[2] <- "F_vals"
F_vals <- as.data.frame(F_vals)
F_vals$F_vals <- as.numeric(F_vals$F_vals)

for (i in 1:length(models)){
  my_pred <- pred[,i]
  F_vals$F_vals[i] <- F_meas(my_pred,test_set$flag, beta = 2)
  rm(my_pred)
}

# Plot F2 Values
p1 <- ggplot(F_vals, aes(x = models, y = F_vals)) +
  geom_point(color = "#00AFBB", size = 3) +
  ggtitle("F2 values for the 6 Models used") +
  scale_y_continuous(name = "F2 values") +
  scale_x_discrete(name = "Models") +
  theme_dark()

# Save the plot

png(file = "./figs/F2s_fig.png", width = 600, height = 400)
p1
dev.off()

# Train knn with cross validation

set.seed(719)
train_knn <- train(flag ~ ., method = "knn", 
                   data = train_set,
                   tuneGrid = data.frame(k = seq(5, 13, 1)))

p1 <- ggplot(train_knn, highlight = TRUE) +   
  ggtitle("Accuracy at Varying K values for KNN Model") +
  scale_y_continuous(name = "Accuracy (Bootstrap)") +
  scale_x_continuous(name = "K Neighbors") +
  theme_dark()

png(file = "./figs/knn_fig.png", width = 600, height = 400)
p1
dev.off()

train_knn$bestTune

confusionMatrix(predict(train_knn, test_set, type = "raw"),
                test_set$flag)

# Train an random forest model

nodesize <- seq(1, 71, 10)
acc <- sapply(nodesize, function(ns){
  train(flag ~ ., method = "rf", data = train_set,
        tuneGrid = data.frame(mtry = 2),
        nodesize = ns)$results$Accuracy
})

node_acc <- as.data.frame(cbind(nodesize,acc))
colnames(node_acc) <- c("Nodesize", "Accuracy")

p1 <- ggplot(node_acc, aes(x = Nodesize, y = Accuracy)) +
  geom_point(color = "#00AFBB", size = 3) +
  ggtitle("Nodesize vs Accuracy in RF Model") +
  scale_y_continuous(name = "Accuracy") +
  scale_x_continuous(name = "Nodesize") +
  theme_dark()

png(file = "./figs/rf_fig.png", width = 600, height = 400)
p1
dev.off()

train_rf_2 <- randomForest(flag ~ ., data=train_set,
                           nodesize = nodesize[which.max(acc)])

rf_preds <- predict(train_rf_2, test_set)

confusionMatrix(rf_preds, test_set$flag)
