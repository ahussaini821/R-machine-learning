#Post-Genomic Bioinformatics assignment 2

#Install and load the required package
install.packages("e1071")
install.packages("snowfall")
install.packages("neldermead")
install.packages("optimbase")
install.packages("R.matlab")

library(e1071)
library(snowfall)
library(neldermead)
library(optimbase)
library(R.matlab)
library(classyfire)

# Clear global environment and close graphics
rm(list=ls())
graphics.off()

#Exercise 1

#Importing the data
bgw <- readMat("data/blood/BWG_BL_CDvCTRL.mat")
attributes(bgw)

#extract the data that will form our input data matrix X (independent variable) and the associated class vector y (dependent or target variable)
X <- bgw$XTIC #extract individual elements by using the $ sign followed by the name of the attribute we want to retrieve
y <- bgw$CLASS
rownames(X) <- as.character(unlist(bgw$SAM))

#Checkthe imported data, and in particular verify the dimensionality
dim(X)
dim(y)

#Class frequencies
table(y)
table(X)

#Plot the total ion chromatogram (TIC) of samples

# Get the retention times
RT <- bgw$RT

# Select a row number to investigate as "i"
i = 2
plot(RT, X[i,], type="l", xlab="retention time (min)", ylab="Intensity",
     main=sprintf("TIC profile for sample: %s", rownames(X)[i]))

#Data pre-processing
Xsc <- scale(X, center = TRUE, scale = TRUE) # each column has a mean of 0 and a variance of 1

#Feature extraction & dimensionality reduction
colour <- c(rep("#ba91e5",18), rep("#b2edae",14)) #first 18 are purple (all the controls), second 14 are green (all the cases)

Xpca <- prcomp((X), scale=TRUE)

Xloadings <- Xpca$rotation
plot(Xloadings, xlab = "PC1", ylab="PC2", pch=21, bg=colour, cex=0.7, cex.lab=0.7, cex.axis=0.7)

#Exercise 2

s <- summary(Xpca)

#explain variance
expVar <- s$importance[2,]*100
expVar

barplot(expVar, xlab = "Principal Components", ylab = "Explained Variance (%)", col = "#e56d9d")

#cumulative variance
cumVar <- s$importance[3,]*100
cumVar

plot(cumVar, type = "o", col = "black", pch=21, bg="skyblue", cex=0.8, ylab="Cumulative vairance", xlab = "Prinicpal Components")

#scored plot for the first 2 PC's
Xpca <- prcomp((X), scale = TRUE)
Xscores <- Xpca$x
plot(Xscores, xlab="PC1", ylab="PC2", pch=21, bg=colour, cex=0.7,
     cex.lab=0.7, cex.axis=0.7)

#classyfire R package
#Building a classification ensemble with classyfire
ens <- cfBuild(inputData = Xscores[,1:9], inputClass = y, bootNum=50, ensNum=50)#contains all the individual classification models that form the ensemble
#Xloadings[,1:9] as it explains over 95% of the variance

attributes(ens)

#view the individual fields
ens$testAcc

#view the overall running time (in sec) of our analysis
ens$totalTime

#get a high level view of how difficult this particular classification problem was
getAvgAcc(ens)$Test

#if we want to know which samples proved particularly challenging to identify
getConfMatr(ens)

#produce a graphical representation of this result
ggClassPred(ens, displayAll = TRUE, fillBrewer = TRUE, showText = TRUE)

#Permutation testing on our model (testing if it works with a permNum of 5)
#The cfPermute function performs permutation testing on a classification emsemble produced by cfBuild.
#Does the cfBuild perform better than random chance
perm <- cfPermute(inputData = Xscores[,1:9], inputClass = y, bootNum=50, ensNum=50, permNum = 5)

#To look at the average test accuracy across all ensembles within each permutation iteration
perm$avgAcc

#Shows the average accuracy of the ensemble and permutest 
enshist <- density(ens$testAcc)
permhist <- density(perm$avgAcc)
plot(enshist, xlim = c(20, 85), ylim = c(0, 0.20), col = "skyblue")
lines(permhist, col = "#bbaaff")
#they overlap means we have a lot of false positives when using our model
#our classifier is not very good

#To look at the overall execution time of permutation testing
perm$totalTime

#To look at the individual execution times for each permutation round
perm$execTime

#trying with a permNum of 20
perm <- cfPermute(inputData = Xscores[,1:9], inputClass = y, bootNum=50, ensNum=50, permNum = 20)

perm$avgAcc

enshist <- density(ens$testAcc)
permhist <- density(perm$avgAcc)
plot(enshist, xlim = c(20, 85), ylim = c(0, 0.10), col = "skyblue")
lines(permhist, col = "#bbaaff")
#they overlap means we have a lot of false positives when using our model
#our classifier is not very good