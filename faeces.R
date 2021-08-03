library(R.matlab)
library(classyfire)
library(ggplot2)

rm(list=ls()) 
graphics.off()

##### Importing data ####
bgw <- readMat("data/faecal/BWG_FA_CDvCTRL.mat")
attributes(bgw)

# Extracting important data; x(independent) and y(dependent)
# X are variables that classification algorithm uses to predict outcome
# Y holds information about the classes (In this case, control vs Crohn's)
X <- bgw$XTIC
y <- bgw$CLASS

rownames(X) <- as.character(unlist(bgw$SAM))

dim(X)

table(y) # Looking at inbalances in X and Y values

# y[,1] <- as.factor(y[,1])
# 
# y=as.data.frame(y)
# 
# ggplot(y,aes(x=V1))+
#   geom_bar()

# Get the retention times
RT <- bgw$RT
i = 1

plot(RT, X[i,], type="l", xlab="retention time(min)", ylab="Intensity",
     main=sprintf("TIC profile for sample: %s", rownames(X)[i]))

#### The one with the transposed data ####
# Scaling something? I dunno
# Xsc <- scale(t(X), scale=TRUE)
# 
# # Create PCA
# Xpca <- prcomp(Xsc, scale=TRUE)
# Xpca2 <- prcomp(t(X), scale=TRUE)
# 
# # Rotations object or something
# xscore <- Xpca$rotation
# 
# ens <- cfBuild(xscore[,1:8], y, bootNum = 50, ensNum=25, cpus=3)


# x2 <- Xpca$x
# sum <- summary(Xpca)
# cumVar <- sum$importance[3,]*100
# plot(cumVar, type="o", col="black", pch=21, bg="blue", cex=0.8, ylab="Cumulative Variance", xlab="PC")

#### Untransposed #####
Xsc <- scale(X, center = TRUE, scale=TRUE) # Apply scaling on input data
Xpca2 <- prcomp((X), scale = TRUE)
attributes(Xpca2)

s = summary(Xpca2)

# Explained variance for each successive PC
expVar <- s$importance[2,] * 100
expVar

barplot(expVar, xlab = "Principal Components", ylab = "Explained Variance (%)") # Barplot for the exlained variance

# The cumulative variance
cumVar <- s$importance[3,] * 100
cumVar
plot(cumVar, type="o" ,col="black", pch=21, bg="red", ylab="Cumulative Variance",xlab="Principal Components")

Xscores <- Xpca2$x
plot(Xscores, xlab="PC1", ylab="PC2", pch=20, col = c( "red","steelblue"))


legend("topright", legend = c("Control","Crohn's Diease"),  pch =20, col = c("red","steelblue"))

Xloadings <- Xpca2$rotation
#ens <- cfBuild(xscore2[,1:8], y, bootNum = 50, ensNum=25, cpus=3)
ens <- cfBuild(Xscores[,1:4], y, bootNum = 50, ensNum=25, cpus=3)
ggEnsTrend(ens, ylim=c(0,100), showText = TRUE)
getAvgAcc(ens)$Test
ggClassPred(ens, displayAll=TRUE, fillBrewer=TRUE, showText=TRUE)

#### Permutation stuff ####
permObj <- cfPermute(inputData = Xscores[,1:4], inputClass = y, bootNum = 50, ensNum = 25, permNum = 50, parallel = TRUE, cpus = 4, type = "SOCK")
permObj$avgAcc
permObj$permList[[1]]

getPerm5Num(permObj)
ggPermHist(permObj)

# Density plot
ggPermHist(permObj, density=TRUE)

# Density plot that highlights additional descriptive statistics
ggPermHist(permObj, density=TRUE, percentiles = TRUE, mean = TRUE)
ggPermHist(permObj, density=TRUE, percentiles = FALSE, median = TRUE, mean = TRUE)

ggFusedHist(ens, permObj)




ensHist <- density(ens$testAcc) # returns the density data
plot(ensHist, ylim=c(0,0.06), xlim=c(30,120), col="blue", xlab="Average accuracy (%)", main="Accuracy density plot for Ensemble and Permutations")

permHist <- density(permObj$avgAcc)
lines(permHist, col="red")
lines()
legend("topright", legend=c("Permutation", "Ensemble"), col=c("red", "blue"), lty=c(1))

t.test(ens$testAcc, permObj$avgAcc) # P = 4.456e-11

