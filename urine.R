#R 3.5.3
#11 Dec, 2019
#Assignment 2

#install the packages and call the library
install.packages('R.matlab')
install.packages('ggplot2')
install.packages('e1071')
install.packages('snowfall')
install.packages('neldermead')
install.packages('optimbase')
install.packages('classyfire') #install it from the file in the directory
install.packages('factoextra')
library(R.matlab)
library("classyfire")
library(ggplot2)
library(factoextra)

#Test on urine                                    
#importing the data
bgw_ur <- readMat('BWG_UR_CDvCTRL.mat')
#view the data elements stored in the bgw object
attributes(bgw_ur)

X <- bgw_ur$XTIC
y <- bgw_ur$CLASS

#assign the rownames to X
rownames(X) <- unlist(bgw_ur$SAM)

#check the dimension 
dim(X)
#22 row and 4349 columns

#class frequencies
table(y)
#there are 14 belonging to 1 (Control)
#there are 8 belonging to 2 (Crohn's disease)


#Exploratory data analysis #####

#plot the total ion chormatogram (TIC) of samples
#Get the retention times
RT <- bgw_ur$RT

#Select a row number to investigate as "i" 
i = 1
plot(RT, X[i,], type = 'l', xlab = 'retention time(min)', ylab = 'intensity',
     main=sprintf("TIC profile for sample: %s", rownames(X)[i]))

#Exercise 2 ######
#PCA
Xpca <- prcomp(X, scale = TRUE)
summ <- summary(Xpca)
summ
#the first 4 pca captured about 91% of the total variance

#cumulative pca
summ$importance[3,]
#the first 4 pca captured about 91% of the total variance

#score plots for the first PCs
plot(Xpca$x[,1:2], col=c('red','green')[y], pch=20)
legend('topright',legend=c('Control', "CROHN'S"), col=c('red', 'green')[y], pch=c(20,20))

#screeplot
fviz_eig(Xpca)
fviz_pca_biplot(Xpca, repel = TRUE, col.var = "blue", col.ind = "red")

#classyfire ######
#build a classification ensemble of SVM's
library(classyfire)
#build a classification of SVM'm using classyfire's cfBuild()
ens <- cfBuild(Xpca$x[, 1:3], y, cpus = 2, bootNum = 50, ensNum = 50)

#retrieve a list of returned fields
attributes(ens)

#investigate the individual predictions of each optimised SVM model in the
#to investigate the individual predictions of each optimised SVM model in the
#classification ensemble
ens$testAcc

#investigating the overall percentage of correctly classified test object
getAvgAcc(ens)$Train
#mean of accuracy = 61.75

#overall running time
ens$totalTime
# user  system elapsed 
#0.120   0.153 139.851 

#to find out which samples proved particularly challenging to identify
getConfMatr(ens)
#the graphical representation of it
ggClassPred(ens, displayAll = TRUE, fillBrewer = TRUE, showText = TRUE)

#to find out how many SVM's are needed in the classification ensemble
ggEnsTrend(ens, ylim=c(50,100))

#bootstrap=50, ensemble=50, accuracy = 82.43, pc=3
#bootstrap=50, ensemble=50, pc=4, accuracy = 82.57
ens2 <- cfBuild(Xpca$x[, 1:4], y, cpus = 2, bootNum = 50, ensNum = 50)
getAvgAcc(ens2)$Train

#bootstrap=100, ensemble=100, pc=3, accuracy = 85.07
ens3 <- cfBuild(Xpca$x[, 1:3], y, cpus = 2, bootNum = 100, ensNum = 100)
getAvgAcc(ens3)$Train
ggClassPred(ens3, displayAll = TRUE, fillBrewer = TRUE, showText = TRUE)

#bootstrap=100, ensemble=50, pc=3, accuracy = 85.43
ens4 <- cfBuild(Xpca$x[, 1:3], y, cpus = 2, bootNum = 100, ensNum = 50)
getAvgAcc(ens4)$Train
ggClassPred(ens4, displayAll = TRUE, fillBrewer = TRUE, showText = TRUE)


#Permutation test ######
permObj <- cfPermute(Xpca$x[, 1:3], y, cpus = 2, bootNum = 100, ensNum = 50, permNum=20)

#Get the average accuracy
permObj$avgAcc

ensHist_ur <- density(ens4$testAcc)
plot(ensHist_ur, xlim=c(0,120), ylim = c(0.00, 0.15),xlab="Average accuracy (%)", col='#85E3FF')

permHist_ur <- density(permObj$avgAcc)
lines(permHist_ur, col='#FF9CEE')
legend('topright', legend=c('Permutation','Ensemble'), col=c('#FF9CEE','#85E3FF'), lty=c(1))
#The histogram shows that there is a lot of overlap between the permutation and our ensemble model
#This indicates that our model is not any better at prediction than by chance.

t.test(permObj$avgAcc, ens4$testAcc)
#there is no significant difference between the accuracy of the permutation test and the accuracy of our model
#p > 0.05