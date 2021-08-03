
#BREATH#####

bgw_Br <- readMat("data/breath/BWG_BR_CDvCTRL.mat")
attributes(bgw_Br)

#Extract required attributes
#The input features,X, are the variables that the classification algorithm uses to predict the outcome. 
#The output label, y, holds the information about the classes (in this case, control vs. Crohn’s disease).

X_Br <- bgw_Br$XTIC  #Independent variable
y_Br <- bgw_Br$CLASS #Dependent/target variable
rownames(X_Br) <- as.character(unlist(bgw_Br$SAM))

#Verify dimensionality and check data has been imported in correct format
dim(X_Br)

#Class frequencies
#Before applying any classification algorithm, a crucial step is to investigate how the associated labels (vector
#y) are distributed. Imbalances in distribution of labels can often lead to poor classification results for the
#“minority” class, even if the classification results for the “majority” class are very good.

table(y_Br) #1:19, 2:16

#Exploratory analysis#
# Get the retention times
RT_Br <- bgw_Br$RT
# Select a row number to investigate as "i"
i = 1
plot(RT_Br, X_Br[i,], type="l", xlab="retention time (min)", ylab="Intensity",
     main=sprintf("TIC profile for sample: %s", rownames(X_Br)[i]))

j = 10
plot(RT_Br, X_Br[j,], type="l", xlab="retention time (min)", ylab="Intensity",
     main=sprintf("TIC profile for sample: %s", rownames(X_Br)[j]))



#PCA#
Xpca_Br <-prcomp(X_Br, scale= TRUE)
attributes(Xpca_Br) #The returned attributes of the PCA function which will be used for further analysis.

summary(Xpca_Br)
s_Br <-summary(Xpca_Br)
str(s_Br) #Summary of PCA results.

#From the summary the percentage variance and cumulative variance for each successive PC can be extracted.
# Explained Variance
expVar_Br <- s_Br$importance[2,] * 100
expVar_Br

#Present results in bar graph
barplot(expVar_Br, xlab = "Principal Components", ylab = "Explained Variance (%)",col ="steelblue")

#Cumulative Variance
cumVar_Br <- s_Br$importance[3,] * 100
cumVar_Br
plot(cumVar_Br, type="o" ,col="black", pch=21, bg="blue", cex=0.8, ylab="Cumulative Variance",xlab="Principal Components")


#Plot of PCA scores can be helpful to assess the results of the PCA.
Xscores_Br   <- Xpca_Br$x
plot(Xscores_Br, xlab="PC1", ylab="PC2", pch=21, bg=c('green', 'red')[y_Br],
     cex=0.7,cex.lab=0.7, cex.axis = 0.7)
legend('topright',legend=c('Control', 'CD'), col=c('green', 'red'))


Xloadings_Br <- Xpca_Br$rotation
plot(Xloadings_Br, xlab="PC1", ylab="PC2", pch=21, bg=c('green', 'red'), cex=0.7,cex.lab=0.7, cex.axis = 0.7)


#SVM- Support vector machines#
#Classy fire

ens_Br<- cfBuild(Xpca_Br$x[, 1:6], y_Br, bootNum=20, ensNum=20)
attributes(ens_Br)

#View individual fields
ens_Br$testAcc
#Overall percentage of correctly classified test objects:
ens_acc_Br<-getAvgAcc(ens_Br)$Test #66.37%
ens_acc_Br
#Confusion matrix:
getConfMatr(ens_Br)

#Graphical
ggClassPred(ens_Br, displayAll = TRUE, fillBrewer = TRUE, showText = TRUE)
ggEnsTrend(ens_Br, ylim=c(50,100))

