#' Function facilitates using partial resampling for regression analysis on methylation datasets
#' @param x getData function output
#' @param response Character indicating the column name your response variable of interest in info file (i.e "Disease")
#' @param numReps Integer indicating umber of boostrap replicates to perform
#' @param method Character indicating which glm family to use ("binomial", "gaussian", "poisson", etc.)
#' @param testFraction Integer between 0 and 1 indicating the fraction of samples to use for each bootstrap replicate
#' @import caret
#' @return Returns a dataframe with each column representing p_values for a cpg and each row indicating the replicate number
#' @export

consensusRegress=function(x, response, predictors=NULL, numReps, method, testFraction=0.8){ 
  
  if(is.null(predictors)){
    #merge pertinent info file columns with data file for analysis
    toRegress = merge( x[[1]][,response, drop=F], x[[2]], by.x = "row.names", by.y = "row.names" )
    #adjust remove extra row names column from merging
    toRegress <- data.frame(toRegress[,-1], row.names = toRegress[,1]) 
    #empty frame to fill
    result<-data.frame()
    for(i in 1:numReps){
      #print status
      print(paste("working on rep ", i, sep=""))
      #new seed for new partition
      set.seed(i)
      #generate specified data partition with caret package createDataPartition function
      inTrain<-createDataPartition(y = toRegress[,response], p = testFraction, list=F )
      toTest<-toRegress[inTrain,]
      #add data partition's result to the overall dataframe
      result<-rbind(result,data.frame(lapply(toTest[,-1], function(s) summary(glm(toTest[,1]~s, family = method))$coefficients[2,4]))) 
    }
  }
  
  if(length(predictors)==1){
    #merge pertinent info file columns with data file for analysis
    toRegress = merge( x[[1]][,c(response,predictors)], x[[2]], by.x = "row.names", by.y = "row.names" )
    #adjust remove extra row names column from merging
    toRegress <- data.frame(toRegress[,-1], row.names = toRegress[,1])
    #empty frame to fill
    result<-data.frame()
    for(i in 1:numReps){
      print(paste("working on rep ", i, sep=""))
      #new seed for new partition
      set.seed(i)
      #generate specified data partition with caret package createDataPartition function
      inTrain<-createDataPartition(y = toRegress[,response], p = testFraction, list=F )
      toTest<-toRegress[inTrain,]
      #add data partition's result to the overall dataframe
      result<-rbind(result,data.frame(lapply(toTest[,-c(1:(length(predictors)+1))], function(s) summary(glm(toTest[,1]~s+toTest[,2], family = method))$coefficients[2,4])))      
    }
  }
  
  if(length(predictors)==2){
    #merge pertinent info file columns with data file for analysis
    toRegress = merge( x[[1]][,c(response,predictors)], x[[2]], by.x = "row.names", by.y = "row.names" )
    #adjust remove extra row names column from merging
    toRegress <- data.frame(toRegress[,-1], row.names = toRegress[,1])
    #empty frame to fill
    result<-data.frame()
    for(i in 1:numReps){
      print(paste("working on rep ", i, sep=""))
      #new seed for new partition
      set.seed(i)
      #generate specified data partition with caret package createDataPartition function
      inTrain<-createDataPartition(y = toRegress[,response], p = testFraction, list=F )
      toTest<-toRegress[inTrain,]
      #add data partition's result to the overall dataframe
      result<-rbind(result,data.frame(lapply(toTest[,-c(1:(length(predictors)+1))], function(s) summary(glm(toTest[,1]~s+toTest[,2]+toTest[,3], family = method))$coefficients[2,4])))      
    }
  }
  
  if(length(predictors)==3){
    #merge pertinent info file columns with data file for analysis
    toRegress = merge( x[[1]][,c(response,predictors)], x[[2]], by.x = "row.names", by.y = "row.names" )
    #adjust remove extra row names column from merging
    toRegress <- data.frame(toRegress[,-1], row.names = toRegress[,1])
    #empty frame to fill
    result<-data.frame()
    for(i in 1:numReps){
      print(paste("working on rep ", i, sep=""))
      #new seed for new partition
      set.seed(i)
      #generate specified data partition with caret package createDataPartition function
      inTrain<-createDataPartition(y = toRegress[,response], p = testFraction, list=F )
      toTest<-toRegress[inTrain,]
      #add data partition's result to the overall dataframe
      result<-rbind(result,data.frame(lapply(toTest[,-c(1:(length(predictors)+1))], function(s) summary(glm(toTest[,1]~s+toTest[,2]+toTest[,3]+toTest[,4], family = method))$coefficients[2,4])))      
    }
  }
  
  if(length(predictors)==4){
    #merge pertinent info file columns with data file for analysis
    toRegress = merge( x[[1]][,c(response,predictors)], x[[2]], by.x = "row.names", by.y = "row.names" )
    #adjust remove extra row names column from merging
    toRegress <- data.frame(toRegress[,-1], row.names = toRegress[,1])
    #empty frame to fill
    result<-data.frame()
    for(i in 1:numReps){
      print(paste("working on rep ", i, sep=""))
      #new seed for new partition
      set.seed(i)
      #generate specified data partition with caret package createDataPartition function
      inTrain<-createDataPartition(y = toRegress[,response], p = testFraction, list=F )
      toTest<-toRegress[inTrain,]
      #add data partition's result to the overall dataframe
      result<-rbind(result,data.frame(lapply(toTest[,-c(1:(length(predictors)+1))], function(s) summary(glm(toTest[,1]~s+toTest[,2]+toTest[,3]+toTest[,4]+toTest[,5], family = method))$coefficients[2,4])))      
    }
  }
  
  return(result)
}