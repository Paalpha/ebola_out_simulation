# Algorithmic framework for MCAR 
Model.comMCAR<-function(rep=1,rd,p,pw,dat){
  #rep is the number of iterations of the model
  #rd is the % missingness generate in the simulated data
  #p is the proportion simulated data used for training
  #pw is the % of the simulated data used for the analysis
  #options(java.parameters="-Xmx12g")
  #Samplesize of the simulated data
  u.sex<-sample(1:nrow(dat), round(pw*nrow(dat)), replace=FALSE)
  dat<-dat[u.sex,]
  CFR_Obs<-1-mean(dat$dead==1,na.rm=TRUE)#"Ground truth"
  #levelofrandomness=rd
  sizeofdata<- nrow(dat)
  #for (i in 1){ # at some point when working well, for (i in 1:100){
  for(j in  c(1:ncol(dat))){
    # intorduce missingness
    # do the following for the outcome and the covariates
    fz <- which(runif(n=sizeofdata,0,1) < rd) 
    dat[fz,j] <- NA   #covariate i or outcome, here data is a datframe with outcome and covariates
    # redo thge analysis below with this new dataset
    # + include analysis with imputed data when possible -> extract CFR for the full dataset (known outcome + imputed outcome)
  }
  #Logistic Regression
  imp<-mlr::impute(dat, target = "dead", cols = list(country= imputeMode(), casedef=imputeMode(),delay=imputeMean(),hc=imputeMode(),mnths=imputeMode(),
                                                     unbld=imputeMode(),confs=imputeMode(),dbrth=imputeMode(),jntp=imputeMode(),
                                                     jnd=imputeMode(),cnjs=imputeMode(),fvr=imputeMode(),ftg=imputeMode(),
                                                     anx=imputeMode(),vmt=imputeMode(),dia=imputeMode(),hed=imputeMode(),
                                                     mus=imputeMode(),chp=imputeMode(),agecat=imputeHist()))
  data2<-imp$data
  
  output<-list()
  for(i in 1:rep){ 
    nonmiss.data<-data2[which(!(is.na(data2$dead))),]
    miss.data<-data2[which((is.na(data2$dead))),]
    #nonmiss.dataC<-nonmiss.data[complete.cases(nonmiss.data),]#because glm doesn't inherently handle missingness in the covariates 
    #data1<-simulated_data
    #Partitioning data into training and validation (testing) set. 
    u.training <-sample(1:nrow(nonmiss.data), round(p*nrow(nonmiss.data)), replace=FALSE)
    #u.training <-createDataPartition(nonmiss.data$dead, p=0.8, list=FALSE) 
    training.data<-nonmiss.data[u.training,]
    #testing.dat$mnths <- factor(testing.dat$mnths, levels = levels(training.dat$mnths))
    #I assuming that the CFR outomces are missing for the testing set. Is this MCAR or MAR?
    testing.data<-nonmiss.data[-u.training,]
    # I create the Modellling task for the alogrithms. By doing this, we dont have to do anymore data processing for each algorithm. 
    #listLearners("classif", check.packages = TRUE, properties = "missings")[c("class","package")]
    #summarizeColumns(training.data)
    #summarizeColumns(testing.data)
    #Recreate training and testing task
    trainTask <- makeRegrTask(data = training.data,target = "dead")
    testTask <- makeRegrTask(data = testing.data, target = "dead")
    removeConstantFeatures(trainTask)
    removeConstantFeatures(testTask)
    #Standardised Predictors
    trainTask <- normalizeFeatures(trainTask,method = "standardize")
    testTask<- normalizeFeatures(testTask,method = "standardize")
    #install.packages("fda.usc")
    #library("fda.usc")
    #Create the learner
    logistic.learner <- makeLearner("regr.glm",predict.type = "response",fix.factors.prediction = TRUE)
    #CV to optimise the the GLM
    #cv.logistic <- mlr::crossval(learner = logistic.learner,task = trainTask,iters = 2,stratify = TRUE)
    #cross validation accuracy
    #cv.logistic$aggr
    #cv.logistic$measures.test
    #Use the best model to for training 
    fmodel <- mlr::train(logistic.learner,trainTask)
    getLearnerModel(fmodel)
    #predict on test data using best model from the CV
    fpmodel <- predict(fmodel, testTask)
    preds<-fpmodel$data$response
    #Model Performance
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_GLM= accuracy$PCC[2]
    sens.Valid_GLM= accuracy$sensitivity[2]
    spec.Valid_GLM= accuracy$specificity[2]
    ROC.Valid_GLM= accuracy$AUC[2]
    #Predict on the data with missing outcomes
    fpmodel<- predict(fmodel, newdata=miss.data)
    preds.impute<-fpmodel$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_GLMPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    #CFR_Obs<-1-mean(dat$dead==1,na.rm=TRUE)#"Ground truth"
    #ANN
    #getParamSet("regr.nnet")
    g.net <- makeLearner("regr.nnet", predict.type = "response",fix.factors.prediction = TRUE)
    #specify tuning method
    rancontrol <- makeTuneControlRandom(maxit = 20L)
    #3 fold cross validation
    set_cv <- makeResampleDesc("CV",iters = 5L,stratify = FALSE)
    #parameters
    net_par<- makeParamSet(
      makeIntegerParam("size", lower = 3, upper = 10), #number of trees
      #makeIntegerParam("maxit", lower = 100, upper = 200), #depth of tree
      makeNumericParam("decay", lower = 1e-08, upper =0.1),
      makeNumericParam("abstol",lower = 0.0001, upper =0.0002),
      makeNumericParam("reltol",lower = 1e-08, upper = 1e-06)
    )
    #tune parameters
    tune_net <- tuneParams(learner = g.net, task = trainTask,resampling = set_cv,par.set = net_par,control = rancontrol)
    #check CV accuracy
    #tune_net$y
    #set parameters
    final_net <- setHyperPars(learner = g.net, par.vals = tune_net$x)
    #train
    to.net <- mlr::train(final_net, trainTask)
    #test 
    pr.net <- predict(to.net, testTask)
    preds<-pr.net$data$response
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_ANN=accuracy$PCC[2]
    sens.Valid_ANN = accuracy$sensitivity[2]
    spec.Valid_ANN = accuracy$specificity[2]
    ROC.Valid_ANN=accuracy$AUC[2]
    pr.net<- predict(to.net, newdata=miss.data)
    preds.impute<-pr.net$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_ANNPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    #########################################################################################################################
    nonmiss.data<-dat[which(!(is.na(dat$dead))),]
    miss.data<-dat[which((is.na(dat$dead))),]
    u.training <-sample(1:nrow(nonmiss.data), round(p*nrow(nonmiss.data)), replace=FALSE)
    #u.training <-createDataPartition(nonmiss.data$dead, p=0.8, list=FALSE) 
    training.data<-nonmiss.data[u.training,]
    #testing.dat$mnths <- factor(testing.dat$mnths, levels = levels(training.dat$mnths))
    #I assuming that the CFR outomces are missing for the testing set. Is this MCAR or MAR?
    testing.data<-nonmiss.data[-u.training,]
    # I create the Modellling task for the alogrithms. By doing this, we dont have to do anymore data processing for each algorithm. 
    #listLearners("classif", check.packages = TRUE, properties = "missings")[c("class","package")]
    #summarizeColumns(training.data)
    #summarizeColumns(testing.data)
    #Recreate training and testing task
    trainTask <- makeRegrTask(data = training.data,target = "dead")
    testTask <- makeRegrTask(data = testing.data, target = "dead")
    removeConstantFeatures(trainTask)
    removeConstantFeatures(testTask)
    #Standardised Predictors
    trainTask <- normalizeFeatures(trainTask,method = "standardize")
    testTask<- normalizeFeatures(testTask,method = "standardize")
    ######RandomForest
    #install.packages("randomForestSRC")
    #library("randomForestSRC")
    getParamSet("regr.randomForestSRC")
    #create a learner
    rf <- makeLearner("regr.randomForestSRC", predict.type = "response", par.vals = list(ntree = 200, mtry = 3),fix.factors.prediction = TRUE)
    #rf$par.vals <- list(importance = TRUE)
    #set tunable parameters
    # random search to find optimal hyperparameters ( manual grid search will be more computationaly expensive )
    rf_param <- makeParamSet(
      makeIntegerParam("ntree",lower = 50, upper = 250),
      makeIntegerParam("mtry", lower = 3, upper = 10),
      makeIntegerParam("nodesize", lower = 10, upper = 30)
    )
    rancontrol <- makeTuneControlRandom(maxit = 20L)
    #set 10 fold cross validation
    set_cv <- makeResampleDesc("CV",iters = 5L,stratify = FALSE)
    #hyperparamter tuning
    rf_tune <- tuneParams(learner = rf, resampling = set_cv, task = trainTask, par.set = rf_param, control = rancontrol)
    #using hyperparameters for modeling
    rf.tree <- setHyperPars(rf, par.vals = rf_tune$x)
    #train a model using the tuned hyperparameters
    rforest <- mlr::train(rf.tree, trainTask)
    getLearnerModel(rforest)
    #make predictions
    #glmmodel$xlevels[["district"]] <- union(glmmodel$xlevels[["district"]], levels(testing.data$district))
    #predict on the validation set
    rfmodel <- predict(rforest, testTask)
    preds<-rfmodel$data$response
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_RF= accuracy$PCC[2]
    sens.Valid_RF= accuracy$sensitivity[2]
    spec.Valid_RF= accuracy$specificity[2]
    ROC.Valid_RF= accuracy$AUC[2]
    #Predict on the data with missing outcomes
    #miss.data<-miss.data[complete.cases(miss.data),]#Random forest cannot handle missing predictors in new data
    #imp = impute(miss.data, classes = list(integer = imputeMean(), factor = imputeMode()),dummy.classes = "integer")
    rfpmodel<- predict(rforest, newdata=miss.data)
    preds.impute<-rfpmodel$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_RFPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    #CFR_GLMObs<-1-mean(data$dead==1,na.rm=TRUE)
    ####################BRT###########################
    #The steps are similar to those for random forest, the only change is that we are tuning a BRT and using that for predicting on the validation data.
    #load GBM
    getParamSet("regr.gbm")
    g.gbm <- makeLearner("regr.gbm", predict.type = "response",fix.factors.prediction = TRUE)
    
    #specify tuning method
    rancontrol <- makeTuneControlRandom(maxit = 20L)
    #3 fold cross validation
    set_cv <- makeResampleDesc("CV",iters = 5L,stratify = FALSE)
    #parameters
    gbm_par<- makeParamSet(
      makeDiscreteParam("distribution", values = "gaussian"),
      makeIntegerParam("n.trees", lower = 100, upper = 3000), #number of trees
      makeIntegerParam("interaction.depth", lower = 2, upper = 10), #depth of tree
      makeIntegerParam("n.minobsinnode", lower = 10, upper = 80),
      makeNumericParam("bag.fraction",lower=0.5, upper = 0.75),
      makeNumericParam("shrinkage",lower = 0.001, upper = 0.05)
    )
    #tune parameters
    tune_gbm <- tuneParams(learner = g.gbm, task = trainTask,resampling = set_cv,par.set = gbm_par,control = rancontrol)
    #check CV accuracy
    #tune_gbm$y
    #set parameters
    final_gbm <- setHyperPars(learner = g.gbm, par.vals = tune_gbm$x)
    #train
    to.gbm <- mlr::train(final_gbm, trainTask)
    #test 
    pr.gbm <- predict(to.gbm, testTask)
    preds<-pr.gbm$data$response
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_BRT=accuracy$PCC[2]
    sens.Valid_BRT = accuracy$sensitivity[2]
    spec.Valid_BRT = accuracy$specificity[2]
    ROC.Valid_BRT=accuracy$AUC[2]
    
    pr.gbm<- predict(to.gbm, newdata=miss.data)
    preds.impute<-pr.gbm$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_BRTPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    #CFR_GLMObs<-1-mean(data$dead==1,na.rm=TRUE)
    
    #load BART
    getParamSet("regr.bartMachine")
    g.bart <- makeLearner("regr.bartMachine", predict.type = "response",fix.factors.prediction = TRUE)
    
    #specify tuning method
    rancontrol <- makeTuneControlRandom(maxit = 20L)
    #3 fold cross validation
    set_cv <- makeResampleDesc("CV",iters = 5L,stratify = FALSE)
    #parameters
    bart_par<- makeParamSet(
      makeIntegerParam("num_trees", lower = 5, upper = 10), #number of trees
      makeIntegerParam("num_burn_in", lower = 20, upper = 30), #depth of tree
      makeIntegerParam("num_iterations_after_burn_in", lower = 70, upper = 100), #depth of tree
      #makeIntegerParam("alpha", lower = 0.95, upper = 3),
      makeNumericParam("beta",lower = 1, upper = 2),
      makeLogicalParam("mem_cache_for_speed",default = FALSE,tunable = FALSE),
      makeLogicalParam("run_in_sample",default = FALSE,tunable = FALSE),
      makeNumericParam("k",lower = 1, upper = 2)
    )
    #tune parameters
    #set_bart_machine_num_cores(5)
    tune_bart <- tuneParams(learner = g.bart, task = trainTask,resampling = set_cv,par.set = bart_par,control = rancontrol)
    #check CV accuracy
    tune_bart$y
    #set parameters
    final_bart <- setHyperPars(learner = g.bart, par.vals = tune_bart$x)
    #train
    to.bart <-  mlr::train(final_bart, trainTask)
    #test 
    pr.bart <- predict(to.bart, testTask)
    preds<-pr.bart$data$response
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_BART=accuracy$PCC[2]
    sens.Valid_BART = accuracy$sensitivity[2]
    spec.Valid_BART = accuracy$specificity[2]
    ROC.Valid_BART=accuracy$AUC[2]
    
    pr.bart<- predict(to.bart, newdata=miss.data)
    preds.impute<-pr.bart$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_BARTPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    
    output[[i]]<- c(p=p,pw=pw,rd=rd,CFR_Obs=c(CFR_Obs),
                    PCC.Valid_GLM=c(PCC.Valid_GLM),
                    sens.Valid_GLM=c(sens.Valid_GLM),
                    spec.Valid_GLM=c(spec.Valid_GLM),
                    ROC.Valid_GLM=c(ROC.Valid_GLM),
                    CFR_GLMPred=c(CFR_GLMPred),
                    PCC.Valid_RF=c(PCC.Valid_RF),
                    sens.Valid_RF=c( sens.Valid_RF),
                    spec.Valid_RF=c(spec.Valid_RF),
                    ROC.Valid_RF=c(ROC.Valid_RF),
                    CFR_RFPred=c(CFR_RFPred),
                    PCC.Valid_BRT=c(PCC.Valid_BRT),
                    sens.Valid_BRT=c(sens.Valid_BRT),
                    spec.Valid_BRT=c(spec.Valid_BRT),
                    ROC.Valid_BRT=c(ROC.Valid_BRT),
                    CFR_BRTPred=c(CFR_BRTPred),
                    PCC.Valid_BART=c(PCC.Valid_BART),
                    sens.Valid_BART=c(sens.Valid_BART),
                    spec.Valid_BART=c( spec.Valid_BART),
                    ROC.Valid_BART=c(ROC.Valid_BART),
                    CFR_BARTPred=c(CFR_BARTPred),
                    PCC.Valid_ANN=c(PCC.Valid_ANN),
                    sens.Valid_ANN=c(sens.Valid_ANN),
                    spec.Valid_ANN=c(spec.Valid_ANN),
                    ROC.Valid_ANN=c( ROC.Valid_ANN),
                    CFR_ANNPred=c(CFR_ANNPred))
  }
  return(output)
}      

# Algorithmic framework for MNAR

Model.comMNAR<-function(rep=1,rd,p,pw,dat){
  #rep is the number of iterations of the model
  #rd is the % missingness generate in the simulated data
  #p is the proportion simulated data used for training
  #pw is the % of the simulated data used for the analysis
  #options(java.parameters="-Xmx12g")
  #Samplesize of the simulated data
  u.sex<-sample(1:nrow(dat), round(pw*nrow(dat)), replace=FALSE)
  dat<-dat[u.sex,]
  CFR_Obs<-1-mean(dat$dead==1,na.rm=TRUE)#"Ground truth"
  levelofrandomness=rd
  sizeofdat<- nrow(dat)
  #for (i in 1){ # at some point when working well, for (i in 1:100){
  #simulateMNAR for Outcome
  dx<-glm(dead~.,family = "binomial", data = dat)
  design_X <-model.matrix(~country+delay+casedef+hc+mnths+unbld+confs+
                            dbrth+jntp+jnd+cnjs+fvr+ftg+anx+vmt+dia+hed+mus+chp+agecat,data=dat)
  coef<-coefficients(dx)
  logit_CFR<- design_X%*%coef
  dat$missing<- 1/( 1 + exp(-logit_CFR))
  #simulateMNAR for Outcome
  #data$missing<-runif(n=sizeofdata,0,1) 
  data.y<-sort(dat$missing,decreasing = TRUE)
  nmar<-data.y[ceiling(rd*length(dat$missing))]
  dat$dead[dat$missing>nmar]<-NA
  #data$dead<-ifelse(data$missing>nmar,NA,data$dead)
  drop <- c("missing")
  dat<-dat[,!(names(dat) %in% drop)]
  for(j in  c(1:ncol(dat[,-21]))){
    # intorduce missingness
    # do the following for the outcome and the covariates
    fz <- which(runif(n=sizeofdat,0,1) < rd) 
    dat[fz,j] <- NA   #covariate i or outcome, here data is a datframe with outcome and covariates
    # redo thge analysis below with this new dataset
    # + include analysis with imputed data when possible -> extract CFR for the full dataset (known outcome + imputed outcome)
  }
  #Logistic Regression
  imp<-mlr::impute(dat, target = "dead", cols = list(country= imputeMode(), casedef=imputeMode(),delay=imputeMean(),hc=imputeMode(),mnths=imputeMode(),
                                                     unbld=imputeMode(),confs=imputeMode(),dbrth=imputeMode(),jntp=imputeMode(),
                                                     jnd=imputeMode(),cnjs=imputeMode(),fvr=imputeMode(),ftg=imputeMode(),
                                                     anx=imputeMode(),vmt=imputeMode(),dia=imputeMode(),hed=imputeMode(),
                                                     mus=imputeMode(),chp=imputeMode(),agecat=imputeHist()))
  data2<-imp$data
  output<-list()
  for(i in 1:rep){ 
    nonmiss.data<-data2[which(!(is.na(data2$dead))),]
    miss.data<-data2[which((is.na(data2$dead))),]
    #nonmiss.dataC<-nonmiss.data[complete.cases(nonmiss.data),]#because glm doesn't inherently handle missingness in the covariates 
    #data1<-simulated_data
    #Partitioning data into training and validation (testing) set. 
    u.training <-sample(1:nrow(nonmiss.data), round(p*nrow(nonmiss.data)), replace=FALSE)
    #u.training <-createDataPartition(nonmiss.data$dead, p=0.8, list=FALSE) 
    training.data<-nonmiss.data[u.training,]
    #testing.dat$mnths <- factor(testing.dat$mnths, levels = levels(training.dat$mnths))
    #I assuming that the CFR outomces are missing for the testing set. Is this MCAR or MAR?
    testing.data<-nonmiss.data[-u.training,]
    # I create the Modellling task for the alogrithms. By doing this, we dont have to do anymore data processing for each algorithm. 
    #listLearners("classif", check.packages = TRUE, properties = "missings")[c("class","package")]
    #summarizeColumns(training.data)
    #summarizeColumns(testing.data)
    #Recreate training and testing task
    trainTask <- makeRegrTask(data = training.data,target = "dead")
    testTask <- makeRegrTask(data = testing.data, target = "dead")
    removeConstantFeatures(trainTask)
    removeConstantFeatures(testTask)
    #Standardised Predictors
    trainTask <- normalizeFeatures(trainTask,method = "standardize")
    testTask<- normalizeFeatures(testTask,method = "standardize")
    #install.packages("fda.usc")
    #library("fda.usc")
    #Create the learner
    logistic.learner <- makeLearner("regr.glm",predict.type = "response",fix.factors.prediction = TRUE)
    #CV to optimise the the GLM
    #cv.logistic <- mlr::crossval(learner = logistic.learner,task = trainTask,iters = 2,stratify = TRUE)
    #cross validation accuracy
    #cv.logistic$aggr
    #cv.logistic$measures.test
    #Use the best model to for training 
    fmodel <- mlr::train(logistic.learner,trainTask)
    getLearnerModel(fmodel)
    #predict on test data using best model from the CV
    fpmodel <- predict(fmodel, testTask)
    preds<-fpmodel$data$response
    #Model Performance
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_GLM= accuracy$PCC[2]
    sens.Valid_GLM= accuracy$sensitivity[2]
    spec.Valid_GLM= accuracy$specificity[2]
    ROC.Valid_GLM= accuracy$AUC[2]
    #Predict on the data with missing outcomes
    fpmodel<- predict(fmodel, newdata=miss.data)
    preds.impute<-fpmodel$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_GLMPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    #CFR_Obs<-1-mean(dat$dead==1,na.rm=TRUE)#"Ground truth"
    #ANN
    #getParamSet("regr.nnet")
    g.net <- makeLearner("regr.nnet", predict.type = "response",fix.factors.prediction = TRUE)
    #specify tuning method
    rancontrol <- makeTuneControlRandom(maxit = 20L)
    #3 fold cross validation
    set_cv <- makeResampleDesc("CV",iters = 5L,stratify = FALSE)
    #parameters
    net_par<- makeParamSet(
      makeIntegerParam("size", lower = 3, upper = 10), #number of trees
      #makeIntegerParam("maxit", lower = 100, upper = 200), #depth of tree
      makeNumericParam("decay", lower = 1e-08, upper =0.1),
      makeNumericParam("abstol",lower = 0.0001, upper =0.0002),
      makeNumericParam("reltol",lower = 1e-08, upper = 1e-06)
    )
    #tune parameters
    tune_net <- tuneParams(learner = g.net, task = trainTask,resampling = set_cv,par.set = net_par,control = rancontrol)
    #check CV accuracy
    #tune_net$y
    #set parameters
    final_net <- setHyperPars(learner = g.net, par.vals = tune_net$x)
    #train
    to.net <- mlr::train(final_net, trainTask)
    #test 
    pr.net <- predict(to.net, testTask)
    preds<-pr.net$data$response
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_ANN=accuracy$PCC[2]
    sens.Valid_ANN = accuracy$sensitivity[2]
    spec.Valid_ANN = accuracy$specificity[2]
    ROC.Valid_ANN=accuracy$AUC[2]
    pr.net<- predict(to.net, newdata=miss.data)
    preds.impute<-pr.net$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_ANNPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    #########################################################################################################################
    nonmiss.data<-dat[which(!(is.na(dat$dead))),]
    miss.data<-dat[which((is.na(dat$dead))),]
    u.training <-sample(1:nrow(nonmiss.data), round(p*nrow(nonmiss.data)), replace=FALSE)
    #u.training <-createDataPartition(nonmiss.data$dead, p=0.8, list=FALSE) 
    training.data<-nonmiss.data[u.training,]
    #testing.dat$mnths <- factor(testing.dat$mnths, levels = levels(training.dat$mnths))
    #I assuming that the CFR outomces are missing for the testing set. Is this MCAR or MAR?
    testing.data<-nonmiss.data[-u.training,]
    # I create the Modellling task for the alogrithms. By doing this, we dont have to do anymore data processing for each algorithm. 
    #listLearners("classif", check.packages = TRUE, properties = "missings")[c("class","package")]
    #summarizeColumns(training.data)
    #summarizeColumns(testing.data)
    #Recreate training and testing task
    trainTask <- makeRegrTask(data = training.data,target = "dead")
    testTask <- makeRegrTask(data = testing.data, target = "dead")
    removeConstantFeatures(trainTask)
    removeConstantFeatures(testTask)
    #Standardised Predictors
    trainTask <- normalizeFeatures(trainTask,method = "standardize")
    testTask<- normalizeFeatures(testTask,method = "standardize")
    ######RandomForest
    #install.packages("randomForestSRC")
    #library("randomForestSRC")
    getParamSet("regr.randomForestSRC")
    #create a learner
    rf <- makeLearner("regr.randomForestSRC", predict.type = "response", par.vals = list(ntree = 200, mtry = 3),fix.factors.prediction = TRUE)
    #rf$par.vals <- list(importance = TRUE)
    #set tunable parameters
    # random search to find optimal hyperparameters ( manual grid search will be more computationaly expensive )
    rf_param <- makeParamSet(
      makeIntegerParam("ntree",lower = 50, upper = 300),
      makeIntegerParam("mtry", lower = 3, upper = 10),
      makeIntegerParam("nodesize", lower = 10, upper = 30)
    )
    rancontrol <- makeTuneControlRandom(maxit = 20L)
    #set 10 fold cross validation
    set_cv <- makeResampleDesc("CV",iters = 5L,stratify = FALSE)
    #hyperparamter tuning
    rf_tune <- tuneParams(learner = rf, resampling = set_cv, task = trainTask, par.set = rf_param, control = rancontrol)
    #using hyperparameters for modeling
    rf.tree <- setHyperPars(rf, par.vals = rf_tune$x)
    #train a model using the tuned hyperparameters
    rforest <- mlr::train(rf.tree, trainTask)
    getLearnerModel(rforest)
    #make predictions
    #glmmodel$xlevels[["district"]] <- union(glmmodel$xlevels[["district"]], levels(testing.data$district))
    #predict on the validation set
    rfmodel <- predict(rforest, testTask)
    preds<-rfmodel$data$response
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_RF= accuracy$PCC[2]
    sens.Valid_RF= accuracy$sensitivity[2]
    spec.Valid_RF= accuracy$specificity[2]
    ROC.Valid_RF= accuracy$AUC[2]
    #Predict on the data with missing outcomes
    #miss.data<-miss.data[complete.cases(miss.data),]#Random forest cannot handle missing predictors in new data
    #imp = impute(miss.data, classes = list(integer = imputeMean(), factor = imputeMode()),dummy.classes = "integer")
    rfpmodel<- predict(rforest, newdata=miss.data)
    preds.impute<-rfpmodel$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_RFPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    #CFR_GLMObs<-1-mean(data$dead==1,na.rm=TRUE)
    ####################BRT###########################
    #The steps are similar to those for random forest, the only change is that we are tuning a BRT and using that for predicting on the validation data.
    #load GBM
    getParamSet("regr.gbm")
    g.gbm <- makeLearner("regr.gbm", predict.type = "response",fix.factors.prediction = TRUE)
    
    #specify tuning method
    rancontrol <- makeTuneControlRandom(maxit = 20L)
    #3 fold cross validation
    set_cv <- makeResampleDesc("CV",iters = 5L,stratify = FALSE)
    #parameters
    gbm_par<- makeParamSet(
      makeDiscreteParam("distribution", values = "gaussian"),
      makeIntegerParam("n.trees", lower = 100, upper = 3000), #number of trees
      makeIntegerParam("interaction.depth", lower = 2, upper = 10), #depth of tree
      makeIntegerParam("n.minobsinnode", lower = 10, upper = 80),
      makeNumericParam("bag.fraction",lower=0.5, upper = 0.75),
      makeNumericParam("shrinkage",lower = 0.001, upper = 0.05)
    )
    #tune parameters
    tune_gbm <- tuneParams(learner = g.gbm, task = trainTask,resampling = set_cv,par.set = gbm_par,control = rancontrol)
    #check CV accuracy
    #tune_gbm$y
    #set parameters
    final_gbm <- setHyperPars(learner = g.gbm, par.vals = tune_gbm$x)
    #train
    to.gbm <- mlr::train(final_gbm, trainTask)
    #test 
    pr.gbm <- predict(to.gbm, testTask)
    preds<-pr.gbm$data$response
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_BRT=accuracy$PCC[2]
    sens.Valid_BRT = accuracy$sensitivity[2]
    spec.Valid_BRT = accuracy$specificity[2]
    ROC.Valid_BRT=accuracy$AUC[2]
    
    pr.gbm<- predict(to.gbm, newdata=miss.data)
    preds.impute<-pr.gbm$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_BRTPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    #CFR_GLMObs<-1-mean(data$dead==1,na.rm=TRUE)
    
    #load BART
    getParamSet("regr.bartMachine")
    g.bart <- makeLearner("regr.bartMachine", predict.type = "response",fix.factors.prediction = TRUE)
    
    #specify tuning method
    rancontrol <- makeTuneControlRandom(maxit = 20L)
    #3 fold cross validation
    set_cv <- makeResampleDesc("CV",iters = 5L,stratify = FALSE)
    #parameters
    bart_par<- makeParamSet(
      makeIntegerParam("num_trees", lower = 5, upper = 10), #number of trees
      makeIntegerParam("num_burn_in", lower = 20, upper = 30), #depth of tree
      makeIntegerParam("num_iterations_after_burn_in", lower = 70, upper = 100), #depth of tree
      #makeIntegerParam("alpha", lower = 0.95, upper = 3),
      makeNumericParam("beta",lower = 1, upper = 2),
      makeLogicalParam("mem_cache_for_speed",default = FALSE,tunable = FALSE),
      makeLogicalParam("run_in_sample",default = FALSE,tunable = FALSE),
      makeNumericParam("k",lower = 1, upper = 2)
    )
    #tune parameters
    #set_bart_machine_num_cores(5)
    tune_bart <- tuneParams(learner = g.bart, task = trainTask,resampling = set_cv,par.set = bart_par,control = rancontrol)
    #check CV accuracy
    tune_bart$y
    #set parameters
    final_bart <- setHyperPars(learner = g.bart, par.vals = tune_bart$x)
    #train
    to.bart <-  mlr::train(final_bart, trainTask)
    #test 
    pr.bart <- predict(to.bart, testTask)
    preds<-pr.bart$data$response
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_BART=accuracy$PCC[2]
    sens.Valid_BART = accuracy$sensitivity[2]
    spec.Valid_BART = accuracy$specificity[2]
    ROC.Valid_BART=accuracy$AUC[2]
    
    pr.bart<- predict(to.bart, newdata=miss.data)
    preds.impute<-pr.bart$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_BARTPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    
    output[[i]]<- c(p=p,pw=pw,rd=rd,CFR_Obs=c(CFR_Obs),
                    PCC.Valid_GLM=c(PCC.Valid_GLM),
                    sens.Valid_GLM=c(sens.Valid_GLM),
                    spec.Valid_GLM=c(spec.Valid_GLM),
                    ROC.Valid_GLM=c(ROC.Valid_GLM),
                    CFR_GLMPred=c(CFR_GLMPred),
                    PCC.Valid_RF=c(PCC.Valid_RF),
                    sens.Valid_RF=c( sens.Valid_RF),
                    spec.Valid_RF=c(spec.Valid_RF),
                    ROC.Valid_RF=c(ROC.Valid_RF),
                    CFR_RFPred=c(CFR_RFPred),
                    PCC.Valid_BRT=c(PCC.Valid_BRT),
                    sens.Valid_BRT=c(sens.Valid_BRT),
                    spec.Valid_BRT=c(spec.Valid_BRT),
                    ROC.Valid_BRT=c(ROC.Valid_BRT),
                    CFR_BRTPred=c(CFR_BRTPred),
                    PCC.Valid_BART=c(PCC.Valid_BART),
                    sens.Valid_BART=c(sens.Valid_BART),
                    spec.Valid_BART=c( spec.Valid_BART),
                    ROC.Valid_BART=c(ROC.Valid_BART),
                    CFR_BARTPred=c(CFR_BARTPred),
                    PCC.Valid_ANN=c(PCC.Valid_ANN),
                    sens.Valid_ANN=c(sens.Valid_ANN),
                    spec.Valid_ANN=c(spec.Valid_ANN),
                    ROC.Valid_ANN=c( ROC.Valid_ANN),
                    CFR_ANNPred=c(CFR_ANNPred))
  }
  return(output)
}     

# Algorithmic framework for MAR 

Model.comMAR<-function(rep=1,rd,p,pw,dat){
  #rep is the number of iterations of the model
  #rd is the % missingness generate in the simulated data
  #p is the proportion simulated data used for training
  #pw is the % of the simulated data used for the analysis
  #options(java.parameters="-Xmx12g")
  #Samplesize of the simulated data
  u.sex<-sample(1:nrow(dat), round(pw*nrow(dat)), replace=FALSE)
  dat<-dat[u.sex,]
  CFR_Obs<-1-mean(dat$dead==1,na.rm=TRUE)#"Ground truth"
  levelofrandomness=rd
  data.y<-sort(dat$delay,decreasing = TRUE)
  nmar<-data.y[ceiling(levelofrandomness*length(dat$delay))]
  dat$dead[dat$delay>nmar]<-NA
  sizeofdata<- nrow(dat)
  for(j in  c(1:ncol(dat[,-21]))){
    # intorduce missingness
    # do the following for the outcome and the covariates
    fz <- which(runif(n=sizeofdata,0,1) < rd) 
    dat[fz,j] <- NA   #covariate i or outcome, here data is a datframe with outcome and covariates
    # redo thge analysis below with this new dataset
    # + include analysis with imputed data when possible -> extract CFR for the full dataset (known outcome + imputed outcome)
  }
  #Logistic Regression
  imp<-mlr::impute(dat, target = "dead", cols = list(country= imputeMode(), casedef=imputeMode(),delay=imputeMean(),hc=imputeMode(),mnths=imputeMode(),
                                                     unbld=imputeMode(),confs=imputeMode(),dbrth=imputeMode(),jntp=imputeMode(),
                                                     jnd=imputeMode(),cnjs=imputeMode(),fvr=imputeMode(),ftg=imputeMode(),
                                                     anx=imputeMode(),vmt=imputeMode(),dia=imputeMode(),hed=imputeMode(),
                                                     mus=imputeMode(),chp=imputeMode(),agecat=imputeHist()))
  data2<-imp$data
  output<-list()
  for(i in 1:rep){ 
    nonmiss.data<-data2[which(!(is.na(data2$dead))),]
    miss.data<-data2[which((is.na(data2$dead))),]
    #nonmiss.dataC<-nonmiss.data[complete.cases(nonmiss.data),]#because glm doesn't inherently handle missingness in the covariates 
    #data1<-simulated_data
    #Partitioning data into training and validation (testing) set. 
    u.training <-sample(1:nrow(nonmiss.data), round(p*nrow(nonmiss.data)), replace=FALSE)
    #u.training <-createDataPartition(nonmiss.data$dead, p=0.8, list=FALSE) 
    training.data<-nonmiss.data[u.training,]
    #testing.dat$mnths <- factor(testing.dat$mnths, levels = levels(training.dat$mnths))
    #I assuming that the CFR outomces are missing for the testing set. Is this MCAR or MAR?
    testing.data<-nonmiss.data[-u.training,]
    # I create the Modellling task for the alogrithms. By doing this, we dont have to do anymore data processing for each algorithm. 
    #listLearners("classif", check.packages = TRUE, properties = "missings")[c("class","package")]
    #summarizeColumns(training.data)
    #summarizeColumns(testing.data)
    #Recreate training and testing task
    trainTask <- makeRegrTask(data = training.data,target = "dead")
    testTask <- makeRegrTask(data = testing.data, target = "dead")
    removeConstantFeatures(trainTask)
    removeConstantFeatures(testTask)
    #Standardised Predictors
    trainTask <- normalizeFeatures(trainTask,method = "standardize")
    testTask<- normalizeFeatures(testTask,method = "standardize")
    #install.packages("fda.usc")
    #library("fda.usc")
    #Create the learner
    logistic.learner <- makeLearner("regr.glm",predict.type = "response",fix.factors.prediction = TRUE)
    #CV to optimise the the GLM
    #cv.logistic <- mlr::crossval(learner = logistic.learner,task = trainTask,iters = 2,stratify = TRUE)
    #cross validation accuracy
    #cv.logistic$aggr
    #cv.logistic$measures.test
    #Use the best model to for training 
    fmodel <- mlr::train(logistic.learner,trainTask)
    getLearnerModel(fmodel)
    #predict on test data using best model from the CV
    fpmodel <- predict(fmodel, testTask)
    preds<-fpmodel$data$response
    #Model Performance
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_GLM= accuracy$PCC[2]
    sens.Valid_GLM= accuracy$sensitivity[2]
    spec.Valid_GLM= accuracy$specificity[2]
    ROC.Valid_GLM= accuracy$AUC[2]
    #Predict on the data with missing outcomes
    fpmodel<- predict(fmodel, newdata=miss.data)
    preds.impute<-fpmodel$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_GLMPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    #CFR_Obs<-1-mean(dat$dead==1,na.rm=TRUE)#"Ground truth"
    #ANN
    #getParamSet("regr.nnet")
    g.net <- makeLearner("regr.nnet", predict.type = "response",fix.factors.prediction = TRUE)
    #specify tuning method
    rancontrol <- makeTuneControlRandom(maxit = 20L)
    #3 fold cross validation
    set_cv <- makeResampleDesc("CV",iters = 5L,stratify = FALSE)
    #parameters
    net_par<- makeParamSet(
      makeIntegerParam("size", lower = 3, upper = 10), #number of trees
      #makeIntegerParam("maxit", lower = 100, upper = 200), #depth of tree
      makeNumericParam("decay", lower = 1e-08, upper =0.1),
      makeNumericParam("abstol",lower = 0.0001, upper =0.0002),
      makeNumericParam("reltol",lower = 1e-08, upper = 1e-06)
    )
    #tune parameters
    tune_net <- tuneParams(learner = g.net, task = trainTask,resampling = set_cv,par.set = net_par,control = rancontrol)
    #check CV accuracy
    #tune_net$y
    #set parameters
    final_net <- setHyperPars(learner = g.net, par.vals = tune_net$x)
    #train
    to.net <- mlr::train(final_net, trainTask)
    #test 
    pr.net <- predict(to.net, testTask)
    preds<-pr.net$data$response
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_ANN=accuracy$PCC[2]
    sens.Valid_ANN = accuracy$sensitivity[2]
    spec.Valid_ANN = accuracy$specificity[2]
    ROC.Valid_ANN=accuracy$AUC[2]
    pr.net<- predict(to.net, newdata=miss.data)
    preds.impute<-pr.net$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_ANNPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    #########################################################################################################################
    nonmiss.data<-dat[which(!(is.na(dat$dead))),]
    miss.data<-dat[which((is.na(dat$dead))),]
    u.training <-sample(1:nrow(nonmiss.data), round(p*nrow(nonmiss.data)), replace=FALSE)
    #u.training <-createDataPartition(nonmiss.data$dead, p=0.8, list=FALSE) 
    training.data<-nonmiss.data[u.training,]
    #testing.dat$mnths <- factor(testing.dat$mnths, levels = levels(training.dat$mnths))
    #I assuming that the CFR outomces are missing for the testing set. Is this MCAR or MAR?
    testing.data<-nonmiss.data[-u.training,]
    # I create the Modellling task for the alogrithms. By doing this, we dont have to do anymore data processing for each algorithm. 
    #listLearners("classif", check.packages = TRUE, properties = "missings")[c("class","package")]
    #summarizeColumns(training.data)
    #summarizeColumns(testing.data)
    #Recreate training and testing task
    trainTask <- makeRegrTask(data = training.data,target = "dead")
    testTask <- makeRegrTask(data = testing.data, target = "dead")
    removeConstantFeatures(trainTask)
    removeConstantFeatures(testTask)
    #Standardised Predictors
    trainTask <- normalizeFeatures(trainTask,method = "standardize")
    testTask<- normalizeFeatures(testTask,method = "standardize")
    ######RandomForest
    #install.packages("randomForestSRC")
    #library("randomForestSRC")
    getParamSet("regr.randomForestSRC")
    #create a learner
    rf <- makeLearner("regr.randomForestSRC", predict.type = "response", par.vals = list(ntree = 200, mtry = 3),fix.factors.prediction = TRUE)
    #rf$par.vals <- list(importance = TRUE)
    #set tunable parameters
    # random search to find optimal hyperparameters ( manual grid search will be more computationaly expensive )
    rf_param <- makeParamSet(
      makeIntegerParam("ntree",lower = 50, upper = 300),
      makeIntegerParam("mtry", lower = 3, upper = 10),
      makeIntegerParam("nodesize", lower = 10, upper = 30)
    )
    rancontrol <- makeTuneControlRandom(maxit = 20L)
    #set 10 fold cross validation
    set_cv <- makeResampleDesc("CV",iters = 5L,stratify = FALSE)
    #hyperparamter tuning
    rf_tune <- tuneParams(learner = rf, resampling = set_cv, task = trainTask, par.set = rf_param, control = rancontrol)
    #using hyperparameters for modeling
    rf.tree <- setHyperPars(rf, par.vals = rf_tune$x)
    #train a model using the tuned hyperparameters
    rforest <- mlr::train(rf.tree, trainTask)
    getLearnerModel(rforest)
    #make predictions
    #glmmodel$xlevels[["district"]] <- union(glmmodel$xlevels[["district"]], levels(testing.data$district))
    #predict on the validation set
    rfmodel <- predict(rforest, testTask)
    preds<-rfmodel$data$response
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_RF= accuracy$PCC[2]
    sens.Valid_RF= accuracy$sensitivity[2]
    spec.Valid_RF= accuracy$specificity[2]
    ROC.Valid_RF= accuracy$AUC[2]
    #Predict on the data with missing outcomes
    #miss.data<-miss.data[complete.cases(miss.data),]#Random forest cannot handle missing predictors in new data
    #imp = impute(miss.data, classes = list(integer = imputeMean(), factor = imputeMode()),dummy.classes = "integer")
    rfpmodel<- predict(rforest, newdata=miss.data)
    preds.impute<-rfpmodel$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_RFPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    #CFR_GLMObs<-1-mean(data$dead==1,na.rm=TRUE)
    ####################BRT###########################
    #The steps are similar to those for random forest, the only change is that we are tuning a BRT and using that for predicting on the validation data.
    #load GBM
    getParamSet("regr.gbm")
    g.gbm <- makeLearner("regr.gbm", predict.type = "response",fix.factors.prediction = TRUE)
    
    #specify tuning method
    rancontrol <- makeTuneControlRandom(maxit = 20L)
    #3 fold cross validation
    set_cv <- makeResampleDesc("CV",iters = 5L,stratify = FALSE)
    #parameters
    gbm_par<- makeParamSet(
      makeDiscreteParam("distribution", values = "gaussian"),
      makeIntegerParam("n.trees", lower = 100, upper = 3000), #number of trees
      makeIntegerParam("interaction.depth", lower = 2, upper = 10), #depth of tree
      makeIntegerParam("n.minobsinnode", lower = 10, upper = 80),
      makeNumericParam("bag.fraction",lower=0.5, upper = 0.75),
      makeNumericParam("shrinkage",lower = 0.001, upper = 0.05)
    )
    #tune parameters
    tune_gbm <- tuneParams(learner = g.gbm, task = trainTask,resampling = set_cv,par.set = gbm_par,control = rancontrol)
    #check CV accuracy
    #tune_gbm$y
    #set parameters
    final_gbm <- setHyperPars(learner = g.gbm, par.vals = tune_gbm$x)
    #train
    to.gbm <- mlr::train(final_gbm, trainTask)
    #test 
    pr.gbm <- predict(to.gbm, testTask)
    preds<-pr.gbm$data$response
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_BRT=accuracy$PCC[2]
    sens.Valid_BRT = accuracy$sensitivity[2]
    spec.Valid_BRT = accuracy$specificity[2]
    ROC.Valid_BRT=accuracy$AUC[2]
    
    pr.gbm<- predict(to.gbm, newdata=miss.data)
    preds.impute<-pr.gbm$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_BRTPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    #CFR_GLMObs<-1-mean(data$dead==1,na.rm=TRUE)
    
    #load BART
    getParamSet("regr.bartMachine")
    g.bart <- makeLearner("regr.bartMachine", predict.type = "response",fix.factors.prediction = TRUE)
    
    #specify tuning method
    rancontrol <- makeTuneControlRandom(maxit = 20L)
    #3 fold cross validation
    set_cv <- makeResampleDesc("CV",iters = 5L,stratify = FALSE)
    #parameters
    bart_par<- makeParamSet(
      makeIntegerParam("num_trees", lower = 5, upper = 10), #number of trees
      makeIntegerParam("num_burn_in", lower = 20, upper = 30), #depth of tree
      makeIntegerParam("num_iterations_after_burn_in", lower = 70, upper = 100), #depth of tree
      #makeIntegerParam("alpha", lower = 0.95, upper = 3),
      makeNumericParam("beta",lower = 1, upper = 2),
      makeLogicalParam("mem_cache_for_speed",default = FALSE,tunable = FALSE),
      makeLogicalParam("run_in_sample",default = FALSE,tunable = FALSE),
      makeNumericParam("k",lower = 1, upper = 2)
    )
    #tune parameters
    #set_bart_machine_num_cores(5)
    tune_bart <- tuneParams(learner = g.bart, task = trainTask,resampling = set_cv,par.set = bart_par,control = rancontrol)
    #check CV accuracy
    tune_bart$y
    #set parameters
    final_bart <- setHyperPars(learner = g.bart, par.vals = tune_bart$x)
    #train
    to.bart <-  mlr::train(final_bart, trainTask)
    #test 
    pr.bart <- predict(to.bart, testTask)
    preds<-pr.bart$data$response
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_BART=accuracy$PCC[2]
    sens.Valid_BART = accuracy$sensitivity[2]
    spec.Valid_BART = accuracy$specificity[2]
    ROC.Valid_BART=accuracy$AUC[2]
    
    pr.bart<- predict(to.bart, newdata=miss.data)
    preds.impute<-pr.bart$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_BARTPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    
    output[[i]]<- c(p=p,pw=pw,rd=rd,CFR_Obs=c(CFR_Obs),
                    PCC.Valid_GLM=c(PCC.Valid_GLM),
                    sens.Valid_GLM=c(sens.Valid_GLM),
                    spec.Valid_GLM=c(spec.Valid_GLM),
                    ROC.Valid_GLM=c(ROC.Valid_GLM),
                    CFR_GLMPred=c(CFR_GLMPred),
                    PCC.Valid_RF=c(PCC.Valid_RF),
                    sens.Valid_RF=c( sens.Valid_RF),
                    spec.Valid_RF=c(spec.Valid_RF),
                    ROC.Valid_RF=c(ROC.Valid_RF),
                    CFR_RFPred=c(CFR_RFPred),
                    PCC.Valid_BRT=c(PCC.Valid_BRT),
                    sens.Valid_BRT=c(sens.Valid_BRT),
                    spec.Valid_BRT=c(spec.Valid_BRT),
                    ROC.Valid_BRT=c(ROC.Valid_BRT),
                    CFR_BRTPred=c(CFR_BRTPred),
                    PCC.Valid_BART=c(PCC.Valid_BART),
                    sens.Valid_BART=c(sens.Valid_BART),
                    spec.Valid_BART=c( spec.Valid_BART),
                    ROC.Valid_BART=c(ROC.Valid_BART),
                    CFR_BARTPred=c(CFR_BARTPred),
                    PCC.Valid_ANN=c(PCC.Valid_ANN),
                    sens.Valid_ANN=c(sens.Valid_ANN),
                    spec.Valid_ANN=c(spec.Valid_ANN),
                    ROC.Valid_ANN=c( ROC.Valid_ANN),
                    CFR_ANNPred=c(CFR_ANNPred))
  }
  return(output)
}      
#########################################################################################################
#Sensitivity analysis for adjusting for model imperfection for MCAR,MAR and MNAR
Model.comMAR4<-function(rep=1,rd,p,pw,dat){
  #rep is the number of iterations of the model
  #rd is the % missingness generate in the simulated data
  #p is the proportion simulated data used for training
  #pw is the % of the simulated data used for the analysis
  #options(java.parameters="-Xmx12g")
  #Samplesize of the simulated data
  u.sex<-sample(1:nrow(dat), round(pw*nrow(dat)), replace=FALSE)
  dat<-dat[u.sex,]
  CFR_Obs<-1-mean(dat$dead==1,na.rm=TRUE)#"Ground truth"
  levelofrandomness=rd
  data.y<-sort(dat$delay,decreasing = TRUE)
  nmar<-data.y[ceiling(levelofrandomness*length(dat$delay))]
  dat$dead[dat$delay>nmar]<-NA
  sizeofdata<- nrow(dat)
  for(j in  c(1:ncol(dat[,-21]))){
    # intorduce missingness
    # do the following for the outcome and the covariates
    fz <- which(runif(n=sizeofdata,0,1) < rd) 
    dat[fz,j] <- NA   #covariate i or outcome, here data is a datframe with outcome and covariates
    # redo thge analysis below with this new dataset
    # + include analysis with imputed data when possible -> extract CFR for the full dataset (known outcome + imputed outcome)
  }
  #Logistic Regression
  imp<-mlr::impute(dat, target = "dead", cols = list(country= imputeMode(), casedef=imputeMode(),delay=imputeMean(),hc=imputeMode(),mnths=imputeMode(),
                                                     unbld=imputeMode(),confs=imputeMode(),dbrth=imputeMode(),jntp=imputeMode(),
                                                     jnd=imputeMode(),cnjs=imputeMode(),fvr=imputeMode(),ftg=imputeMode(),
                                                     anx=imputeMode(),vmt=imputeMode(),dia=imputeMode(),hed=imputeMode(),
                                                     mus=imputeMode(),chp=imputeMode(),agecat=imputeHist()))
  data2<-imp$data
  output<-list()
  for(i in 1:rep){ 
    nonmiss.data<-data2[which(!(is.na(data2$dead))),]
    miss.data<-data2[which((is.na(data2$dead))),]
    #nonmiss.dataC<-nonmiss.data[complete.cases(nonmiss.data),]#because glm doesn't inherently handle missingness in the covariates 
    #data1<-simulated_data
    #Partitioning data into training and validation (testing) set. 
    u.training <-sample(1:nrow(nonmiss.data), round(p*nrow(nonmiss.data)), replace=FALSE)
    #u.training <-createDataPartition(nonmiss.data$dead, p=0.8, list=FALSE) 
    training.data<-nonmiss.data[u.training,]
    #testing.dat$mnths <- factor(testing.dat$mnths, levels = levels(training.dat$mnths))
    #I assuming that the CFR outomces are missing for the testing set. Is this MCAR or MAR?
    testing.data<-nonmiss.data[-u.training,]
    # I create the Modellling task for the alogrithms. By doing this, we dont have to do anymore data processing for each algorithm. 
    #listLearners("classif", check.packages = TRUE, properties = "missings")[c("class","package")]
    #summarizeColumns(training.data)
    #summarizeColumns(testing.data)
    #Recreate training and testing task
    trainTask <- makeRegrTask(data = training.data,target = "dead")
    testTask <- makeRegrTask(data = testing.data, target = "dead")
    removeConstantFeatures(trainTask)
    removeConstantFeatures(testTask)
    #Standardised Predictors
    trainTask <- normalizeFeatures(trainTask,method = "standardize")
    testTask<- normalizeFeatures(testTask,method = "standardize")
    #install.packages("fda.usc")
    #library("fda.usc")
    #Create the learner
    logistic.learner <- makeLearner("regr.glm",predict.type = "response",fix.factors.prediction = TRUE)
    #CV to optimise the the GLM
    #cv.logistic <- mlr::crossval(learner = logistic.learner,task = trainTask,iters = 2,stratify = TRUE)
    #cross validation accuracy
    #cv.logistic$aggr
    #cv.logistic$measures.test
    #Use the best model to for training 
    fmodel <- mlr::train(logistic.learner,trainTask)
    getLearnerModel(fmodel)
    #predict on test data using best model from the CV
    fpmodel <- predict(fmodel, testTask)
    preds<-fpmodel$data$response
    #Model Performance
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_GLM= accuracy$PCC[2]
    sens.Valid_GLM= accuracy$sensitivity[2]
    spec.Valid_GLM= accuracy$specificity[2]
    ROC.Valid_GLM= accuracy$AUC[2]
    #Predict on the data with missing outcomes
    fpmodel<- predict(fmodel, newdata=miss.data)
    preds.impute<-fpmodel$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_GLMPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    
    ds<- sum(new.data$dead==0,na.rm = TRUE)
    if (sum(!is.na(new.data$dead))>0) {
      ans<-epi.prev(pos=ds,tested=sum(!is.na(new.data$dead)),
                    se=accuracy$sensitivity[2], sp=accuracy$specificity[2],units = 1)
      fs<-ifelse(ans$tp[[1]]>=0,min(1,ans$tp[[1]])*sum(!is.na(new.data$dead)),
                 max(0,ans$tp[[1]])*sum(!is.na(new.data$dead)))
      tpCFR_GLM<-(fs+sum(nonmiss.data$dead==0,na.rm = TRUE))/sum(!is.na(final.data$dead))
    } else if (sum(!is.na(new.data$dead))==0)
    {tpCFR_GLM=CFR_GLMPred} 
    
    #CFR_Obs<-1-mean(dat$dead==1,na.rm=TRUE)#"Ground truth"
    #ANN
    #getParamSet("regr.nnet")
    g.net <- makeLearner("regr.nnet", predict.type = "response",fix.factors.prediction = TRUE)
    #specify tuning method
    rancontrol <- makeTuneControlRandom(maxit = 20L)
    #3 fold cross validation
    set_cv <- makeResampleDesc("CV",iters = 5L,stratify = FALSE)
    #parameters
    net_par<- makeParamSet(
      makeIntegerParam("size", lower = 3, upper = 10), #number of trees
      #makeIntegerParam("maxit", lower = 100, upper = 200), #depth of tree
      makeNumericParam("decay", lower = 1e-08, upper =0.1),
      makeNumericParam("abstol",lower = 0.0001, upper =0.0002),
      makeNumericParam("reltol",lower = 1e-08, upper = 1e-06)
    )
    #tune parameters
    tune_net <- tuneParams(learner = g.net, task = trainTask,resampling = set_cv,par.set = net_par,control = rancontrol)
    #check CV accuracy
    #tune_net$y
    #set parameters
    final_net <- setHyperPars(learner = g.net, par.vals = tune_net$x)
    #train
    to.net <- mlr::train(final_net, trainTask)
    #test 
    pr.net <- predict(to.net, testTask)
    preds<-pr.net$data$response
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_ANN=accuracy$PCC[2]
    sens.Valid_ANN = accuracy$sensitivity[2]
    spec.Valid_ANN = accuracy$specificity[2]
    ROC.Valid_ANN=accuracy$AUC[2]
    pr.net<- predict(to.net, newdata=miss.data)
    preds.impute<-pr.net$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_ANNPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    
    ds<-sum(new.data$dead==0,na.rm = TRUE)
    if (sum(!is.na(new.data$dead))>0) {
      ans<-epi.prev(pos=ds,tested=sum(!is.na(new.data$dead)),
                    se=accuracy$sensitivity[2], sp=accuracy$specificity[2],units = 1)
      fs<-ifelse(ans$tp[[1]]>=0,min(1,ans$tp[[1]])*sum(!is.na(new.data$dead)),
                 max(0,ans$tp[[1]])*sum(!is.na(new.data$dead)))
      tpCFR_ANN<-(fs+sum(nonmiss.data$dead==0,na.rm = TRUE))/sum(!is.na(final.data$dead))
    } else if (sum(!is.na(new.data$dead))==0)
    {tpCFR_ANN=CFR_ANNPred} 
    #########################################################################################################################
    nonmiss.data<-dat[which(!(is.na(dat$dead))),]
    miss.data<-dat[which((is.na(dat$dead))),]
    u.training <-sample(1:nrow(nonmiss.data), round(p*nrow(nonmiss.data)), replace=FALSE)
    #u.training <-createDataPartition(nonmiss.data$dead, p=0.8, list=FALSE) 
    training.data<-nonmiss.data[u.training,]
    #testing.dat$mnths <- factor(testing.dat$mnths, levels = levels(training.dat$mnths))
    #I assuming that the CFR outomces are missing for the testing set. Is this MCAR or MAR?
    testing.data<-nonmiss.data[-u.training,]
    # I create the Modellling task for the alogrithms. By doing this, we dont have to do anymore data processing for each algorithm. 
    #listLearners("classif", check.packages = TRUE, properties = "missings")[c("class","package")]
    #summarizeColumns(training.data)
    #summarizeColumns(testing.data)
    #Recreate training and testing task
    trainTask <- makeRegrTask(data = training.data,target = "dead")
    testTask <- makeRegrTask(data = testing.data, target = "dead")
    removeConstantFeatures(trainTask)
    removeConstantFeatures(testTask)
    #Standardised Predictors
    trainTask <- normalizeFeatures(trainTask,method = "standardize")
    testTask<- normalizeFeatures(testTask,method = "standardize")
    ######RandomForest
    #install.packages("randomForestSRC")
    #library("randomForestSRC")
    getParamSet("regr.randomForestSRC")
    #create a learner
    rf <- makeLearner("regr.randomForestSRC", predict.type = "response", par.vals = list(ntree = 200, mtry = 3),fix.factors.prediction = TRUE)
    #rf$par.vals <- list(importance = TRUE)
    #set tunable parameters
    # random search to find optimal hyperparameters ( manual grid search will be more computationaly expensive )
    rf_param <- makeParamSet(
      makeIntegerParam("ntree",lower = 50, upper = 300),
      makeIntegerParam("mtry", lower = 3, upper = 10),
      makeIntegerParam("nodesize", lower = 10, upper = 30)
    )
    rancontrol <- makeTuneControlRandom(maxit = 20L)
    #set 10 fold cross validation
    set_cv <- makeResampleDesc("CV",iters = 5L,stratify = FALSE)
    #hyperparamter tuning
    rf_tune <- tuneParams(learner = rf, resampling = set_cv, task = trainTask, par.set = rf_param, control = rancontrol)
    #using hyperparameters for modeling
    rf.tree <- setHyperPars(rf, par.vals = rf_tune$x)
    #train a model using the tuned hyperparameters
    rforest <- mlr::train(rf.tree, trainTask)
    getLearnerModel(rforest)
    #make predictions
    #glmmodel$xlevels[["district"]] <- union(glmmodel$xlevels[["district"]], levels(testing.data$district))
    #predict on the validation set
    rfmodel <- predict(rforest, testTask)
    preds<-rfmodel$data$response
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_RF= accuracy$PCC[2]
    sens.Valid_RF= accuracy$sensitivity[2]
    spec.Valid_RF= accuracy$specificity[2]
    ROC.Valid_RF= accuracy$AUC[2]
    #Predict on the data with missing outcomes
    #miss.data<-miss.data[complete.cases(miss.data),]#Random forest cannot handle missing predictors in new data
    #imp = impute(miss.data, classes = list(integer = imputeMean(), factor = imputeMode()),dummy.classes = "integer")
    rfpmodel<- predict(rforest, newdata=miss.data)
    preds.impute<-rfpmodel$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_RFPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    
    ds<-sum(new.data$dead==0,na.rm = TRUE)
    if (sum(!is.na(new.data$dead))>0) {
      ans<-epi.prev(pos=ds,tested=sum(!is.na(new.data$dead)),
                    se=accuracy$sensitivity[2], sp=accuracy$specificity[2],units = 1)
      fs<-ifelse(ans$tp[[1]]>=0,min(1,ans$tp[[1]])*sum(!is.na(new.data$dead)),
                 max(0,ans$tp[[1]])*sum(!is.na(new.data$dead)))
      tpCFR_RF<-(fs+sum(nonmiss.data$dead==0,na.rm = TRUE))/sum(!is.na(final.data$dead))
    } else if (sum(!is.na(new.data$dead))==0)
    {tpCFR_RF=CFR_RFPred} 
    #CFR_GLMObs<-1-mean(data$dead==1,na.rm=TRUE)
    ####################BRT###########################
    #The steps are similar to those for random forest, the only change is that we are tuning a BRT and using that for predicting on the validation data.
    #load GBM
    getParamSet("regr.gbm")
    g.gbm <- makeLearner("regr.gbm", predict.type = "response",fix.factors.prediction = TRUE)
    
    #specify tuning method
    rancontrol <- makeTuneControlRandom(maxit = 20L)
    #3 fold cross validation
    set_cv <- makeResampleDesc("CV",iters = 5L,stratify = FALSE)
    #parameters
    gbm_par<- makeParamSet(
      makeDiscreteParam("distribution", values = "gaussian"),
      makeIntegerParam("n.trees", lower = 100, upper = 3000), #number of trees
      makeIntegerParam("interaction.depth", lower = 2, upper = 10), #depth of tree
      makeIntegerParam("n.minobsinnode", lower = 10, upper = 80),
      makeNumericParam("bag.fraction",lower=0.5, upper = 0.75),
      makeNumericParam("shrinkage",lower = 0.001, upper = 0.05)
    )
    #tune parameters
    tune_gbm <- tuneParams(learner = g.gbm, task = trainTask,resampling = set_cv,par.set = gbm_par,control = rancontrol)
    #check CV accuracy
    #tune_gbm$y
    #set parameters
    final_gbm <- setHyperPars(learner = g.gbm, par.vals = tune_gbm$x)
    #train
    to.gbm <- mlr::train(final_gbm, trainTask)
    #test 
    pr.gbm <- predict(to.gbm, testTask)
    preds<-pr.gbm$data$response
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_BRT=accuracy$PCC[2]
    sens.Valid_BRT = accuracy$sensitivity[2]
    spec.Valid_BRT = accuracy$specificity[2]
    ROC.Valid_BRT=accuracy$AUC[2]
    
    pr.gbm<- predict(to.gbm, newdata=miss.data)
    preds.impute<-pr.gbm$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_BRTPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    #CFR_GLMObs<-1-mean(data$dead==1,na.rm=TRUE)
    
    ds<-sum(new.data$dead==0,na.rm = TRUE)
    if (sum(!is.na(new.data$dead))>0) {
      ans<-epi.prev(pos=ds,tested=sum(!is.na(new.data$dead)),
                    se=accuracy$sensitivity[2], sp=accuracy$specificity[2],units = 1)
      fs<-ifelse(ans$tp[[1]]>=0,min(1,ans$tp[[1]])*sum(!is.na(new.data$dead)),
                 max(0,ans$tp[[1]])*sum(!is.na(new.data$dead)))
      tpCFR_BRT<-(fs+sum(nonmiss.data$dead==0,na.rm = TRUE))/sum(!is.na(final.data$dead))
    } else if (sum(!is.na(new.data$dead))==0)
    {tpCFR_BRT=CFR_BRTPred} 
    #load BART
    getParamSet("regr.bartMachine")
    g.bart <- makeLearner("regr.bartMachine", predict.type = "response",fix.factors.prediction = TRUE)
    
    #specify tuning method
    rancontrol <- makeTuneControlRandom(maxit = 20L)
    #3 fold cross validation
    set_cv <- makeResampleDesc("CV",iters = 5L,stratify = FALSE)
    #parameters
    bart_par<- makeParamSet(
      makeIntegerParam("num_trees", lower = 5, upper = 10), #number of trees
      makeIntegerParam("num_burn_in", lower = 20, upper = 30), #depth of tree
      makeIntegerParam("num_iterations_after_burn_in", lower = 70, upper = 100), #depth of tree
      #makeIntegerParam("alpha", lower = 0.95, upper = 3),
      makeNumericParam("beta",lower = 1, upper = 2),
      makeLogicalParam("mem_cache_for_speed",default = FALSE,tunable = FALSE),
      makeLogicalParam("run_in_sample",default = FALSE,tunable = FALSE),
      makeNumericParam("k",lower = 1, upper = 2)
    )
    #tune parameters
    #set_bart_machine_num_cores(5)
    tune_bart <- tuneParams(learner = g.bart, task = trainTask,resampling = set_cv,par.set = bart_par,control = rancontrol)
    #check CV accuracy
    tune_bart$y
    #set parameters
    final_bart <- setHyperPars(learner = g.bart, par.vals = tune_bart$x)
    #train
    to.bart <-  mlr::train(final_bart, trainTask)
    #test 
    pr.bart <- predict(to.bart, testTask)
    preds<-pr.bart$data$response
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_BART=accuracy$PCC[2]
    sens.Valid_BART = accuracy$sensitivity[2]
    spec.Valid_BART = accuracy$specificity[2]
    ROC.Valid_BART=accuracy$AUC[2]
    
    pr.bart<- predict(to.bart, newdata=miss.data)
    preds.impute<-pr.bart$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_BARTPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    
    ds<-sum(new.data$dead==0,na.rm = TRUE)
    if (sum(!is.na(new.data$dead))>0) {
      ans<-epi.prev(pos=ds,tested=sum(!is.na(new.data$dead)),
                    se=accuracy$sensitivity[2], sp=accuracy$specificity[2],units = 1)
      fs<-ifelse(ans$tp[[1]]>=0,min(1,ans$tp[[1]])*sum(!is.na(new.data$dead)),
                 max(0,ans$tp[[1]])*sum(!is.na(new.data$dead)))
      tpCFR_BART<-(fs+sum(nonmiss.data$dead==0,na.rm = TRUE))/sum(!is.na(final.data$dead))
    } else if (sum(!is.na(new.data$dead))==0)
    {tpCFR_BART=CFR_BARTPred} 
    
    output[[i]]<- c(p=p,pw=pw,rd=rd,CFR_Obs=c(CFR_Obs),
                    PCC.Valid_GLM=c(PCC.Valid_GLM),
                    sens.Valid_GLM=c(sens.Valid_GLM),
                    spec.Valid_GLM=c(spec.Valid_GLM),
                    ROC.Valid_GLM=c(ROC.Valid_GLM),
                    CFR_GLMPred=c(CFR_GLMPred),
                    PCC.Valid_RF=c(PCC.Valid_RF),
                    sens.Valid_RF=c( sens.Valid_RF),
                    spec.Valid_RF=c(spec.Valid_RF),
                    ROC.Valid_RF=c(ROC.Valid_RF),
                    CFR_RFPred=c(CFR_RFPred),
                    PCC.Valid_BRT=c(PCC.Valid_BRT),
                    sens.Valid_BRT=c(sens.Valid_BRT),
                    spec.Valid_BRT=c(spec.Valid_BRT),
                    ROC.Valid_BRT=c(ROC.Valid_BRT),
                    CFR_BRTPred=c(CFR_BRTPred),
                    PCC.Valid_BART=c(PCC.Valid_BART),
                    sens.Valid_BART=c(sens.Valid_BART),
                    spec.Valid_BART=c( spec.Valid_BART),
                    ROC.Valid_BART=c(ROC.Valid_BART),
                    CFR_BARTPred=c(CFR_BARTPred),
                    PCC.Valid_ANN=c(PCC.Valid_ANN),
                    sens.Valid_ANN=c(sens.Valid_ANN),
                    spec.Valid_ANN=c(spec.Valid_ANN),
                    ROC.Valid_ANN=c( ROC.Valid_ANN),
                    CFR_ANNPred=c(CFR_ANNPred),
                    tpCFR_GLM=c(tpCFR_GLM),
                    tpCFR_ANN=c(tpCFR_ANN),
                    tpCFR_RF=c(tpCFR_RF),
                    tpCFR_BRT=c(tpCFR_BRT),
                    tpCFR_BART=c(tpCFR_BART))
  }
  return(output)
}      
###########################################################################################################################
Model.comMCAR4<-function(rep=1,rd,p,pw,dat){
  #rep is the number of iterations of the model
  #rd is the % missingness generate in the simulated data
  #p is the proportion simulated data used for training
  #pw is the % of the simulated data used for the analysis
  #options(java.parameters="-Xmx12g")
  #Samplesize of the simulated data
  u.sex<-sample(1:nrow(dat), round(pw*nrow(dat)), replace=FALSE)
  dat<-dat[u.sex,]
  CFR_Obs<-1-mean(dat$dead==1,na.rm=TRUE)#"Ground truth"
  #levelofrandomness=rd
  sizeofdata<- nrow(dat)
  #for (i in 1){ # at some point when working well, for (i in 1:100){
  for(j in  c(1:ncol(dat))){
    # intorduce missingness
    # do the following for the outcome and the covariates
    fz <- which(runif(n=sizeofdata,0,1) < rd) 
    dat[fz,j] <- NA   #covariate i or outcome, here data is a datframe with outcome and covariates
    # redo thge analysis below with this new dataset
    # + include analysis with imputed data when possible -> extract CFR for the full dataset (known outcome + imputed outcome)
  }
  #Logistic Regression
  imp<-mlr::impute(dat, target = "dead", cols = list(country= imputeMode(), casedef=imputeMode(),delay=imputeMean(),hc=imputeMode(),mnths=imputeMode(),
                                                     unbld=imputeMode(),confs=imputeMode(),dbrth=imputeMode(),jntp=imputeMode(),
                                                     jnd=imputeMode(),cnjs=imputeMode(),fvr=imputeMode(),ftg=imputeMode(),
                                                     anx=imputeMode(),vmt=imputeMode(),dia=imputeMode(),hed=imputeMode(),
                                                     mus=imputeMode(),chp=imputeMode(),agecat=imputeHist()))
  data2<-imp$data
  output<-list()
  for(i in 1:rep){ 
    nonmiss.data<-data2[which(!(is.na(data2$dead))),]
    miss.data<-data2[which((is.na(data2$dead))),]
    #nonmiss.dataC<-nonmiss.data[complete.cases(nonmiss.data),]#because glm doesn't inherently handle missingness in the covariates 
    #data1<-simulated_data
    #Partitioning data into training and validation (testing) set. 
    u.training <-sample(1:nrow(nonmiss.data), round(p*nrow(nonmiss.data)), replace=FALSE)
    #u.training <-createDataPartition(nonmiss.data$dead, p=0.8, list=FALSE) 
    training.data<-nonmiss.data[u.training,]
    #testing.dat$mnths <- factor(testing.dat$mnths, levels = levels(training.dat$mnths))
    #I assuming that the CFR outomces are missing for the testing set. Is this MCAR or MAR?
    testing.data<-nonmiss.data[-u.training,]
    # I create the Modellling task for the alogrithms. By doing this, we dont have to do anymore data processing for each algorithm. 
    #listLearners("classif", check.packages = TRUE, properties = "missings")[c("class","package")]
    #summarizeColumns(training.data)
    #summarizeColumns(testing.data)
    #Recreate training and testing task
    trainTask <- makeRegrTask(data = training.data,target = "dead")
    testTask <- makeRegrTask(data = testing.data, target = "dead")
    removeConstantFeatures(trainTask)
    removeConstantFeatures(testTask)
    #Standardised Predictors
    trainTask <- normalizeFeatures(trainTask,method = "standardize")
    testTask<- normalizeFeatures(testTask,method = "standardize")
    #install.packages("fda.usc")
    #library("fda.usc")
    #Create the learner
    logistic.learner <- makeLearner("regr.glm",predict.type = "response",fix.factors.prediction = TRUE)
    #CV to optimise the the GLM
    #cv.logistic <- mlr::crossval(learner = logistic.learner,task = trainTask,iters = 2,stratify = TRUE)
    #cross validation accuracy
    #cv.logistic$aggr
    #cv.logistic$measures.test
    #Use the best model to for training 
    fmodel <- mlr::train(logistic.learner,trainTask)
    getLearnerModel(fmodel)
    #predict on test data using best model from the CV
    fpmodel <- predict(fmodel, testTask)
    preds<-fpmodel$data$response
    #Model Performance
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_GLM= accuracy$PCC[2]
    sens.Valid_GLM= accuracy$sensitivity[2]
    spec.Valid_GLM= accuracy$specificity[2]
    ROC.Valid_GLM= accuracy$AUC[2]
    #Predict on the data with missing outcomes
    fpmodel<- predict(fmodel, newdata=miss.data)
    preds.impute<-fpmodel$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_GLMPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    
    ds<- sum(new.data$dead==0,na.rm = TRUE)
    if (sum(!is.na(new.data$dead))>0) {
      ans<-epi.prev(pos=ds,tested=sum(!is.na(new.data$dead)),
                    se=accuracy$sensitivity[2], sp=accuracy$specificity[2],units = 1)
      fs<-ifelse(ans$tp[[1]]>=0,min(1,ans$tp[[1]])*sum(!is.na(new.data$dead)),
                 max(0,ans$tp[[1]])*sum(!is.na(new.data$dead)))
      tpCFR_GLM<-(fs+sum(nonmiss.data$dead==0,na.rm = TRUE))/sum(!is.na(final.data$dead))
    } else if (sum(!is.na(new.data$dead))==0)
    {tpCFR_GLM=CFR_GLMPred} 
    
    #CFR_Obs<-1-mean(dat$dead==1,na.rm=TRUE)#"Ground truth"
    #ANN
    #getParamSet("regr.nnet")
    g.net <- makeLearner("regr.nnet", predict.type = "response",fix.factors.prediction = TRUE)
    #specify tuning method
    rancontrol <- makeTuneControlRandom(maxit = 20L)
    #3 fold cross validation
    set_cv <- makeResampleDesc("CV",iters = 5L,stratify = FALSE)
    #parameters
    net_par<- makeParamSet(
      makeIntegerParam("size", lower = 3, upper = 10), #number of trees
      #makeIntegerParam("maxit", lower = 100, upper = 200), #depth of tree
      makeNumericParam("decay", lower = 1e-08, upper =0.1),
      makeNumericParam("abstol",lower = 0.0001, upper =0.0002),
      makeNumericParam("reltol",lower = 1e-08, upper = 1e-06)
    )
    #tune parameters
    tune_net <- tuneParams(learner = g.net, task = trainTask,resampling = set_cv,par.set = net_par,control = rancontrol)
    #check CV accuracy
    #tune_net$y
    #set parameters
    final_net <- setHyperPars(learner = g.net, par.vals = tune_net$x)
    #train
    to.net <- mlr::train(final_net, trainTask)
    #test 
    pr.net <- predict(to.net, testTask)
    preds<-pr.net$data$response
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_ANN=accuracy$PCC[2]
    sens.Valid_ANN = accuracy$sensitivity[2]
    spec.Valid_ANN = accuracy$specificity[2]
    ROC.Valid_ANN=accuracy$AUC[2]
    pr.net<- predict(to.net, newdata=miss.data)
    preds.impute<-pr.net$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_ANNPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    
    ds<-sum(new.data$dead==0,na.rm = TRUE)
    if (sum(!is.na(new.data$dead))>0) {
      ans<-epi.prev(pos=ds,tested=sum(!is.na(new.data$dead)),
                    se=accuracy$sensitivity[2], sp=accuracy$specificity[2],units = 1)
      fs<-ifelse(ans$tp[[1]]>=0,min(1,ans$tp[[1]])*sum(!is.na(new.data$dead)),
                 max(0,ans$tp[[1]])*sum(!is.na(new.data$dead)))
      tpCFR_ANN<-(fs+sum(nonmiss.data$dead==0,na.rm = TRUE))/sum(!is.na(final.data$dead))
    } else if (sum(!is.na(new.data$dead))==0)
    {tpCFR_ANN=CFR_ANNPred} 
    #########################################################################################################################
    nonmiss.data<-dat[which(!(is.na(dat$dead))),]
    miss.data<-dat[which((is.na(dat$dead))),]
    u.training <-sample(1:nrow(nonmiss.data), round(p*nrow(nonmiss.data)), replace=FALSE)
    #u.training <-createDataPartition(nonmiss.data$dead, p=0.8, list=FALSE) 
    training.data<-nonmiss.data[u.training,]
    #testing.dat$mnths <- factor(testing.dat$mnths, levels = levels(training.dat$mnths))
    #I assuming that the CFR outomces are missing for the testing set. Is this MCAR or MAR?
    testing.data<-nonmiss.data[-u.training,]
    # I create the Modellling task for the alogrithms. By doing this, we dont have to do anymore data processing for each algorithm. 
    #listLearners("classif", check.packages = TRUE, properties = "missings")[c("class","package")]
    #summarizeColumns(training.data)
    #summarizeColumns(testing.data)
    #Recreate training and testing task
    trainTask <- makeRegrTask(data = training.data,target = "dead")
    testTask <- makeRegrTask(data = testing.data, target = "dead")
    removeConstantFeatures(trainTask)
    removeConstantFeatures(testTask)
    #Standardised Predictors
    trainTask <- normalizeFeatures(trainTask,method = "standardize")
    testTask<- normalizeFeatures(testTask,method = "standardize")
    ######RandomForest
    #install.packages("randomForestSRC")
    #library("randomForestSRC")
    getParamSet("regr.randomForestSRC")
    #create a learner
    rf <- makeLearner("regr.randomForestSRC", predict.type = "response", par.vals = list(ntree = 200, mtry = 3),fix.factors.prediction = TRUE)
    #rf$par.vals <- list(importance = TRUE)
    #set tunable parameters
    # random search to find optimal hyperparameters ( manual grid search will be more computationaly expensive )
    rf_param <- makeParamSet(
      makeIntegerParam("ntree",lower = 50, upper = 300),
      makeIntegerParam("mtry", lower = 3, upper = 10),
      makeIntegerParam("nodesize", lower = 10, upper = 30)
    )
    rancontrol <- makeTuneControlRandom(maxit = 20L)
    #set 10 fold cross validation
    set_cv <- makeResampleDesc("CV",iters = 5L,stratify = FALSE)
    #hyperparamter tuning
    rf_tune <- tuneParams(learner = rf, resampling = set_cv, task = trainTask, par.set = rf_param, control = rancontrol)
    #using hyperparameters for modeling
    rf.tree <- setHyperPars(rf, par.vals = rf_tune$x)
    #train a model using the tuned hyperparameters
    rforest <- mlr::train(rf.tree, trainTask)
    getLearnerModel(rforest)
    #make predictions
    #glmmodel$xlevels[["district"]] <- union(glmmodel$xlevels[["district"]], levels(testing.data$district))
    #predict on the validation set
    rfmodel <- predict(rforest, testTask)
    preds<-rfmodel$data$response
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_RF= accuracy$PCC[2]
    sens.Valid_RF= accuracy$sensitivity[2]
    spec.Valid_RF= accuracy$specificity[2]
    ROC.Valid_RF= accuracy$AUC[2]
    #Predict on the data with missing outcomes
    #miss.data<-miss.data[complete.cases(miss.data),]#Random forest cannot handle missing predictors in new data
    #imp = impute(miss.data, classes = list(integer = imputeMean(), factor = imputeMode()),dummy.classes = "integer")
    rfpmodel<- predict(rforest, newdata=miss.data)
    preds.impute<-rfpmodel$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_RFPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    
    ds<-sum(new.data$dead==0,na.rm = TRUE)
    if (sum(!is.na(new.data$dead))>0) {
      ans<-epi.prev(pos=ds,tested=sum(!is.na(new.data$dead)),
                    se=accuracy$sensitivity[2], sp=accuracy$specificity[2],units = 1)
      fs<-ifelse(ans$tp[[1]]>=0,min(1,ans$tp[[1]])*sum(!is.na(new.data$dead)),
                 max(0,ans$tp[[1]])*sum(!is.na(new.data$dead)))
      tpCFR_RF<-(fs+sum(nonmiss.data$dead==0,na.rm = TRUE))/sum(!is.na(final.data$dead))
    } else if (sum(!is.na(new.data$dead))==0)
    {tpCFR_RF=CFR_RFPred} 
    #CFR_GLMObs<-1-mean(data$dead==1,na.rm=TRUE)
    ####################BRT###########################
    #The steps are similar to those for random forest, the only change is that we are tuning a BRT and using that for predicting on the validation data.
    #load GBM
    getParamSet("regr.gbm")
    g.gbm <- makeLearner("regr.gbm", predict.type = "response",fix.factors.prediction = TRUE)
    
    #specify tuning method
    rancontrol <- makeTuneControlRandom(maxit = 20L)
    #3 fold cross validation
    set_cv <- makeResampleDesc("CV",iters = 5L,stratify = FALSE)
    #parameters
    gbm_par<- makeParamSet(
      makeDiscreteParam("distribution", values = "gaussian"),
      makeIntegerParam("n.trees", lower = 100, upper = 3000), #number of trees
      makeIntegerParam("interaction.depth", lower = 2, upper = 10), #depth of tree
      makeIntegerParam("n.minobsinnode", lower = 10, upper = 80),
      makeNumericParam("bag.fraction",lower=0.5, upper = 0.75),
      makeNumericParam("shrinkage",lower = 0.001, upper = 0.05)
    )
    #tune parameters
    tune_gbm <- tuneParams(learner = g.gbm, task = trainTask,resampling = set_cv,par.set = gbm_par,control = rancontrol)
    #check CV accuracy
    #tune_gbm$y
    #set parameters
    final_gbm <- setHyperPars(learner = g.gbm, par.vals = tune_gbm$x)
    #train
    to.gbm <- mlr::train(final_gbm, trainTask)
    #test 
    pr.gbm <- predict(to.gbm, testTask)
    preds<-pr.gbm$data$response
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_BRT=accuracy$PCC[2]
    sens.Valid_BRT = accuracy$sensitivity[2]
    spec.Valid_BRT = accuracy$specificity[2]
    ROC.Valid_BRT=accuracy$AUC[2]
    
    pr.gbm<- predict(to.gbm, newdata=miss.data)
    preds.impute<-pr.gbm$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_BRTPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    #CFR_GLMObs<-1-mean(data$dead==1,na.rm=TRUE)
    
    ds<-sum(new.data$dead==0,na.rm = TRUE)
    if (sum(!is.na(new.data$dead))>0) {
      ans<-epi.prev(pos=ds,tested=sum(!is.na(new.data$dead)),
                    se=accuracy$sensitivity[2], sp=accuracy$specificity[2],units = 1)
      fs<-ifelse(ans$tp[[1]]>=0,min(1,ans$tp[[1]])*sum(!is.na(new.data$dead)),
                 max(0,ans$tp[[1]])*sum(!is.na(new.data$dead)))
      tpCFR_BRT<-(fs+sum(nonmiss.data$dead==0,na.rm = TRUE))/sum(!is.na(final.data$dead))
    } else if (sum(!is.na(new.data$dead))==0)
    {tpCFR_BRT=CFR_BRTPred} 
    #load BART
    getParamSet("regr.bartMachine")
    g.bart <- makeLearner("regr.bartMachine", predict.type = "response",fix.factors.prediction = TRUE)
    
    #specify tuning method
    rancontrol <- makeTuneControlRandom(maxit = 20L)
    #3 fold cross validation
    set_cv <- makeResampleDesc("CV",iters = 5L,stratify = FALSE)
    #parameters
    bart_par<- makeParamSet(
      makeIntegerParam("num_trees", lower = 5, upper = 10), #number of trees
      makeIntegerParam("num_burn_in", lower = 20, upper = 30), #depth of tree
      makeIntegerParam("num_iterations_after_burn_in", lower = 70, upper = 100), #depth of tree
      #makeIntegerParam("alpha", lower = 0.95, upper = 3),
      makeNumericParam("beta",lower = 1, upper = 2),
      makeLogicalParam("mem_cache_for_speed",default = FALSE,tunable = FALSE),
      makeLogicalParam("run_in_sample",default = FALSE,tunable = FALSE),
      makeNumericParam("k",lower = 1, upper = 2)
    )
    #tune parameters
    #set_bart_machine_num_cores(5)
    tune_bart <- tuneParams(learner = g.bart, task = trainTask,resampling = set_cv,par.set = bart_par,control = rancontrol)
    #check CV accuracy
    tune_bart$y
    #set parameters
    final_bart <- setHyperPars(learner = g.bart, par.vals = tune_bart$x)
    #train
    to.bart <-  mlr::train(final_bart, trainTask)
    #test 
    pr.bart <- predict(to.bart, testTask)
    preds<-pr.bart$data$response
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_BART=accuracy$PCC[2]
    sens.Valid_BART = accuracy$sensitivity[2]
    spec.Valid_BART = accuracy$specificity[2]
    ROC.Valid_BART=accuracy$AUC[2]
    
    pr.bart<- predict(to.bart, newdata=miss.data)
    preds.impute<-pr.bart$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_BARTPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    
    ds<-sum(new.data$dead==0,na.rm = TRUE)
    if (sum(!is.na(new.data$dead))>0) {
      ans<-epi.prev(pos=ds,tested=sum(!is.na(new.data$dead)),
                    se=accuracy$sensitivity[2], sp=accuracy$specificity[2],units = 1)
      fs<-ifelse(ans$tp[[1]]>=0,min(1,ans$tp[[1]])*sum(!is.na(new.data$dead)),
                 max(0,ans$tp[[1]])*sum(!is.na(new.data$dead)))
      tpCFR_BART<-(fs+sum(nonmiss.data$dead==0,na.rm = TRUE))/sum(!is.na(final.data$dead))
    } else if (sum(!is.na(new.data$dead))==0)
    {tpCFR_BART=CFR_BARTPred} 
    
    output[[i]]<- c(p=p,pw=pw,rd=rd,CFR_Obs=c(CFR_Obs),
                    PCC.Valid_GLM=c(PCC.Valid_GLM),
                    sens.Valid_GLM=c(sens.Valid_GLM),
                    spec.Valid_GLM=c(spec.Valid_GLM),
                    ROC.Valid_GLM=c(ROC.Valid_GLM),
                    CFR_GLMPred=c(CFR_GLMPred),
                    PCC.Valid_RF=c(PCC.Valid_RF),
                    sens.Valid_RF=c( sens.Valid_RF),
                    spec.Valid_RF=c(spec.Valid_RF),
                    ROC.Valid_RF=c(ROC.Valid_RF),
                    CFR_RFPred=c(CFR_RFPred),
                    PCC.Valid_BRT=c(PCC.Valid_BRT),
                    sens.Valid_BRT=c(sens.Valid_BRT),
                    spec.Valid_BRT=c(spec.Valid_BRT),
                    ROC.Valid_BRT=c(ROC.Valid_BRT),
                    CFR_BRTPred=c(CFR_BRTPred),
                    PCC.Valid_BART=c(PCC.Valid_BART),
                    sens.Valid_BART=c(sens.Valid_BART),
                    spec.Valid_BART=c( spec.Valid_BART),
                    ROC.Valid_BART=c(ROC.Valid_BART),
                    CFR_BARTPred=c(CFR_BARTPred),
                    PCC.Valid_ANN=c(PCC.Valid_ANN),
                    sens.Valid_ANN=c(sens.Valid_ANN),
                    spec.Valid_ANN=c(spec.Valid_ANN),
                    ROC.Valid_ANN=c( ROC.Valid_ANN),
                    CFR_ANNPred=c(CFR_ANNPred),
                    tpCFR_GLM=c(tpCFR_GLM),
                    tpCFR_ANN=c(tpCFR_ANN),
                    tpCFR_RF=c(tpCFR_RF),
                    tpCFR_BRT=c(tpCFR_BRT),
                    tpCFR_BART=c(tpCFR_BART))
  }
  return(output)
}  
#########################################################################################################################
Model.comMNAR4<-function(rep=1,rd,p,pw,dat){
  #rep is the number of iterations of the model
  #rd is the % missingness generate in the simulated data
  #p is the proportion simulated data used for training
  #pw is the % of the simulated data used for the analysis
  #options(java.parameters="-Xmx12g")
  #Samplesize of the simulated data
  u.sex<-sample(1:nrow(dat), round(pw*nrow(dat)), replace=FALSE)
  dat<-dat[u.sex,]
  CFR_Obs<-1-mean(dat$dead==1,na.rm=TRUE)#"Ground truth"
  levelofrandomness=rd
  sizeofdat<- nrow(dat)
  #for (i in 1){ # at some point when working well, for (i in 1:100){
  #simulateMNAR for Outcome
  dx<-glm(dead~.,family = "binomial", data = dat)
  design_X <-model.matrix(~country+delay+casedef+hc+mnths+unbld+confs+
                            dbrth+jntp+jnd+cnjs+fvr+ftg+anx+vmt+dia+hed+mus+chp+agecat,data=dat)
  coef<-coefficients(dx)
  logit_CFR<- design_X%*%coef
  dat$missing<- 1/( 1 + exp(-logit_CFR))
  #simulateMNAR for Outcome
  #data$missing<-runif(n=sizeofdata,0,1) 
  data.y<-sort(dat$missing,decreasing = TRUE)
  nmar<-data.y[ceiling(rd*length(dat$missing))]
  dat$dead[dat$missing>nmar]<-NA
  #data$dead<-ifelse(data$missing>nmar,NA,data$dead)
  drop <- c("missing")
  dat<-dat[,!(names(dat) %in% drop)]
  for(j in  c(1:ncol(dat[,-21]))){
    # intorduce missingness
    # do the following for the outcome and the covariates
    fz <- which(runif(n=sizeofdat,0,1) < rd) 
    dat[fz,j] <- NA   #covariate i or outcome, here data is a datframe with outcome and covariates
    # redo thge analysis below with this new dataset
    # + include analysis with imputed data when possible -> extract CFR for the full dataset (known outcome + imputed outcome)
  }
  #Logistic Regression
  imp<-mlr::impute(dat, target = "dead", cols = list(country= imputeMode(), casedef=imputeMode(),delay=imputeMean(),hc=imputeMode(),mnths=imputeMode(),
                                                     unbld=imputeMode(),confs=imputeMode(),dbrth=imputeMode(),jntp=imputeMode(),
                                                     jnd=imputeMode(),cnjs=imputeMode(),fvr=imputeMode(),ftg=imputeMode(),
                                                     anx=imputeMode(),vmt=imputeMode(),dia=imputeMode(),hed=imputeMode(),
                                                     mus=imputeMode(),chp=imputeMode(),agecat=imputeHist()))
  data2<-imp$data
  output<-list()
  for(i in 1:rep){ 
    nonmiss.data<-data2[which(!(is.na(data2$dead))),]
    miss.data<-data2[which((is.na(data2$dead))),]
    #nonmiss.dataC<-nonmiss.data[complete.cases(nonmiss.data),]#because glm doesn't inherently handle missingness in the covariates 
    #data1<-simulated_data
    #Partitioning data into training and validation (testing) set. 
    u.training <-sample(1:nrow(nonmiss.data), round(p*nrow(nonmiss.data)), replace=FALSE)
    #u.training <-createDataPartition(nonmiss.data$dead, p=0.8, list=FALSE) 
    training.data<-nonmiss.data[u.training,]
    #testing.dat$mnths <- factor(testing.dat$mnths, levels = levels(training.dat$mnths))
    #I assuming that the CFR outomces are missing for the testing set. Is this MCAR or MAR?
    testing.data<-nonmiss.data[-u.training,]
    # I create the Modellling task for the alogrithms. By doing this, we dont have to do anymore data processing for each algorithm. 
    #listLearners("classif", check.packages = TRUE, properties = "missings")[c("class","package")]
    #summarizeColumns(training.data)
    #summarizeColumns(testing.data)
    #Recreate training and testing task
    trainTask <- makeRegrTask(data = training.data,target = "dead")
    testTask <- makeRegrTask(data = testing.data, target = "dead")
    removeConstantFeatures(trainTask)
    removeConstantFeatures(testTask)
    #Standardised Predictors
    trainTask <- normalizeFeatures(trainTask,method = "standardize")
    testTask<- normalizeFeatures(testTask,method = "standardize")
    #install.packages("fda.usc")
    #library("fda.usc")
    #Create the learner
    logistic.learner <- makeLearner("regr.glm",predict.type = "response",fix.factors.prediction = TRUE)
    #CV to optimise the the GLM
    #cv.logistic <- mlr::crossval(learner = logistic.learner,task = trainTask,iters = 2,stratify = TRUE)
    #cross validation accuracy
    #cv.logistic$aggr
    #cv.logistic$measures.test
    #Use the best model to for training 
    fmodel <- mlr::train(logistic.learner,trainTask)
    getLearnerModel(fmodel)
    #predict on test data using best model from the CV
    fpmodel <- predict(fmodel, testTask)
    preds<-fpmodel$data$response
    #Model Performance
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_GLM= accuracy$PCC[2]
    sens.Valid_GLM= accuracy$sensitivity[2]
    spec.Valid_GLM= accuracy$specificity[2]
    ROC.Valid_GLM= accuracy$AUC[2]
    #Predict on the data with missing outcomes
    fpmodel<- predict(fmodel, newdata=miss.data)
    preds.impute<-fpmodel$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_GLMPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    
    ds<- sum(new.data$dead==0,na.rm = TRUE)
    if (sum(!is.na(new.data$dead))>0) {
      ans<-epi.prev(pos=ds,tested=sum(!is.na(new.data$dead)),
                    se=accuracy$sensitivity[2], sp=accuracy$specificity[2],units = 1)
      fs<-ifelse(ans$tp[[1]]>=0,min(1,ans$tp[[1]])*sum(!is.na(new.data$dead)),
                 max(0,ans$tp[[1]])*sum(!is.na(new.data$dead)))
      tpCFR_GLM<-(fs+sum(nonmiss.data$dead==0,na.rm = TRUE))/sum(!is.na(final.data$dead))
    } else if (sum(!is.na(new.data$dead))==0)
    {tpCFR_GLM=CFR_GLMPred} 
    
    #CFR_Obs<-1-mean(dat$dead==1,na.rm=TRUE)#"Ground truth"
    #ANN
    #getParamSet("regr.nnet")
    g.net <- makeLearner("regr.nnet", predict.type = "response",fix.factors.prediction = TRUE)
    #specify tuning method
    rancontrol <- makeTuneControlRandom(maxit = 20L)
    #3 fold cross validation
    set_cv <- makeResampleDesc("CV",iters = 5L,stratify = FALSE)
    #parameters
    net_par<- makeParamSet(
      makeIntegerParam("size", lower = 3, upper = 10), #number of trees
      #makeIntegerParam("maxit", lower = 100, upper = 200), #depth of tree
      makeNumericParam("decay", lower = 1e-08, upper =0.1),
      makeNumericParam("abstol",lower = 0.0001, upper =0.0002),
      makeNumericParam("reltol",lower = 1e-08, upper = 1e-06)
    )
    #tune parameters
    tune_net <- tuneParams(learner = g.net, task = trainTask,resampling = set_cv,par.set = net_par,control = rancontrol)
    #check CV accuracy
    #tune_net$y
    #set parameters
    final_net <- setHyperPars(learner = g.net, par.vals = tune_net$x)
    #train
    to.net <- mlr::train(final_net, trainTask)
    #test 
    pr.net <- predict(to.net, testTask)
    preds<-pr.net$data$response
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_ANN=accuracy$PCC[2]
    sens.Valid_ANN = accuracy$sensitivity[2]
    spec.Valid_ANN = accuracy$specificity[2]
    ROC.Valid_ANN=accuracy$AUC[2]
    pr.net<- predict(to.net, newdata=miss.data)
    preds.impute<-pr.net$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_ANNPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    
    ds<-sum(new.data$dead==0,na.rm = TRUE)
    if (sum(!is.na(new.data$dead))>0) {
      ans<-epi.prev(pos=ds,tested=sum(!is.na(new.data$dead)),
                    se=accuracy$sensitivity[2], sp=accuracy$specificity[2],units = 1)
      fs<-ifelse(ans$tp[[1]]>=0,min(1,ans$tp[[1]])*sum(!is.na(new.data$dead)),
                 max(0,ans$tp[[1]])*sum(!is.na(new.data$dead)))
      tpCFR_ANN<-(fs+sum(nonmiss.data$dead==0,na.rm = TRUE))/sum(!is.na(final.data$dead))
    } else if (sum(!is.na(new.data$dead))==0)
    {tpCFR_ANN=CFR_ANNPred} 
    #########################################################################################################################
    nonmiss.data<-dat[which(!(is.na(dat$dead))),]
    miss.data<-dat[which((is.na(dat$dead))),]
    u.training <-sample(1:nrow(nonmiss.data), round(p*nrow(nonmiss.data)), replace=FALSE)
    #u.training <-createDataPartition(nonmiss.data$dead, p=0.8, list=FALSE) 
    training.data<-nonmiss.data[u.training,]
    #testing.dat$mnths <- factor(testing.dat$mnths, levels = levels(training.dat$mnths))
    #I assuming that the CFR outomces are missing for the testing set. Is this MCAR or MAR?
    testing.data<-nonmiss.data[-u.training,]
    # I create the Modellling task for the alogrithms. By doing this, we dont have to do anymore data processing for each algorithm. 
    #listLearners("classif", check.packages = TRUE, properties = "missings")[c("class","package")]
    #summarizeColumns(training.data)
    #summarizeColumns(testing.data)
    #Recreate training and testing task
    trainTask <- makeRegrTask(data = training.data,target = "dead")
    testTask <- makeRegrTask(data = testing.data, target = "dead")
    removeConstantFeatures(trainTask)
    removeConstantFeatures(testTask)
    #Standardised Predictors
    trainTask <- normalizeFeatures(trainTask,method = "standardize")
    testTask<- normalizeFeatures(testTask,method = "standardize")
    ######RandomForest
    #install.packages("randomForestSRC")
    #library("randomForestSRC")
    getParamSet("regr.randomForestSRC")
    #create a learner
    rf <- makeLearner("regr.randomForestSRC", predict.type = "response", par.vals = list(ntree = 200, mtry = 3),fix.factors.prediction = TRUE)
    #rf$par.vals <- list(importance = TRUE)
    #set tunable parameters
    # random search to find optimal hyperparameters ( manual grid search will be more computationaly expensive )
    rf_param <- makeParamSet(
      makeIntegerParam("ntree",lower = 50, upper = 300),
      makeIntegerParam("mtry", lower = 3, upper = 10),
      makeIntegerParam("nodesize", lower = 10, upper = 30)
    )
    rancontrol <- makeTuneControlRandom(maxit = 20L)
    #set 10 fold cross validation
    set_cv <- makeResampleDesc("CV",iters = 5L,stratify = FALSE)
    #hyperparamter tuning
    rf_tune <- tuneParams(learner = rf, resampling = set_cv, task = trainTask, par.set = rf_param, control = rancontrol)
    #using hyperparameters for modeling
    rf.tree <- setHyperPars(rf, par.vals = rf_tune$x)
    #train a model using the tuned hyperparameters
    rforest <- mlr::train(rf.tree, trainTask)
    getLearnerModel(rforest)
    #make predictions
    #glmmodel$xlevels[["district"]] <- union(glmmodel$xlevels[["district"]], levels(testing.data$district))
    #predict on the validation set
    rfmodel <- predict(rforest, testTask)
    preds<-rfmodel$data$response
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_RF= accuracy$PCC[2]
    sens.Valid_RF= accuracy$sensitivity[2]
    spec.Valid_RF= accuracy$specificity[2]
    ROC.Valid_RF= accuracy$AUC[2]
    #Predict on the data with missing outcomes
    #miss.data<-miss.data[complete.cases(miss.data),]#Random forest cannot handle missing predictors in new data
    #imp = impute(miss.data, classes = list(integer = imputeMean(), factor = imputeMode()),dummy.classes = "integer")
    rfpmodel<- predict(rforest, newdata=miss.data)
    preds.impute<-rfpmodel$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_RFPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    
    ds<-sum(new.data$dead==0,na.rm = TRUE)
    if (sum(!is.na(new.data$dead))>0) {
      ans<-epi.prev(pos=ds,tested=sum(!is.na(new.data$dead)),
                    se=accuracy$sensitivity[2], sp=accuracy$specificity[2],units = 1)
      fs<-ifelse(ans$tp[[1]]>=0,min(1,ans$tp[[1]])*sum(!is.na(new.data$dead)),
                 max(0,ans$tp[[1]])*sum(!is.na(new.data$dead)))
      tpCFR_RF<-(fs+sum(nonmiss.data$dead==0,na.rm = TRUE))/sum(!is.na(final.data$dead))
    } else if (sum(!is.na(new.data$dead))==0)
    {tpCFR_RF=CFR_RFPred} 
    #CFR_GLMObs<-1-mean(data$dead==1,na.rm=TRUE)
    ####################BRT###########################
    #The steps are similar to those for random forest, the only change is that we are tuning a BRT and using that for predicting on the validation data.
    #load GBM
    getParamSet("regr.gbm")
    g.gbm <- makeLearner("regr.gbm", predict.type = "response",fix.factors.prediction = TRUE)
    
    #specify tuning method
    rancontrol <- makeTuneControlRandom(maxit = 20L)
    #3 fold cross validation
    set_cv <- makeResampleDesc("CV",iters = 5L,stratify = FALSE)
    #parameters
    gbm_par<- makeParamSet(
      makeDiscreteParam("distribution", values = "gaussian"),
      makeIntegerParam("n.trees", lower = 100, upper = 3000), #number of trees
      makeIntegerParam("interaction.depth", lower = 2, upper = 10), #depth of tree
      makeIntegerParam("n.minobsinnode", lower = 10, upper = 80),
      makeNumericParam("bag.fraction",lower=0.5, upper = 0.75),
      makeNumericParam("shrinkage",lower = 0.001, upper = 0.05)
    )
    #tune parameters
    tune_gbm <- tuneParams(learner = g.gbm, task = trainTask,resampling = set_cv,par.set = gbm_par,control = rancontrol)
    #check CV accuracy
    #tune_gbm$y
    #set parameters
    final_gbm <- setHyperPars(learner = g.gbm, par.vals = tune_gbm$x)
    #train
    to.gbm <- mlr::train(final_gbm, trainTask)
    #test 
    pr.gbm <- predict(to.gbm, testTask)
    preds<-pr.gbm$data$response
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_BRT=accuracy$PCC[2]
    sens.Valid_BRT = accuracy$sensitivity[2]
    spec.Valid_BRT = accuracy$specificity[2]
    ROC.Valid_BRT=accuracy$AUC[2]
    
    pr.gbm<- predict(to.gbm, newdata=miss.data)
    preds.impute<-pr.gbm$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_BRTPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    #CFR_GLMObs<-1-mean(data$dead==1,na.rm=TRUE)
    
    ds<-sum(new.data$dead==0,na.rm = TRUE)
    if (sum(!is.na(new.data$dead))>0) {
      ans<-epi.prev(pos=ds,tested=sum(!is.na(new.data$dead)),
                    se=accuracy$sensitivity[2], sp=accuracy$specificity[2],units = 1)
      fs<-ifelse(ans$tp[[1]]>=0,min(1,ans$tp[[1]])*sum(!is.na(new.data$dead)),
                 max(0,ans$tp[[1]])*sum(!is.na(new.data$dead)))
      tpCFR_BRT<-(fs+sum(nonmiss.data$dead==0,na.rm = TRUE))/sum(!is.na(final.data$dead))
    } else if (sum(!is.na(new.data$dead))==0)
    {tpCFR_BRT=CFR_BRTPred} 
    #load BART
    getParamSet("regr.bartMachine")
    g.bart <- makeLearner("regr.bartMachine", predict.type = "response",fix.factors.prediction = TRUE)
    
    #specify tuning method
    rancontrol <- makeTuneControlRandom(maxit = 20L)
    #3 fold cross validation
    set_cv <- makeResampleDesc("CV",iters = 5L,stratify = FALSE)
    #parameters
    bart_par<- makeParamSet(
      makeIntegerParam("num_trees", lower = 5, upper = 10), #number of trees
      makeIntegerParam("num_burn_in", lower = 20, upper = 30), #depth of tree
      makeIntegerParam("num_iterations_after_burn_in", lower = 70, upper = 100), #depth of tree
      #makeIntegerParam("alpha", lower = 0.95, upper = 3),
      makeNumericParam("beta",lower = 1, upper = 2),
      makeLogicalParam("mem_cache_for_speed",default = FALSE,tunable = FALSE),
      makeLogicalParam("run_in_sample",default = FALSE,tunable = FALSE),
      makeNumericParam("k",lower = 1, upper = 2)
    )
    #tune parameters
    #set_bart_machine_num_cores(5)
    tune_bart <- tuneParams(learner = g.bart, task = trainTask,resampling = set_cv,par.set = bart_par,control = rancontrol)
    #check CV accuracy
    tune_bart$y
    #set parameters
    final_bart <- setHyperPars(learner = g.bart, par.vals = tune_bart$x)
    #train
    to.bart <-  mlr::train(final_bart, trainTask)
    #test 
    pr.bart <- predict(to.bart, testTask)
    preds<-pr.bart$data$response
    outcome.col <-which(colnames(training.data) %in% "dead")
    data.auc <- cbind(seq(1, dim(testing.data)[1], 1), testing.data[,outcome.col], preds)
    opt.thresh <- optimal.thresholds(data.auc, opt.methods=c("Sens=Spec"))
    accuracy<- presence.absence.accuracy(data.auc, threshold=as.numeric(opt.thresh))
    PCC.Valid_BART=accuracy$PCC[2]
    sens.Valid_BART = accuracy$sensitivity[2]
    spec.Valid_BART = accuracy$specificity[2]
    ROC.Valid_BART=accuracy$AUC[2]
    
    pr.bart<- predict(to.bart, newdata=miss.data)
    preds.impute<-pr.bart$data$response
    preds.outcome<- preds.impute
    preds.outcome[preds.impute >= as.numeric(opt.thresh[2])] <- 1 
    preds.outcome[preds.impute< as.numeric(opt.thresh[2])] <- 0
    # return the results
    new.data<-cbind(miss.data,preds.outcome)
    drop <- c("dead")
    new.data<-new.data[,!(names(new.data) %in% drop)]
    names(new.data)[names(new.data) == "preds.outcome"] <- "dead"
    #nonmiss.dat1<-nonmiss.dat[sample (1:nrow(nonmiss.dat) ,nrow(nonmiss.dat) , replace=TRUE),]
    #CFR with the predictions 
    final.data<-rbind(nonmiss.data,new.data)
    CFR_BARTPred<-1-mean(final.data$dead==1,na.rm=TRUE)
    
    ds<-sum(new.data$dead==0,na.rm = TRUE)
    if (sum(!is.na(new.data$dead))>0) {
      ans<-epi.prev(pos=ds,tested=sum(!is.na(new.data$dead)),
                    se=accuracy$sensitivity[2], sp=accuracy$specificity[2],units = 1)
      fs<-ifelse(ans$tp[[1]]>=0,min(1,ans$tp[[1]])*sum(!is.na(new.data$dead)),
                 max(0,ans$tp[[1]])*sum(!is.na(new.data$dead)))
      tpCFR_BART<-(fs+sum(nonmiss.data$dead==0,na.rm = TRUE))/sum(!is.na(final.data$dead))
    } else if (sum(!is.na(new.data$dead))==0)
    {tpCFR_BART=CFR_BARTPred} 
    
    output[[i]]<- c(p=p,pw=pw,rd=rd,CFR_Obs=c(CFR_Obs),
                    PCC.Valid_GLM=c(PCC.Valid_GLM),
                    sens.Valid_GLM=c(sens.Valid_GLM),
                    spec.Valid_GLM=c(spec.Valid_GLM),
                    ROC.Valid_GLM=c(ROC.Valid_GLM),
                    CFR_GLMPred=c(CFR_GLMPred),
                    PCC.Valid_RF=c(PCC.Valid_RF),
                    sens.Valid_RF=c( sens.Valid_RF),
                    spec.Valid_RF=c(spec.Valid_RF),
                    ROC.Valid_RF=c(ROC.Valid_RF),
                    CFR_RFPred=c(CFR_RFPred),
                    PCC.Valid_BRT=c(PCC.Valid_BRT),
                    sens.Valid_BRT=c(sens.Valid_BRT),
                    spec.Valid_BRT=c(spec.Valid_BRT),
                    ROC.Valid_BRT=c(ROC.Valid_BRT),
                    CFR_BRTPred=c(CFR_BRTPred),
                    PCC.Valid_BART=c(PCC.Valid_BART),
                    sens.Valid_BART=c(sens.Valid_BART),
                    spec.Valid_BART=c( spec.Valid_BART),
                    ROC.Valid_BART=c(ROC.Valid_BART),
                    CFR_BARTPred=c(CFR_BARTPred),
                    PCC.Valid_ANN=c(PCC.Valid_ANN),
                    sens.Valid_ANN=c(sens.Valid_ANN),
                    spec.Valid_ANN=c(spec.Valid_ANN),
                    ROC.Valid_ANN=c( ROC.Valid_ANN),
                    CFR_ANNPred=c(CFR_ANNPred),
                    tpCFR_GLM=c(tpCFR_GLM),
                    tpCFR_ANN=c(tpCFR_ANN),
                    tpCFR_RF=c(tpCFR_RF),
                    tpCFR_BRT=c(tpCFR_BRT),
                    tpCFR_BART=c(tpCFR_BART))
  }
  return(output)
}      











































