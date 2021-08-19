SVDRdgReg = function(X,Y, nCV, CVShuffle, VRb = FALSE, Print = TRUE){
  #' X - Covariate Matrix
  #' Y - Response vector
  #' nCV - Number of Cross-Validation Groups
  #' CVShuffle - Do we Want to shuffle the CV groups?
  #' VRb - Verbose mode
  #' Print - Print out graphs (set as false for bootstrapping)
  
  #Generate Cross-Validation Index set
  CVIndexSet = CVIndex(n = nrow(X), nGp = nCV, shuffle = CVShuffle)
  n = nrow(X) #Calculate the number of observations
  
  
  lambda_seq = 10^seq(20, -20, by = -.05) # Declare set of ridge parameters, lambda
  MSE = rep(0, length(lambda_seq)) #Declare MSE vector
  
  #Perform Svd on complete data
  
  #Center X data
  cX = X - matrix(colMeans(X), nrow=nrow(X), ncol=length(colMeans(X)), byrow = TRUE) 
  
  #Perform SVD
  decompo = svd(cX,nu=50, nv = 50)
  Lxo = decompo$d
  Uxo = decompo$u
  Vxo = decompo$v
  L = lambda_seq*sum(Lxo^2) #Scale lambdas
  PredictedValues = matrix(ncol = length(L), nrow=n) #Declare predicted values matrix(holds predicted values for each lambda)
  
  
  if (VRb == TRUE){print("Running Cross Validation")}
  #Perform cross-validation
  for (i in 1:nCV){
    #Generate data sets from cross-validation index set
    CVTrainX = X[CVIndexSet!=i,]
    CVTrainY = Y[CVIndexSet!=i]
    CVTestX = X[CVIndexSet==i,]
    CVTestY = Y[CVIndexSet==i]
    
    
    #Centre data wrt mean
    CentCVTrainX = CVTrainX - matrix(colMeans(CVTrainX), nrow=nrow(CVTrainX), ncol=length(colMeans(CVTrainX)), byrow = TRUE)
    CentCVTrainY = CVTrainY - mean(CVTrainY)
    CentCVTestX = CVTestX - matrix(colMeans(CVTrainX), nrow=nrow(CVTestX), ncol=length(colMeans(CVTrainX)), byrow = TRUE)
    CentCVTestY = CVTestY - mean(CVTrainY)
    
    
    
    #Singular value decomposition on cross-validated data
    decomp = svd(CentCVTrainX, nu=ncol(CentCVTrainX))
    tmp = decomp$d > 0.0000001*max(decomp$d)
    Ux = decomp$u[,tmp==TRUE]      
    Lx = decomp$d[tmp==TRUE]
    Vx = decomp$v[,tmp==TRUE]
    
    tmp1 = CentCVTestX%*%Vx #used to calculate Yhat
    tmp3 = t(Ux)%*%CentCVTrainY #used to calculate Yhat
    
    #Loop over lambdas
    for (j in 1:length(L)){
      
      bottom = Lx^2 + rep(1, length(Lx))*L[j] # used to calculate tmp2
      tmp2 =Lx/bottom #used to calculate Yhat
      
      Yhat =  tmp1%*%diag(tmp2)%*%tmp3 #Calculate predicted values
      MSE[j] = MSE[j] + sum((Yhat - CentCVTestY)^2) #Calculate Mean Square Error(MSE)
      PredictedValues[CVIndexSet==i,j] = Yhat + rep(1,length=length(as.vector(Yhat)))%*%t(mean(CVTrainY)) #Calculate the predicted values of the test set
    }
  }
  if (VRb == TRUE){print("Calculating minMSE")}
  
  MSE = MSE/nrow(X) #Calculate the final MSE
  minMSE = which.min(MSE) #Find the optimum lambda (Minimum MSE)
  
  if(minMSE == 1){ #Check if the ridge parameter is one of the edge cases (do we need to increase the range of the lambdas?)
    print("WARNING: Optimum lambda is the minimum shrinkage parameter!")
  } else if (minMSE == length(MSE)){
    print("WARNING: Optimum lambda is the maximum shrinkage parameter!")
  } else {
    if (VRb == TRUE){print("Good shrinkage parameter")}
  }
  Yhat = PredictedValues[,minMSE] #Predicted values of the optimum lambda cross-validation set
  if (VRb == TRUE){print("Plotting MSE")}
  if (Print == TRUE){
    plot(log(L,10), MSE, main = "A plot of log lambdas and its corresponding Mean Square Error")
  }
  
  
  if (VRb == TRUE){print("Calculate Final Parameters")}
  
  #Centre full data wrt mean
  
  CentX = X - matrix(colMeans(X), nrow=nrow(X), ncol=length(colMeans(X)), byrow = TRUE)
  CentY = Y - mean(Y)
  
  #Full singular value decomposition
  decomp = 0
  decomp = svd(CentX, nu=ncol(CentX))
  tmp = decomp$d > 0.0000001*max(decomp$d)
  Ux = decomp$u[,tmp==TRUE] 
  Lx = decomp$d[tmp==TRUE] 
  Vx = decomp$v[,tmp==TRUE] 
  
  tmp1 = Vx #used to calculate Yhat
  tmp3 = t(Ux)%*%CentY #used to calculate Yhat
  bottom = Lx^2 + rep(1, length(Lx))*L[minMSE] #used to calculate tmp2
  tmp2 =Lx/bottom #used to calculate Yhat
  
  if (VRb == TRUE){print("Returning Final Parameters")}
  
  Betahat = tmp1%*%diag(tmp2)%*%tmp3 #Calculate optimum betas
  
  return(return(list(L[minMSE], Betahat, colMeans(X), mean(Y),Yhat,MSE,L))) #return
}


#This function allows the user to create a cross validation index set that can be used to select rows of a dataset.
CVIndex = function(n = 0, nGp = 0, shuffle = FALSE){
  #Make an index set for traditional k-fold cross validation
  #'
  #' n - Number of observations in the original dataset
  #' nGp - number of groups in k means cross validation i.e. nGp = k
  #' shuffle - select TRUE if you want the output of the function to be randomised(with the same number of groups nGp)
  #'
  
  if (n == 0){ #Check n is provided
    stop("Please provide n")
  }
  if (nGp == 0){ #Check nGp is provided
    stop("Please provide a number of groups nGp")
  }
  CV = vector(length =  n) #Prepare output vector
  if (n == nGp){ #Group length is 1 so output a sequence
    CV = seq(1:n)
  } else {
    #Check that nGp is a divisor of n
    if (n%%nGp != 0){
      print(n) #Print group size 
      stop("The number of groups is not a divisor of n. Please pick a suitable group size")
    }
    Gpsize = floor(n/nGp) #Find the size of each group
    j=1
    a=1
    #create  groups
    for (i in 1:n){
      if (a<Gpsize){
        CV[i] = j
        a=a+1
      } else {
        CV[i] = j
        a=1
        j = j+1
      }
    }
  }
  
  if (shuffle == TRUE){ #Shuffle the output
    CV = sample(CV,replace = FALSE)
  }
  return(CV) #Return index set
  
}
