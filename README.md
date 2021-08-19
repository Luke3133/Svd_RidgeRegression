# Svd_RidgeRegression
Perform Ridge Regression using spectral value decomposition.

main.r contains 2 functions: SVDRdgReg and CVIndex.

CVIndex is used to generate an index set which can be used for cross validation. This function outputs numbers from 1 to nGp (Number of groups) which can be used to split the data in the SVDRdgReg into groups. For example, CVIndex(n=12, nGp = 4, Shuffle=False) outputs the vector [1,1,1,2,2,2,3,3,3,4,4,4] which reperesents the first 3 rows of data going to group 1, the next 3 to group 2 etc... We can set Shuffle=True to output [1,4,2,1,1,3,2,4,3,2,4,3] which can be used to reduce dependence between data.

The second function SVDRdgReg is what actually does the ridge regression. The idea is to centre the X data, perform svd on the full data and use the diagonal of the sigma matrix to scale the test ridge parameters to an appropriate range. 
Performing SVD on the cross-validated sets allows us to find yhat for each lambda as well as calculating the mean square error. Using the mean square error, we find the optimum ridge parameter lambda and re-perform svd on the full data and calculate the optimum parameters Betahats using the optimmum lambda.

To understand further how this works, look into Principle Component Analysis.
