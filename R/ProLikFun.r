#################################################
#                                               #
# R function to fit the logistic Box-Cox model  #
# using profile likelihood                      #
#                                               #
#################################################
#
# Jan 18th, 2014
# Revised April 9, 2018
# Revised April 1, 2021

##' myfit.prof
##'
##' this function is aiming at estimates the LBC model parameter
##'
##' @param Predictor predictor in the model
##' @param Outcome response in the model
##'
##' @return beta0, beta1, lambda, iteration
##' @return calculate the beta0, beta1, lambda, iteration in LBC model
##' @examples
##'xx <- runif(100,min=0,max=1)
##'yy <- rbind(matrix(0,50,1),matrix(1,50,1))
##'myprofit <- myfit.prof(xx, yy)
##'init <- myprofit[1:3]
##'#mymle <- maxLik(logLik=LogLikeFun, grad = ScoreFun, start=init, ixx=xx, iyy=yy )
##'#beta0 <- mymle$estimate[1]
##'#beta1 <- mymle$estimate[2]
##'#lambda <- mymle$estimate[3]
##'#yy.gen <- beta0 + beta1*(xx^lambda-1)/lambda
##'
myfit.prof <- function(Predictor,Outcome)
   # Predictor = predictor in the model
   # Outcome = response in the model
   # Predictor = xx; Outcome = yy
{
   # Check if data are appropriate

   N1 <- length(Predictor)
   N2 <- length(Outcome)
   if (N1 == 0) {
      stop(" Null Predictor ")
   }
   if (N2 == 0) {
      stop(" Null Outcome ")
   }
   if (N1 != N2) {
      stop(" unequal length of response and predictor......")
   }

   # Set initial value of lambda and the upper boundary

   UpperBound = 2
   # when lambda is bigger than 2, it is associated with an unrealistic, rapid
   # change in risk
   ilambda = 0

   NoOFIni <- 10000
   myseq <- seq(0, UpperBound, length = NoOFIni) # form the sequence
   mylik <-
      rep(NA, NoOFIni)  # form the sequence of NA with same length

   QQ <- quantile(Predictor, 0.90) # get the last 10% biggest data
   Predictor[which(Predictor > QQ)] <- QQ

   for (i in 1:NoOFIni)
   {
      lambda0 = myseq[i]  # lambda determines the type of transformation

      if (lambda0 == 0) {
         myV0 <- log(Predictor)
      }
      # different transformation depends on lambda
      if (lambda0 != 0) {
         myV0 <- (Predictor ^ lambda0 - 1) / lambda0
      }
      # box-cox transformation

      mylogit <-
         glm(Outcome ~ myV0, family = "binomial") # use the logistic regression
      mycoef <- mylogit$coef # get the coefficients
      mybeta0 <- mycoef[1]  # get the beta0
      mybeta1 <- mycoef[2]  # get the beta1

      myPstar <- mybeta0 + mybeta1 * myV0  #
      mylik[i] <- sum(Outcome * myPstar - log(1 + exp(myPstar)))
   }

   myID <- order(mylik)[NoOFIni] # get the order
   ilambda <- myseq[myID]   # get the lambda

   if (ilambda == 0) {
      myV0 <- log(Predictor)
   }   # the same process as above
   if (ilambda != 0) {
      myV0 <- (Predictor ^ lambda0 - 1) / lambda0
   }
   mylogit <- glm(Outcome ~ myV0, family = "binomial")
   mycoef <- mylogit$coef
   mybeta0 <- mycoef[1] # the model parameter estimates
   mybeta1 <- mycoef[2]
   myresult <- c(mybeta0, mybeta1, ilambda, UpperBound)
   names(myresult) <-
      c("beta0", "beta1", "lambda", "iteration") # get the result estimated

   return(myresult)
}
