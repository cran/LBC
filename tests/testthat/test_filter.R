# 1 November, 2021
# revised 7 November, 2021
# author：Haoxue Wang (haoxwang@student.ethz.ch)


library(LBC)
library(testthat)


test_that('myfit.prof, the key function.', {

xx <- runif(100,min=0,max=1)
yy <- rbind(matrix(0,50,1),matrix(1,50,1))
myprofit <- myfit.prof(xx, yy)
init <- myprofit[1:3]

expect_equal(length(init), 3)
})

