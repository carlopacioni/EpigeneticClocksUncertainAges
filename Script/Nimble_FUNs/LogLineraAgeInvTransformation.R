llin3.inv <- nimbleFunction(
  run = function(LLin3Age=double(0), maturity=double(0), trans=logical(0)) {
    if(trans == TRUE) {
      if(LLin3Age < 0) {
        age_yr <- (maturity + 1.5) * exp(LLin3Age) - 1.5
      } else {
        age_yr <- (maturity + 1.5) * LLin3Age + maturity
      }
    } else {
      age_yr <- LLin3Age
    }
    returnType(double(0))
  return(age_yr)
})

# Testing 
# fun_llin3.inv(c(0, -1, 1), 1)
# llin3.inv(LLin3Age = 0, maturity = 1, trans=TRUE)
# llin3.inv(LLin3Age = -1, maturity = 1, trans=TRUE)
# llin3.inv(LLin3Age = 1, maturity = 1, trans=TRUE)
# llin3.inv(LLin3Age = 1, maturity = 1, trans=FALSE)
# cllin3.inv <- compileNimble(llin3.inv)
# cllin3.inv(1, 1, TRUE)
