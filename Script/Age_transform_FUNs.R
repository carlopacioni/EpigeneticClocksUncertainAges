# Clock 1 - Log transform
trans_Clock1 <- function(Age) {
  tAge <- log(Age)
  tAge
}

revtrans_Clock1 <- function(tAge) {
  Age <- exp(tAge)
  Age
}

# Clock 2 - Loglog transformation of relative age, using maxlifespan and gestation time
trans_Clock2 <- function(Age, GestationTimeInYears, HighmaxAge) {
  RelativeAge <- (Age + GestationTimeInYears)/(HighmaxAge + GestationTimeInYears)
  tAge <- -log(-log(RelativeAge))
  tAge
}

revtrans_Clock2 <- function(tAge, GestationTimeInYears, HighmaxAge) {
  Age <- ((exp(-exp(-1*tAge)))*(HighmaxAge + GestationTimeInYears)) - GestationTimeInYears
  Age
} 


# Clock 3 
revtrans_Clock3 <- function(tAge, averagedMaturity.yrs, GestationTimeInYears) {
  Age <- tAge * (averagedMaturity.yrs + GestationTimeInYears) - GestationTimeInYears
  Age
} 

F1_logli <- function(age1, m1, m2 = m1, c1=1){
  ifelse(age1 >= m1, (age1-m1)/m2 , c1*log((age1-m1)/m2/c1 +1) )
}
#RelativeAdultAge
F2_revtrsf_clock3 <- function(y.pred, m1, m2 = m1, c1=1){
  ifelse(y.pred<0, (exp(y.pred/c1)-1)*m2*c1 + m1, y.pred*m2+m1 )
}

############# 
# The `loglifn` function shows how to calculate m1 for the transformation
# It is the `a_Logli` in the function

F3_loglifn = function(dat1,b1=1,max_tage = 4,
                      c1=5, c2 = 0.38, c0=0){
  n=nrow(dat1)
  
  age1 = (dat1$maxAge+dat1$GestationTimeInYears)/(dat1$averagedMaturity.yrs+dat1$GestationTimeInYears)
  
  a1 = age1/(1+max_tage)
  dat1$a1_Logli = a1 #x/m1 in manuscript
  
  a2 = (dat1$GestationTimeInYears + c0)/(dat1$averagedMaturity.yrs) 
  dat1$a_Logli = a_Logli = c1*a2^c2
  #m=5*(G/ASM)^0.38 from regression analysis/formula(7)
  
  
  x = dat1$Age + dat1$GestationTimeInYears
  t2 = dat1$averagedMaturity.yrs*b1 + dat1$GestationTimeInYears
  x2 = x/t2 #### log(x/t2)
  y = F1_logli(x2, a_Logli, a_Logli)
  
  dat1$LogliAge <- y
  return(dat1)
}


