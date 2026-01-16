# Gsens_ST-MH
Example codes for exploring genetic confounding of the associations between screen time and depressive symptoms

Code for exploring the association between screen time at age 16 and mental health at age 26 is used as an example:

## Linear regression models 

```r
ST16_S<-c("ST16S", "TV16S","COMP16S","TEXT16S","TALK16S","TV16WS","COMP16WS","TEXT16WS","TALK16WS")
regdta_16<-data.frame()
for (j in 1:9){
  ST <- ST16_S[j]
  print(ST)
#choose the model needed
#formula <- as.formula(paste("MH ~", ST))
#formula <- as.formula(paste("MH ~", ST,"+SEXF+MARRF+PEDUF+PCLASSF+NEET16F"))
#formula <- as.formula(paste("MH ~", ST,"+SEXF+MARRF+PEDUF+PCLASSF+NEET16F+MDD"))
print(formula)
REGM <- lm(formula, data = imp_sem[[i]])
  cr_ci<-coefficientr(REGM)
    regdta_16<-rbind(regdta_16, cr_ci) 
  #AIC<-AIC(REGM)
  #rsq<-summary(REGM)$r.squared
  #fit<-c(AIC,rsq)
  #fits<-rbind(fits, fit) 
 }
print(regdta_16)
```

## Gsens model and parameters used in the model
Gsens model (https://github.com/JBPG/Gsens) 

```r
ST16_S<-c("ST16S", "TV16S","COMP16S","TEXT16S","TALK16S","TV16WS","COMP16WS","TEXT16WS","TALK16WS")

r_data_16_26<-data.frame(col1="ST", col2="rgx", col3="rgy", col4="rxy")
regdta_16_26<-regdta_16
t1<-cor.test(imp_sem[[i]]$MH, imp_sem[[i]]$MDD, method = "pearson")
mh_pgs_p1_r<-t1$estimate

for (j in 1:9)
{print(ST16_S[j])
r_data_16_26[j,1]<-ST16_S[j]
#rgy path
r_data_16_26[j,3]<-mh_pgs_p1_r
#rxy path
r_data_16_26[j,4]<-regdta_16_26[j,1]
#rgx path
for (i in 1:50) {
  imp_sem[[i]]$ST <- imp_sem[[i]][[ST16_S[j]]]}
t2<-cor.test(imp_sem[[i]]$ST, imp_sem[[i]]$MDD, method = "pearson")
r_data_16_26[j,2]<-t2$estimate}

r_data_16_26[ , 2:4] <- lapply(r_data_16_26[ , 2:4], function(x) as.numeric(as.character(x)))
sapply(r_data_16_26, class)

gsensdta_16_26<-data.frame()
for (j in 1:9){
   M1_p1<-gsensY_prop(rxy=r_data_16_26[j,4],
                   rgx =r_data_16_26[j,2],
                   rgy =r_data_16_26[j,3],
                   n=N,
                   h2=0.05 #SNP-based heritability estimates)
results<-BCI(M1_p1)
results_select<-gsensout(results)
gsensdta_16_26<-rbind(gsensdta_16_26, as.data.frame(t(results_select)))
}

```

## Functions used in the codes

```r
#function to draw coefficients (standardized mean difference)
coefficientr<-function(M){
  data.frame(
    Estimate=coef(M)[2],
    LowerCI=confint(M)[2,1],
    UpperCI=confint(M)[2,2],
    ALLCI=sprintf("%.2f (%.2f, %.2f)", round(coef(M)[2], 3),round(confint(M)[2,1],3),round(confint(M)[2,2],3))
  )}

#function to convert smd to correlation, and get 95% CI
#references link: https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/es-calc.html
correlationr<-function(smd,N,n1,n2){
  r=smd / sqrt(smd^2 + ((N^2-2*N) / (n1 * n2)))
  z = 0.5 * (log(1 + r) - log(1 - r))
  se = 1 / sqrt(N - 3)
  z_ci_lower = z - qnorm(0.975) * se
  z_ci_upper = z + qnorm(0.975) * se
  r_ci_lower <- z_to_r(z_ci_lower)
  r_ci_upper <- z_to_r(z_ci_upper)
  ci95=sprintf("%.2f (%.2f, %.2f)", round(r, 3),round(r_ci_lower,3),round(r_ci_upper,3))
  r_ci<-cbind(r,r_ci_lower,r_ci_upper,ci95 )
  return(r_ci)
}
#function to convert z to r
z_to_r <- function(z) {
  r = (exp(2 * z) - 1) / (exp(2 * z) + 1)
  return(r)
}

#For the functions gsensY_prop and BCI, please refer to https://github.com/jr-baldwin/ACEs_mental_health_RR
```
