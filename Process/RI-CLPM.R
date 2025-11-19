# RI-CLPM
library(dplyr)
rm(list=ls())
getwd()
load("Input/Base_multinivel.Rdata")

#Keep individuals with at least 3 waves
df <- df %>%
  group_by(folio) %>%
  filter(n_distinct(ola) >= 3)


# Reshape data

library(panelr)
df_long <- panel_data(df, id = folio, wave = ola)
df_wide <- widen_panel(df_long, separator = "_")

# promediar cantidad y frecuencia
df_wide$quant_1<-(df_wide$cont_inter_1+df_wide$frec_inter_1)/2
df_wide$quant_2<-(df_wide$cont_inter_2+df_wide$frec_inter_2)/2
df_wide$quant_3<-(df_wide$cont_inter_3+df_wide$frec_inter_3)/2
df_wide$quant_4<-(df_wide$cont_inter_4+df_wide$frec_inter_4)/2
#seleccionar variables
df_wide<-select(df_wide, folio, conf_inter_1, conf_inter_2, conf_inter_3, conf_inter_4,
                calidad_inter_1, calidad_inter_2, calidad_inter_3, calidad_inter_4,
                neg_inter_1, neg_inter_2, neg_inter_3, neg_inter_4, 
                presn_1, presn_2, presn_3, presn_4, 
                quant_1, quant_2, quant_3, quant_4, indigena_es)

# RI-CLPM

##########################################
# Complete model (no constraints)
##########################################

library(lavaan)

M6 <- '
             # Definition of random intercepts (between-person component)
             RI_trust      =~  1*conf_inter_1  + 1*conf_inter_2  +  1*conf_inter_3 + 1*conf_inter_4
             RI_quant     =~  1*quant_1     + 1*quant_2     +  1*quant_3    + 1*quant_4
             RI_qual     =~  1*calidad_inter_1     + 1*calidad_inter_2     +  1*calidad_inter_3    + 1*calidad_inter_4
             RI_neg     =~  1*neg_inter_1     + 1*neg_inter_2     +  1*neg_inter_3    + 1*neg_inter_4
             RI_presn     =~  1*presn_1     + 1*presn_2     +  1*presn_3    + 1*presn_4

             # Definition of within-person components
             w_trust_1       =~  1*conf_inter_1
             w_trust_2       =~  1*conf_inter_2
             w_trust_3       =~  1*conf_inter_3
             w_trust_4       =~  1*conf_inter_4
             
             w_quant_1      =~  1*quant_1
             w_quant_2      =~  1*quant_2
             w_quant_3      =~  1*quant_3
             w_quant_4      =~  1*quant_4
             
             w_qual_1      =~  1*calidad_inter_1
             w_qual_2      =~  1*calidad_inter_2
             w_qual_3      =~  1*calidad_inter_3
             w_qual_4      =~  1*calidad_inter_4             
             
             w_neg_1      =~  1*neg_inter_1
             w_neg_2      =~  1*neg_inter_2
             w_neg_3      =~  1*neg_inter_3
             w_neg_4      =~  1*neg_inter_4
             
             w_presn_1      =~  1*presn_1
             w_presn_2      =~  1*presn_2
             w_presn_3      =~  1*presn_3
             w_presn_4      =~  1*presn_4

                          
             # Set residual variances of observed variables to zero
             conf_inter_1   ~~  0*conf_inter_1
             conf_inter_2   ~~  0*conf_inter_2
             conf_inter_3   ~~  0*conf_inter_3
             conf_inter_4   ~~  0*conf_inter_4
             
             quant_1      ~~  0*quant_1
             quant_2      ~~  0*quant_2
             quant_3      ~~  0*quant_3
             quant_4      ~~  0*quant_4
             
             calidad_inter_1 ~~  0*calidad_inter_1
             calidad_inter_2 ~~  0*calidad_inter_2
             calidad_inter_3 ~~  0*calidad_inter_3
             calidad_inter_4 ~~  0*calidad_inter_4

             neg_inter_1 ~~  0*neg_inter_1
             neg_inter_2 ~~  0*neg_inter_2
             neg_inter_3 ~~  0*neg_inter_3
             neg_inter_4 ~~  0*neg_inter_4

             presn_1 ~~  0*presn_1
             presn_2 ~~  0*presn_2
             presn_3 ~~  0*presn_3
             presn_4 ~~  0*presn_4

            
             # Within-person autoregressive paths
             w_trust_2 ~ w_trust_1
             w_trust_3 ~ w_trust_2
             w_trust_4 ~ w_trust_3
             
             w_quant_2 ~ w_trust_1 + w_quant_1
             w_quant_3 ~ w_trust_2 + w_quant_2
             w_quant_4 ~ w_trust_3 + w_quant_3
             
             w_qual_2 ~ w_trust_1 + w_qual_1
             w_qual_3 ~ w_trust_2 + w_qual_2
             w_qual_4 ~ w_trust_3 + w_qual_3

             w_neg_2 ~ w_trust_1 + w_neg_1
             w_neg_3 ~ w_trust_2 + w_neg_2
             w_neg_4 ~ w_trust_3 + w_neg_3

             w_presn_2 ~ w_trust_1 + w_presn_1
             w_presn_3 ~ w_trust_2 + w_presn_2
             w_presn_4 ~ w_trust_3 + w_presn_3

             # Cross-lagged predictors of trust
             w_trust_2 ~ w_quant_1 + w_qual_1 + w_neg_1  + w_presn_1
             w_trust_3 ~ w_quant_2 + w_qual_2 + w_neg_2  + w_presn_2
             w_trust_4 ~ w_quant_3 + w_qual_3 + w_neg_3  + w_presn_3
             
             # Estimate the covariance of the random intercepts
             RI_trust  ~~ RI_quant
             RI_trust  ~~ RI_neg
             RI_trust  ~~ RI_qual
             RI_trust  ~~ RI_presn

             RI_quant  ~~ RI_neg
             RI_quant  ~~ RI_qual
             RI_quant  ~~ RI_presn

             RI_neg    ~~ RI_qual
             RI_neg    ~~ RI_presn

             RI_qual   ~~ RI_presn


             # Estimate the (residual) covariance between the within-person centered variables
             w_trust_1 ~~ w_quant_1
             w_trust_1 ~~ w_qual_1
             w_trust_1 ~~ w_neg_1
             w_trust_1 ~~ w_presn_1

             w_quant_1 ~~ w_qual_1
             w_quant_1 ~~ w_neg_1
             w_quant_1 ~~ w_presn_1

             w_qual_1  ~~ w_neg_1
             w_qual_1  ~~ w_presn_1

             w_neg_1   ~~ w_presn_1

             
             w_trust_2 ~~ w_quant_2
             w_trust_2 ~~ w_qual_2
             w_trust_2 ~~ w_neg_2
             w_trust_2 ~~ w_presn_2

             w_quant_2 ~~ w_qual_2
             w_quant_2 ~~ w_neg_2
             w_quant_2 ~~ w_presn_2

             w_qual_2  ~~ w_neg_2
             w_qual_2  ~~ w_presn_2

             w_neg_2   ~~ w_presn_2

             
             w_trust_3 ~~ w_quant_3
             w_trust_3 ~~ w_qual_3
             w_trust_3 ~~ w_neg_3
             w_trust_3 ~~ w_presn_3

             w_quant_3 ~~ w_qual_3
             w_quant_3 ~~ w_neg_3
             w_quant_3 ~~ w_presn_3

             w_qual_3  ~~ w_neg_3
             w_qual_3  ~~ w_presn_3

             w_neg_3   ~~ w_presn_3


             w_trust_4 ~~ w_quant_4
             w_trust_4 ~~ w_qual_4
             w_trust_4 ~~ w_neg_4
             w_trust_4 ~~ w_presn_4

             w_quant_4 ~~ w_qual_4
             w_quant_4 ~~ w_neg_4
             w_quant_4 ~~ w_presn_4

             w_qual_4  ~~ w_neg_4
             w_qual_4  ~~ w_presn_4

             w_neg_4   ~~ w_presn_4

'

fit_M6 <- sem(M6, data = df_wide, missing = "fiml")

summary(fit_M6, fit.measures=TRUE, standardized=TRUE)


##########################################
# Complete model (time constrained)
##########################################

M7 <- '
             # Definition of random intercepts (between-person component)
             RI_trust      =~  1*conf_inter_1  + 1*conf_inter_2  +  1*conf_inter_3 + 1*conf_inter_4
             RI_quant     =~  1*quant_1     + 1*quant_2     +  1*quant_3    + 1*quant_4
             RI_qual     =~  1*calidad_inter_1     + 1*calidad_inter_2     +  1*calidad_inter_3    + 1*calidad_inter_4
             RI_neg     =~  1*neg_inter_1     + 1*neg_inter_2     +  1*neg_inter_3    + 1*neg_inter_4
             RI_presn     =~  1*presn_1     + 1*presn_2     +  1*presn_3    + 1*presn_4

             # Definition of within-person components
             w_trust_1       =~  1*conf_inter_1
             w_trust_2       =~  1*conf_inter_2
             w_trust_3       =~  1*conf_inter_3
             w_trust_4       =~  1*conf_inter_4
             
             w_quant_1      =~  1*quant_1
             w_quant_2      =~  1*quant_2
             w_quant_3      =~  1*quant_3
             w_quant_4      =~  1*quant_4
             
             w_qual_1      =~  1*calidad_inter_1
             w_qual_2      =~  1*calidad_inter_2
             w_qual_3      =~  1*calidad_inter_3
             w_qual_4      =~  1*calidad_inter_4             
             
             w_neg_1      =~  1*neg_inter_1
             w_neg_2      =~  1*neg_inter_2
             w_neg_3      =~  1*neg_inter_3
             w_neg_4      =~  1*neg_inter_4
             
             w_presn_1      =~  1*presn_1
             w_presn_2      =~  1*presn_2
             w_presn_3      =~  1*presn_3
             w_presn_4      =~  1*presn_4
                          
             # Set residual variances of observed variables to zero
             conf_inter_1   ~~  0*conf_inter_1
             conf_inter_2   ~~  0*conf_inter_2
             conf_inter_3   ~~  0*conf_inter_3
             conf_inter_4   ~~  0*conf_inter_4
             
             quant_1      ~~  0*quant_1
             quant_2      ~~  0*quant_2
             quant_3      ~~  0*quant_3
             quant_4      ~~  0*quant_4
             
             calidad_inter_1 ~~  0*calidad_inter_1
             calidad_inter_2 ~~  0*calidad_inter_2
             calidad_inter_3 ~~  0*calidad_inter_3
             calidad_inter_4 ~~  0*calidad_inter_4

             neg_inter_1 ~~  0*neg_inter_1
             neg_inter_2 ~~  0*neg_inter_2
             neg_inter_3 ~~  0*neg_inter_3
             neg_inter_4 ~~  0*neg_inter_4

             presn_1 ~~  0*presn_1
             presn_2 ~~  0*presn_2
             presn_3 ~~  0*presn_3
             presn_4 ~~  0*presn_4
            
             # Within-person autoregressive paths
             w_trust_2 ~ a*w_trust_1
             w_trust_3 ~ a*w_trust_2
             w_trust_4 ~ a*w_trust_3
             
             w_quant_2 ~ g*w_trust_1 + b*w_quant_1
             w_quant_3 ~ g*w_trust_2 + b*w_quant_2
             w_quant_4 ~ g*w_trust_3 + b*w_quant_3
             
             w_qual_2 ~ h*w_trust_1 + c*w_qual_1
             w_qual_3 ~ h*w_trust_2 + c*w_qual_2
             w_qual_4 ~ h*w_trust_3 + c*w_qual_3

             w_neg_2 ~ i*w_trust_1 + d*w_neg_1
             w_neg_3 ~ i*w_trust_2 + d*w_neg_2
             w_neg_4 ~ i*w_trust_3 + d*w_neg_3

             w_presn_2 ~ k*w_trust_1 + f*w_presn_1
             w_presn_3 ~ k*w_trust_2 + f*w_presn_2
             w_presn_4 ~ k*w_trust_3 + f*w_presn_3

             # Cross-lagged predictors of trust
             w_trust_2 ~ l*w_quant_1 + m*w_qual_1 + n*w_neg_1 + p*w_presn_1
             w_trust_3 ~ l*w_quant_2 + m*w_qual_2 + n*w_neg_2 + p*w_presn_2
             w_trust_4 ~ l*w_quant_3 + m*w_qual_3 + n*w_neg_3 + p*w_presn_3
             
             # Estimate the covariance of the random intercepts
             RI_trust  ~~ RI_quant
             RI_trust  ~~ RI_neg
             RI_trust  ~~ RI_qual
             RI_trust  ~~ RI_presn

             RI_quant  ~~ RI_neg
             RI_quant  ~~ RI_qual
             RI_quant  ~~ RI_presn

             RI_neg    ~~ RI_qual
             RI_neg    ~~ RI_presn

             RI_qual   ~~ RI_presn


             # Estimate the (residual) covariance between the within-person centered variables
             w_trust_1 ~~ w_quant_1
             w_trust_1 ~~ w_qual_1
             w_trust_1 ~~ w_neg_1
             w_trust_1 ~~ w_presn_1

             w_quant_1 ~~ w_qual_1
             w_quant_1 ~~ w_neg_1
             w_quant_1 ~~ w_presn_1

             w_qual_1  ~~ w_neg_1
             w_qual_1  ~~ w_presn_1

             w_neg_1   ~~ w_presn_1

             
             w_trust_2 ~~ w_quant_2
             w_trust_2 ~~ w_qual_2
             w_trust_2 ~~ w_neg_2
             w_trust_2 ~~ w_presn_2

             w_quant_2 ~~ w_qual_2
             w_quant_2 ~~ w_neg_2
             w_quant_2 ~~ w_presn_2

             w_qual_2  ~~ w_neg_2
             w_qual_2  ~~ w_presn_2

             w_neg_2   ~~ w_presn_2

             
             w_trust_3 ~~ w_quant_3
             w_trust_3 ~~ w_qual_3
             w_trust_3 ~~ w_neg_3
             w_trust_3 ~~ w_presn_3

             w_quant_3 ~~ w_qual_3
             w_quant_3 ~~ w_neg_3
             w_quant_3 ~~ w_presn_3

             w_qual_3  ~~ w_neg_3
             w_qual_3  ~~ w_presn_3

             w_neg_3   ~~ w_presn_3


             w_trust_4 ~~ w_quant_4
             w_trust_4 ~~ w_qual_4
             w_trust_4 ~~ w_neg_4
             w_trust_4 ~~ w_presn_4

             w_quant_4 ~~ w_qual_4
             w_quant_4 ~~ w_neg_4
             w_quant_4 ~~ w_presn_4

             w_qual_4  ~~ w_neg_4
             w_qual_4  ~~ w_presn_4

             w_neg_4   ~~ w_presn_4

'

fit_M7 <- sem(M7, data = df_wide, missing = "fiml")

summary(fit_M7, fit.measures=TRUE, standardized=TRUE)



##########################################
# Complete model (time constrained, multigroup)
##########################################

library(lavaan)
M8 <- '
             # Definition of random intercepts (between-person component)
             RI_trust      =~  1*conf_inter_1  + 1*conf_inter_2  +  1*conf_inter_3 + 1*conf_inter_4
             RI_quant     =~  1*quant_1     + 1*quant_2     +  1*quant_3    + 1*quant_4
             RI_qual     =~  1*calidad_inter_1     + 1*calidad_inter_2     +  1*calidad_inter_3    + 1*calidad_inter_4
             RI_neg     =~  1*neg_inter_1     + 1*neg_inter_2     +  1*neg_inter_3    + 1*neg_inter_4
             RI_presn     =~  1*presn_1     + 1*presn_2     +  1*presn_3    + 1*presn_4

             # Definition of within-person components
             w_trust_1       =~  1*conf_inter_1
             w_trust_2       =~  1*conf_inter_2
             w_trust_3       =~  1*conf_inter_3
             w_trust_4       =~  1*conf_inter_4
             
             w_quant_1      =~  1*quant_1
             w_quant_2      =~  1*quant_2
             w_quant_3      =~  1*quant_3
             w_quant_4      =~  1*quant_4
             
             w_qual_1      =~  1*calidad_inter_1
             w_qual_2      =~  1*calidad_inter_2
             w_qual_3      =~  1*calidad_inter_3
             w_qual_4      =~  1*calidad_inter_4             
             
             w_neg_1      =~  1*neg_inter_1
             w_neg_2      =~  1*neg_inter_2
             w_neg_3      =~  1*neg_inter_3
             w_neg_4      =~  1*neg_inter_4
             
             w_presn_1      =~  1*presn_1
             w_presn_2      =~  1*presn_2
             w_presn_3      =~  1*presn_3
             w_presn_4      =~  1*presn_4
             
                          
             # Set residual variances of observed variables to zero
             conf_inter_1   ~~  0*conf_inter_1
             conf_inter_2   ~~  0*conf_inter_2
             conf_inter_3   ~~  0*conf_inter_3
             conf_inter_4   ~~  0*conf_inter_4
             
             quant_1      ~~  0*quant_1
             quant_2      ~~  0*quant_2
             quant_3      ~~  0*quant_3
             quant_4      ~~  0*quant_4
             
             calidad_inter_1 ~~  0*calidad_inter_1
             calidad_inter_2 ~~  0*calidad_inter_2
             calidad_inter_3 ~~  0*calidad_inter_3
             calidad_inter_4 ~~  0*calidad_inter_4

             neg_inter_1 ~~  0*neg_inter_1
             neg_inter_2 ~~  0*neg_inter_2
             neg_inter_3 ~~  0*neg_inter_3
             neg_inter_4 ~~  0*neg_inter_4

             presn_1 ~~  0*presn_1
             presn_2 ~~  0*presn_2
             presn_3 ~~  0*presn_3
             presn_4 ~~  0*presn_4

            
             # Within-person autoregressive paths
             w_trust_2 ~ c(yy1, yy2)*w_trust_1
             w_trust_3 ~ c(yy1, yy2)*w_trust_2
             w_trust_4 ~ c(yy1, yy2)*w_trust_3
             
             w_quant_2 ~ c(yu1,yu2)*w_trust_1 + c(uu1,uu2)*w_quant_1
             w_quant_3 ~ c(yu1,yu2)*w_trust_2 + c(uu1,uu2)*w_quant_2
             w_quant_4 ~ c(yu1,yu2)*w_trust_3 + c(uu1,uu2)*w_quant_3
             
             w_qual_2 ~ c(yv1,yv2)*w_trust_1 + c(vv1,vv2)*w_qual_1
             w_qual_3 ~ c(yv1,yv2)*w_trust_2 + c(vv1,vv2)*w_qual_2
             w_qual_4 ~ c(yv1,yv2)*w_trust_3 + c(vv1,vv2)*w_qual_3

             w_neg_2 ~ c(yw1,yw2)*w_trust_1 + c(ww1,ww2)*w_neg_1
             w_neg_3 ~ c(yw1,yw2)*w_trust_2 + c(ww1,ww2)*w_neg_2
             w_neg_4 ~ c(yw1,yw2)*w_trust_3 + c(ww1,ww2)*w_neg_3


             w_presn_2 ~ c(yz1,yz2)*w_trust_1 + c(zz1,zz2)*w_presn_1
             w_presn_3 ~ c(yz1,yz2)*w_trust_2 + c(zz1,zz2)*w_presn_2
             w_presn_4 ~ c(yz1,yz2)*w_trust_3 + c(zz1,zz2)*w_presn_3

             # Cross-lagged predictors of trust
             w_trust_2 ~ c(uy1,uy2)*w_quant_1 + c(vy1,vy2)*w_qual_1 + c(wy1,wy2)*w_neg_1 + c(zy1,zy2)*w_presn_1
             w_trust_3 ~ c(uy1,uy2)*w_quant_2 + c(vy1,vy2)*w_qual_2 + c(wy1,wy2)*w_neg_2 + c(zy1,zy2)*w_presn_2
             w_trust_4 ~ c(uy1,uy2)*w_quant_3 + c(vy1,vy2)*w_qual_3 + c(wy1,wy2)*w_neg_3 + c(zy1,zy2)*w_presn_3
             
             # Estimate the covariance of the random intercepts
             RI_trust  ~~ RI_quant
             RI_trust  ~~ RI_neg
             RI_trust  ~~ RI_qual
             RI_trust  ~~ RI_presn

             RI_quant  ~~ RI_neg
             RI_quant  ~~ RI_qual
             RI_quant  ~~ RI_presn

             RI_neg    ~~ RI_qual
             RI_neg    ~~ RI_presn

             RI_qual   ~~ RI_presn


             # Estimate the (residual) covariance between the within-person centered variables
             w_trust_1 ~~ w_quant_1
             w_trust_1 ~~ w_qual_1
             w_trust_1 ~~ w_neg_1
             w_trust_1 ~~ w_presn_1

             w_quant_1 ~~ w_qual_1
             w_quant_1 ~~ w_neg_1
             w_quant_1 ~~ w_presn_1

             w_qual_1  ~~ w_neg_1
             w_qual_1  ~~ w_presn_1

             w_neg_1   ~~ w_presn_1


             w_trust_2 ~~ w_quant_2
             w_trust_2 ~~ w_qual_2
             w_trust_2 ~~ w_neg_2
             w_trust_2 ~~ w_presn_2

             w_quant_2 ~~ w_qual_2
             w_quant_2 ~~ w_neg_2
             w_quant_2 ~~ w_presn_2

             w_qual_2  ~~ w_neg_2
             w_qual_2  ~~ w_presn_2

             w_neg_2   ~~ w_presn_2


             w_trust_3 ~~ w_quant_3
             w_trust_3 ~~ w_qual_3
             w_trust_3 ~~ w_neg_3
             w_trust_3 ~~ w_presn_3

             w_quant_3 ~~ w_qual_3
             w_quant_3 ~~ w_neg_3
             w_quant_3 ~~ w_presn_3

             w_qual_3  ~~ w_neg_3
             w_qual_3  ~~ w_presn_3

             w_neg_3   ~~ w_presn_3


             w_trust_4 ~~ w_quant_4
             w_trust_4 ~~ w_qual_4
             w_trust_4 ~~ w_neg_4
             w_trust_4 ~~ w_presn_4

             w_quant_4 ~~ w_qual_4
             w_quant_4 ~~ w_neg_4
             w_quant_4 ~~ w_presn_4

             w_qual_4  ~~ w_neg_4
             w_qual_4  ~~ w_presn_4

             w_neg_4   ~~ w_presn_4

'

fit_M8 <- sem(M8, data = df_wide, missing = "fiml", group = "indigena_es")

summary(fit_M8, fit.measures=TRUE, standardized=TRUE)


##################################
# Fitness comparison
##################################

# Load required libraries
library(lavaan)
library(flextable)
library(officer)
library(dplyr)

# Extract fit measures for each model
fit_M6_measures <- fitMeasures(fit_M6, c("chisq", "df", "pvalue", "cfi", "rmsea"))
fit_M7_measures <- fitMeasures(fit_M7, c("chisq", "df", "pvalue", "cfi", "rmsea"))
fit_M8_measures <- fitMeasures(fit_M8, c("chisq", "df", "pvalue", "cfi", "rmsea"))

# Create a data frame with the fit indices
fit_table <- data.frame(
  Model = c("Unconstrained", "Time constrained", "Multigroup"),
  χ2 = c(fit_M6_measures["chisq"], fit_M7_measures["chisq"], fit_M8_measures["chisq"]),
  df = c(fit_M6_measures["df"], fit_M7_measures["df"], fit_M8_measures["df"]),
  p = c(fit_M6_measures["pvalue"], fit_M7_measures["pvalue"], fit_M8_measures["pvalue"]),
  CFI = c(fit_M6_measures["cfi"], fit_M7_measures["cfi"], fit_M8_measures["cfi"]),
  RMSEA = c(fit_M6_measures["rmsea"], fit_M7_measures["rmsea"], fit_M8_measures["rmsea"])
)

# Format the values for publication
fit_table_formatted <- fit_table %>%
  mutate(
    χ2 = round(χ2, 2),
    p = ifelse(p < 0.001, "< .001", sprintf("%.3f", p)),
    CFI = round(CFI, 3),
    RMSEA = round(RMSEA, 3)
  )

# Create publication-ready table
word_table <- fit_table_formatted %>%
  flextable() %>%
  set_caption("Table 6: Model Fit Indices for Random intercept cross-lagged panel models") %>%
  set_header_labels(
    Model = "Model",
    ChiSq = "χ²",
    df = "df",
    p = "p-value",
    CFI = "CFI",
    RMSEA = "RMSEA"
  ) %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left") %>%  # Left align model names
  bold(part = "header") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 11, part = "all") %>%
  border_outer(part = "all") %>%
  border_inner_h(part = "body") %>%
  padding(padding = 4, part = "all") %>%
  add_footer_lines("Note: CFI = Comparative Fit Index; RMSEA = Root Mean Square Error of Approximation") %>%
  fontsize(size = 10, part = "footer")

doc <- read_docx() %>%  # Create a new Word document
  body_add_flextable(word_table)  # Add the table

# Save the Word document
print(doc, target = "Output/model_fit_indices_table.docx")

# Also display the table in R for verification
word_table

# Report:

#Final model included the theoretically-necessary time constraints 
#and multigroup parameters, demonstrated acceptable fit to the data: 
#χ²(286) = 439.50, CFI = .945, RMSEA = .030. While alternative models 
#showed slightly better fit, this model was retained as it directly addressed 
#the research questions."


# Regression parameters table

# Load required libraries
library(lavaan)
library(flextable)
library(officer)
library(dplyr)
library(tidyr)

# Extract parameter estimates for Model 6
model8_estimates <- parameterEstimates(fit_M8, standardized = TRUE)

# Filter for regression coefficients and extract time information
time_constrained_estimates <- model8_estimates %>%
  filter(op == "~") %>%
  # Extract time information from variable names
  mutate(
    lhs_time = as.numeric(gsub(".*_(\\d+)$", "\\1", lhs)),
    rhs_time = as.numeric(gsub(".*_(\\d+)$", "\\1", rhs)),
    # Get base variable name (without time suffix)
    lhs_var = gsub("_(\\d+)$", "", lhs),
    rhs_var = gsub("_(\\d+)$", "", rhs),
    group = ifelse(group == 1, "No indígena", "Indígena")
  ) %>%
  # Filter for time-constrained paths (consecutive time points)
  filter(lhs_time == rhs_time + 1) %>%
  # Create a unique identifier for each type of effect
  mutate(
    effect_type = case_when(
      lhs_var == rhs_var ~ "Autoregressive",
      lhs_var != rhs_var ~ "Cross-lagged"
    ),
    # Replace variable names with meaningful labels
    lhs_var_clean = case_when(
      lhs_var == "w_neg" ~ "Negative contact",
      lhs_var == "w_trust" ~ "Outgroup trust",
      lhs_var == "w_qual" ~ "Contact quality",
      lhs_var == "w_presn" ~ "Injunctive norms",
      lhs_var == "w_quant" ~ "Contact quantity",
      TRUE ~ lhs_var
    ),
    rhs_var_clean = case_when(
      rhs_var == "w_neg" ~ "Negative contact",
      rhs_var == "w_trust" ~ "Outgroup trust",
      rhs_var == "w_qual" ~ "Contact quality",
      rhs_var == "w_presn" ~ "Injunctive norms",
      rhs_var == "w_quant" ~ "Contact quantity",
      TRUE ~ rhs_var
    ),
    effect_label = paste(rhs_var_clean, "→", lhs_var_clean)
  ) %>%
  # Keep only one estimate per effect type (since they're time-constrained)
  group_by(group, effect_label, effect_type) %>%
  slice(1) %>%  # Take the first estimate for each effect type
  ungroup() %>%
  # Format the estimates
  mutate(
    Significance = case_when(
      pvalue < 0.001 ~ "***",
      pvalue < 0.01 ~ "**",
      pvalue < 0.05 ~ "*",
      TRUE ~ ""
    ),
    Estimate = sprintf("%.3f%s", est, Significance),
    SE = sprintf("(%.3f)", se),
    Std_Estimate = sprintf("%.3f%s", std.all, Significance),
    p = ifelse(pvalue < 0.001, "< .001", 
               ifelse(pvalue < 0.01, "< .010",
                      ifelse(pvalue < 0.05, "< .050", sprintf("%.3f", pvalue))))
  )

# Create the main table
main_table <- time_constrained_estimates %>%
  select(group, effect_type, effect_label, Estimate, SE, Std_Estimate, p) %>%
  arrange(group, effect_type, effect_label)

main_table$group[main_table$group=="Indígena"]<-"Indigenous"
main_table$group[main_table$group=="No indígena"]<-"Non-indigenous"
# Create a publication-ready flextable
word_table <- main_table %>%
  flextable() %>%
  set_caption("Time-Constrained Regression Estimates for Random Intercept Cross-Lagged Panel Model (Multigroup model)") %>%
  set_header_labels(
    group = "Group",
    effect_type = "Effect Type",
    effect_label = "Path",
    Estimate = "Unstd. Estimate",
    SE = "(SE)",
    Std_Estimate = "Std. Estimate",
    p = "p-value"
  ) %>%
  merge_v(j = c("group", "effect_type")) %>%
  theme_booktabs() %>%
  align(align = "center", part = "all") %>%
  align(j = c("effect_type", "effect_label"), align = "left") %>%
  bold(part = "header") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  padding(padding = 3, part = "all") %>%
  autofit()

word_table
# Export to Word
doc <- read_docx() %>%
  body_add_flextable(word_table) %>%
  body_add_par("Note: *** p < .001, ** p < .01, * p < .05. Estimates represent time-constrained paths (e.g., T1 → T2, T2 → T3, T3 → T4).", 
               style = "Normal")

print(doc, target = "Output/model8_time_constrained__multi_estimates.docx")





