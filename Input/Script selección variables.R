##################################
# Confianza Intergrupal
# Modelos Multinivel
##################################

# Load df
rm(list=ls())
set.seed(202578)
getwd()
library(haven)
library(misty)
df<-read_dta("Input/BBDD_ELRI_LONG.dta")

library(dplyr)

# Create variable to identify consistent identifiers
df <- df %>%
  group_by(folio) %>%
  mutate(same_a1 = all(a1 == first(a1))) %>%
  ungroup()

df<-subset(df, same_a1==TRUE)
df <- df %>% 
  mutate(indigena_es = case_when(
    a1 >= 12 ~ "No indígena",
    a1 < 12   ~ "Indígena"))
table(df$indigena_es)
df <- as.na(df, na = c(66, 88, 99, 8888, 9999))

# 1. VARIABLE DEPENDIENTE

#######################
# Trust in Outgroup
#######################
library(Hmisc)
library(dplyr)
df$conf_inter<-NA
df <- df %>%
  mutate(conf_inter = ifelse(indigena_es=="Indígena", c2, conf_inter))

df <- df %>%
  mutate(conf_inter = ifelse(indigena_es=="No indígena", c5, conf_inter))

label(df$conf_inter)<-"Trust in Outgroup"

######################
# Ingroup trust
######################
df$conf_ingroup<-NA
df <- df %>%
  mutate(conf_ingroup = ifelse(indigena_es=="No indígena", c2, conf_ingroup))

df <- df %>%
  mutate(conf_ingroup = ifelse(indigena_es=="Indígena", c5, conf_ingroup))

label(df$conf_ingroup)<-"Ingroup trust"


# 2. PREDICTORS

# 2.2 COMPONENTE DE CONTACTO

#######################
# Cantidad de conocidos
#######################

table(df$c7_1)
table(df$c12)

df$cont_inter<-NA

df <- df %>%
  mutate(cont_inter = ifelse(indigena_es=="Indígena", c7_1, cont_inter))

df <- df %>%
  mutate(cont_inter = ifelse(indigena_es=="No indígena", c12, cont_inter))

label(df$cont_inter)<-"Number of acquintances"

########################
# Frecuencia de contacto
########################
table(df$c7_2)
table(df$c13)
df$frec_inter<-NA

df <- df %>%
  mutate(frec_inter = ifelse(indigena_es=="Indígena", c7_2, frec_inter))

df <- df %>%
  mutate(frec_inter = ifelse(indigena_es=="No indígena", c13, frec_inter))

label(df$frec_inter)<-"Intergroup contact frequency"

########################
# Cantidad de contacto
#######################
library(polycor)
polychor(df$cont_inter, df$frec_inter) #0.66
df$quant_inter=(df$frec_inter + df$cont_inter)/2

#######################
# Calidad de contacto
#######################
table(df$c7_3)
table(df$c14)
df$calidad_inter<-NA

df <- df %>%
  mutate(calidad_inter = ifelse(indigena_es=="Indígena", c7_3, calidad_inter))

df <- df %>%
  mutate(calidad_inter = ifelse(indigena_es=="No indígena", c14, calidad_inter))

label(df$calidad_inter)<-"Intergroup contact quality"

#######################
# Contacto negativo
#######################
table(df$c8)
table(df$c15)
df$neg_inter<-NA

df <- df %>%
  mutate(neg_inter = ifelse(indigena_es=="Indígena", c8, neg_inter))

df <- df %>%
  mutate(neg_inter = ifelse(indigena_es=="No indígena", c15, neg_inter))

label(df$neg_inter)<-"Frequency of negative intergroup contact"


########################################
# Componente de normas
########################################

# Normas prescriptivas hacia exogrupo (c26_1 para muestra chilenos, c26_2 para muestra indígena)

df$presn<-NA

df <- df %>%
  mutate(presn = ifelse(indigena_es=="Indígena", c26_2, presn))

df <- df %>%
  mutate(presn = ifelse(indigena_es=="No indígena", c26_1, presn))

label(df$presn)<-"Normas prescriptivas hacia exogrupo"


# Fijar sexo, edad y educación en t1
df <- df %>%
  group_by(folio) %>%
  mutate(
    g2 = ifelse(!is.na(first(g2[ola == 1])), first(g2[ola == 1]), g2),
    g18 = ifelse(!is.na(first(g18[ola == 1])), first(g18[ola == 1]), g18),
    g32_1 = ifelse(!is.na(first(g32_1[ola == 3])), first(g32_1[ola == 3]), g32_1)
  ) %>%
  ungroup()

#####################################
# Inverse probability weights (IPW)
#####################################

# Create variable to identify people with at least 3 observations
df <- df %>%
  group_by(folio) %>%
  mutate(three_obs = all(n_distinct(ola) >= 3)) |> ungroup()

df_w1 <- df |> filter(ola == 1)
df_w1 <- df_w1 |> 
  mutate(
    retained = three_obs,
    retained_label = ifelse(retained, 1, 0)
  )

library(mice)
df_w1<-dplyr::select(df_w1, folio, retained_label, indigena_es, g2, g18, g32_1, conf_inter, conf_ingroup, cont_inter, frec_inter, calidad_inter, neg_inter, presn)
imp <- mice(df_w1[2:13], m = 20, method = "pmm")  # or appropriate methods per variable type
fit_list <- with(imp, glm(retained_label ~ indigena_es + g2 + g18 + g32_1 + conf_inter + 
                            conf_ingroup + cont_inter + frec_inter + calidad_inter + 
                            neg_inter + presn, 
                          family = binomial))
pooled <- pool(fit_list)
summary(pooled)
pooled_summary <- summary(pooled, conf.int = TRUE)

# Build a clean results table with OR and 95% CI
results_df <- pooled_summary %>%
  mutate(
    OR       = exp(estimate),
    CI_lower = exp(`2.5 %`),
    CI_upper = exp(`97.5 %`),
    p.value  = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
  ) %>%
  transmute(
    Variable = term,
    OR       = sprintf("%.2f", OR),
    `95% CI` = sprintf("%.2f–%.2f", CI_lower, CI_upper),
    `p-value` = p.value
  )
#Export table to excel
library(openxlsx)
write.xlsx(results_df, "Output/sensitivity analysis/results_table.xlsx", rowNames = FALSE)

#Calculate IPW
pred_list <- lapply(1:imp$m, function(i) {
  d <- complete(imp, i)
  predict(fit_list$analyses[[i]], newdata = d, type = "response")
})

pred_matrix <- do.call(cbind, pred_list)   # N x m matrix
p_retained  <- rowMeans(pred_matrix) 

p_marginal <- mean(df_w1$retained_label, na.rm = TRUE)

sipw <- ifelse(
  df_w1$retained_label == 1,
  p_marginal / p_retained,
  (1 - p_marginal) / (1 - p_retained)
)

df_w1$sipw <- sipw
df_sipw <- df_w1 |> select(folio, sipw)
df<-merge(df, df_sipw, by = "folio", all.x = TRUE)
df<-dplyr::select(df, folio, ola, urbano_rural, g2, g18, g32_1, a1, 
                  indigena_es, conf_inter, conf_ingroup, cont_inter, frec_inter, quant_inter, calidad_inter, 
                  neg_inter, presn, same_a1, pond, sipw, three_obs)

save(df, file = "Input/Base_multinivel.rdata")

