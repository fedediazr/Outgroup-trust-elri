##################################
# Confianza Intergrupal
# Modelos Multinivel
##################################

# Load df
rm(list=ls())
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

df<-dplyr::select(df, folio, ola, urbano_rural, g2, g18, g32_1, a1, 
                  indigena_es, conf_inter, conf_ingroup, cont_inter, frec_inter, quant_inter, calidad_inter, 
                  neg_inter, presn, same_a1, pond)


save(df, file = "Input/Base_multinivel.rdata")

