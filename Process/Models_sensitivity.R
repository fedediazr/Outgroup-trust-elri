#Multilevel models--------------------------------------------------------------
rm(list=ls())
library(dplyr)
getwd()
load("Input/Base_multinivel.Rdata")

#Keep individuals with at least 3 waves
df<-subset(df, three_obs==TRUE)
#df <- df %>%
#group_by(folio) %>%
#filter(n_distinct(ola) >= 3)
df$pond_sens<-df$pond*df$sipw

## Variable centering ----------------------------------------------------------
df<- df %>% group_by(folio) %>% mutate(conf_inter_group=mean(conf_inter,na.rm=T))
df$conf_cwc<- df$conf_inter-df$conf_inter_group 
summary(df$conf_cwc)
summary(df$conf_inter_group)

df<- df %>% group_by(folio) %>% mutate(conf_ingroup_group=mean(conf_ingroup,na.rm=T))
df$conf_ingroup_cwc<- df$conf_ingroup-df$conf_ingroup_group 
summary(df$conf_ingroup_cwc)
summary(df$conf_ingroup_group)

df<- df %>% group_by(folio) %>% mutate(quant_inter_group=mean(quant_inter,na.rm=T))
df$quant_inter_cwc<- df$quant_inter-df$quant_inter_group 
summary(df$quant_inter_cwc)
summary(df$quant_inter_group)

df<- df %>% group_by(folio) %>% mutate(calidad_inter_group=mean(calidad_inter,na.rm=T))
df$calidad_inter_cwc<- df$calidad_inter-df$calidad_inter_group 
summary(df$calidad_inter_cwc)
summary(df$calidad_inter_group)

df<- df %>% group_by(folio) %>% mutate(neg_inter_group=mean(neg_inter,na.rm=T))
df$neg_inter_cwc<- df$neg_inter-df$neg_inter_group 
summary(df$neg_inter_cwc)
summary(df$neg_inter_group)

df<- df %>% group_by(folio) %>% mutate(presn_group=mean(presn,na.rm=T))
df$presn_cwc<- df$presn-df$presn_group 
summary(df$presn_cwc)
summary(df$presn_group)

# Descriptive statistic --------------------------------------------------------
library(dplyr)
## Attrition
df %>%
  group_by(ola) %>%
  summarise(n_people = n_distinct(folio)) %>%
  arrange(ola) %>%
  mutate(attrition_from_wave1 = 1 - n_people / first(n_people))

## Intraclass correlation-------------------------------------------------------
df_cl<-subset(df, indigena_es=="No indígena")
df_ind<-subset(df, indigena_es=="Indígena")
library(lme4)
m0 <- lmer(conf_inter ~ 1 + (1 |folio), weights = pond_sens, data=df_cl)
library(sjPlot)
sjPlot::tab_model(m0, title = "Modelo nulo", dv.labels = "Trust in Outgroup")

m0_ind <- lmer(conf_inter ~ 1 + (1 |folio), weights = pond_sens, data=df_ind)
library(sjPlot)
sjPlot::tab_model(m0_ind, title = "Modelo nulo", dv.labels = "Trust in Outgroup")

## Time specification-----------------------------------------------------------
#load packages
library(texreg)
library(GLMMadaptive)
library(lme4)
# non indigenous
#null model
m01c<-lmer(conf_inter~ (1|folio), weights = pond_sens, data=df_cl)


# Year Fixed model
m02c<-lmer(conf_inter~ as.numeric(ola) + (1|folio), weights = pond_sens, data=df_cl)


#Year as categorical
m04c<-lmer(conf_inter~ as.factor(ola) + (1|folio), weights = pond_sens, data=df_cl)


wordreg(c(m01c,m02c,m04c), file = "Output/sensitivity analysis/table1.docx", single.row = F, 
        custom.coef.names = c("Intercept", "Wave (Numeric)", 
                              "Wave 2 (ref = Wave 1)", "Wave 3 (ref = Wave 1)", 
                              "Wave 4 (ref = Wave 1)"), custom.model.names = 
          c("Null model", "Linear time model", "Categorical time model"))

# Indigenous

m01c_ind<-lmer(conf_inter~ (1|folio), weights = pond_sens, data=df_ind)


# Year Fixed model
m02c_ind<-lmer(conf_inter~ as.numeric(ola) + (1|folio), weights = pond_sens, data=df_ind)


#Year as categorical
m04c_ind<-lmer(conf_inter~ as.factor(ola) + (1|folio), weights = pond_sens, data=df_ind)


wordreg(c(m01c_ind,m02c_ind,m04c_ind), file = "Output/sensitivity analysis/table1_ind.docx", single.row = F, 
        custom.coef.names = c("Intercept", "Wave (Numeric)", 
                              "Wave 2 (ref = Wave 1)", "Wave 3 (ref = Wave 1)", 
                              "Wave 4 (ref = Wave 1)"), custom.model.names = 
          c("Null model", "Linear time model", "Categorical time model"))



## Non-Indigenous sample model--------------------------------------------------------

library(lme4)
library(texreg)
library(GLMMadaptive)

## Between variables Model
m6d<- lmer(conf_inter ~as.factor(ola) + g2 + g18 + g32_1 + quant_inter_group +
             calidad_inter_group+neg_inter_group+presn_group+(1|folio), 
           weights = pond_sens, data=df_cl)


## Within variables Model
m3d<- lmer(conf_inter ~as.factor(ola) + g2 + g18 + g32_1 + quant_inter_cwc+
             calidad_inter_cwc+neg_inter_cwc+presn_cwc+(1|folio), 
           weights = pond_sens, data=df_cl)


## Full Model
m9d<- lmer(conf_inter ~as.factor(ola) + g2 + g18 + g32_1 +
             quant_inter_group+
             calidad_inter_group+neg_inter_group+presn_group + quant_inter_cwc+
             calidad_inter_cwc+neg_inter_cwc+presn_cwc + (1|folio), 
           weights = pond_sens, data=df_cl)

## Interaction model
m9di<- lmer(conf_inter ~as.factor(ola) + g2 + g18 + g32_1 +
              quant_inter_group+
              calidad_inter_group+neg_inter_group+presn_group + quant_inter_cwc+
              calidad_inter_cwc+neg_inter_cwc+presn_cwc + calidad_inter_cwc*neg_inter_cwc + calidad_inter_group*neg_inter_group + (1|folio), weights = pond_sens, data=df_cl)


library(texreg)
wordreg(list(m6d, m3d, m9d, m9di), file = "Output/sensitivity analysis/models.docx", digits = 3, single.row = T, 
        custom.model.names = c("Between model", "Within model", "Ful model", "Interaction model"), 
        custom.coef.names = c("Intercept", "Wave 2 (ref = Wave 1)", "Wave 3 (ref = Wave 1)", 
                              "Wave 4 (ref = Wave 1)", "Age", 
                              "Sex","Educational level",  "Contact quantity (between)", 
                              "Contact quality (between)", 
                              "Negative contact frequency (between)", "Prescriptive norms (between)",
                              "Contact quantity (within)",
                              "Contact quality (within)", "Negative contact frequency (within)", 
                              "Prescriptive norms (within)", "Negative contact x Contact quality (within)", 
                              "Negative contact x Contact quality (between)"))

## Indigenous sample model--------------------------------------------------------

library(lme4)
library(texreg)
library(GLMMadaptive)

## Between variables Model
m6d_ind<- lmer(conf_inter ~as.factor(ola) + g2 + g18 + g32_1 + quant_inter_group +
                 calidad_inter_group+neg_inter_group+presn_group+(1|folio), 
               weights = pond_sens, data=df_ind)


## Within variables Model
m3d_ind<- lmer(conf_inter ~as.factor(ola) + g2 + g18 + g32_1 + quant_inter_cwc+
                 calidad_inter_cwc+neg_inter_cwc+presn_cwc+(1|folio), 
               weights = pond_sens, data=df_ind)


## Full Model
m9d_ind<- lmer(conf_inter ~as.factor(ola) + g2 + g18 + g32_1 +
                 quant_inter_group+
                 calidad_inter_group+neg_inter_group+presn_group + quant_inter_cwc+
                 calidad_inter_cwc+neg_inter_cwc+presn_cwc + (1|folio), 
               weights = pond_sens, data=df_ind)

## Interaction model
m9di_ind<- lmer(conf_inter ~as.factor(ola) + g2 + g18 + g32_1 +
                  quant_inter_group+
                  calidad_inter_group+neg_inter_group+presn_group + quant_inter_cwc+
                  calidad_inter_cwc+neg_inter_cwc+presn_cwc + calidad_inter_cwc*neg_inter_cwc + calidad_inter_group*neg_inter_group + (1|folio), weights = pond_sens, data=df_ind)


library(texreg)
wordreg(list(m6d_ind, m3d_ind, m9d_ind, m9di_ind), file = "Output/sensitivity analysis/models_ind.docx", digits = 3, single.row = T, 
        custom.model.names = c("Between model", "Within model", "Ful model", "Interaction model"), 
        custom.coef.names = c("Intercept", "Wave 2 (ref = Wave 1)", "Wave 3 (ref = Wave 1)", 
                              "Wave 4 (ref = Wave 1)", "Age", 
                              "Sex","Educational level",  "Contact quantity (between)", 
                              "Contact quality (between)", 
                              "Negative contact frequency (between)", "Prescriptive norms (between)",
                              "Contact quantity (within)",
                              "Contact quality (within)", "Negative contact frequency (within)", 
                              "Prescriptive norms (within)", "Negative contact x Contact quality (within)", 
                              "Negative contact x Contact quality (between)"))
## Plots -----------------------------------------------------------------------

library(interplot)
plot1<-interplot(m = m9di, var1 = "neg_inter_group", var2 = "calidad_inter_group") +  # Add labels for X and Y axes
  xlab("Contact quality levels (average)") +
  ylab("Estimated Coefficient for negative contact frequency") +
  # Change the background
  theme_bw() +
  # Add the title
  ggtitle("Estimated Coefficient of negative contact frequency by contact quality levels") +
  theme(plot.title = element_text(face="bold")) +
  # Add a horizontal line at y = 0
  geom_hline(yintercept = 0, linetype = "dashed")

plot2<-interplot(m = m9di, var1 = "calidad_inter_cwc", var2 = "neg_inter_cwc") +  # Add labels for X and Y axes
  xlab("negative contact frequency change (within)") +
  ylab("Estimated Coefficient for contact quality levels change (within)") +
  #Change the background
  theme_bw() +
  #Add the title
  ggtitle("Estimated Coefficient of contact quality increase over time by negative contact frequency change over time") +
  theme(plot.title = element_text(face="bold")) +
  #Add a horizontal line at y = 0
  geom_hline(yintercept = 0, linetype = "dashed")

plot3<-interplot(m = m9di_ind, var1 = "neg_inter_group", var2 = "calidad_inter_group") +  # Add labels for X and Y axes
  xlab("Contact quality levels (average)") +
  ylab("Estimated Coefficient for negative contact frequency") +
  # Change the background
  theme_bw() +
  # Add the title
  ggtitle("Estimated Coefficient of negative contact frequency by contact quality levels") +
  theme(plot.title = element_text(face="bold")) +
  # Add a horizontal line at y = 0
  geom_hline(yintercept = 0, linetype = "dashed")

graficos <- list(plot1, plot2, plot3)
# Nombres de archivos correspondientes
nombres_archivos <- c("Output/sensitivity analysis/Plot1.png", "Output/sensitivity analysis/Plot2.png", "Output/sensitivity analysis/Plot3_ind.png") # Añade los nombres correspondientes

# Dimensiones
width <- 38.09
height <- 18.95
units <- "cm"

# Guardar cada gráfico
for (i in 1:length(graficos)) {
  ggsave(nombres_archivos[i], plot = graficos[[i]], width = width, height = height, units = units)
}
