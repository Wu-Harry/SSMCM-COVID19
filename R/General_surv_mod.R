---
title: "General_surv_mod"
date: '2023-03-14'
---

# This script is to fit general paramatric survival model use 'flexsurv' package.


# general survival models for PSM data (alpha & delta)
data_all <- matchdata[,c("Age", "Gender", "Vaccination.Status", "Patient.Status", "FLU", "Symptom.Fever", "Fever.History", "Symptom.Sore.Throat", "Symptom.Cough", "Symptom.Diarrehea", "Symptom.Breathing.issue", "Headache", "Cardiovascular.disease.including.hypertension", "chronic.lung.disease", "Diabetes", "Pregnancy", "Is.Home.Quarantine.", "ICU.Admission", "Put.On.Ventilator", "T", "Period")] %>% .[.$T!="",]
for (i in 5:12) {
  data_all[which(data_all[,i]==""), i] <- "No"
}
# data <- data_all[data_all$ICU.Admission=="Yes",]

# Alpha / Delta
library(flexsur)
data <- na.omit(data_all[data_all$Period=="alpha",]) # delta
# SOD
data <- as.data.frame(data[data$Patient.Status!="Recovered",])
data$Status = ""
for (i in 1:length(data$Patient.Status)) {
  if(data$Patient.Status[i] == "Died")
    data$Status[i] = 2
  else
    data$Status[i] = 1
}
data$T <- as.numeric(data$T)
data$Status <- as.numeric(data$Status)
data$T[which(data$T==0)] = 0.5
res <- flexsurvreg(formula = Surv(T, Status) ~ 1, data = data, dist = "gamma")
res <- flexsurvreg(formula = Surv(T, Status) ~ 1, data = data, dist = "weibull")
res <- flexsurvreg(formula = Surv(T, Status) ~ 1, data = data, dist = "llogis")
res <- flexsurvreg(formula = Surv(T, Status) ~ 1, data = data, dist = "lnorm")
mean_llogis(shape = 1.7282, scale = 7.9754) # alpha
mean_llogis(shape = 1.6539, scale = 7.5756)
mean_llogis(shape = 1.8059, scale = 8.3964)
mean_llogis(shape = 1.8432, scale = 8.2644) # delta
mean_llogis(shape = 1.7657, scale = 7.8882)
mean_llogis(shape = 1.9242, scale = 8.6586)
# For alpha, Log-logistic distribution provides the smallest AIC (9838.21), the median SOD is 14.95196.
# For delta, Log-logistic distribution provides the smallest AIC (9518.959), the median SOD is 14.21274.

# SOR
data <- data_all[data_all$Period=="alpha",] # delta
data <- as.data.frame(data[data$Patient.Status!="Died",])
for (i in 1:length(data$Patient.Status)) {
  if(data$Patient.Status[i] == "Recovered")
    data$Status[i] = 2
  else
    data$Status[i] = 1
}
data$T <- as.numeric(data$T)
data$Status <- as.numeric(data$Status)
data$T[which(data$T==0)] = 0.5
res <- flexsurvreg(formula = Surv(T, Status) ~ 1, data = data, dist = "gamma")
res <- flexsurvreg(formula = Surv(T, Status) ~ 1, data = data, dist = "weibull")
res <- flexsurvreg(formula = Surv(T, Status) ~ 1, data = data, dist = "llogis")
res <- flexsurvreg(formula = Surv(T, Status) ~ 1, data = data, dist = "lnorm")
mean_llogis(shape = 3.3151, scale = 15.3793) # alpha
mean_llogis(shape = 3.8336, scale = 14.2480) # delta
# For alpha, Log-logistic distribution provides smallest AIC (271939.7), the median SOR is 17.94756.
# For delta, Log-logistic distribution provides smallest AIC (243788.7), the median SOR is 15.97718.


# SOD
data_all <- dat[,c("Age", "Gender", "Vaccination.Status", "Patient.Status", "FLU", "Symptom.Fever", "Fever.History", "Symptom.Sore.Throat", "Symptom.Cough", "Symptom.Diarrehea", "Symptom.Breathing.issue", "Headache", "Cardiovascular.disease.including.hypertension", "chronic.lung.disease", "Diabetes", "Pregnancy", "Is.Home.Quarantine.", "ICU.Admission", "Put.On.Ventilator", "T", "Period")] %>% .[.$T!="",]
data <- data_all[data_all$ICU.Admission=="No",] # Yes
data <- as.data.frame(data[data$Patient.Status!="Recovered",])
data$Status = ""
for (i in 1:length(data$Patient.Status)) {
  if(data$Patient.Status[i] == "Died")
    data$Status[i] = 2
  else
    data$Status[i] = 1
}
data$T <- as.numeric(data$T)
data$Status <- as.numeric(data$Status)
data$T[which(data$T==0)] = 0.5
res <- flexsurvreg(formula = Surv(T, Status) ~ 1, data = data, dist = "gamma")
res <- flexsurvreg(formula = Surv(T, Status) ~ 1, data = data, dist = "weibull")
res <- flexsurvreg(formula = Surv(T, Status) ~ 1, data = data, dist = "llogis")
res <- flexsurvreg(formula = Surv(T, Status) ~ 1, data = data, dist = "lnorm")
mean_llogis(shape = 1.7375, scale = 8.1082) # non-ICU SOD
mean_llogis(shape = 1.7031, scale = 6.3297) # ICU SOD
mean_llogis(shape = 1.732, scale = 7.977) # All
# For non-ICU SOD, Log-logistic distribution provides smallest AIC (37672.66), the median SOD is 15.08327.
# For ICU SOD, Log-logistic distribution provides smallest AIC (2573.261), the median SOD is 12.12784.
# For all, Log-logistic distribution provides smallest AIC (40264.94), the median SOD is 14.90727.

# SOR
data <- data_all[data_all$ICU.Admission=="No",] # Yes
data <- as.data.frame(data[data$Patient.Status!="Died",])
data$Status = ""
for (i in 1:length(data$Patient.Status)) {
  if(data$Patient.Status[i] == "Recovered")
    data$Status[i] = 2
  else
    data$Status[i] = 1
}
data$T <- as.numeric(data$T)
data$Status <- as.numeric(data$Status)
data$T[which(data$T==0)] = 0.5
res <- flexsurvreg(formula = Surv(T, Status) ~ 1, data = data, dist = "gamma")
res <- flexsurvreg(formula = Surv(T, Status) ~ 1, data = data, dist = "weibull")
res <- flexsurvreg(formula = Surv(T, Status) ~ 1, data = data, dist = "llogis")
res <- flexsurvreg(formula = Surv(T, Status) ~ 1, data = data, dist = "lnorm")
mean_llogis(shape = 2.42186, scale = 16.46294) # non-ICU SOR
mean_llogis(shape = 1.770, scale = 27.016) # ICU SOR
mean_llogis(shape = 2.42130, scale = 16.46549) # All
# For non-ICU SOR, Log-logistic distribution provides smallest AIC (1377163), the median SOR is 22.18053.
# For ICU SOR, Log-logistic distribution provides smallest AIC (569.5771), the median SOR is 48.96752.
# For all SOR, Log-logistic distribution provides smallest AIC (1377761), the median SOR is 22.18723.



# visualization
library(ggplot2)
data_filled$T <- as.numeric(data_filled$T)
plot_SSLR_Male <- ggplot(data.frame(x = c(1:1000)), aes(x)) + 
  geom_histogram(data = data_filled[data_filled$Gendermale == 1,], 
                 aes(x = T, y = ..density..), 
                 color = "black", fill = "grey") + 
  geom_line(aes(y = dllogis(x, shape = shape, scale = scale)), 
            color = "red", size = 2) + 
  xlim(0, 400) + 
  ylim(0, 0.015) + 
  theme_bw() +  # turn background to white
  labs(x = "Symptoms onset to death (days)", y = "Density", title = "Male cases", caption = "caption") + 
  theme(panel.grid.major = element_blank(),     # hide grid line
        panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust = 0.5)) # title centered

plot_SSLR_Female <- ggplot(data.frame(x = c(1:1000)), aes(x)) + 
  geom_histogram(data = data_filled[data_filled$Gendermale == 0,], 
                 aes(x = T, y = ..density..), 
                 color = "black", fill = "grey") + 
  geom_line(aes(y = dllogis(x, shape = shape, scale = scale)), 
            color = "red", size = 2) + 
  xlim(0, 400) + 
  ylim(0, 0.015) + 
  theme_bw() +  # turn background to white
  labs(x = "Symptoms onset to death (days)", y = "Density", title = "Female cases", caption = "caption") + 
  theme(panel.grid.major = element_blank(),     # hide grid line
        panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust = 0.5)) # title centered
plot_SSLR_Male + plot_SSLR_Female

