---
title: "New_suvr_mod"
date: '2023-05-11'
---

# This script is to fit parametric mixture cure model in SSMCM.

# multiple propensity score matching (PSM) for vaccination
data_filled <- dat_filled
data_filled$Vaccination.Status <- ifelse(data_filled$Vaccination.Status0==1, 0, ifelse(data_filled$Vaccination.Status1==1, 1, ifelse(data_filled$Vaccination.Status2==1, 2, -1)))
data_filled <- data_filled[data_filled$Vaccination.Status!=-1,]


catVars <- c("Gendermale", "Vaccination.Status", "Symptom.FeverYes", "Cardiovascular.disease.including.hypertensionYes", "Is.Home.Quarantine.Yes")
dat_var <- data_filled[data_filled$Vaccination.Status %in% c(0,1),]
dat_var[,catVars] <- lapply(dat_var[,catVars], as.factor)
match_mod <- matchit(Vaccination.Status ~ Gendermale + Symptom.FeverYes + Cardiovascular.disease.including.hypertensionYes + Is.Home.Quarantine.Yes, data = dat_var, method = "nearest")
matchdata01 <- match.data(match_mod)

dat_var <- data_filled[data_filled$Vaccination.Status %in% c(1,2),]
# dat_var$Symptom.Fever[which(dat_var$Symptom.Fever!="Yes")] <- "No"
dat_var[,catVars] <- lapply(dat_var[,catVars], as.factor)
match_mod <- matchit(Vaccination.Status ~ Gendermale + Symptom.FeverYes + Cardiovascular.disease.including.hypertensionYes + Is.Home.Quarantine.Yes, data = dat_var, method = "nearest")
matchdata12 <- match.data(match_mod)

matchdata <- rbind(matchdata01, matchdata12[matchdata12$Vaccination.Status==2,])


Vars <- c("Age", "Gender", "Vaccination.Status", "Patient.Status", "Symptom.Fever", "Cardiovascular.disease.including.hypertension", "Is.Home.Quarantine.", "Put.On.Ventilator", "Period")
catVars <- c("Gender", "Vaccination.Status", "Patient.Status", "Symptom.Fever", "Cardiovascular.disease.including.hypertension", "Is.Home.Quarantine.", "Put.On.Ventilator", "Period")
nonVars <- "Age"
matchdata$Age <- as.numeric(matchdata$Age)
table1_gender <- CreateTableOne(vars = Vars, 
                                factorVars = catVars,
                                strata = "Gender",
                                data = matchdata,
                                addOverall = T)
table1_gender_p <- print(table1_gender, showAllLevels = T,
                         catDigits = 1, contDigits = 0, pDigits = 4,
                         nonnormal = nonVars)
# table1_p <- cbind(table1_ICU_p, table1_gender_p)
# table1_p <- table1_p[,-c(6,12)]
write.csv(table1_gender_p, "/Users/wuhanyu/Desktop/Research/LBS/Covid-19/table1_gender.csv", quote = F)


# multiple propensity score matching (PSM) for gender
catVars <- c("Gendermale", "Vaccination.Status", "Symptom.Fever", "Cardiovascular.disease.including.hypertension", "Is.Home.Quarantine.")
Vars <- c("Age", "Gender", "Vaccination.Status", "Patient.Status", "Symptom.Fever", "Cardiovascular.disease.including.hypertension", "Is.Home.Quarantine.", "Put.On.Ventilator", "Period")
catVars <- c("Gender", "Vaccination.Status", "Patient.Status", "Symptom.Fever", "Cardiovascular.disease.including.hypertension", "Is.Home.Quarantine.", "Put.On.Ventilator", "Period")
nonVars <- "Age"
dat_var <- data_all
dat_var[,catVars] <- lapply(dat_var[,catVars], as.factor)
match_mod <- matchit(Gender ~ Vaccination.Status + Symptom.Fever + Cardiovascular.disease.including.hypertension + Is.Home.Quarantine., data = dat_var, method = "nearest")
matchdata <- match.data(match_mod, distance = "Score", weights = "Weight", subclass = "Subclass")


# Weibull distribution
# MLE
nll_weibull <- function(lambda, k){
  eta <- as.numeric(data_filled$eta)
  t <- as.numeric(data_filled$T)
  B <- as.numeric(data_filled$B)
  pi <- as.numeric(data_filled$.pred)
  -sum( (1-B)*(log(pi) - lambda*(t)^k + eta*(log(k)+log(lambda)+(k-1)*log(t))) + B*log(1-pi) )
}
nll_log_logis <- function(lambda, k){
  eta <- as.numeric(data_filled$eta)
  t <- as.numeric(data_filled$T)
  B <- as.numeric(data_filled$B)
  pi <- as.numeric(data_filled$.pred)
  -sum( (1-B)*(log(pi) - (1+eta)*log(1+(lambda*t)^k) + eta*(log(k)+log(lambda)+(k-1)*log(lambda*t))) + B*log(1-pi) )
}
nll_log_norm <- function(mu, sigma){
  eta <- as.numeric(data_filled$eta)
  t <- as.numeric(data_filled$T)
  B <- as.numeric(data_filled$B)
  pi <- as.numeric(data_filled$.pred)
  Z <- (log(t)-mu)/sigma
  -sum( (1-B)*(log(pi) + (1-eta)*log(1-dnorm(Z)) + eta*(log(pnorm(Z))-log(sigma*t))) + B*log(1-pi) )
}
nll_gamma <- function(beta, gamma){
  eta <- as.numeric(data_filled$eta)
  t <- as.numeric(data_filled$T)
  B <- as.numeric(data_filled$B)
  pi <- as.numeric(data_filled$.pred)
  Z <- t/beta
  ga <- expint::gammainc(gamma, t)
  -sum( (1-B)*(log(pi) + (1-eta)*log(ga/gamma(gamma)) + eta*((gamma-1)*log(Z)-Z-log(beta)-log(gamma(gamma))) ) + B*log(1-pi) )
}


# optimizing
result <- data_frame()
SOD=vector()
AIC=vector()
method=vector()
low=vector()
sup=vector()
par1=vector()
par2=vector()
for (i in 0:3) {
  if(i %in% matchdata$Vaccination.Status)
    data_filled <- matchdata[matchdata$Vaccination.Status==i,]
  else
    data_filled <- matchdata
  # data_filled <- matchdata[matchdata$Gendermale==i,]
  # make sure the nll is always finite
  data_filled[data_filled$.pred==1, ".pred"] <- 0.999
  data_filled[data_filled$.pred==0, ".pred"] <- 0.001
  data_filled <- data_filled[data_filled$T!="", ]
  data_filled[data_filled$T==0, "T"] <- 0.5
  
  est <- bbmle::mle2(minuslogl = nll_weibull, start = list(lambda=1, k=1), method = "Nelder-Mead")
  lambda = est@coef[1]
  k = est@coef[2]
  SOD[i+1] <- lambda^(-1/k)*gamma(1+1/k)
  par1[i+1] <- lambda
  par2[i+1] <- k
  lambda = bbmle::confint(est)[1,1]
  k = bbmle::confint(est)[2,1]
  low[i+1] <- lambda^(-1/k)*gamma(1+1/k)
  lambda = bbmle::confint(est)[1,2]
  k = bbmle::confint(est)[2,2]
  sup[i+1] <- lambda^(-1/k)*gamma(1+1/k)
  AIC[i+1] = AIC(est)
  method[i+1] = "weibull"
  
  est <- bbmle::mle2(minuslogl = nll_log_logis, start = list(lambda=1,k=1), method = "Nelder-Mead")
  if(AIC(est) < AIC[i+1]){
    shape = as.numeric(est@coef[2])
    scale = 1/as.numeric(est@coef[1])
    SOD[i+1] = mean_llogis(shape = shape, scale = scale)
    par1[i+1] <- shape
    par2[i+1] <- scale
    shape = bbmle::confint(est)[2,1]
    scale = 1/bbmle::confint(est)[1,1]
    low[i+1] = mean_llogis(shape = shape, scale = scale)
    shape = bbmle::confint(est)[2,2]
    scale = 1/bbmle::confint(est)[1,2]
    sup[i+1] = mean_llogis(shape = shape, scale = scale)
    AIC[i+1] = AIC(est)
    method[i+1] = "llogis"
  }
  
  # est <- bbmle::mle2(minuslogl = nll_log_norm, start = list(mu=0, sigma=1), method = "Nelder-Mead")
  # if(AIC(est) < AIC[i+1]){
  #   mu = as.numeric(est@coef[1])
  #   sigma = as.numeric(est@coef[2])
  #   SOD[i+1] = exp(mu + sigma^2/2)
  #   AIC[i+1] = AIC(est)
  #   method[i+1] = "lnorm"
  # }
  
  est <- bbmle::mle2(minuslogl = nll_gamma, start = list(beta=1, gamma=1), method = "Nelder-Mead")
  if(AIC(est) < AIC[i+1]){
    beta = as.numeric(est@coef[1])
    gamma = as.numeric(est@coef[2])
    SOD[i+1] = beta*gamma
    par1[i+1] <- beta
    par2[i+1] <- gamma
    beta = bbmle::confint(est)[1,1]
    gamma = bbmle::confint(est)[2,1]
    low[i+1] = beta*gamma
    beta = bbmle::confint(est)[1,2]
    gamma = bbmle::confint(est)[2,2]
    sup[i+1] = beta*gamma
    AIC[i+1] = AIC(est)
    method[i+1] = "gamma"
  }
}
result <- data_frame(Vaccine=0:3, Median=SOD, Model=method, lower=sup, upper=low, AIC=AIC, par1=par1, par2=par2)
write.csv(result, "/Users/wuhanyu/Desktop/Research/LBS/Covid-19/Result/Vaccine_table.csv", quote = F, row.names = F)


# Visualization for gender
shapeM=2.637973
scaleM=16.45674
data_M <- matchdata[matchdata$Gendermale==1,]
plot_SSLR_Male <- ggplot(data.frame(x = c(1:1000)), aes(x)) + 
  geom_histogram(data = data_M, 
                 aes(x = T, y = ..density..), 
                 color = "black", fill = "grey") + 
  geom_line(aes(y = dllogis(x, shape = shapeM, scale = scaleM)), 
            color = "red", size = 2) + 
  xlim(0, 300) + 
  ylim(0, 0.05) + 
  theme_bw() +  # turn background to white
  labs(x = "Symptoms onset to death (days)", y = "Density", title = "Male cases") + 
  theme(panel.grid.major = element_blank(),     # hide grid line
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour ="black"),
        plot.title = element_text(hjust = 0.5)) # title centered

shapeF=3.46373 
scaleF=14.8509 
data_F <- matchdata[matchdata$Gendermale==0,]
plot_SSLR_Female <- ggplot(data.frame(x = c(1:1000)), aes(x)) + 
  geom_histogram(data = data_F, 
                 aes(x = T, y = ..density..), 
                 color = "black", fill = "grey") + 
  geom_line(aes(y = dllogis(x, shape = shapeF, scale = scaleF)), 
            color = "red", size = 2) + 
  xlim(0, 300) + 
  ylim(0, 0.065) + 
  theme_bw() +  # turn background to white
  labs(x = "Symptoms onset to death (days)", y = "Density", title = "Female cases") + 
  theme(panel.grid.major = element_blank(),     # hide grid line
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour ="black"),
        plot.title = element_text(hjust = 0.5)) # title centered

shapeA=3.11827
scaleA=15.50752
data_filled_A <- matchdata
plot_SSLR_all <- ggplot(data.frame(x = c(1:1000)), aes(x)) + 
  geom_histogram(data = data_filled_A, 
                 aes(x = T, y = ..density..), 
                 color = "black", fill = "grey") + 
  geom_line(aes(y = dllogis(x, shape = shapeA, scale = scaleA)), 
            color = "red", size = 2) + 
  xlim(0, 300) + 
  ylim(0, 0.06) + 
  theme_bw() +  # turn background to white
  labs(x = "Symptoms onset to death (days)", y = "Density", title = "All cases") + 
  theme(panel.grid.major = element_blank(),     # hide grid line
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour ="black"),
        plot.title = element_text(hjust = 0.5)) # title centered

# Visualization for vaccination
shape0=3.405978
scale0=15.03074
data_0 = matchdata[matchdata$Vaccination.Status == 0,]
plot_SSLR_Vacc0 <- ggplot(data.frame(x = c(1:1000)), aes(x)) + 
  geom_histogram(data_0, 
                 mapping=aes(x = T, y = ..density..), 
                 color = "black", fill = "grey") + 
  geom_line(aes(y = dllogis(x, shape = shape0, scale = scale0)), 
            color = "red", size = 2) + 
  xlim(0, 300) + 
  ylim(0, 0.065) + 
  theme_bw() +  # turn background to white
  labs(x = "Symptoms onset to death (days)", y = "Density", title = "Vaccine_0 cases") + 
  theme(panel.grid.major = element_blank(),     # hide grid line
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour ="black"),
        plot.title = element_text(hjust = 0.5)) # title centered

shape1=3.288829
scale1=15.23077
data_1 = matchdata[matchdata$Vaccination.Status == 1,]
plot_SSLR_Vacc1 <- ggplot(data.frame(x = c(1:1000)), aes(x)) + 
  geom_histogram(data_1, 
                 mapping=aes(x = T, y = ..density..), 
                 color = "black", fill = "grey") + 
  geom_line(aes(y = dllogis(x, shape = shape1, scale = scale1)), 
            color = "red", size = 2) + 
  xlim(0, 300) + 
  ylim(0, 0.06) + 
  theme_bw() +  # turn background to white
  labs(x = "Symptoms onset to death (days)", y = "Density", title = "Vaccine_1 cases") + 
  theme(panel.grid.major = element_blank(),     # hide grid line
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour ="black"),
        plot.title = element_text(hjust = 0.5)) # title centered

shape2=2.74632
scale2=16.37457
data_2 = matchdata[matchdata$Vaccination.Status == 2,]
plot_SSLR_Vacc2 <- ggplot(data.frame(x = c(1:1000)), aes(x)) + 
  geom_histogram(data_2, 
                 mapping=aes(x = T, y = ..density..), 
                 color = "black", fill = "grey") + 
  geom_line(aes(y = dllogis(x, shape = shape2, scale = scale2)), 
            color = "red", size = 2) + 
  xlim(0, 300) + 
  ylim(0, 0.05) + 
  theme_bw() +  # turn background to white
  labs(x = "Symptoms onset to death (days)", y = "Density", title = "Vaccine_2 cases") + 
  theme(panel.grid.major = element_blank(),     # hide grid line
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour ="black"),
        plot.title = element_text(hjust = 0.5)) # title centered

(plot_SSLR_Male + plot_SSLR_Female + plot_SSLR_all) / (plot_SSLR_Vacc0 + plot_SSLR_Vacc1 + plot_SSLR_Vacc2) + plot_annotation(tag_levels = "a")

