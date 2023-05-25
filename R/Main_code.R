---
title: "Main_code"
---

# Main code of the whole project, includes data cleaning, PSM, SSL-RF, mixture cure model, Cox model, etc.


# Read data
library(xlsx)
library(dplyr)
dat <- read.csv("/Users/wuhanyu/Desktop/Research/LBS/Covid-19/Data/DGHS_Covid 19-2020-2021.csv", header = T, row.names = 1) %>% 
  na.omit() %>% 
  .[.$Date.of.status!="",] %>% 
  .[.$Date.of.Illness!="",] %>% 
  .[.$Age!="Nil",] %>%
  .[.$Gender!="other",]
dat[which(dat$Vaccination.Status==c("partially", "yes")),"Vaccination.Status"] = 1
dat[which(dat$Vaccination.Status=="fully"),"Vaccination.Status"] = 2
dat[which(dat$Vaccination.Status=="no"),"Vaccination.Status"] = 0
dat[which(dat$Patient.Status=="Expired"),"Patient.Status"] = "Died"
dat[which(dat$Vaccination.Status %in% c("partially","yes")),"Vaccination.Status"] = 1

# Split dataset by patient status
dat_dead <- dat[dat$Patient.Status=="Died",]
dat_reco <- dat[dat$Patient.Status=="Recovered",]
dat_cens <- dat[dat$Patient.Status=="Active",]


# calculate SOD/SOR/censored time
library(tidyr)
dat$ID <- rownames(dat)
for (i in 1:length(dat$Date.of.Recover)) {
  if(dat$Date.of.Recover[i]!="")
    dat$Date[i] <- dat$Date.of.Recover[i]
}
for (i in 1:length(dat$Date.of.Death)) {
  if(dat$Date.of.Death[i]!="")
    dat$Date[i] <- dat$Date.of.Death[i]
}
for (i in 1:length(dat$Date.of.Death)) {
  if(dat$Patient.Status[i]=="Active")
    dat$Date[i] <- dat$Date.of.status[i]
}
date <- dat[, c("Date.of.Illness","Date", "ID")]
date <- na.omit(date)
date <- separate(date, Date.of.Illness, into = c("day1", "month1", "year1"), sep = "/", remove = F) %>%
  separate(Date, into = c("day2", "month2", "year2"), sep = "/", remove = F)
date$date1 <- paste(date$year1, "-", date$month1, "-", date$day1, sep = "")
date$date2 <- paste(date$year2, "-", date$month2, "-", date$day2, sep = "")
date$time <- difftime(date$date2, date$date1, units = "day")
duration <- separate(date, time, into = c("time", "unit")) %>% .[,colnames(.)=="time"] %>% cbind(., date$ID)
colnames(duration) <- c("T", "ID")
dat <- merge(dat, duration, by = "ID")


# identify alpha & delta period
date <- as.data.frame(lapply(date, as.numeric))
date$period <- NA
date[which((date$year1==2021)&
             ((date$month1==3)|
              (date$month1==4)|
              ((date$month1==2)&(date$day1>25))|
              ((date$month1==5)&(date$day1<16)))), 12] <- "alpha"
date[which((date$year1==2021)&
             ((date$month1==8)|
              (date$month1==9)|
              (date$month1==10)|
              ((date$month1==7)&(date$day1>6))|
              ((date$month1==11)&(date$day1<31)))), 12] <- "delta"
dat <- cbind(dat, date$period) %>% rename("Period" = "date$period")


# baseline table (Table 1)
library(tableone)
dat$Age <- as.numeric(dat$Age)
Vars <- c("Age", "Gender", "Vaccination.Status", "Patient.Status", "Symptom.Fever", "Cardiovascular.disease.including.hypertension", "Is.Home.Quarantine.", "Put.On.Ventilator", "Period")
catVars <- c("Gender", "Vaccination.Status", "Patient.Status", "Symptom.Fever", "Cardiovascular.disease.including.hypertension", "Is.Home.Quarantine.", "Put.On.Ventilator", "Period")
nonVars <- "Age"
table1_ICU <- CreateTableOne(vars = Vars, 
                         factorVars = catVars,
                         strata = "ICU.Admission",
                         data = dat,
                         addOverall = T)
table1_ICU_p <- print(table1_ICU, showAllLevels = T,
      catDigits = 1, contDigits = 0, pDigits = 4,
      nonnormal = nonVars)

# propensity score matching (PSM) for 2 periods
library(MatchIt)
dat_var <- na.omit(dat)
dat_var[,catVars] <- lapply(dat_var[,catVars], as.factor)
# match_mod <- matchit(Period ~ Age + Gender + Vaccination.Status + Patient.Status + Symptom.Fever + Cardiovascular.disease.including.hypertension + Is.Home.Quarantine. + Put.On.Ventilator + ICU.Admission, data = dat_var, method = "nearest")
match_mod <- matchit(Period ~ Age + Gender, data = dat_var, method = "nearest", caliper = 0.01)
matchdata <- match.data(match_mod)
table1_period <- CreateTableOne(vars = Vars, 
                                factorVars = catVars,
                                strata = "Period",
                                data = matchdata,
                                addOverall = T)
table1_period_p <- print(table1_period, showAllLevels = T,
                         catDigits = 1, contDigits = 0, pDigits = 4,
                         nonnormal = nonVars)
table1_p <- cbind(table1_ICU_p, table1_period_p)
table1_p <- table1_p[,-c(6,12)]
write.csv(table1_p, "/Users/wuhanyu/Desktop/Research/LBS/Covid-19/table1.csv", quote = F)


# SSL-RF
library(SSLR)
library(caret)
# Data used to estimate \pi(x)
data_all <- dat[,c("Age", "Gender", "Vaccination.Status", "Patient.Status", "FLU", "Symptom.Fever", "Fever.History", "Symptom.Sore.Throat", "Symptom.Cough", "Symptom.Diarrehea", "Symptom.Breathing.issue", "Headache", "Cardiovascular.disease.including.hypertension", "chronic.lung.disease", "Diabetes", "Pregnancy", "Is.Home.Quarantine.", "ICU.Admission", "Put.On.Ventilator", "T", "Period")] %>% .[.$T!="",]
for (i in 5:12) {
  data_all[which(data_all[,i]==""), i] <- "No"
}
data <- data_all
data <- data_all[data_all$Period=="alpha",]
data <- data_all[data_all$ICU.Admission=="Yes",]

# data <- as.data.frame(data.matrix(data))

# data.matrix <- model.matrix(~ Gender + Vaccination.Status + FLU + Symptom.Fever + Fever.History + Symptom.Sore.Throat + Symptom.Cough + Headache + Diabetes + Is.Home.Quarantine. + ICU.Admission + Put.On.Ventilator, data = data)
# data.matrix <- as.data.frame(data.matrix)
# reg <- as.data.frame(model.matrix(~ Patient.Status, data = data))

# active=1, recover=3, dead=2
dat_labeled <- data[data$Patient.Status!="Active",]
time_labeled <- dat_labeled$T
dat_unlabeled <- data[data$Patient.Status=="Active",]
time_unlabeled <- dat_unlabeled$T
patient_status <- which(colnames(data) == "Patient.Status")
# data_unlabeled[,patient_status] <- NA
age_labeled <- dat_labeled[,1]
age_unlabeled <- dat_unlabeled[,1]
debug_contrast_error(dat_labeled) # check which column has only one class
debug_contrast_error(dat_unlabeled)

# Create labeled data for fitting
data.matrix <- model.matrix(
  ~ Patient.Status + Gender + Vaccination.Status + FLU + Symptom.Fever + Fever.History + Symptom.Sore.Throat + Symptom.Cough + Headache + Cardiovascular.disease.including.hypertension + Diabetes + Is.Home.Quarantine., 
  data = dat_labeled)
data_labeled <- as.data.frame(data.matrix) %>% cbind(.,age_labeled)
names(data_labeled)[names(data_labeled)=="age_labeled"] <- "age"

# Create unlabeled data for fitting
data.matrix <- model.matrix(
  ~ Gender + Vaccination.Status + Symptom.Fever + FLU + Fever.History + Symptom.Sore.Throat + Symptom.Cough + Headache + Cardiovascular.disease.including.hypertension + Diabetes + Is.Home.Quarantine., 
  data = dat_unlabeled)
data_unlabeled <- as.data.frame(data.matrix) %>% cbind(.,age_unlabeled)
data_unlabeled$Patient.StatusRecovered = rep(NA, length(data_unlabeled[,1])) # add patient status to the data
names(data_unlabeled)[names(data_unlabeled)=="age_unlabeled"] <- "age"

# Split train & test data
set.seed(1)
data_died <- data_labeled[data_labeled$Patient.StatusRecovered==0,] # died samples all in
data_recovered <- data_labeled[data_labeled$Patient.StatusRecovered==1,] # recovered samples partially in
train_index <- createDataPartition(y = data_labeled$Patient.StatusRecovered, p = 0.2, list = F)
data_train <- merge(data_recovered[train_index,], data_died, all = T)
data_train <- merge(data_train, data_unlabeled, all = T)
data_test <- data_recovered[-train_index,]
# data_train1 <- data_train[c(1:2000,119000:121947),]

# Semi-supervised RF
# m <- SSLRDecisionTree(min_samples_split = round(length(train_index) * 0.25), w = 0.3)
m <- SSLRRandomForest(trees = 100, w = 0.9) %>% fit(Patient.StatusRecovered ~ ., data = data_train)

# predict for unlabeled
predict_unlabeled <- predict(m, data_unlabeled)
unlabeled <- cbind(data_unlabeled, predict_unlabeled, time_unlabeled)
unlabeled$Patient.StatusRecovered <- as.numeric(unlabeled$.pred>0.5) # fill up cure information
unlabeled <- rename(unlabeled, "T" = "time_unlabeled")
unlabeled$censored <- 0

# predict for labeled
predict_labeled <- predict(m, data_labeled)
labeled <- cbind(data_labeled, predict_labeled, time_labeled) %>% 
  rename("T" = "time_labeled")
labeled$censored <- 1
# labeled$.pred[labeled$Patient.StatusRecovered==1] <- 1 # use 1/0 instead of predicted .pred for labeled
# labeled$.pred[labeled$Patient.StatusRecovered==0] <- 0
data_filled <- merge(unlabeled, labeled, all = T)

# mod_RF <- SSLRDecisionTree(max_depth = 20, w = 0.3)
# fit(Patient.Status ~ ., data = data_train,object = mod_RF)
# fit_decision_tree(object = mod_RF, y = Patient.Status, x = data_train)

# prediction accuracy (这部分意义不大，用ROC替换)
predict_labeled$status <- (predict_labeled$.pred>0.5)
predict_labeled_comb <- as.data.frame(cbind(predict_labeled$status, !!(data_labeled$Patient.StatusRecovered)))
predict_labeled_comb[[2]] <- as.factor(predict_labeled_comb[[2]])
predict_labeled_comb[[1]] <- as.factor(predict_labeled_comb[[1]])
yardstick::metrics(predict_labeled_comb, truth = "V2", estimate = "V1", na.rm = F)

# ROC curve and AUC
label_test <- as.data.frame(cbind(data_labeled$Patient.StatusRecovered, predict_labeled))
colnames(label_test) <- c("true", "pred")
ROCres <- pROC::roc(label_test$true, label_test$pred)
pROC::coords(ROCres, "best") # find best criteria
pROC::plot.roc(ROCres) # ROC curve
pROC::auc(ROCres) # AUC

# Platt Calibration for the model
calib_res <- eRic::prCalibrate(r.calib = label_test$true, p.calib = label_test$pred, nbins = 10)
ROCres <- pROC::roc(calib_res$responses, calib_res$cal.probs, auc = TRUE)
pROC::coords(ROCres, "best", ret="threshold", best.method="closest.topleft") # find best criteria by closest.topleft
pROC::coords(ROCres, "best", ret="threshold", best.method="youden") # find best threshold by Youden's Index
# threshold specificity sensitivity: 0.9813315	0.8543024	0.8249631		
pROC::plot.roc(ROCres)
plot(ROCres, legacy.axes = TRUE, main="ROC curve best threshold", 
     thresholds="best", # 基于youden指数选择roc曲线最佳阈值点
     print.thres="best") # 在roc曲线上显示最佳阈值点


# Compare different w (0.1~0.9) 服务器跑好后存到本地，直接读结果然后计算AUC
j=1 
thres=vector()
auc=vector()
speci=vector()
sensi=vector()
for (i in seq(0.1:0.9, by=0.2)) {
  path <- paste("/Users/wuhanyu/Desktop/Research/LBS/Covid-19/Semi-supRF/Semi-supRF-",i, "/predict_labeled.csv", sep = "")
  predict_labeled <- read.csv(path)[,-1]
  label_test <- as.data.frame(cbind(data_labeled$Patient.StatusRecovered, predict_labeled))
  colnames(label_test) <- c("true", "pred")
  calib_res <- eRic::prCalibrate(r.calib = label_test$true, p.calib = label_test$pred, nbins = 10)
  ROCres <- pROC::roc(calib_res$responses, calib_res$cal.probs, auc = TRUE)
  thres[j] <- pROC::coords(ROCres, "best", ret="threshold", best.method="youden")
  auc[j] <- pROC::auc(ROCres)
  speci[j] <- pROC::coords(ROCres, "best", ret="specificity", best.method="youden")
  sensi[j] <- pROC::coords(ROCres, "best", ret="sensitivity", best.method="youden")
  j=j+1
}
roc_res <- cbind(w=seq(0.1:0.9, by=0.2),thres,auc,speci,sensi)
write.csv(roc_res, "/Users/wuhanyu/Desktop/Research/LBS/Covid-19/Result/ROC_table.csv", row.names = F, quote = F)
rm(thres,auc,speci,sensi,j)
plot(ROCres, legacy.axes = TRUE, main="ROC curve best threshold", 
     thresholds="best", # 基于youden指数选择roc曲线最佳阈值点
     print.thres="best",
     ylim=c(0,1)) # 在roc曲线上显示最佳阈值点

predict_labeled$status <- (predict_labeled$.pred>0.9813)


# Compare different w (0.9~1) 由于服务器崩坏，使用小数据集在本地跑
set.seed(1)
unl_index <- createDataPartition(y = data_unlabeled$Gendermale, p = 0.01, list = F)
unl <- data_unlabeled[unl_index,]
test_index <- createDataPartition(y = data_labeled$Patient.StatusRecovered, p = 0.001, list = F)
test <- data_labeled[test_index,]

roc_res = data.frame()
for (part in seq(0.0005, 0.002, by=0.0005)) {
  rec_index <- createDataPartition(y = data_recovered$Gendermale, p = part, list = F)
  rec <- data_recovered[rec_index,]
  die_index <- createDataPartition(y = data_died$Gendermale, p = part, list = F)
  die <- data_died[die_index,]
  data_train <- rbind(rec,die,unl)
  j=1 
  thres=vector()
  auc=vector()
  speci=vector()
  sensi=vector()
  for (w in seq(0.9,1,by=0.02)) {
    m <- SSLR::SSLRRandomForest(trees = 50, w=w) %>% fit(Patient.StatusRecovered ~ ., data = data_train)
    predict_labeled <- predict(m, test)
    label_test <- as.data.frame(cbind(test$Patient.StatusRecovered, predict_labeled))
    colnames(label_test) <- c("true", "pred")
    calib_res <- eRic::prCalibrate(r.calib = label_test$true, p.calib = label_test$pred, nbins = 10)
    ROCres <- pROC::roc(calib_res$responses, calib_res$cal.probs, auc = TRUE)
    thres[j] <- pROC::coords(ROCres, "best", ret="threshold", best.method="youden")
    auc[j] <- pROC::auc(ROCres)
    speci[j] <- pROC::coords(ROCres, "best", ret="specificity", best.method="youden")
    sensi[j] <- pROC::coords(ROCres, "best", ret="sensitivity", best.method="youden")
    j=j+1
    cat("finished", part, w)
  }
  roc_res <- rbind(roc_res, cbind(w=seq(0.9:1, by=0.02),thres,auc,speci,sensi))
}
rm(rec_index, rec, die_index, die, test_index, test, thres,auc,speci,sensi,j)
roc_res$p <- as.character(c(rep(0.93,6), rep(0.95,6), rep(0.97,6), rep(0.99,6)))
roc_res$w <- as.numeric(roc_res$w)
roc_res$auc <- as.numeric(roc_res$auc)
ggplot(roc_res, mapping = aes(x=w,y=auc,group=p,color=p)) + geom_point() + geom_line() + 
  theme_bw() +theme(panel.grid.major = element_blank(), panel.border = element_blank(), axis.line = element_line(colour ="black"),plot.title = element_text(hjust = 0.5))


# The SOD mean of all dead samples
mean(as.numeric(data[data$Patient.Status=="Died",]$T))
# = 13.17697


# multiple propensity score matching (PSM) for vaccination
data_filled$Vaccination.Status <- ifelse(data_filled$Vaccination.Status0==1, 0, ifelse(data_filled$Vaccination.Status1==1, 1, ifelse(data_filled$Vaccination.Status2==1, 2, -1)))
data_filled <- data_filled[data_filled$Vaccination.Status!=-1,]

catVars <- c("Gendermale", "Vaccination.Status", "Symptom.FeverYes", "Cardiovascular.disease.including.hypertensionYes", "Is.Home.Quarantine.Yes")
dat_var <- data_filled[data_filled$Vaccination.Status %in% c(0,1),]
dat_var[,catVars] <- lapply(dat_var[,catVars], as.factor)
match_mod <- matchit(Vaccination.Status ~ Gendermale + Symptom.FeverYes + Cardiovascular.disease.including.hypertensionYes + Is.Home.Quarantine.Yes, data = dat_var, method = "nearest")
matchdata01 <- match.data(match_mod)

dat_var <- data_filled[data_filled$Vaccination.Status %in% c(1,2),]
dat_var$Symptom.Fever[which(dat_var$Symptom.Fever!="Yes")] <- "No"
dat_var[,catVars] <- lapply(dat_var[,catVars], as.factor)
match_mod <- matchit(Vaccination.Status ~ Gendermale + Symptom.FeverYes + Cardiovascular.disease.including.hypertensionYes + Is.Home.Quarantine.Yes, data = dat_var, method = "nearest")
matchdata12 <- match.data(match_mod)

matchdata <- rbind(matchdata01, matchdata12[matchdata12$Vaccination.Status==2,])


# Weibull distribution
# MLE
# data_filled <- matchdata[matchdata$Vaccination.Status==0,]
# data_filled <- data_filled[data_filled$Vaccination.Status0==1,]
# data_filled <- data_filled[data_filled$Gendermale==0,]
# data_filled <- as.numeric(data_filled)
# data_filled$.pred[data_filled$Patient.StatusRecovered==1] <- 0.999
# data_filled$.pred[data_filled$Patient.StatusRecovered==0] <- 0.001
# data_filled$censored <- 0
nll_weibull <- function(lambda, k){
  eta <- as.numeric(data_filled$censored)
  t <- as.numeric(data_filled$T)
  B <- as.numeric(data_filled$Patient.StatusRecovered)
  pi <- as.numeric(data_filled$.pred)
  -sum( (1-B)*(log(pi) - lambda*(t)^k + eta*(log(k)+log(lambda)+(k-1)*log(t))) + B*log(1-pi) )
}
# make sure the nll is always finite
data_filled[data_filled$.pred==1, ".pred"] <- 0.999
data_filled[data_filled$.pred==0, ".pred"] <- 0.001
data_filled <- data_filled[data_filled$T!="", ]
data_filled[data_filled$T==0, "T"] <- 0.5
# optimizing
est_weibull <- bbmle::mle2(minuslogl = nll_weibull, start = list(lambda=1, k=1), method = "Nelder-Mead")
est_weibull
# mean SOD
lambda = 0.02245693
k = 1.60654771
lambda^(-1/k)*gamma(1+1/k)
0.07714061^(-1/0.98212835)*gamma(1+1/0.98212835) # mean = 13.68839
# confident interval
bbmle::confint(est_weibull)
0.07286052^(-1/0.96632667)*gamma(1+1/0.96632667) # = 15.26558
0.08154168^(-1/0.99845329)*gamma(1+1/0.99845329) # = 12.31946
bbmle::AIC(est_weibull) # = 419084.7


# Log-logistic distribution
# MLE
nll_log_logis <- function(lambda, k){
  eta <- as.numeric(data_filled$censored)
  t <- as.numeric(data_filled$T)
  B <- as.numeric(data_filled$Patient.StatusRecovered)
  pi <- as.numeric(data_filled$.pred)
  -sum( (1-B)*(log(pi) - (1+eta)*log(1+(lambda*t)^k) + eta*(log(k)+log(lambda)+(k-1)*log(lambda*t))) + B*log(1-pi) )
}
est_log_logis <-bbmle::mle2(minuslogl = nll_log_logis, start = list(lambda=1,k=1), method = "Nelder-Mead")
est_log_logis
shape <- 1.7874553
scale <- 1/0.1295311
mean_llogis(shape = shape, scale = scale)
pi/(0.1108574*1.8664299*sin(pi/1.8664299)) # mean = 15.28
bbmle::confint(est_log_logis)
pi/(0.1083358*1.8256083*sin(pi/1.8256083)) # = 16.06489
pi/(0.1134774*1.9059008*sin(pi/1.9059008)) # = 14.5696
bbmle::AIC(est_log_logis) # = 417884.3


# Log-normal distribution
# MLE
nll_log_norm <- function(mu, sigma){
  eta <- as.numeric(data_filled$censored)
  t <- as.numeric(data_filled$T)
  B <- as.numeric(data_filled$Patient.StatusRecovered)
  pi <- as.numeric(data_filled$.pred)
  Z <- (log(t)-mu)/sigma
  -sum( (1-B)*(log(pi) + (1-eta)*log(1-dnorm(Z)) + eta*(log(pnorm(Z))-log(sigma*t))) + B*log(1-pi) )
}
est_log_norm <- bbmle::mle2(minuslogl = nll_log_norm, start = list(mu=0, sigma=1), method = "Nelder-Mead")
est_log_norm
mu <- -8.039871e-01
sigma <- 2.360352e-34
exp(mu + sigma^2/2)
exp(-7.441506e-01) # mean = mu
bbmle::confint(est_log_norm)
AIC(est_log_norm)


# Gamma distribution
# MLE
nll_gamma <- function(beta, gamma){
  eta <- as.numeric(data_filled$censored)
  t <- as.numeric(data_filled$T)
  B <- as.numeric(data_filled$Patient.StatusRecovered)
  pi <- as.numeric(data_filled$.pred)
  Z <- t/beta
  ga <- expint::gammainc(gamma, t)
  -sum( (1-B)*(log(pi) + (1-eta)*log(ga/gamma(gamma)) + eta*((gamma-1)*log(Z)-Z-log(beta)-log(gamma(gamma))) ) + B*log(1-pi) )
}
est_gamma <- bbmle::mle2(minuslogl = nll_gamma, start = list(beta=1, gamma=1), method = "Nelder-Mead")
est_gamma
gamma <- 1.666221
beta <- 6.812172
gamma*beta
# mean = gamma*beta = 13.20122
bbmle::confint(est_gamma)
9.167379*1.345994 # = 12.33924
9.873567*1.429604 # = 14.11529
bbmle::AIC(est_gamma) # = 424315.8


# AFT Weibull model
data_filled <- data_filled[data_filled$age!="",]
rownames(data_filled) <- c(1:nrow(data_filled))
nll_AFT_1 <- function(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, sigma){
  beta <- c(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17)
  eta <- as.numeric(data_filled$censored)
  t <- as.numeric(data_filled$time)
  B <- as.numeric(data_filled$Patient.StatusRecovered)
  pi <- as.numeric(data_filled$.pred)
  X <- data_filled[,1:17]
  for (i in 1:17) {
  X[[i]] <- as.numeric(X[[i]])
  }
  X <- data.matrix(X)
  eps <- (log(t) - X %*% beta)/sigma
  -sum( (1-B)*(log(pi) - exp(eps) + eta*(eps) ) + B*log(1-pi) )
}
beta0 <- as.matrix(rep(1,17))
est_AFT_1 <- bbmle::mle2(minuslogl = nll_AFT_1, start = list(a1=1,a2=1,a3=1,a4=1,a5=1,a6=1,a7=1,a8=1,a9=1,a10=1,a11=1,a12=1,a13=1,a14=1,a15=1,a16=1,a17=1,sigma=1), method = "BFGS")



# --------- SSMC Cox model --------- #

## Data preparation
cut <- 0.981
# unlabeled
predict_unlabeled <- read.csv("/Users/wuhanyu/Desktop/Research/LBS/Covid-19/Semi-supRF/Semi-supRF-0.9/predict_unlabeled.csv")[,-1]
unlabeled <- cbind(data_unlabeled, predict_unlabeled, time_unlabeled)
unlabeled <- rename(unlabeled, "T" = "time_unlabeled", ".pred" = "predict_unlabeled", "B" = "Patient.StatusRecovered,")
unlabeled$B <- as.numeric(!unlabeled$.pred>cut) # fill up cure information
unlabeled$eta <- 0 # eta = 0 (censored) or 1 (uncensored)
# labeled
predict_labeled <- read.csv("/Users/wuhanyu/Desktop/Research/LBS/Covid-19/Semi-supRF/Semi-supRF-0.9/predict_labeled.csv")[,-1]
labeled <- cbind(data_labeled, predict_labeled, time_labeled) %>% 
  rename("T" = "time_labeled", ".pred" = "predict_labeled", "B" = "Patient.StatusRecovered,")
labeled$B <- as.numeric(!labeled$B)
labeled$eta <- 1
# labeled$.pred[labeled$Patient.StatusRecovered==1] <- 1 # use 1/0 instead of predicted .pred for labeled
# labeled$.pred[labeled$Patient.StatusRecovered==0] <- 0
data_filled <- merge(unlabeled, labeled, all = T)
write.csv(data_filled, "/Users/wuhanyu/Desktop/Research/LBS/Covid-19/Semi-supRF/Semi-supRF-0.9/data_filled", row.names = F, quote = F)

## Covariates used in regression
covar <- ""
for (i in 2:17) {
  covar <- paste(covar, colnames(data_filled)[i])
  if(i < 17) covar <- paste(covar, "+")
}
print(covar)

## Initial beta
library(smcure)
data_filled$T <- as.numeric(data_filled$T)
data_filled$age <- as.numeric(data_filled$age)
data_filled <- na.omit(data_filled)
c <- coxph(Surv(T, eta)~Gendermale + Vaccination.Status0 + Vaccination.Status1 + Vaccination.Status2 + Vaccination.Statuspartially + Vaccination.Statusyes + Symptom.FeverYes + FLUYes + Fever.HistoryYes + Symptom.Sore.ThroatYes + Symptom.CoughYes + HeadacheYes + Cardiovascular.disease.including.hypertensionYes + DiabetesYes + Is.Home.Quarantine.Yes + age, data = data_filled, method="breslow")
beta <- c$coefficients
beta[is.na(beta)] <- 0

## Baseline survival function
s <- smsurv(Time = data_filled$T, 
            Status = data_filled$B,
            X = data_filled[,1:17],
            beta = beta,
            w = data_filled$B,
            model = "ph")$survival
s[which(s==0)] <- 1e-10 # in case that s = 0 and hazard is undefined
data_baseline <- cbind(s, "T"=data_filled$T) %>% as.data.frame() %>% unique() %>% arrange(T)
h <- -diff(log(data_baseline$s))/diff(data_baseline$T)
h[length(h)+1] <- h[length(h)]
data_baseline <- data_baseline %>% cbind("h"=h)
for (i in 1:length(data_filled$T)) {
  data_filled$h[i] <- data_baseline$h[which(data_baseline$T==data_filled$T[i])]
  data_filled$S0[i] <- data_baseline$s[which(data_baseline$T==data_filled$T[i])]
}
data_filled$h[which(data_filled$h<=0)] = 0.0001 # hazard must positive

## Survival function estimation
# likelihood
Likelihood <- function(beta){
  X = as.matrix(data_filled[,1:17])
  y = data_filled$T
  B = data_filled$B
  eta = data_filled$eta
  beta = as.matrix(beta, ncol=1)
  mu = exp(as.matrix(data_filled[,1:17]) %*% beta)
  S0 = data_filled$S0
  h = data_filled$h
  -sum( (B*eta)*log(h*mu) + mu*B*log(S0) )
}
beta_init <- rep(1,17)
# optimize
est <- bbmle::mle2(minuslogl = Likelihood, start = list(beta_init), method = "Nelder-Mead")
est <- optim(par = beta_init, fn = Likelihood, method = "BFGS")
# variance, CI, p-Value of estimated beta
hess <- numDeriv::hessian(f = Likelihood, x = est$par)
vcov <- solve(hess) # variance-covariance matrix
na.index <- which(diag(vcov)<0) # variables with inf var set to 0
vcov[,na.index] <- 0
vcov[na.index,] <- 0
variances <- diag(vcov) # variance
se <- sqrt(variances) # standard error
tvals <- est$par / se # t-statistics
pvals <- 2 * pt(abs(tvals), df = length(data_filled$T) - 2, lower.tail = FALSE) # p-values
pvals <- p.adjust(pvals, method = "BH") # Benjamini-Hochberg (BH) adjust
beta_est <- data.frame(est=est$par, var=variances, lower=est$par + qnorm(0.025) * se, upper=est$par + qnorm(0.975) * se, pValue=pvals, HR=exp(est$par))
rownames(beta_est) <- colnames(data_filled)[1:17]

## visualization
# estimated survival and its variance
for (i in 1:length(data_filled$T)) {
  x = data_filled[i,1:17]
  data_filled$S[i] = data_filled$S0[i]^(exp(beta_est$est %*% t(x)))
  grad = as.matrix(-data_filled$S0[i] * exp(beta_est$est %*% t(x)) * x)
  data_filled$varS[i] = grad %*% vcov %*% t(grad) # delta method
}
for (i in 1:length(data_baseline$s)) {
  data_baseline$S[i] = mean(data_filled$S[data_filled$T==data_baseline$T[i]])
  data_baseline$varS[i] = mean(data_filled$varS[data_filled$T==data_baseline$T[i]])
  data_baseline$lowerS[i] = data_baseline$S[i] + qnorm(0.025) * sqrt(data_baseline$varS[i])
  data_baseline$upperS[i] = data_baseline$S[i] + qnorm(0.975) * sqrt(data_baseline$varS[i])
}
ggplot(data_baseline, mapping = aes(x=T,y=S)) +
  geom_line(color = "#2E9FDF") + 
  geom_ribbon(aes(ymin=lowerS, ymax=upperS), alpha=0.3, fill = "#2E9FDF") + 
  xlim(0,40) + 
  labs(x="Time (days)",y="Survival probability") + 
  theme_bw() + theme(panel.border = element_blank())


# --------- General Cox model --------- #

# data_cox <- data_filled[data_filled$B==1,]
# data_cox$eta <- data_cox$eta+1
data_cox <- dat_filled[!(dat_filled$B==0)&(dat_filled$eta==1),]
fit_cox <- coxph(Surv(T, eta) ~ Gendermale + Vaccination.Status0 + Vaccination.Status1 + Vaccination.Status2 + Vaccination.Statuspartially + Vaccination.Statusyes + Symptom.FeverYes + FLUYes + Fever.HistoryYes + Symptom.Sore.ThroatYes + Symptom.CoughYes + HeadacheYes + Cardiovascular.disease.including.hypertensionYes + DiabetesYes + Is.Home.Quarantine.Yes + age, data = data_cox)
beta_cox <- fit_cox$coefficients
beta_cox[is.na(beta_cox)] <- 0
# predicted survival and CI
fit_cox$coefficients[is.na(fit_cox$coefficients)] <- 0
fit_cox_sum <- summary(fit_cox)
beta_est_cox <- data.frame(est=exp(fit_cox$coefficients), pValue=fit_cox_sum$coefficients[,5])
data_baseline$index <- data_baseline$T %in% fit$time
for (i in 1:length(data_baseline$T)) {
  
  if(data_baseline$index[i]){
    data_baseline$Scox[i] <- mean(fit$surv[fit$time==data_baseline$T[i]])
    data_baseline$lowerScox[i] <- summary(fit, times=data_baseline$T[i])$lower
    data_baseline$upperScox[i] <- summary(fit, times=data_baseline$T[i])$upper
  }
  else{
    j=i 
    k=i
    while (!data_baseline$index[j]) {
      j=j-1
    }
    while (!data_baseline$index[k] && k<length(data_baseline$T)) {
      k=k+1
    }
    if(k<length(data_baseline$T)){
      data_baseline$Scox[i] <- (mean(fit$surv[fit$time==data_baseline$T[j]])*(i-j) + mean(fit$surv[fit$time==data_baseline$T[k]])*(k-i))/(k-j)
    data_baseline$lowerScox[i] <- (summary(fit, times=data_baseline$T[j])$lower*(i-j) + summary(fit, times=data_baseline$T[k])$lower*(k-i))/(k-j)
    data_baseline$upperScox[i] <- (summary(fit, times=data_baseline$T[j])$upper*(i-j) + summary(fit, times=data_baseline$T[k])$upper*(k-i))/(k-j) # 插值法补充t不在样本中的情况（无t时刻死亡）
    }
    else{
      data_baseline$Scox[i] <- (mean(fit$surv[fit$time==data_baseline$T[j]]))*j/i
      data_baseline$lowerScox[i] <- (summary(fit, times=data_baseline$T[j]))$lower*j/i
      data_baseline$upperScox[i] <- (summary(fit, times=data_baseline$T[j]))$upper*j/i
    }
  }
}
vcov_cox <- as.matrix(vcov(fit_cox))
var_cox <- diag(vcov_cox)
ggplot(data_baseline) +
  geom_line( mapping = aes(x=T,y=S), color = "#2E9FDF") + 
  geom_line( mapping = aes(x=T,y=Scox), color = "orange") + 
  geom_ribbon(aes(x=T, y=S, ymin=lowerS, ymax=upperS), alpha=0.3, fill = "#2E9FDF") + 
  geom_ribbon(aes(x=T, y=Scox, ymin=lowerScox, ymax=upperScox), alpha=0.3, fill = "orange") +
  xlim(0,40) + 
  labs(x="Time (days)",y="Survival probability") + 
  theme_bw() + theme(panel.border = element_blank())



fit.p <- signif(as.matrix(fit$coefficients)[,5],2)
fit.OR <- signif(as.matrix(fit$coefficients)[,2],2)
fit.lower <- signif(fit$conf.int[,3],2)
fit.upper <- signif(p.adjust(fit$conf.int[,4], method = "BH"),2)
fit_table <- data.frame(OR = fit.OR, 
                        CI = paste(" (",fit.lower,", ",fit.upper,")",sep=""),
                        p.value = fit.p,
                        stringsAsFactors = F
) %>% rename("95%CI" = "CI", "p-Value*" = "p.value") %>% na.omit()
fit <- survfit(fit_cox, newdata = data_cox)
ggsurvplot(fit, color = "#2E9FDF",pval = TRUE,  conf.int = TRUE, 
           ggtheme = theme_minimal(), 
           data = data_cox, 
           xlim = c(0, 100),
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = theme_cleantable())


# --------- General Mixture Cure model --------- #
# 可观测治愈视作删失
library(smcure)
data_mxc <- data_filled
data_mxc$eta[data_mxc$B==0] = 0
for (i in 1:length(colnames(data_mxc))) {
  data_mxc[,i] <- as.numeric(data_mxc[,i])
}
data_mxc <- na.omit(data_mxc)
fit_mxc <- smcure(Surv(T, eta) ~ Gendermale + Vaccination.Status0 + Vaccination.Status1 + Vaccination.Status2 + Vaccination.Statuspartially + Vaccination.Statusyes + Symptom.FeverYes + FLUYes + Fever.HistoryYes + Symptom.Sore.ThroatYes + Symptom.CoughYes + HeadacheYes + Cardiovascular.disease.including.hypertensionYes + DiabetesYes + Is.Home.Quarantine.Yes + age, 
                  cureform = ~ Gendermale + Vaccination.Status0 + Vaccination.Status1 + Vaccination.Status2 + Vaccination.Statuspartially + Vaccination.Statusyes + Symptom.FeverYes + FLUYes + Fever.HistoryYes + Symptom.Sore.ThroatYes + Symptom.CoughYes + HeadacheYes + Cardiovascular.disease.including.hypertensionYes + DiabetesYes + Is.Home.Quarantine.Yes + age, 
                  data = data_mxc, 
                  model = "ph")
