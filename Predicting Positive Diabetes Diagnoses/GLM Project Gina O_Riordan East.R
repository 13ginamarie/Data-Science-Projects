#load libraries
library(ggplot2)

#load the data
library(faraway)
data(diabetes)
help(diabetes)

#---------------
#Step 1 - Exploratory Data

#quick look at dataset
head(diabetes)
#review for missing data
summary(diabetes)
#remove bs.2s and bp.2d for missing data 
diabetes$bp.2s = NULL
diabetes$bp.2d = NULL
#remove bp.1s, bp.1d, and time.ppn for irrelevance
diabetes$bp.1s = NULL
diabetes$bp.1d = NULL
diabetes$time.ppn = NULL
head(diabetes)
#collapsing glyhb into a binary response, 1 = diabetes
diabetes$diagnosis <- ifelse(diabetes$glyhb > 7, 1, 0)
#adding variables for easy of graphing
diabetes$hdl <- as.numeric(diabetes$hdl)
diabetes$diagchar <- as.character(diabetes$diagnosis)
#remove missing values from the dataset
diabetes <- na.omit(diabetes)
#review
head(diabetes)
summary(diabetes)

##predictor 1 - weight
scplot_weight <- ggplot(data=diabetes,aes(x=diagnosis, y=weight)) +
  geom_point() + geom_jitter() +
  labs(title = "Scatterplot Predictor 1 - Jittered Weight", x = "Diagnosis", y="Weight")
scplot_weight

hist_weight <- ggplot(data = diabetes, aes(x=weight,fill=diagchar)) + 
  geom_histogram(position = "dodge", binwidth = 10) + 
  labs(title = "Histogram Predictor 1 - Weight", x="Weight", y="Count") + 
  scale_fill_discrete(name="Diagnosis", breaks=c(0,1), labels=c("Non-Diabetic", "Diabetic"))
hist_weight

##predictor 2 - Stabilized Glucose
scplot_glu <- ggplot(data=diabetes, aes(x=diagnosis, y=stab.glu)) +
  geom_point() + geom_jitter() +
  labs(title = "Scatterplot Predictor 2 - Jittered Stabilized Glucose", x="Diagnosis", y="Stabilized Glucose")
scplot_glu

hist_glu <- ggplot(data = diabetes, aes(x=stab.glu, fill=diagchar)) + 
  geom_histogram(position = "dodge", binwidth = 10) + 
  labs(title = "Histogram Predictor 2 - Stabilized Glucose", x="Stabilized Glucose", y="Count") + 
  scale_fill_discrete(name="Diagnosis", breaks=c(0,1), labels=c("Non-Diabetic","Diabetic"))
hist_glu


#---------------
#Step 2 - Variable Selection

##build full model
#removed glyhb (basis of response variable)
mod_full <- glm(diagnosis~id+chol+stab.glu+hdl+ratio+location+age+gender+height+weight+frame+waist+hip, 
                family=binomial, data=diabetes)
summary(mod_full)

##Run step() function to perform AIC selection
aic_step <- step(mod_full) 

##Obtain formula for reduced model
formula(aic_step)
#reduced model includes chol, stab.glu, hdl, age, time.ppn

##Subtract number of terms in reduced model from number 
##of terms in full model for # predictors removed
length(attr(terms(formula(mod_full)), "term.labels")) - 
  length(attr(terms(formula(aic_step)), "term.labels"))
#9 predictors removed
mod_red <- glm(diagnosis ~ chol+stab.glu+age+waist, family=binomial, data=diabetes)
summary(mod_red)


#---------------
#Step 3 - Assess Model Fit

##lowess plot for chol
#Obtain y_bar, smoothed responses using the loess function
y_smooth_chol <- predict(loess(diagnosis~chol, data=diabetes))
#Which values of y_bar are less than zero or greater than one?
zero_one_chol <- which(y_smooth_chol>0 & y_smooth_chol<1)
#Make plot of X, y_bar for observations with y_bar in (0,1)
plot(jitter(diabetes$chol)[zero_one_chol],
     log(y_smooth_chol[zero_one_chol]/(1-y_smooth_chol[zero_one_chol])),
     xlab = "Total Cholesterol", ylab = "Diagnosis")
#linear = no spline needed

##lowess plot for stab.glu
#Obtain y_bar, smoothed responses using the loess function
y_smooth_glu <- predict(loess(diagnosis~stab.glu, data=diabetes))
#Which values of y_bar are less than zero or greater than one?
zero_one_glu <- which(y_smooth_glu>0 & y_smooth_glu<1)
#Make plot of X, y_bar for observations with y_bar in (0,1)
plot(jitter(diabetes$stab.glu)[zero_one_glu],
     log(y_smooth_glu[zero_one_glu]/(1-y_smooth_glu[zero_one_glu])),
     xlab = "Stabilized Glucose", ylab = "Diagnosis")
#not linear = spline needed

##lowess plot for age
#Obtain y_bar, smoothed responses using the loess function
y_smooth_age <- predict(loess(diagnosis~age, data=diabetes))
#Which values of y_bar are less than zero or greater than one?
zero_one_age <- which(y_smooth_age>0 & y_smooth_age<1)
#Make plot of X, y_bar for observations with y_bar in (0,1)
plot(jitter(diabetes$age)[zero_one_age],
     log(y_smooth_age[zero_one_age]/(1-y_smooth_age[zero_one_age])),
     xlab = "Age", ylab = "Diagnosis")
#not linear = spline needed

##lowess plot for waist
#Obtain y_bar, smoothed responses using the loess function
y_smooth_waist <- predict(loess(diagnosis~waist, data=diabetes))
#Which values of y_bar are less than zero or greater than one?
zero_one_waist <- which(y_smooth_waist>0 & y_smooth_waist<1)
#Make plot of X, y_bar for observations with y_bar in (0,1)
plot(jitter(diabetes$waist)[zero_one_waist],
     log(y_smooth_waist[zero_one_waist]/(1-y_smooth_waist[zero_one_waist])),
     xlab = "Waist Size", ylab = "Diagnosis")
#linear = no spline needed


##function to create splines, x is predictor, k is vector of knots. Can only do k=3 knots
my.4splines <- function(x,k){
  k1 <- k[1]
  k2 <- k[2]
  k3 <- k[3]
  x.l1 <- NULL
  x.l2 <- NULL
  x.l3 <- NULL
  x.l4 <- NULL
  for (i in 1:length(x)){
    x.l1[i] <- min(x[i],k1)
    x.l2[i] <- max(min(x[i],k2),k1)-k1
    x.l3[i] <- max(min(x[i],k3),k2)-k2
    x.l4[i] <- max(x[i],k3)-k3
  }
  x.s <- cbind(x.l1,x.l2,x.l3,x.l4)
}

#Find knots for stab.glu, age, time.ppn
knots_glu <- quantile(diabetes$stab.glu,c(.10,.50,.90))
knots_age <- quantile(diabetes$age,c(.10,.50,.90))
#Create splines for stab.glu, age, time.ppn
glu_spline <- my.4splines(diabetes$stab.glu,knots_glu) 
age_spline <- my.4splines(diabetes$age,knots_age) 

#remake model with splines
mod_final <- glm(diagnosis ~ chol + glu_spline + age_spline + waist, 
                 family=binomial, data=diabetes)
summary(mod_final)

##Hosmer-Lemeshow (HL) goodness of fit test
library(generalhoslem)
#for binomial logistic regression model
hoslem.test(diabetes$diagnosis, mod_final$fitted.values)
#model fit is adequate (not rejecting the null, high p-value)
AIC(mod_final)

#---------------
#Step 4 - Model Inferences

##p-values for significant predictors - not including splines
summary(mod_final)
##confidence intervals of significant predictors
beta_chol <- mod_final$coefficients[2]
se_chol <- 0.005481
z_stat <- 1.96
chol_upper <- beta_chol + z_stat*se_chol
chol_lower <- beta_chol - z_stat*se_chol
chol_lower
chol_upper

beta_waist <- mod_final$coefficients[11]
se_waist <- 0.03781
waist_upper <- beta_waist + z_stat*se_waist
waist_lower <- beta_waist - z_stat*se_waist
waist_lower
waist_upper

  
#influential observations
n_final <- dim(model.matrix(mod_final))[1]
p_final <- dim(model.matrix(mod_final))[2]
num_df_final <- p_final
den_df_final<- n_final-p_final
F_threshold_final <- qf(0.5,num_df_final, den_df_final)
#compare cooks distance to F threshold
which(cooks.distance(mod_final)>F_threshold_final)
#no influential observations found

#---------------
#Step 5 - Assess Predictive Power
##ROC curve
library(ROCR)
##arguments for pred are model predicted probabilities and response y
pred <- prediction(mod_final$fitted.values, diabetes$diagnosis)
#Now, can just copy and paste rest of code
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(perf, col=rainbow(10))
#Obtain area under the curve (AUC)
auc.tmp <- performance(pred,"auc"); 
auc <- as.numeric(auc.tmp@y.values)
auc
#outstanding model discrimination

#plots to summarize model discrimination
par(mfrow=c(2,2))
plot(mod_final$fitted.values, jitter(diabetes$diagnosis, factor=0.1), 
     main = "Model Fitted Values", xlab = "Fitted Values", ylab = "Jittered Diagnosis")
plot(perf, col=rainbow(10), main = "ROC Curve")
hist(mod_final$fitted.values[diabetes$diagnosis==0], 
     main="Negative Diagnosis (Diagnosis = 0)", xlab = "Non-Diabetic Diagnosis", breaks = 20)
hist(mod_final$fitted.values[diabetes$diagnosis==1], 
     main="Positive Diagnosis (Diagnosis = 1)", xlab = "Diabetic Diagnosis", breaks = 20)
