# Heart Disease Dataset


# This data set dates from 1988 and consists of four databases: Cleveland, Hungary, 
# Switzerland, and Long Beach V. It contains 76 attributes, including the predicted 
# attribute, but all published experiments refer to using a subset of 14 of them. 
# The "target" field refers to the presence of heart disease in the patient. 
# It is integer valued 0 = no disease and 1 = disease.

# This is multivariate type of dataset which means providing or involving a variety of 
# separate mathematical or statistical variables, multivariate numerical data analysis. 
# It is composed of 14 attributes.
# 
# One of the major tasks on this dataset is to predict based on the given attributes of a 
# patient that whether that particular person has a heart disease or not and other is the 
# experimental task to diagnose and find out various insights from this dataset which could 
# help in understanding the problem more.


# There are various factors associated in the process of determining whether a person will 
# have a heart disease or not. In this project, we will do the analysis on some hypothesis and
# will come up with some conclusions for the same.

# Dataset columns:

# age: The person’s age in years

# sex: The person’s sex (1 = male, 0 = female)

# cp: chest pain type
# — Value 0: asymptomatic
# — Value 1: atypical angina
# — Value 2: non-anginal pain
# — Value 3: typical angina

# trestbps: The person’s resting blood pressure (mm Hg on admission to the hospital)

# chol: The person’s cholesterol measurement in mg/dl

# fbs: The person’s fasting blood sugar (> 120 mg/dl, 1 = true; 0 = false)

# restecg: resting electrocardiographic results
#   Value 0: showing probable or definite left ventricular hypertrophy by Estes’ criteria
#   Value 1: normal
#   Value 2: having ST-T wave abnormality (T wave inversions and/or ST elevation or depression of > 0.05 mV)

# thalach: The person’s maximum heart rate achieved

# exang: Exercise induced angina (1 = yes; 0 = no)

# oldpeak: ST depression induced by exercise relative to rest (‘ST’ relates to positions on the ECG plot. See more here)

# slope: the slope of the peak exercise ST segment — 0: downsloping; 1: flat; 2: upsloping
# 0: downsloping; 1: flat; 2: upsloping

# ca: The number of major vessels (0–3)

# thal: A blood disorder called thalassemia Value 0: NULL (dropped from the dataset previously
#   Value 1: fixed defect (no blood flow in some part of the heart)
#   Value 2: normal blood flow
#   Value 3: reversible defect (a blood flow is observed but it is not normal)

# target: Heart disease (1 = no, 0= yes)


library(magrittr) 
library(dplyr) 
library(ggplot2)
library(forcats)
library(rsample)
library(tidyverse)
library(tidymodels)
library(gridExtra)
library(pROC)
library(readr)
library(corrplot)

heart <- read_csv("/Users/sarju/Desktop/MITA Sem 2/MVA/Individual_Project/heart.csv")
heart

dim(heart) # 1025 rows & 14 columns
colnames(heart)
str(heart)
summary(heart)


# str function says that target column is a num type... so will make it a factor of 2 groups
# heart$target <- as.factor(heart$target)

# also converting cp column into factor of 4 groups
# heart$cp <- as.factor(heart$cp)
# 
# heart$sex <- as.character(heart$sex)
str(heart)
# 
# heart[heart$sex==1, 'sex'] <- "male"
# heart[heart$sex==0, 'sex'] <- "female"
# heart$sex


heart_x = subset(heart, select = c(age,cp, trestbps, chol, fbs, restecg, thalach, exang, oldpeak, slope, ca, thal))
heart_x   # heart_x contains all the columns except the target column...

heart_y = heart$target
heart_y  

heart_cm <- colMeans(heart_x)
heart_cm

heart_S <- cov(heart_x)
heart_S

d <- apply(heart_x, MARGIN = 1, function(heart_x)t(heart_x - heart_cm) %*% solve(heart_S) %*% (heart_x - heart_cm)) # taking the distance... (distance formula)
d

# mahalanobis
library(stats)

heart_MD <- mahalanobis(heart_x, heart_cm, heart_S)
heart_MD
heart$pvalues <- pchisq(heart_MD, df=3, lower.tail=FALSE)
heart


# Create a matrix of categorical variables
# heart_cat = subset(heart, select = c(sex,cp, fbs, restecg, exang, slope, ca, thal))
# 
# heart_cat_matrix <- as.matrix(heart_cat)
# heart_cat_matrix
# 
# dim(heart_cat_matrix)
# length(heart$target)


# The warning message suggests that the chi-squared approximation may be incorrect, 
# which can happen when the expected frequency counts are too low. One way to address this 
# is to merge categories with low frequencies into a single category.
# For example, let's say that the categories "0" and "1" have low frequencies in the "thal" column. 
# We can merge them into a single category "0/1" using the following code:

heart$thal
heart$thal <- as.character(heart$thal) # convert to character
heart$thal[heart$thal %in% c("0", "1")] <- "0/1" # merge categories
heart$thal <- factor(heart$thal) # convert back to factor


# Convert categorical variables to factors
heart$sex <- as.factor(heart$sex)
heart$cp <- as.factor(heart$cp)
heart$fbs <- as.factor(heart$fbs)
heart$restecg <- as.factor(heart$restecg)
heart$exang <- as.factor(heart$exang)
heart$slope <- as.factor(heart$slope)
heart$ca <- as.factor(heart$ca)
heart$thal <- as.factor(heart$thal)

# Create contingency tables
sex_table <- table(heart$sex, heart$target)
cp_table <- table(heart$cp, heart$target)
fbs_table <- table(heart$fbs, heart$target)
restecg_table <- table(heart$restecg, heart$target)
exang_table <- table(heart$exang, heart$target)
slope_table <- table(heart$slope, heart$target)
ca_table <- table(heart$ca, heart$target)
thal_table <- table(heart$thal, heart$target)

# Perform chi-squared test on each table
sex_chi <- chisq.test(sex_table)
cp_chi <- chisq.test(cp_table)
fbs_chi <- chisq.test(fbs_table)
restecg_chi <- chisq.test(restecg_table)
exang_chi <- chisq.test(exang_table)
slope_chi <- chisq.test(slope_table)
ca_chi <- chisq.test(ca_table)
thal_chi <- chisq.test(thal_table)


# Print p-values for each variable
print(paste0("Sex p-value: ", sex_chi$p.value))
print(paste0("Chest pain type p-value: ", cp_chi$p.value))
print(paste0("Fasting blood sugar > 120 mg/dl p-value: ", fbs_chi$p.value))
print(paste0("Resting electrocardiographic results p-value: ", restecg_chi$p.value))
print(paste0("Exercise induced angina p-value: ", exang_chi$p.value))
print(paste0("ST segment slope p-value: ", slope_chi$p.value))
print(paste0("Number of major vessels colored by flourosopy p-value: ", ca_chi$p.value))
print(paste0("Thalassemia p-value: ", thal_chi$p.value))


# Create contingency tables for each feature
heart_cat <- heart[, c("sex", "cp", "fbs", "restecg", "exang", "slope", "ca", "thal")]
cont_tables <- lapply(heart_cat, function(x) table(x, heart$target))
cont_tables

# Perform chi-squared tests for each contingency table
chi_results <- lapply(cont_tables, function(x) chisq.test(x))
chi_results

# Extract p-values from each test
p_values <- sapply(chi_results, function(x) x$p.value)
p_values

# Choose a significance threshold
sig_threshold <- 0.05

# Select significant features based on p-values
sig_features <- names(p_values[p_values < sig_threshold])

# Print significant features
print(sig_features)

# So, we got the significant features:
# "sex", "cp", "restecg", "exang", "slope", "ca", "thal"   



# chi square plot
plot(qchisq((1:nrow(heart_x) - 1/2) / nrow(heart_x), df = 3), sort(d),
     xlab = expression(paste(chi[3]^2, " Quantile")),
     ylab = "Ordered distances")
abline(a = 0, b = 1)

# From the chi square plot, we see that the distances appear to be roughly linearly related 
# to the quantiles of the chi-squared distribution, and the points fall close to the diagonal line. 
# This suggests that the Mahalanobis distance metric may be appropriate for this dataset. 
# However, further analysis and interpretation may be needed to draw more definitive conclusions.

# Correlations
pairs(heart_x)

# Now, let us see the distribution of Heart Disease

# First, we will be looking at how many data points do we have where the client has or does not 
# have heart disease. We would like to make sure that the dataset is has a balance of 
# people with and without heart disease in order to create an accurate model.
heart %>%
  count(target)

# As we can see from the result the dataset is also balanced as we have 526 datapoints 
# which dont have heart disease and 499 data points with heart disease. 
# While the dataset might seem balanced on the surface when we look further we will see that 
# the dataset is not balanced based on gender

# 1) Hypothesis 1
# Sex: Studies have shown that men are more likely to develop heart disease than women. 
# In the Heart Disease dataset, we can use a contingency table and chi-square test to 
# investigate the association between sex and the presence of heart disease.

# Which population is diagnosed more with a heart disease? Male Population or Female population?

# The signs of a woman having a heart attack are much less noticeable than the signs of a male. 
# In women, heart attacks may feel uncomfortable squeezing, pressure, fullness, or pain in the 
# center of the chest. It may also cause pain in one or both arms, the back, neck, jaw or stomach, 
# shortness of breath, nausea and other symptoms. Men experience typical symptoms of heart attack, 
# such as chest pain , discomfort, and stress. They may also experience pain in other areas, 
# such as arms, neck , back, and jaw, and shortness of breath, sweating, and discomfort that 
# mimics heartburn.

# source: healthblog.uofmhealth

# The proportion of females and males in the dataset.

heart %>% 
  group_by(sex) %>% 
  summarise(percent = 100 * n() / nrow(heart))


# sex   percent
# <fct>   <dbl>
#  0       30.4
#  1       69.6

# There are 30.4 % females and 69.6% males in the dataset


# Distribution of Heart Disease among males and females
# Create a contingency table of sex and target (heart disease)

Sub_female <- table(heart[heart$sex==0,]$target)
Sub_male <- table(heart[heart$sex==1,]$target)
FM_combine <- rbind(Sub_female,Sub_male)

#Rename columns names and rows names.
colnames(FM_combine) <- c("Has heart disease", "Does not have heart disease")
rownames(FM_combine) <- c("Females", "Males")

#Display the table
FM_combine


# There are 86 females out of 312 who have diagnosed with heart disease and 413 males out of 713 
# were diagnosed with heart disease.
# This indicates that 58% of males in this dataset are diagnosed with heart disease 
# where is only 27% of females are diagnosed with heart disease.


# Perform a chi-square test on the contingency table
chi_test <- chisq.test(FM_combine)
chi_test


# The output from performing Pearson's chi-squared test with Yates' continuity correction 
# on a contingency table:

# Pearson's Chi-squared test with Yates' continuity correction
# data:  cont_table
# X-squared = 78.863, df = 1, p-value < 2.2e-16


# The chi-squared test tests the null hypothesis that there is no association between the two 
# variables, and the alternative hypothesis that there is a significant association.
# The test statistic X-squared is 78.863 with 1 degree of freedom (df) and a p-value less than 
# 2.2e-16, which is smaller than any significance level commonly used (0.05).
# This indicates strong evidence against the null hypothesis and suggests that there is a 
# significant association between the "sex" and "target" variables in the Heart Disease dataset. 
# Specifically, the test suggests that men are more likely to have heart disease than women 
# since the contingency table shows that a higher proportion of men have heart disease compared 
# to women, and so we can infer that men are more likely to have heart disease than women.



# This is a mosaic plot, helps to visualize the statistical association between two variables.
library(mosaic)
mosaicplot(heart$sex ~ heart$target,
           main="Heart disease outcome by Gender", shade=FALSE,color=TRUE,
           xlab="Gender", ylab="Heart disease")


# Ans:  We can conclude that males are diagnosed more with a heart disease than females



# 2) Hypothesis 
# Age: As age increases, the likelihood of developing heart disease also increases. 
# This is supported by previous studies on cardiovascular disease. In the Heart Disease dataset, 
# we can use a scatter plot to visualize the relationship between age and the presence of heart 
# disease, and a correlation test to measure the strength of the relationship.

ggplot(heart, aes(x = age, y = factor(target), color = target)) +
  geom_point() +
  labs(x = "Age", y = "Presence of Heart Disease", color = "Heart Disease")

# heart$target <- as.numeric(heart$target)
cor.test(heart$age, heart$target, method = "pearson")

# true correlation is not equal to 0 means that there is a significant linear relationship 
# between the two variables being tested.

# Pearson's product-moment correlation
# 
# data:  heart$age and heart$target
# t = -7.5356, df = 1023, p-value = 1.068e-13
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.2865322 -0.1704854
# sample estimates:
#        cor 
# -0.2293236 


# If the p-value is less than the chosen significance level (0.05), we can conclude 
# that there is a significant correlation between age and the presence of heart disease.

# The test statistic is the t-value, which in this case is -7.5356. 
# The degrees of freedom (df) are 1023, and the p-value is 1.068e-13, which is very small 
# and indicates that there is a significant correlation between age and target.

# "alternative hypothesis: true correlation is not equal to 0": 
  # This line specifies the alternative hypothesis, which is that there is a non-zero 
  # correlation between age and target.

# "95 percent confidence interval: -0.2865322 -0.1704854": 
  # This line provides the 95% confidence interval for the correlation coefficient. 
  # In this case, it ranges from -0.2865322 to -0.1704854.

# "sample estimates: cor -0.2293236": 
  # This line provides the sample estimate of the correlation coefficient, which is -0.2293236. 
  # This indicates a moderate negative correlation between age and target.

# Based on the output of the Pearson's correlation test, we can conclude that there is a 
# statistically significant negative correlation between age and the presence of heart disease 
# (target variable) in the Heart Disease dataset.
# The correlation coefficient is -0.2293, which indicates a moderate negative correlation 
# between the two variables. As age increases, the likelihood of having heart disease decreases. 
# However, it's important to note that correlation does not necessarily imply causation, 
# and further analysis would be needed to establish causality.

# Let us create a A faceted histogram which shows the relationship between heart disease and age 

age.plot <- ggplot(heart, mapping = aes(x = age, fill = target)) +
  geom_histogram() +
  facet_wrap(vars(target)) +
  labs(title = "Prevelance of Heart Disease Across Age", x = "Age (years)", y = "Count", fill = "Heart Disease")

age.plot


# The histograms for age faceted on the presence and absence of heart disease have 
# different distribution shapes, suggesting that age does have a relationship with heart disease. 
# The distribution of the presence of heart disease is left skewed while the distribution 
# of the absence of heart disease appears more normally distributed. 
# These graphics suggests that there are more older people with heart disease than younger 
# people with heart disease.

# This is a boxplot to displays the age distribution of heart diagnosis.
boxplot(heart$age ~ heart$target,
        main="Heart disease diagnosis distribution by Age",
        ylab="Age",xlab="Heart disease diagnosed")

# Ans: We can conclude that people with a higher age are more likely diagnised with a heart disease



# 3) Hypothesis 3
# Chest pain type: The type of chest pain a patient experiences can provide important 
# diagnostic information about heart disease. 
# cp: chest pain type
# — Value 0: asymptomatic
# — Value 1: atypical angina
# — Value 2: non-anginal pain
# — Value 3: typical angina


#The proportion of patients with chest pain types.

heart %>% group_by(cp) %>% 
  summarise( percent = 100 * n() / nrow( heart ))

#  cp    percent
# <fct>   <dbl>
#  0       48.5 
#  1       16.3 
#  2       27.7 
#  3       7.51

# There are 7% of patient with typical angina chest pain, 16 % of patients with atypical angina, 
# 27% with non-angina pain, and 48% with asymptomatic chest pain.


# In the Heart Disease dataset, we can use a box plot to visualize the distribution of chest 
# pain type among patients with and without heart disease, and a t-test to compare the 
# means of the two groups.

# create a box plot to visualize the distribution of chest pain type

ggplot(heart, aes(x = target, y = factor(cp))) + 
  geom_boxplot() +
  xlab("Presence of Heart Disease") +
  ylab("Chest Pain Type")

# The boxplot suggests that chest pain type may be a significant predictor of heart disease in 
# patients.


# This plot to visualize the Heart disease diagnosis Distributions by Chest pain. 
# There are four types of chest pain:(0) asymptomatic, (1) atypical angina, (2) non-anginal pain, 
# and (3) typical angina.

ggplot(data = heart, aes(x = target, fill = cp)) + 
  geom_bar(position = "fill") +
  labs(title = "Heart disease diagnosis Distributions by Chest pain",
       x = "Heart disease diagnosis",
       y = "chest pain") +
  theme_test()


str(heart)
heart$target <- as.numeric(heart$target)
heart$cp <- as.numeric(heart$cp)

# perform a t-test to compare the means of chest pain type between patients with and without heart disease
t.test(cp ~ target, data = heart)

# Welch Two Sample t-test
# 
# data:  cp by target
# t = -15.462, df = 1022.9, p-value < 2.2e-16
# alternative hypothesis: true difference in means between group 1 and group 2 is not equal to 0
# 95 percent confidence interval:
#   -1.0089917 -0.7817304
# sample estimates:
#   mean in group 1 mean in group 2 
# 1.482966        2.378327 

# This is the output of a Welch two-sample t-test which is used to test if there is a significant 
# difference between the means of two groups. In this case, the groups are defined by the presence 
# or absence of heart disease (target variable) and the chest pain type (cp variable).

# The output shows the t-value, degrees of freedom (df), and the p-value. 
# The null hypothesis is that the means of the two groups are equal, 
# and the alternative hypothesis is that they are not equal. 
# A p-value less than the significance level (typically 0.05) indicates that there is 
# evidence to reject the null hypothesis.
# In this case, the p-value is extremely small (less than 2.2e-16), 
# so we can conclude that there is a significant difference between the means of the two groups. 
# The 95% confidence interval shows the range within which we can be 95% confident that the 
# true difference between the means of the two groups lies. 
# The sample estimates show the mean value of chest pain type (cp variable) for the two groups, 
# with group 0 indicating those without heart disease and group 1 indicating those with 
# heart disease.


# A a bar chart is shown for the relationship between heart disease and chest pain type. 

cp.plot <- ggplot(heart, mapping = aes(x=target, fill = factor(cp))) +
  geom_bar(position = "dodge") +
  labs(title = "Prevelance of Heart Disease for Different Chest Pain Types", x = "Heart Disease", y = "Count", fill = "Chest Pain Type")

cp.plot

# There does appear to be a relationship between type of chest pain and heart disease. 
# Interestingly, asymptomatic chest pain type (Value 0) has the highest count for the presence of 
# heart disease, while typical angina pain has the lowest count. There is a higher count 
# of people without heart disease that have atypical or typical angina chest pain compared 
# to people with heart disease. Angina is listed as one of the most common symptoms of 
# heart attack and so this result is skeptical and needs further investigation, but we will 
# assume it is correct for the current analysis.



# Ans: Therefore, chest pain type is a significant predictor of heart disease in patients.




attach(heart)


# Visualizations of each variable are created in order to familiarize ourselves with 
# the predictors. 
# Proportional bar charts are created for our 3 binary variables: sex, fasting blood sugar, 
# and exercise induced angina. Box plots are shown for our 3 numerical variables: 
#   resting blood pressure, serum cholesterol, and maximum heart rate.


heart$target <- as.factor(heart$target)

# sex: The person’s sex (1 = male, 0 = female)
sex.plot <- ggplot(heart, mapping = aes(x = sex, fill = target)) +
  geom_bar(position = "fill") +
  labs(x = "Sex", y = "Proportion", fill = "Heart Disease") +
  theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12))

# fbs: The person’s fasting blood sugar (> 120 mg/dl, 1 = true; 0 = false)
fbs.plot <- ggplot(heart, mapping = aes(x=fbs, fill=target)) +
  geom_bar(position = "fill") +
  labs(x = "Fasting Blood Sugar", y = "Proportion", fill = "Heart Disease") +
  scale_x_discrete(labels = c("low", "high"))+
  theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12))

# exang: Exercise induced angina (1 = yes; 0 = no)
exang.plot <- ggplot(heart, mapping = aes(x = exang, fill = target)) +
  geom_bar(position = "fill") +
  labs(x = "Exercise induced angina", y = "Proportion", fill = "Heart Disease") +
  theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 12))


grid.arrange(sex.plot, fbs.plot, exang.plot, nrow=2)

# The bar plot on the top left shows a higher proportion of males with heart disease than 
# females with heart disease, suggesting that there is a relationship between sex and heart disease. 
# There is an even larger distinction for heart disease as it relates to exercise induced angina, 
# for a much higher proportion of people with exercise induced angina have heart disease compared 
# to people without exercise induced angina (bottom left plot). 
# Fasting blood sugar levels do not appear to have a correlation with heart disease, as there 
# appears to be a similar proportion of presence and absence of heart disease for people with 
# high and low fasting blood sugar levels (top right plot).

# More visualizations

# trestbps: The person’s resting blood pressure (mm Hg on admission to the hospital)
trestbps.plot <- ggplot(heart, mapping = aes(x=trestbps, y=target)) +
  geom_boxplot() +
  labs(x = "Resting Blood Pressure (mm Hg)", y = "Heart Disease") +
  theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12))

# chol: The person’s cholesterol measurement in mg/dl
chol.plot <- ggplot(heart, mapping = aes(x=chol, y=target)) +
  geom_boxplot() +
  labs(x = "Serum Cholestoral (mg/dl)", y = "Heart Disease") +
  theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12))

# thalach: The person’s maximum heart rate achieved
maxhr.plot <- ggplot(heart, mapping = aes(x = thalach, y = target)) +
  geom_boxplot() +
  labs(x = "Maximum Heart Rate (bpm)", y = "Heart Disease") +
  theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12))

grid.arrange(trestbps.plot, chol.plot, maxhr.plot, nrow=2)

# The box plots for resting blood pressure (top left) appear similar for the presence and 
# absence of heart disease; the boxes (middle 50% of data) overlap and have similar centers 
# and boundaries. These similar distributions suggest that there is likely not a strong 
# relationship between resting blood pressure and heart disease. 
# No strong relationship between serum cholesterol and heart disease is visible either, 
# for these distributions are close in center and spread (top right). 
# However, the last box plot does show a relationship between maximum heart rate achieved 
# during exercise and heart disease. The trend of maximum heart rate appears to be higher for 
# the absence of heart disease compared to the presence of heart disease, and there is 
# little overlap between the two boxes.

# And, to visualize the relationship between restecg and target, we can create a stacked bar plot. 
# This plot shows the count of each value of restecg for each value of target.

# restecg: resting electrocardiographic results
#   Value 0: showing probable or definite left ventricular hypertrophy by Estes’ criteria
#   Value 1: normal
#   Value 2: having ST-T wave abnormality (T wave inversions and/or ST elevation or depression of > 0.05 mV)
ggplot(heart, aes(x = target, fill = factor(restecg))) +
  geom_bar()

# From the above plot, we can infer that restecg value 0 and value 2 is found more in common with 
# people having a heart disease. 
# So, restecg is also a important feature.


# From the above observations, we got the significant features which is correct:
# "sex", "cp", "restecg", "exang", "slope", "ca", "thal"   



# More useful visualizations 


# Another plot to visualize heart disease diagnosis Distributions by Number of major vessels.
# ca: The number of major vessels (0–3)

ggplot(data = heart, aes(x = target, fill = ca)) + 
  geom_bar(position = "fill") +
  labs(title = "Heart disease diagnosis Distributions by Number of major vessels ",
       x = "Heart disease diagnosis",
       y = "thal") +
  theme_test()



# Here is the scatterplot showing the relationship between age and maximum heart rate achieved 
# in the heart dataset:
  
ggplot(heart, aes(x=age, y=thalach, color=target)) + 
  geom_point(alpha=0.5) +
  labs(x="Age", y="Max heart rate", color="Target") +
  theme_minimal()

# The downward trend indicates that as age increases, maximum heart rate achieved tends to decrease. 
# This relationship appears to be somewhat linear, although there is some variability in the data.


# Histogram of patient’s age and gender

meanAge <- data.frame(sex = c(0, 1), age = c(mean(heart[heart$sex==0,]$age),mean(heart[heart$sex==1,]$age)))

#ggplot of age of the patients categorized by sex
Plot <- ggplot(heart, aes(x=age, fill=as.factor(sex))) +
  geom_histogram(alpha=0.5, position="identity")+ 
  geom_vline(aes(xintercept = age), meanAge)+
  facet_wrap(~as.factor(sex))+
  labs(title="Histogram of patients's age by gender", 
       x="Age of patients", y="Count", fill="Sex")+
  geom_text(meanAge, mapping=aes(x=age, y=8.5, label=paste("Mean=", signif(age,4))),
            size=4, angle=90, vjust=-0.4, hjust=0)+
  scale_fill_discrete(breaks=c("0", "1"),
                      labels=c("0 - Female", "1 - Male"))

#display the Plot
Plot


# Heart disease diagnosis frequency by Resting electrocardiographic results and sex
heart$restecg <- as.factor(heart$restecg)
heart %>%
  ggplot(aes(x = target, fill=restecg)) + 
  geom_bar(position = "dodge") +
  facet_grid(~sex) +
  scale_fill_brewer(palette = "Dark2") +
  labs(title="Heart disease diagnosis frequency by restecg and sex")



heart$age <- as.numeric(heart$age)
heart$sex <- as.numeric(heart$sex)
heart$cp <- as.numeric(heart$cp)
heart$trestbps <- as.numeric(heart$trestbps)
heart$chol <- as.numeric(heart$chol)
heart$fbs <- as.numeric(heart$fbs)
heart$restecg <- as.numeric(heart$restecg)
heart$thalach <- as.numeric(heart$thalach)
heart$exang <- as.numeric(heart$exang)
heart$oldpeak <- as.numeric(heart$oldpeak)
heart$slope <- as.numeric(heart$slope)
heart$ca <- as.numeric(heart$ca)
heart$thal <- as.numeric(heart$thal)

correlations <- cor(heart[,1:13])
corrplot(correlations, method="circle")


# A dot-representation was used where blue represents positive correlation and red negative. 
# The larger the dot the larger the correlation.

# -> As age increases, cholesterol, stress of heart during exercise and resting bp also increase. 
#    On the other hand, maximum heart rate falls with old age.
# -> As cholesterol increases, stress of heart during exercise and resting bp increase, 
     # while maximum heart rate falls.
# -> As ST depression rises, i.e. stress of the heart falls, resting bp rises.
# -> Resting bp also has a negative relation with maximum heart rate.
# -> The degree of correlation is very small between all variables. However, age and maximum 
     # heart rate show a slightly higher correlation.
# -> St depression and maximum heart rate also show similar results.




# Feature Importance
# There are different ways to identify the important features in the data.
# 
# 1- Correlation
# 
# 2- Random Forest: Gini Importance or Mean Decrease in Impurity (MDI) calculates each feature 
# importance as the sum over the number of splits (across all tress) that include the feature, 
# proportionally to the number of samples it splits.

corelations = data.frame(cor(heart[,1:13], use = "complete.obs"))
corelations

# Splitting the data for females and males in order to find the most important factor 
# leading to heart disease in each gender.


#Create a subset for males only.
Male_Data <- subset(heart, sex==1)
Male_Data
#Create another subset for females only.
Female_Data <- subset(heart, sex != 1)
Female_Data

# converting the target column into a factor of 2 groups
Male_Data$target <- as.factor(Male_Data$target)
Female_Data$target <- as.factor(Female_Data$target)

library(randomForest)
library(caret)
#Feature selection using random forest technique
Feature_Importance_Males = randomForest(target~., data=Male_Data)
# Create an importance based on mean decreasing gini
importance(Feature_Importance_Males)

varImp(Feature_Importance_Males)

varImpPlot(Feature_Importance_Males, col= "red", pch= 20)

Feature_Importance_Females = randomForest(target~., data=Female_Data)
# Create an importance based on mean decreasing gini
importance(Feature_Importance_Females)

varImp(Feature_Importance_Females)

varImpPlot(Feature_Importance_Females, col= "red", pch= 20)



# 4) Hypothesis

# H0 = There is no association between chest pain and heart disease diagnosis

# HA = There is association between chest pain and heart disease diagnosis

qqnorm(heart$age)
qqline(heart$age)


# Scatterplot
library(ggplot2)
gg <- ggplot(heart, aes(x=chol, y=trestbps)) + 
  geom_point(aes(col=target, size=oldpeak)) + 
  geom_smooth(method="loess", se=F) + 
  xlim(c(100, 430)) + 
  ylim(c(75, 200)) + 
  labs(subtitle="trestbps Vs chol", 
       y="trestbps", 
       x="chol", 
       title="Scatterplot", 
       caption = "Source: midwest", 
       bins = 30)
plot(gg)

# We will try logistic regression to predict which patients have heart disease.

# Logistic Regression Prediction


library(caTools)
library(MASS)
library(ROCR)
library(Information)
library(OptimalCutpoints)
library(pROC)
library(caret)
library(ggplot2)
library(lattice)

#Split the Data to training and testing data to conduct a logistic regression model
set.seed(123)
split=sample.split(heart$target, SplitRatio = 0.75)
Train_Data=subset(heart,split == TRUE)
Test_Data=subset(heart,split == FALSE)


# Fit a logistic regression model: We will fit a logistic regression model using the "glm" 
# function in R.

Log_model <- glm(target ~., data=Train_Data, family = "binomial")
summary(Log_model)


# Model evaluation: We will evaluate the performance of the model using various metrics such as 
# accuracy, sensitivity, specificity, and ROC curve.


# Make predictions on test data
predictions <- predict(Log_model, newdata = Test_Data, type = "response")

# Measure the accuracy of prediction in the test data
# The common practice is to take the probability cutoff as 0.5.
# If the probability of Y is > 0.5, then it can be classified an event (presence of heart disease).
# So if pred is greater than 0.5, it is positive(heart disease =yes) else it is negative

y_pred_num <- ifelse(predictions > 0.5, 1, 0)
y_pred <- factor(y_pred_num, levels=c(0, 1))
y_act <- Test_Data$target

# Result : Prediction Accuracy (Proportion of predicted target that matches with actual target)
mean(y_pred == y_act)
# 0.8754864

# Calculate the optimal threshold
pred <- prediction(predictions, Test_Data$target)
perf <- performance(pred, "tpr", "fpr")
thresholds <- perf@alpha.values[[1]]

# Plot the ROC curve
plot(perf)

# Add the optimal threshold point to the ROC curve
abline(a = 0, b = 1)
points(thresholds, thresholds, col = "red", pch = 19)


# Get p-values for all features
summary(Log_model)$coef[, "Pr(>|z|)"]

# Suppose you have a vector of p-values in scientific notation
p_values <- c(601880e-03, 6.525117e-01, 1.709865e-09, 3.400166e-13, 2.303392e-03, 3.193722e-02, 
              2.370040e-01, 4.173028e-02, 4.954128e-05, 1.587542e-04, 1.448314e-05, 8.981884e-02, 
              4.013419e-10, 5.346086e-08 )

# Use format() to convert to normal form
formatted_p_values <- format(p_values, scientific = FALSE)

# View the result
formatted_p_values


# The p-value obtained for the coefficient of the chest pain variable (cp) in the logistic 
# regression model indicates the statistical significance of the association between 
# chest pain and heart disease diagnosis. 
# The p-value for cp is given as "3.40e-13" which is very small (less than 0.05), 
# indicating strong evidence against the null hypothesis that there is no association 
# between chest pain and heart disease diagnosis.

# Therefore, based on the p-value obtained, we can reject the null hypothesis and conclude 
# that there is a significant association between chest pain and heart disease diagnosis.


predictTrain = predict(Log_model, type='response')
#Confusion matrix using threshold of 0.5
table(Train_Data$target, predictTrain>0.5)


#    FALSE TRUE
# 0   301   73
# 1    35  359


# The table shows the number of instances that fall into each combination of actual and predicted classes.
# In this case, there were 301 instances where the actual target was 0 and the model 
# predicted 0 (true negative), 73 instances where the actual target was 0 and the 
# model predicted 1 (false positive), 35 instances where the actual target was 1 and the 
# model predicted 0 (false negative), and 359 instances where the actual target was 1 and 
# the model predicted 1 (true positive).

# True Negative (TN): 301 - The number of instances where the model correctly predicted the 
# absence of heart disease (actual target = 0)
# False Positive (FP): 73 - The number of instances where the model predicted the presence 
# of heart disease (actual target = 0)
# False Negative (FN): 35 - The number of instances where the model predicted the absence 
# of heart disease (actual target = 1)
# True Positive (TP): 359 - The number of instances where the model correctly predicted the 
# presence of heart disease (actual target = 1)


#Calculate the accuracy on the training set
(301+359)/nrow(Train_Data)  # 0.859375

#Predictions on Test set
predictTest = predict(Log_model, newdata=Test_Data, type='response')
#Confusion matrix using threshold of 0.5
table(Test_Data$target, predictTest>0.5)

# FALSE TRUE
# 0   105   20
# 1    12  120


# Anova test
anova(Log_model, test="Chisq")

# ANOVA (Analysis of Variance) test is used to determine whether there are significant differences 
# between the means of two or more groups. In the context of logistic regression, an ANOVA test 
# is used to determine if there is a significant difference in the deviance between two or more 
# nested models.
# 
# In the output, we can see the ANOVA table for the logistic regression model Log_model. 
# The table shows the change in deviance for each variable added to the model, and the 
# p-value associated with that change in deviance. The null hypothesis is that the model 
# without the variable is not significantly different from the model with the variable. 
# The alternative hypothesis is that the model with the variable is significantly better 
# than the model without the variable.
# 
# In this case, we can see that all the variables are statistically significant except 
# for slope, restecg, and fbs (which has a p-value of 0.65). This suggests that the 
# model is a good fit for the data and that all the variables are important in predicting 
# the target variable.


# t-test
# t-test is statistical method used to determine the significant difference between the means of two groups.

# t-test to confirm the association between chest pain and heart disease

ttest_cp <- t.test(heart$cp ~ heart$target, var.equal= TRUE)  
ttest_cp


# The output shows the test statistic t, which is -15.445, the degrees of freedom (df) which 
# is 1023, and the p-value which is less than 2.2e-16, indicating strong evidence against 
# the null hypothesis. The null hypothesis is that there is no significant difference in 
# the means of chest pain between individuals with and without heart disease. 
# The alternative hypothesis is that there is a significant difference in the means of 
# chest pain between the two groups.
# 
# The 95% confidence interval is also shown, which is -1.009114 to -0.781608. The sample 
# means for each group are given as well, which are 1.482966 for group 0 (no heart disease) 
# and 2.378327 for group 1 (heart disease).
# 
# Therefore, the results of the t-test indicate a significant difference in the means of 
# chest pain between individuals with and without heart disease. Individuals with heart 
# disease tend to have a higher mean chest pain score than those without heart disease.

# Chi-test
# let's perform chi-squared test of independence to investigate the association between two 
# categorical variables, chest pain (heart$cp) and heart disease diagnosis (heart$target).
CHI_cp <- chisq.test(heart$cp, heart$target) 

# Print the results to see if p<0.05.
print(CHI_cp)

# The p-value is less than 2.2e-16, which indicates strong evidence against the null hypothesis 
# of independence. Therefore, we can conclude that there is a significant association between 
# chest pain and heart disease diagnosis.



# Calculate the confusion matrix
conf_matrix <- confusionMatrix(y_pred, y_act)
print(conf_matrix)

# Calculate the accuracy
accuracy <- conf_matrix$overall[["Accuracy"]]
print(paste0("Accuracy: ", accuracy))

# Compute metrics
sensitivity <- conf_matrix$byClass["Sensitivity"]
print(sensitivity) #  0.84 


specificity <- conf_matrix$byClass["Specificity"]
print(specificity) # 0.9090909 

# ROC curve
library(pROC)
roc_data <- roc(Test_Data$target, predictions)
plot(roc_data, col = "blue", print.auc = TRUE, auc.polygon = TRUE, max.auc.polygon = TRUE,
     grid=c(0.1,0.2), grid.col=c("green", "yellow"), lwd=2, main="ROC Curve")
legend("bottomright", legend = paste("AUC =", round(auc(roc_data),2)), col = "blue", lty=1, cex=1.5)


# The logistic regression model has a sensitivity of 0.84, which means that the model 
# correctly identifies 84% of the people who have heart disease. 

# The specificity of the model is 0.909, which means that the model correctly identifies 
# 90.9% of the people who do not have heart disease.

# The ROC curve plot shows the tradeoff between sensitivity and specificity for 
# different probability thresholds of the model. 
# The closer the curve is to the top left corner, the better the model performance. 
# The AUC (Area Under the Curve) value is a measure of the model's overall performance, 
# where an AUC of 1 represents a perfect model and an AUC of 0.5 represents a random classifier. 
# In this case, the AUC is 0.91, which indicates a good performance of the model.

# Overall, the logistic regression model seems to perform well in predicting the 
# presence of heart disease in the test data.




# 5) Hypothesis 5
# Is there a significant difference in the mean cholesterol levels between individuals 
# with and without heart disease?

# chol: The person’s cholesterol measurement in mg/dl

# Two-sample t-test to compare mean cholesterol levels between individuals with and without heart disease
ttest_chol <- t.test(heart$chol ~ heart$target, var.equal=TRUE)

# Print the results
print(ttest_chol)

# p-value = 0.001353

# Since this value is less than 0.05, we can reject the null hypothesis and conclude that there 
# is a significant difference in the mean cholesterol levels between individuals with and 
# without heart disease.


# Factor Analysis


library(cluster)
library(factoextra)
library(magrittr)
library(NbClust)
library(data.table)
library(dplyr)
library(psych)

heart_data <- heart
heart_data

colnames(heart_data)
str(heart_data)


attach(heart_data)
#heart_data[1:6]
heart_data

heart_data$target <- as.numeric(heart_data$target)

fit.pc <- principal(heart_data, nfactors=4, rotate="varimax")
fit.pc

round(fit.pc$values, 3)
fit.pc$loadings

# Loadings with more digits
for (i in c(1,3,2,4)) { print(fit.pc$loadings[[1,i]])}

# Communalities
fit.pc$communality

# Rotated factor scores, Notice the columns ordering: RC1, RC3, RC2 and RC4
# fit.pc$scores
head(fit.pc$scores)

# Play with FA utilities
fa.parallel(heart_data) # See factor recommendation
fa.plot(fit.pc) # See Correlations within Factors
fa.diagram(fit.pc) # Visualize the relationship
vss(heart_data) # See Factor recommendations for a simple structure


# Computing Correlation Matrix
corrm.emp <- cor(heart_data)
corrm.emp
plot(corrm.emp)

heart_data_pca <- prcomp(heart_data, scale=TRUE)
summary(heart_data_pca)
plot(heart_data_pca)

# A table containing eigenvalues and %'s accounted, follows. Eigenvalues are the sdev^2
(eigen_heart_data <- round(heart_data_pca$sdev^2,3))
round(fit.pc$values, 3)
names(eigen_heart_data) <- paste("PC",1:14,sep="")
eigen_heart_data

sumlambdas <- sum(eigen_heart_data)
sumlambdas

propvar <- round(eigen_heart_data/sumlambdas,2)
propvar

cumvar_heart_data <- cumsum(propvar)
cumvar_heart_data

matlambdas <- rbind(eigen_heart_data,propvar,cumvar_heart_data)
matlambdas

rownames(matlambdas) <- c("Eigenvalues","Prop. variance","Cum. prop. variance")
rownames(matlambdas)

eigvec.emp <- heart_data_pca$rotation
print(heart_data_pca)

# Taking the first four PCs to generate linear combinations for all the variables with four factors
pcafactors.emp <- eigvec.emp[,1:4]
pcafactors.emp

# Multiplying each column of the eigenvector’s matrix by the square-root of the corresponding eigenvalue in order to get the factor loadings
unrot.fact.emp <- sweep(pcafactors.emp,MARGIN=2,heart_data_pca$sdev[1:4],`*`)
unrot.fact.emp

# Computing communalities
communalities.emp <- rowSums(unrot.fact.emp^2)
communalities.emp

# Performing the varimax rotation. The default in the varimax function is norm=TRUE thus, Kaiser normalization is carried out
rot.fact.emp <- varimax(unrot.fact.emp)

#View(unrot.fact.emp)
rot.fact.emp

# The print method of varimax omits loadings less than abs(0.1). In order to display all the loadings, it is necessary to ask explicitly the contents of the object $loadings
fact.load.emp <- rot.fact.emp$loadings[1:14,1:4]
fact.load.emp

# Computing the rotated factor scores for the 30 European Countries. Notice that signs are reversed for factors F2 (PC2), F3 (PC3) and F4 (PC4)
scale.emp <- scale(heart_data[1:14])
# scale.emp
head(scale.emp)
head(as.matrix(scale.emp)%*%fact.load.emp%*%solve(t(fact.load.emp)%*%fact.load.emp))




# Different way of PCA Method


library(factoextra)
library(FactoMineR)
library(ggfortify)
library(psych)
library(corrplot)
library(devtools)

res.pca <- PCA(heart[,1:13], graph = FALSE)
print(res.pca)

# Visualize and Interpret PCA using these functions 

#get_eigenvalue(res.pca): Extract the eigenvalues/variances of principal components
#fviz_eig(res.pca): Visualize the eigenvalues
#get_pca_ind(res.pca), get_pca_var(res.pca): Extract the results for individuals and variables, respectively.
#fviz_pca_ind(res.pca), fviz_pca_var(res.pca): Visualize the results individuals and variables, respectively.
#fviz_pca_biplot(res.pca): Make a biplot of individuals and variables.

eig.val <- get_eigenvalue(res.pca)
eig.val

fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

var <- get_pca_var(res.pca)
#var$coord: coordinates of variables to create a scatter plot
#var$cos2: represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
#var$contrib: contains the contributions (in percentage) of the variables to the principal components. 
#The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
var

# Coordinates
head(var$coord)
# Cos2: quality on the factore map
head(var$cos2)
# Contributions to the principal components
head(var$contrib)

#The plot Below is also known as variable correlation plots. It shows the relationships between all variables. It can be interpreted as follow:

#Positively correlated variables are grouped together.
#Negatively correlated variables are positioned on opposite sides of the plot origin (opposed quadrants).
#The distance between variables and the origin measures the quality of the variables on the factor map. 
#Variables that are away from the origin are well represented on the factor map.

# Correlation circle
fviz_pca_var(res.pca, col.var = "black")

# Quality of representation


corrplot(var$cos2, is.corr=FALSE)
# Total cos2 of variables on Dim.1 and Dim.2
#A high cos2 indicates a good representation of the variable on the principal component. 
#In this case the variable is positioned close to the circumference of the correlation circle.
#A low cos2 indicates that the variable is not perfectly represented by the PCs. 
#In this case the variable is close to the center of the circle.

fviz_cos2(res.pca, choice = "var", axes = 1:2)
fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)
# Change the transparency by cos2 values
fviz_pca_var(res.pca, alpha.var = "cos2")
corrplot(var$contrib, is.corr=FALSE)
# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)
fviz_pca_var(res.pca, alpha.var = "contrib")

fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = heart$target, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)


# Description of PC

res.desc <- dimdesc(res.pca, axes = c(1,2,3,4,5), proba = 0.05)
# Description of dimension 1
res.desc$Dim.1
res.desc$Dim.2
res.desc$Dim.3
res.desc$Dim.4
res.desc$Dim.5

# Graph of Indiviuals
ind <- get_pca_ind(res.pca)
ind

## Principal Component Analysis Results for individuals
##  ===================================================
##   Name       Description                       
## 1 "$coord"   "Coordinates for the individuals" 
## 2 "$cos2"    "Cos2 for the individuals"        
## 3 "$contrib" "contributions of the individuals"
#To get access to the different components, use this:

# Coordinates of individuals
head(ind$coord)
# Quality of individuals
head(ind$cos2)
# Contributions of individuals
head(ind$contrib)

fviz_pca_ind(res.pca)

fviz_pca_ind(res.pca, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)

fviz_pca_ind(res.pca, pointsize = "cos2", 
             pointshape = 21, fill = "#E7B800",
             repel = TRUE # Avoid text overlapping (slow if many points)
)

fviz_pca_ind(res.pca, col.ind = "cos2", pointsize = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)

fviz_cos2(res.pca, choice = "ind")

# Total contribution on PC1 and PC2
fviz_contrib(res.pca, choice = "ind", axes = 1:2)

# Create a random continuous variable of length 23,
# Same length as the number of active individuals in the PCA
set.seed(123)
my.cont.var <- rnorm(1025)
# Color individuals by the continuous variable
fviz_pca_ind(res.pca, col.ind = my.cont.var,
             gradient.cols = c("blue", "yellow", "red"),
             legend.title = "Cont.Var")

heart$target <- as.factor(heart$target)
fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = heart$target, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)

fviz_pca_ind(res.pca, geom.ind = "point", col.ind = heart$target, 
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, ellipse.type = "confidence",
             legend.title = "Groups"
)
fviz_pca_ind(res.pca,
             label = "none", # hide individual labels
             habillage = heart$target, # color by groups
             addEllipses = TRUE, # Concentration ellipses
             palette = "jco"
)
fviz_pca_var(res.pca, geom.var = c("point", "text"))
# Show individuals text labels only
fviz_pca_ind(res.pca, geom.ind =  "text")
# Change the size of arrows an labels
fviz_pca_var(res.pca, arrowsize = 1, labelsize = 5, 
             repel = TRUE)
# Change points size, shape and fill color
# Change labelsize
fviz_pca_ind(res.pca, 
             pointsize = 3, pointshape = 21, fill = "lightblue",
             labelsize = 5, repel = TRUE)

fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (but not "text")
             group.ind = heart$target, # color by groups
             legend.title = "Groups",
             mean.point = FALSE)
fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (but not "text")
             group.ind = heart$target, # color by groups
             legend.title = "Groups",
             mean.point = TRUE)
fviz_pca_var(res.pca, axes.linetype = "blank")



ind.p <- fviz_pca_ind(res.pca, geom = "point", col.ind = heart$target)
ggpubr::ggpar(ind.p,
              title = "Principal Component Analysis",
              subtitle = "Heart disease data",
              caption = "Source: factoextra",
              xlab = "PC1", ylab = "PC2",
              legend.title = "Disease or No Disease", legend.position = "top",
              ggtheme = theme_gray(), palette = "jco"
)

fviz_pca_biplot(res.pca, repel = TRUE,col.ind = heart$target,
                col.var = "#2E9FDF", # Variables color
)

fviz_pca_biplot(res.pca, 
                col.ind = heart$target, palette = "jco", 
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Target") 

fviz_pca_biplot(res.pca, 
                # Fill individuals by groups
                geom.ind = "point",
                pointshape = 21,
                pointsize = 2.5,
                fill.ind = heart$target,
                col.ind = "black",
                # Color variable by groups
                legend.title = list(fill = "Target", color = "Clusters"),
                repel = TRUE        # Avoid label overplotting
)+
  ggpubr::fill_palette("jco")+      # Indiviual fill color
  ggpubr::color_palette("npg")      # Variable colors

fviz_pca_biplot(res.pca, 
                # Individuals
                geom.ind = "point",
                fill.ind = heart$target, col.ind = "black",
                pointshape = 21, pointsize = 2,
                palette = "jco",
                addEllipses = TRUE,
                # Variables
                alpha.var ="contrib", col.var = "contrib",
                gradient.cols = "RdYlBu",
                
                legend.title = list(fill = "Survivorship", color = "Contrib",
                                    alpha = "Contrib")
)









# Conclusion 

# From the above observations, we got the significant features which is correct:
# "sex", "cp", "restecg", "exang", "slope", "ca", "thal"  


# 1- Males are more vulnerable to be diagnosed with heart disease than females.
# 
# 2- Chest Pain is most common factor that leads to heart disease for males and females.
# 
# 3- Maximum heart rate achieved is the highest cause factor to cause heart disease for females where is Thalassemia is the highest to cause heart disease for males.
# 
# 4- There is a high association between chest pain and heart disease diagnosis.

# Limitation
# The dataset is missing some useful information such as smoking, obesity or family history 
# that can help in predicting.


# The following conditions are associated with increased prevalence of heart disease:

# Asymptomatic angina chest pain (relative to typical angina chest pain, atypical angina pain, or non-angina pain)
# Presence of exercise induced angina
# Lower fasting blood sugar
# Flat or down-sloaping peak exercise ST segment
# Presence of left ventricle hypertrophy
# Male
# Higher thelassemia score
# Higher age
# Lower max heart rate achieved
# Higher resting blood pressure
# Higher cholesterol
# Higher ST depression induced by exercise relative to rest

