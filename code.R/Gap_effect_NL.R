rm(list=ls())

library(mgcv)
library("readxl")
library(ggeffects)
library(ggpubr)
library(ggplot2)
library(ggrepel)
library(vegan)
library(dplyr)
library(nlme)
library(visreg)
library(gmodels)
library(metafor)
setwd("D:/Arbeit/Side Projects/Data gap effect")
source("MMK.R")

#import and clean data####
mydata <- as.data.frame(read_excel("Netherland.xlsx"))
#take out multiple species entries
unique(mydata$taxon)
mydata<-mydata[!grepl("/", mydata$taxon),]
mydata<-mydata[!grepl("sp.", mydata$taxon),]
mydata<-mydata[!grepl("spp.", mydata$taxon),]
#write.csv2(mydata,"crosstable.csv",sep=";")

#check species per year
Annual<-mydata %>%
  group_by(year) %>%  
  summarize(c(Species = length(taxon),sites=length(site_id)))
Annual

#fit without anything missing
set.seed(123)
mydata2 <- as.data.frame(read_excel("crosstable_NL.xlsx"))

#estimate richness and tot_abundance
# select columns 2 to 130
mydata22 <- mydata2 %>% select(2:141)
mydata22$Richness <- rowSums(mydata22 > 0)
mydata23 <- mydata2 %>% select(2:141)
mydata23$Total_Abundance <- apply(mydata23, 1, function(x) sum(x[which(x > 0)]))
df_no_gap<-cbind(mydata2$year,mydata2$iteration,mydata22$Richness, mydata23$Total_Abundance)
colnames(df_no_gap)<-c("year","richness","tot_abundance")
df_no_gap<-as.data.frame(df_no_gap)

#calculate annual species richness +SD
df_no_gap <- filter(df_no_gap, !is.na(richness))
Richness_0missing<-df_no_gap %>%
  dplyr::group_by(year) %>%  
  dplyr::summarize(mean_richness=mean(richness),
                   SD_richness=sd(richness))
Richness_0missing
colnames(Richness_0missing)<-c("year","mean_richness","sd_rich")

plot0miss<-ggplot(Richness_0missing, aes(x=year, y=mean_richness))+
  geom_point(color = "blue", size = 2)+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits = c(10, 50), breaks=seq(10, 50,5), labels = seq(10, 50,5))+
  geom_smooth()+
  theme_classic2()+
  theme_cleveland()
plot0miss
svg("NL.svg")
plot0miss
dev.off()


gam_0miss = gam(mean_richness ~ s(year), 
                data = Richness_0missing, method="REML")
fit<-as.data.frame(gam_0miss$fitted.values)
summary(gam_0miss)
#confint(gam_0miss)
AIC(gam_0miss)
colnames(fit)[1]<-"annual_fit"
years <- seq(1991, 2016)
years <- data.frame(year = years)
original_full_fit<-cbind(years,fit)

residuals <- residuals(gam_0miss)
mean_squared_error <- mean(residuals^2)
rmse_orig_gam <- sqrt(mean_squared_error)

lm_0miss = lm(mean_richness ~ year, 
              data = Richness_0missing)
deviance(lm_0miss) / (lm_0miss$df.residual)
confint(lm_0miss, level=0.95)
summary(lm_0miss)
fit<-as.data.frame(lm_0miss$fitted.values)
colnames(fit)[1]<-"annual_fit"
AIC(lm_0miss)
years <- seq(1991, 2016)
years <- data.frame(year = years)
original_full_fit_lm<-cbind(years,fit)

residuals <- residuals(lm_0miss)
mean_squared_error <- mean(residuals^2)
rmse_orig_lm <- sqrt(mean_squared_error)

#1 gap
#deleting rows randommly with 100 iterations#####
set.seed(123)
mydata2 <- as.data.frame(read_excel("crosstable_NL.xlsx"))

final_df_1_missing <- data.frame() 
for(i in 1:100){
  # randomly select a row and store its first column value
  random_row <- sample(2:(nrow(mydata2)-1), 1)
  first_col_value <- mydata2[random_row,1]
  # delete the randomly selected row
  mydata3 <- mydata2[-random_row,]
  # re-add the first column value of the deleted row
  mydata3 <- rbind(mydata3, c(first_col_value,rep(NA,ncol(mydata3)-1)))
  # sort dataframe by first column
  mydata3 <- mydata3[order(mydata3[,1]),]
  # add a new column representing the iteration number
  mydata3$iteration <- i
  #add the new dataframe to the previous one
  final_df_1_missing <- rbind(final_df_1_missing, mydata3)
}

#estimate richness and tot_abundance
# select columns 2 to 130
final_df_1_missing2 <- final_df_1_missing %>% select(2:141)
final_df_1_missing2$Richness <- rowSums(final_df_1_missing2 > 0)
final_df_1_missing3 <- final_df_1_missing %>% select(2:141)
final_df_1_missing3$Total_Abundance <- apply(final_df_1_missing3, 1, function(x) sum(x[which(x > 0)]))
final_df_1_gap<-cbind(final_df_1_missing$year,final_df_1_missing$iteration,final_df_1_missing2$Richness, final_df_1_missing3$Total_Abundance)
colnames(final_df_1_gap)<-c("year","iteration","richness","tot_abundance")
final_df_1_gap<-as.data.frame(final_df_1_gap)

#calculate annual species richness +SD
final_df_1_gap <- filter(final_df_1_gap, !is.na(richness))
Richness_1missing<-final_df_1_gap %>%
  dplyr::group_by(year) %>%  
  dplyr::summarize(mean_richness=mean(richness),
                   SD_richness=sd(richness))
Richness_1missing

colnames(Richness_1missing)<-c("year","mean_richness","sd_rich")

plot1miss<-ggplot(Richness_1missing, aes(x=year, y=mean_richness))+
  geom_point(color = "blue", size = 2)+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits = c(0, 50), breaks=seq(0, 50,5), labels = seq(0, 50,5))+
  geom_smooth()+
  theme_classic2()+
  theme_cleveland()
plot1miss

gam_1miss = gam(mean_richness ~ s(year), 
                data = Richness_1missing, method="REML")
fit<-as.data.frame(gam_1miss$fitted.values)
colnames(fit)[1]<-"annual_fit"
years <- seq(1991, 2016)
years <- data.frame(year = years)
fit1miss<-cbind(years,fit)

#loop to fit every iteration individually and then averaging the fits
#MK trend
MK_1_miss <- as.data.frame(matrix(nrow=0,ncol=11))
columns=c("Corrected Zc","new P-value","N/N*","Original Z","old P.value","Tau","Sen s slope","old.variance","new.variance","S statistic","n")
colnames(MK_1_miss)<-columns

iteration<-unique(final_df_1_gap$iteration)
for(n in iteration){
  subset_df<- final_df_1_gap[final_df_1_gap$iteration == n,]
  MK<-My.mmkh(subset_df$richness)
  MK_1_miss<-rbind(MK_1_miss,MK)
  names(MK_1_miss) <- c("Corrected Zc","new P-value","N/N*","Original Z","old P.value","Tau","Sen s slope","old.variance","new.variance","S statistic","n")
  cat(".")
}

#gam
AIC_1_miss_gam <- as.data.frame(matrix(nrow=0,ncol=1))
columns="AIC"
colnames(AIC_1_miss_gam)<-columns

DEV_1_miss_gam <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(DEV_1_miss_gam)<-columns

RMSE_1_miss_gam <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(RMSE_1_miss_gam)<-columns

loop_1_missing <- as.data.frame(matrix(nrow=0,ncol=1)) 
columns="fit"
colnames(loop_1_missing)<-columns
iteration<-unique(final_df_1_gap$iteration)
for(n in iteration){
  subset_df<- final_df_1_gap[final_df_1_gap$iteration == n,]
  gam_1miss2<-gam(richness~ s(year), data =subset_df)
  fit<-gam_1miss2$fitted.values
  AIC<-AIC(gam_1miss2)
  AIC<-as.data.frame(AIC)
  colnames(AIC)[1]<-"AIC"
  AIC_1_miss_gam<-rbind(AIC_1_miss_gam,AIC)
  residuals <- residuals(gam_1miss2)
  mean_squared_error <- mean(residuals^2)
  rmse <- sqrt(mean_squared_error)
  RMSE_1_miss_gam<-rbind(RMSE_1_miss_gam,rmse)
  DEV<-((gam_1miss2$null.deviance - gam_1miss2$deviance)/gam_1miss2$null.deviance) * 100
  DEV<-as.data.frame(DEV)
  colnames(DEV)[1]<-"DEV"
  DEV_1_miss_gam<-rbind(DEV_1_miss_gam,DEV)
  fit<-as.data.frame(fit)
  colnames(fit)[1]<-"fit"
  loop_1_missing<-rbind(loop_1_missing,fit)
  #fit1miss_lm<-cbind(final_df_1_gap$year,fit2)
  cat(".")
}
fit1miss_gam<-cbind(final_df_1_gap$year,final_df_1_gap$iteration,loop_1_missing)
head(fit1miss_gam)
colnames(fit1miss_gam)<-c("year","iteration","fit")
head(fit1miss_gam)

mean(DEV_1_miss_gam$DEV)
sd(DEV_1_miss_gam$DEV)

#lm
AIC_1_miss_lm <- as.data.frame(matrix(nrow=0,ncol=1))
columns="AIC"
colnames(AIC_1_miss_lm)<-columns

DEV_1_miss_lm <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(DEV_1_miss_lm)<-columns

RMSE_1_miss_lm <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(RMSE_1_miss_gam)<-columns

loop_1_missing <- as.data.frame(matrix(nrow=0,ncol=1)) 
columns="fit"
colnames(loop_1_missing)<-columns
iteration<-unique(final_df_1_gap$iteration)
for(n in iteration){
  subset_df<- final_df_1_gap[final_df_1_gap$iteration == n,]
  lm_1miss2<-lm(richness~ year, data =subset_df)
  fit<-lm_1miss2$fitted.values
  AIC<-AIC(lm_1miss2)
  AIC<-as.data.frame(AIC)
  colnames(AIC)[1]<-"AIC"
  AIC_1_miss_lm<-rbind(AIC_1_miss_lm,AIC)
  residuals <- residuals(lm_1miss2)
  mean_squared_error <- mean(residuals^2)
  rmse <- sqrt(mean_squared_error)
  RMSE_1_miss_lm<-rbind(RMSE_1_miss_lm,rmse)
  DEV<-deviance(lm_1miss2)/lm_1miss2$df.residual
  DEV<-as.data.frame(DEV)
  colnames(DEV)[1]<-"DEV"
  DEV_1_miss_lm<-rbind(DEV_1_miss_lm,DEV)
  fit<-as.data.frame(fit)
  colnames(fit)[1]<-"fit"
  loop_1_missing<-rbind(loop_1_missing,fit)
  #fit1miss_lm<-cbind(final_df_1_gap$year,fit2)
  cat(".")
}
fit1miss_lm<-cbind(final_df_1_gap$year,final_df_1_gap$iteration,loop_1_missing)
head(fit1miss_lm)
colnames(fit1miss_lm)<-c("year","iteration","fit")
head(fit1miss_lm)

mean(DEV_1_miss_lm$DEV)
sd(DEV_1_miss_lm$DEV)

Annual<-fit1miss_gam %>%
  group_by(year) %>%  
  summarize(c(mean = mean(fit),SD=sd(fit)))
Annual
table(fit1miss_gam$fit)


#using fitted values with gaps
#lm
plot1miss_lm<-ggplot()+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits = c(20, 50), breaks=seq(20, 50,5), labels = seq(20, 50,5))+
  geom_smooth(data=fit1miss_lm, aes(x=year, y=fit, group=iteration,colour =iteration),method = "loess", se=FALSE,size=0.25, alpha=0.5)+
  geom_smooth(data=original_full_fit_lm,aes(x=year, y=annual_fit),method = "loess",colour="red")+
  #geom_point(data=fit1miss_lm,aes(x=year, y=fit, group=iteration),color = "blue", size = 3)+
  #geom_point(data=original_full_fit_lm,aes(x=year, y=annual_fit),color = "red", size = 3)+
  theme_classic2()+
  theme_cleveland()
plot1miss_lm

#gam
plot1miss_gam<-ggplot()+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits =c(25, 50), breaks=seq(25, 50,5), labels = seq(25, 50,5))+
  geom_smooth(data=fit1miss_gam, aes(x=year, y=fit, group=iteration,colour =iteration),method = "loess", se=FALSE,size=0.25, alpha=0.5)+
  geom_smooth(data=original_full_fit,aes(x=year, y=annual_fit),method = "loess",se=FALSE,colour="red")+
  theme_classic2()+
  theme_cleveland()
plot1miss_gam

#raw data with loess and gaps
plot1miss_raw<-ggplot()+
  geom_smooth(data=final_df_1_gap, aes(x=year, y=richness, group=iteration,colour =iteration),method = "loess",se=FALSE, size=0.25, alpha=0.5)+
  geom_smooth(data=df_no_gap,aes(x=year, y=richness),method = "loess", colour="red", se=FALSE)+
  geom_point(data=df_no_gap,aes(x=year, y=richness),color = "red", size = 1)+
  theme_classic2()+
  theme_cleveland()
plot1miss_raw

#2 gaps
mydata2 <- as.data.frame(read_excel("crosstable_NL.xlsx"))
final_df_2_missing <- data.frame() 
for(i in 1:100){
  # randomly select a row and store its first column value
  random_row <- sample(2:(nrow(mydata2)-1), 2)
  first_col_value <- mydata2[random_row,1]
  # delete the randomly selected row
  mydata3 <- mydata2[-random_row,]
  # re-add the first column values of the deleted rows
  for (j in 1:length(first_col_value)) {
    mydata3 <- rbind(mydata3, c(first_col_value[j],rep(NA,ncol(mydata3)-1)))
  }
  # sort dataframe by first column
  mydata3 <- mydata3[order(mydata3[,1]),]
  # add a new column representing the iteration number
  mydata3$iteration <- i
  #add the new dataframe to the previous one
  final_df_2_missing <- rbind(final_df_2_missing, mydata3)
}

#estimate richness and tot_abundance
# select columns 2 to 130
final_df_2_missing2 <- final_df_2_missing %>% select(2:141)
final_df_2_missing2$Richness <- rowSums(final_df_2_missing2 > 0)
final_df_2_missing3 <- final_df_2_missing %>% select(2:141)
final_df_2_missing3$Total_Abundance <- apply(final_df_2_missing3, 1, function(x) sum(x[which(x > 0)]))
final_df_2_gap<-cbind(final_df_2_missing$year,final_df_2_missing$iteration,final_df_2_missing2$Richness, final_df_2_missing3$Total_Abundance)
colnames(final_df_2_gap)<-c("year","iteration","richness","tot_abundance")
final_df_2_gap<-as.data.frame(final_df_2_gap)

#calculate annual species richness +SD
final_df_2_gap <- filter(final_df_2_gap, !is.na(richness))
Richness_2missing<-final_df_2_gap %>%
  dplyr::group_by(year) %>%  
  dplyr::summarize(mean_richness=mean(richness),
                   SD_richness=sd(richness))
Richness_2missing
colnames(Richness_2missing)<-c("year","mean_richness","sd_rich")

plot2miss<-ggplot(Richness_2missing, aes(x=year, y=mean_richness))+
  geom_point(color = "blue", size = 2)+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits = c(20, 50), breaks=seq(20, 50,5), labels = seq(20, 50,5))+
  geom_smooth()+
  theme_classic2()+
  theme_cleveland()
plot2miss

gam_2miss = gam(mean_richness ~ s(year), 
                data = Richness_2missing, method="REML")
fit<-as.data.frame(gam_2miss$fitted.values)
colnames(fit)[1]<-"annual_fit"
years <- seq(1991, 2016)
years <- data.frame(year = years)
fit2miss<-cbind(years,fit)

#loop to fit every iteration individually and then averaging the fits
#MK trend
MK_2_miss <- as.data.frame(matrix(nrow=0,ncol=11))
columns=c("Corrected Zc","new P-value","N/N*","Original Z","old P.value","Tau","Sen s slope","old.variance","new.variance","S statistic","n")
colnames(MK_2_miss)<-columns

iteration<-unique(final_df_2_gap$iteration)
for(n in iteration){
  subset_df<- final_df_2_gap[final_df_2_gap$iteration == n,]
  MK<-My.mmkh(subset_df$richness)
  MK_2_miss<-rbind(MK_2_miss,MK)
  names(MK_2_miss) <- c("Corrected Zc","new P-value","N/N*","Original Z","old P.value","Tau","Sen s slope","old.variance","new.variance","S statistic","n")
  cat(".")
}

#gam
AIC_2_miss_gam <- as.data.frame(matrix(nrow=0,ncol=1))
columns="AIC"
colnames(AIC_2_miss_gam)<-columns

DEV_2_miss_gam <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(DEV_2_miss_gam)<-columns

RMSE_2_miss_gam <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(RMSE_2_miss_gam)<-columns

loop_2_missing <- as.data.frame(matrix(nrow=0,ncol=1)) 
columns="fit"
colnames(loop_2_missing)<-columns
iteration<-unique(final_df_2_gap$iteration)
for(n in iteration){
  subset_df<- final_df_2_gap[final_df_2_gap$iteration == n,]
  gam_2miss2<-gam(richness~ s(year), data =subset_df)
  fit<-gam_2miss2$fitted.values
  AIC<-AIC(gam_2miss2)
  AIC<-as.data.frame(AIC)
  colnames(AIC)[1]<-"AIC"
  AIC_2_miss_gam<-rbind(AIC_2_miss_gam,AIC)
  residuals <- residuals(gam_2miss2)
  mean_squared_error <- mean(residuals^2)
  rmse <- sqrt(mean_squared_error)
  RMSE_2_miss_gam<-rbind(RMSE_2_miss_gam,rmse)
  DEV<-((gam_2miss2$null.deviance - gam_2miss2$deviance)/gam_2miss2$null.deviance) * 100
  DEV<-as.data.frame(DEV)
  colnames(DEV)[1]<-"DEV"
  DEV_2_miss_gam<-rbind(DEV_2_miss_gam,DEV)
  fit<-as.data.frame(fit)
  colnames(fit)[1]<-"fit"
  loop_2_missing<-rbind(loop_2_missing,fit)
  #fit2miss_lm<-cbind(final_df_2_gap$year,fit2)
  cat(".")
}
fit2miss_gam<-cbind(final_df_2_gap$year,final_df_2_gap$iteration,loop_2_missing)
head(fit2miss_gam)
colnames(fit2miss_gam)<-c("year","iteration","fit")
head(fit2miss_gam)
mean(DEV_2_miss_gam$DEV)
sd(DEV_2_miss_gam$DEV)

#lm
AIC_2_miss_lm <- as.data.frame(matrix(nrow=0,ncol=1))
columns="AIC"
colnames(AIC_1_miss_lm)<-columns

DEV_2_miss_lm <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(DEV_2_miss_lm)<-columns

RMSE_2_miss_lm <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(RMSE_2_miss_lm)<-columns



loop_2_missing <- as.data.frame(matrix(nrow=0,ncol=1)) 
columns="fit"
colnames(loop_2_missing)<-columns
iteration<-unique(final_df_2_gap$iteration)
for(n in iteration){
  subset_df<- final_df_2_gap[final_df_2_gap$iteration == n,]
  lm_2miss2<-lm(richness~ year, data =subset_df)
  fit<-lm_2miss2$fitted.values
  AIC<-AIC(lm_2miss2)
  AIC<-as.data.frame(AIC)
  colnames(AIC)[1]<-"AIC"
  AIC_2_miss_lm<-rbind(AIC_2_miss_lm,AIC)
  residuals <- residuals(lm_2miss2)
  mean_squared_error <- mean(residuals^2)
  rmse <- sqrt(mean_squared_error)
  RMSE_2_miss_lm<-rbind(RMSE_2_miss_lm,rmse)
  DEV<-deviance(lm_2miss2)/lm_2miss2$df.residual
  DEV<-as.data.frame(DEV)
  colnames(DEV)[1]<-"DEV"
  DEV_2_miss_lm<-rbind(DEV_2_miss_lm,DEV)
  fit<-as.data.frame(fit)
  colnames(fit)[1]<-"fit"
  loop_2_missing<-rbind(loop_2_missing,fit)
  #fit2miss_lm<-cbind(final_df_2_gap$year,fit2)
  cat(".")
}
fit2miss_lm<-cbind(final_df_2_gap$year,final_df_2_gap$iteration,loop_2_missing)
head(fit2miss_lm)
colnames(fit2miss_lm)<-c("year","iteration","fit")
head(fit2miss_lm)

mean(DEV_2_miss_lm$DEV)
sd(DEV_2_miss_lm$DEV)

Annual<-fit2miss_lm %>%
  group_by(year) %>%  
  summarize(c(mean = mean(fit),SD=sd(fit)))
Annual
table(fit2miss_lm$fit)


#using fitted values with gaps
#lm
plot2miss_lm<-ggplot()+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits = c(20, 50), breaks=seq(20, 50,5), labels = seq(20, 50,5))+
  geom_smooth(data=fit2miss_lm, aes(x=year, y=fit, group=iteration,colour =iteration),method = "loess", se=FALSE,size=0.25, alpha=0.5)+
  geom_smooth(data=original_full_fit_lm,aes(x=year, y=annual_fit),method = "loess",colour="red")+
  #geom_point(data=fit2miss_lm,aes(x=year, y=fit, group=iteration),color = "blue", size = 3)+
  #geom_point(data=original_full_fit_lm,aes(x=year, y=annual_fit),color = "red", size = 3)+
  theme_classic2()+
  theme_cleveland()
plot2miss_lm

#gam
plot2miss_gam<-ggplot()+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits =c(25, 50), breaks=seq(25, 50,5), labels = seq(25, 50,5))+
  geom_smooth(data=fit2miss_gam, aes(x=year, y=fit, group=iteration,colour =iteration),method = "loess", se=FALSE,size=0.25, alpha=0.5)+
  geom_smooth(data=original_full_fit,aes(x=year, y=annual_fit),method = "loess",se=FALSE,colour="red")+
  theme_classic2()+
  theme_cleveland()
plot2miss_gam

#raw data with loess and gaps
plot2miss_raw<-ggplot()+
  geom_smooth(data=final_df_2_gap, aes(x=year, y=richness, group=iteration,colour =iteration),method = "loess",se=FALSE, size=0.25, alpha=0.5)+
  geom_smooth(data=df_no_gap,aes(x=year, y=richness),method = "loess", colour="red", se=FALSE)+
  geom_point(data=df_no_gap,aes(x=year, y=richness),color = "red", size = 1)+
  theme_classic2()+
  theme_cleveland()
plot2miss_raw


#5 gaps
mydata2 <- as.data.frame(read_excel("crosstable_NL.xlsx"))
final_df_5_missing <- data.frame() 
for(i in 1:100){
  # randomly select a row and store its first column value
  random_row <- sample(2:(nrow(mydata2)-1), 5)
  first_col_value <- mydata2[random_row,1]
  # delete the randomly selected row
  mydata3 <- mydata2[-random_row,]
  # re-add the first column values of the deleted rows
  for (j in 1:length(first_col_value)) {
    mydata3 <- rbind(mydata3, c(first_col_value[j],rep(NA,ncol(mydata3)-1)))
  }
  # sort dataframe by first column
  mydata3 <- mydata3[order(mydata3[,1]),]
  # add a new column representing the iteration number
  mydata3$iteration <- i
  #add the new dataframe to the previous one
  final_df_5_missing <- rbind(final_df_5_missing, mydata3)
}

#estimate richness and tot_abundance
# select columns 2 to 130
final_df_5_missing2 <- final_df_5_missing %>% select(2:141)
final_df_5_missing2$Richness <- rowSums(final_df_5_missing2 > 0)
final_df_5_missing3 <- final_df_5_missing %>% select(2:141)
final_df_5_missing3$Total_Abundance <- apply(final_df_5_missing3, 1, function(x) sum(x[which(x > 0)]))
final_df_5_gap<-cbind(final_df_5_missing$year,final_df_5_missing$iteration,final_df_5_missing2$Richness, final_df_5_missing3$Total_Abundance)
colnames(final_df_5_gap)<-c("year","iteration","richness","tot_abundance")
final_df_5_gap<-as.data.frame(final_df_5_gap)

#calculate annual species richness +SD
final_df_5_gap <- filter(final_df_5_gap, !is.na(richness))
Richness_5missing<-final_df_5_gap %>%
  dplyr::group_by(year) %>%  
  dplyr::summarize(mean_richness=mean(richness),
                   SD_richness=sd(richness))
Richness_5missing
colnames(Richness_5missing)<-c("year","mean_richness","sd_rich")

plot5miss<-ggplot(Richness_5missing, aes(x=year, y=mean_richness))+
  geom_point(color = "blue", size = 2)+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits = c(20, 50), breaks=seq(20, 50,5), labels = seq(20, 50,5))+
  geom_smooth()+
  theme_classic2()+
  theme_cleveland()
plot5miss

gam_5miss = gam(mean_richness ~ s(year), 
                data = Richness_5missing, method="REML")
fit<-as.data.frame(gam_5miss$fitted.values)
colnames(fit)[1]<-"annual_fit"
years <- seq(1991, 2016)
years <- data.frame(year = years)
fit5miss<-cbind(years,fit)

#loop to fit every iteration individually and then averaging the fits
#MK trend
MK_5_miss <- as.data.frame(matrix(nrow=0,ncol=11))
columns=c("Corrected Zc","new P-value","N/N*","Original Z","old P.value","Tau","Sen s slope","old.variance","new.variance","S statistic","n")
colnames(MK_5_miss)<-columns

iteration<-unique(final_df_5_gap$iteration)
for(n in iteration){
  subset_df<- final_df_5_gap[final_df_5_gap$iteration == n,]
  MK<-My.mmkh(subset_df$richness)
  MK_5_miss<-rbind(MK_5_miss,MK)
  names(MK_5_miss) <- c("Corrected Zc","new P-value","N/N*","Original Z","old P.value","Tau","Sen s slope","old.variance","new.variance","S statistic","n")
  cat(".")
}

#gam
AIC_5_miss_gam <- as.data.frame(matrix(nrow=0,ncol=1))
columns="AIC"
colnames(AIC_5_miss_gam)<-columns

DEV_5_miss_gam <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(DEV_5_miss_gam)<-columns

RMSE_5_miss_gam <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(RMSE_5_miss_gam)<-columns



loop_5_missing <- as.data.frame(matrix(nrow=0,ncol=1)) 
columns="fit"
colnames(loop_5_missing)<-columns
iteration<-unique(final_df_5_gap$iteration)
for(n in iteration){
  subset_df<- final_df_5_gap[final_df_5_gap$iteration == n,]
  gam_5miss2<-gam(richness~ s(year), data =subset_df)
  fit<-gam_5miss2$fitted.values
  AIC<-AIC(gam_5miss2)
  AIC<-as.data.frame(AIC)
  colnames(AIC)[1]<-"AIC"
  AIC_5_miss_gam<-rbind(AIC_5_miss_gam,AIC)
  residuals <- residuals(gam_5miss2)
  mean_squared_error <- mean(residuals^2)
  rmse <- sqrt(mean_squared_error)
  RMSE_5_miss_gam<-rbind(RMSE_5_miss_gam,rmse)
  DEV<-((gam_5miss2$null.deviance - gam_5miss2$deviance)/gam_5miss2$null.deviance) * 100
  DEV<-as.data.frame(DEV)
  colnames(DEV)[1]<-"DEV"
  DEV_5_miss_gam<-rbind(DEV_5_miss_gam,DEV)
  fit<-as.data.frame(fit)
  colnames(fit)[1]<-"fit"
  loop_5_missing<-rbind(loop_5_missing,fit)
  #fit5miss_lm<-cbind(final_df_5_gap$year,fit2)
  cat(".")
}
fit5miss_gam<-cbind(final_df_5_gap$year,final_df_5_gap$iteration,loop_5_missing)
head(fit5miss_gam)
colnames(fit5miss_gam)<-c("year","iteration","fit")
head(fit5miss_gam)
mean(DEV_5_miss_gam$DEV)
sd(DEV_5_miss_gam$DEV)
#lm
AIC_5_miss_lm <- as.data.frame(matrix(nrow=0,ncol=1))
columns="AIC"
colnames(AIC_5_miss_lm)<-columns

DEV_5_miss_lm <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(DEV_5_miss_lm)<-columns

RMSE_5_miss_lm <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(RMSE_5_miss_lm)<-columns

loop_5_missing <- as.data.frame(matrix(nrow=0,ncol=1)) 
columns="fit"
colnames(loop_5_missing)<-columns
iteration<-unique(final_df_5_gap$iteration)
for(n in iteration){
  subset_df<- final_df_5_gap[final_df_5_gap$iteration == n,]
  lm_5miss2<-lm(richness~ year, data =subset_df)
  fit<-lm_5miss2$fitted.values
  AIC<-AIC(lm_5miss2)
  AIC<-as.data.frame(AIC)
  colnames(AIC)[1]<-"AIC"
  AIC_5_miss_lm<-rbind(AIC_5_miss_lm,AIC)
  residuals <- residuals(lm_5miss2)
  mean_squared_error <- mean(residuals^2)
  rmse <- sqrt(mean_squared_error)
  RMSE_5_miss_lm<-rbind(RMSE_5_miss_lm,rmse)
  DEV<-deviance(lm_5miss2)/lm_5miss2$df.residual
  DEV<-as.data.frame(DEV)
  colnames(DEV)[1]<-"DEV"
  DEV_5_miss_lm<-rbind(DEV_5_miss_lm,DEV)
  fit<-as.data.frame(fit)
  colnames(fit)[1]<-"fit"
  loop_5_missing<-rbind(loop_5_missing,fit)
  #fit5miss_lm<-cbind(final_df_5_gap$year,fit2)
  cat(".")
}
fit5miss_lm<-cbind(final_df_5_gap$year,final_df_5_gap$iteration,loop_5_missing)
head(fit5miss_lm)
colnames(fit5miss_lm)<-c("year","iteration","fit")
head(fit5miss_lm)
mean(DEV_5_miss_lm$DEV)
sd(DEV_5_miss_lm$DEV)

Annual<-fit5miss_lm %>%
  group_by(year) %>%  
  summarize(c(mean = mean(fit),SD=sd(fit)))
Annual
table(fit5miss_lm$fit)


#using fitted values with gaps
#lm
plot5miss_lm<-ggplot()+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits = c(20, 50), breaks=seq(20, 50,5), labels = seq(20, 50,5))+
  geom_smooth(data=fit5miss_lm, aes(x=year, y=fit, group=iteration,colour =iteration),method = "loess", se=FALSE,size=0.25, alpha=0.5)+
  geom_smooth(data=original_full_fit_lm,aes(x=year, y=annual_fit),method = "loess",colour="red")+
  #geom_point(data=fit5miss_lm,aes(x=year, y=fit, group=iteration),color = "blue", size = 3)+
  #geom_point(data=original_full_fit_lm,aes(x=year, y=annual_fit),color = "red", size = 3)+
  theme_classic2()+
  theme_cleveland()
plot5miss_lm

#gam
plot5miss_gam<-ggplot()+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits =c(25, 50), breaks=seq(25, 50,5), labels = seq(25, 50,5))+
  geom_smooth(data=fit5miss_gam, aes(x=year, y=fit, group=iteration,colour =iteration),method = "loess", se=FALSE,size=0.25, alpha=0.5)+
  geom_smooth(data=original_full_fit,aes(x=year, y=annual_fit),method = "loess",se=FALSE,colour="red")+
  theme_classic2()+
  theme_cleveland()
plot5miss_gam

#raw data with loess and gaps
plot5miss_raw<-ggplot()+
  geom_smooth(data=final_df_5_gap, aes(x=year, y=richness, group=iteration,colour =iteration),method = "loess",se=FALSE, size=0.25, alpha=0.5)+
  geom_smooth(data=df_no_gap,aes(x=year, y=richness),method = "loess", colour="red", se=FALSE)+
  geom_point(data=df_no_gap,aes(x=year, y=richness),color = "red", size = 1)+
  theme_classic2()+
  theme_cleveland()
plot5miss_raw


#10 gaps
mydata2 <- as.data.frame(read_excel("crosstable_NL.xlsx"))
final_df_10_missing <- data.frame() 
for(i in 1:100){
  # randomly select a row and store its first column value
  random_row <- sample(2:(nrow(mydata2)-1), 10)
  first_col_value <- mydata2[random_row,1]
  # delete the randomly selected row
  mydata3 <- mydata2[-random_row,]
  # re-add the first column values of the deleted rows
  for (j in 1:length(first_col_value)) {
    mydata3 <- rbind(mydata3, c(first_col_value[j],rep(NA,ncol(mydata3)-1)))
  }
  # sort dataframe by first column
  mydata3 <- mydata3[order(mydata3[,1]),]
  # add a new column representing the iteration number
  mydata3$iteration <- i
  #add the new dataframe to the previous one
  final_df_10_missing <- rbind(final_df_10_missing, mydata3)
}

#estimate richness and tot_abundance
# select columns 2 to 130
final_df_10_missing2 <- final_df_10_missing %>% select(2:141)
final_df_10_missing2$Richness <- rowSums(final_df_10_missing2 > 0)
final_df_10_missing3 <- final_df_10_missing %>% select(2:141)
final_df_10_missing3$Total_Abundance <- apply(final_df_10_missing3, 1, function(x) sum(x[which(x > 0)]))
final_df_10_gap<-cbind(final_df_10_missing$year,final_df_10_missing$iteration,final_df_10_missing2$Richness, final_df_10_missing3$Total_Abundance)
colnames(final_df_10_gap)<-c("year","iteration","richness","tot_abundance")
final_df_10_gap<-as.data.frame(final_df_10_gap)

#calculate annual species richness +SD
final_df_10_gap <- filter(final_df_10_gap, !is.na(richness))
Richness_10missing<-final_df_10_gap %>%
  dplyr::group_by(year) %>%  
  dplyr::summarize(mean_richness=mean(richness),
                   SD_richness=sd(richness))
Richness_10missing
colnames(Richness_10missing)<-c("year","mean_richness","sd_rich")

plot10miss<-ggplot(Richness_10missing, aes(x=year, y=mean_richness))+
  geom_point(color = "blue", size = 2)+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits = c(20, 50), breaks=seq(20, 50,5), labels = seq(20, 50,5))+
  geom_smooth()+
  theme_classic2()+
  theme_cleveland()
plot10miss

gam_10miss = gam(mean_richness ~ s(year), 
                 data = Richness_10missing, method="REML")
fit<-as.data.frame(gam_10miss$fitted.values)
colnames(fit)[1]<-"annual_fit"
years <- seq(1991, 2016)
years <- data.frame(year = years)
fit10miss<-cbind(years,fit)

#loop to fit every iteration individually and then averaging the fits
#MK trend
MK_10_miss <- as.data.frame(matrix(nrow=0,ncol=11))
columns=c("Corrected Zc","new P-value","N/N*","Original Z","old P.value","Tau","Sen s slope","old.variance","new.variance","S statistic","n")
colnames(MK_10_miss)<-columns

iteration<-unique(final_df_10_gap$iteration)
for(n in iteration){
  subset_df<- final_df_10_gap[final_df_10_gap$iteration == n,]
  MK<-My.mmkh(subset_df$richness)
  MK_10_miss<-rbind(MK_10_miss,MK)
  names(MK_10_miss) <- c("Corrected Zc","new P-value","N/N*","Original Z","old P.value","Tau","Sen s slope","old.variance","new.variance","S statistic","n")
  cat(".")
}
#gam
AIC_10_miss_gam <- as.data.frame(matrix(nrow=0,ncol=1))
columns="AIC"
colnames(AIC_10_miss_gam)<-columns

DEV_10_miss_gam <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(DEV_10_miss_gam)<-columns

RMSE_10_miss_gam <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(RMSE_10_miss_gam)<-columns



loop_10_missing <- as.data.frame(matrix(nrow=0,ncol=1)) 
columns="fit"
colnames(loop_10_missing)<-columns
iteration<-unique(final_df_10_gap$iteration)
for(n in iteration){
  subset_df<- final_df_10_gap[final_df_10_gap$iteration == n,]
  gam_10miss2<-gam(richness~ s(year), data =subset_df)
  fit<-gam_10miss2$fitted.values
  AIC<-AIC(gam_10miss2)
  AIC<-as.data.frame(AIC)
  colnames(AIC)[1]<-"AIC"
  AIC_10_miss_gam<-rbind(AIC_10_miss_gam,AIC)
  residuals <- residuals(gam_10miss2)
  mean_squared_error <- mean(residuals^2)
  rmse <- sqrt(mean_squared_error)
  RMSE_10_miss_gam<-rbind(RMSE_10_miss_gam,rmse)
  DEV<-((gam_10miss2$null.deviance - gam_10miss2$deviance)/gam_10miss2$null.deviance) * 100
  DEV<-as.data.frame(DEV)
  colnames(DEV)[1]<-"DEV"
  DEV_10_miss_gam<-rbind(DEV_10_miss_gam,DEV)
  fit<-as.data.frame(fit)
  colnames(fit)[1]<-"fit"
  loop_10_missing<-rbind(loop_10_missing,fit)
  #fit10miss_lm<-cbind(final_df_10_gap$year,fit2)
  cat(".")
}
fit10miss_gam<-cbind(final_df_10_gap$year,final_df_10_gap$iteration,loop_10_missing)
head(fit10miss_gam)
colnames(fit10miss_gam)<-c("year","iteration","fit")
head(fit10miss_gam)
mean(DEV_10_miss_gam$DEV)
sd(DEV_10_miss_gam$DEV)

#lm
AIC_10_miss_lm <- as.data.frame(matrix(nrow=0,ncol=1))
columns="AIC"
colnames(AIC_10_miss_lm)<-columns

DEV_10_miss_lm <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(DEV_10_miss_lm)<-columns

RMSE_10_miss_lm <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(RMSE_10_miss_lm)<-columns



loop_10_missing <- as.data.frame(matrix(nrow=0,ncol=1)) 
columns="fit"
colnames(loop_10_missing)<-columns
iteration<-unique(final_df_10_gap$iteration)
for(n in iteration){
  subset_df<- final_df_10_gap[final_df_10_gap$iteration == n,]
  lm_10miss2<-lm(richness~ year, data =subset_df)
  fit<-lm_10miss2$fitted.values
  AIC<-AIC(lm_10miss2)
  AIC<-as.data.frame(AIC)
  colnames(AIC)[1]<-"AIC"
  AIC_10_miss_lm<-rbind(AIC_10_miss_lm,AIC)
  residuals <- residuals(lm_10miss2)
  mean_squared_error <- mean(residuals^2)
  rmse <- sqrt(mean_squared_error)
  RMSE_10_miss_lm<-rbind(RMSE_10_miss_lm,rmse)
  DEV<-deviance(lm_10miss2)/lm_10miss2$df.residual
  DEV<-as.data.frame(DEV)
  colnames(DEV)[1]<-"DEV"
  DEV_10_miss_lm<-rbind(DEV_10_miss_lm,DEV)
  fit<-as.data.frame(fit)
  colnames(fit)[1]<-"fit"
  loop_10_missing<-rbind(loop_10_missing,fit)
  #fit10miss_lm<-cbind(final_df_10_gap$year,fit2)
  cat(".")
}
fit10miss_lm<-cbind(final_df_10_gap$year,final_df_10_gap$iteration,loop_10_missing)
head(fit10miss_lm)
colnames(fit10miss_lm)<-c("year","iteration","fit")
head(fit10miss_lm)
mean(DEV_10_miss_lm$DEV)
sd(DEV_10_miss_lm$DEV)

Annual<-fit10miss_lm %>%
  group_by(year) %>%  
  summarize(c(mean = mean(fit),SD=sd(fit)))
Annual
table(fit10miss_lm$fit)
#using fitted values with gaps
#lm
plot10miss_lm<-ggplot()+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits = c(20, 50), breaks=seq(20, 50,5), labels = seq(20, 50,5))+
  geom_smooth(data=fit10miss_lm, aes(x=year, y=fit, group=iteration,colour =iteration),method = "loess", se=FALSE,size=0.25, alpha=0.5)+
  geom_smooth(data=original_full_fit_lm,aes(x=year, y=annual_fit),method = "loess",colour="red")+
  #geom_point(data=fit10miss_lm,aes(x=year, y=fit, group=iteration),color = "blue", size = 3)+
  #geom_point(data=original_full_fit_lm,aes(x=year, y=annual_fit),color = "red", size = 3)+
  theme_classic2()+
  theme_cleveland()
plot10miss_lm

#gam
plot10miss_gam<-ggplot()+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits =c(25, 50), breaks=seq(25, 50,5), labels = seq(25, 50,5))+
  geom_smooth(data=fit10miss_gam, aes(x=year, y=fit, group=iteration,colour =iteration),method = "gam", se=FALSE,size=0.25, alpha=0.5)+
  geom_smooth(data=original_full_fit,aes(x=year, y=annual_fit),method = "gam",se=FALSE,colour="red")+
  theme_classic2()+
  theme_cleveland()
plot10miss_gam

#raw data with loess and gaps
plot10miss_raw<-ggplot()+
  geom_smooth(data=final_df_10_gap, aes(x=year, y=richness, group=iteration),method = "loess",se=FALSE, size=0.25, alpha=0.5)+
  geom_point(data=final_df_10_gap,aes(x=year, y=richness, group=iteration),color = "blue", size = 3)+
  geom_smooth(data=df_no_gap,aes(x=year, y=richness),method = "loess", colour="red", se=FALSE)+
  geom_point(data=df_no_gap,aes(x=year, y=richness),color = "red", size = 1)+
  theme_classic2()+
  theme_cleveland()
plot10miss_raw



#overall plot

#lm
plot_miss_overall_LM<-ggplot()+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits =c(25, 50), breaks=seq(25, 50,5), labels = seq(25, 50,5))+
  geom_smooth(data=fit10miss_lm, aes(x=year, y=fit, group=iteration,colour =iteration),colour ="red",se=FALSE, size=0.25, alpha=0.5)+
  geom_smooth(data=fit5miss_lm,aes(x=year, y=fit, group=iteration,colour =iteration),colour ="orange",se=FALSE, size=0.25, alpha=0.5)+
  geom_smooth(data=fit2miss_lm, aes(x=year, y=fit, group=iteration,colour =iteration),colour ="yellow",se=FALSE, size=0.25, alpha=0.5)+
  geom_smooth(data=fit1miss_lm, aes(x=year, y=fit, group=iteration,colour =iteration),colour ="green",se=FALSE, size=0.25, alpha=0.5)+
  geom_smooth(data=original_full_fit_lm,aes(x=year, y=annual_fit),colour="black",size=1)+
  theme_classic2()+
  theme_cleveland()
plot_miss_overall_LM

#gam
plot_miss_overall_GAM<-ggplot()+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits = c(25, 50), breaks=seq(25, 50,5), labels = seq(25, 50,5))+
  geom_smooth(data=fit10miss_gam, aes(x=year, y=fit, group=iteration),colour ="red",se=FALSE, size=0.25, alpha=0.5)+
  geom_smooth(data=fit5miss_gam, aes(x=year, y=fit, group=iteration),colour ="orange",se=FALSE, size=0.25, alpha=0.5)+
  geom_smooth(data=fit2miss_gam, aes(x=year, y=fit, group=iteration),colour ="yellow",se=FALSE, size=0.25, alpha=0.5)+
  geom_smooth(data=fit1miss_gam, aes(x=year, y=fit, group=iteration),colour ="green",se=FALSE, size=0.25, alpha=0.5)+
  geom_line(data=original_full_fit,aes(x=year, y=annual_fit),colour="black",size=1)+
  theme_classic2()+
  theme_cleveland()
plot_miss_overall_GAM

#raw data with loess and gaps
plot_miss_overall_raw<-ggplot()+
  geom_smooth(data=final_df_10_gap, aes(x=year, y=richness, group=iteration), colour = "red",method = "loess",se=FALSE, size=0.25, alpha=0.5)+  
  geom_smooth(data=final_df_5_gap, aes(x=year, y=richness, group=iteration), colour = "orange",method = "loess",se=FALSE, size=0.25, alpha=0.5)+  
  geom_smooth(data=final_df_2_gap, aes(x=year, y=richness, group=iteration), colour ="yellow",method = "loess",se=FALSE, size=0.25, alpha=0.5)+  
  geom_smooth(data=final_df_1_gap, aes(x=year, y=richness, group=iteration), colour="green",method = "loess",se=FALSE, size=0.25, alpha=0.5)+
  geom_smooth(data=df_no_gap,aes(x=year, y=richness),method = "loess", colour="black", se=FALSE)+
  geom_point(data=df_no_gap,aes(x=year, y=richness),color = "black", size = 1)+
  theme_classic2()+
  theme_cleveland()
plot_miss_overall_raw

#detach("package:sjPlot", unload=TRUE)
library(cowplot)
plot1<-plot_grid(plot1miss_lm, plot1miss_gam,
                 plot2miss_lm, plot2miss_gam,
                 plot5miss_lm, plot5miss_gam,
                 plot10miss_lm, plot10miss_gam,
                 #labels = "AUTO", 
                 ncol = 2)
plot1
svg("plot_ALL_missing_NL_sep.svg")
plot1
dev.off()

plot2<-plot_grid(plot_miss_overall_LM,
                 plot_miss_overall_GAM,
                 plot_miss_overall_raw,
                 #labels = "AUTO", 
                 ncol = 1)
plot2
svg("plot_ALL_missing_NL.svg")
plot2
dev.off()


#Exporting Avg + SD of AICs
mean(AIC_1_miss_lm$AIC)
sd(AIC_1_miss_lm$AIC)
mean(AIC_2_miss_lm$AIC)
sd(AIC_2_miss_lm$AIC)
mean(AIC_5_miss_lm$AIC)
sd(AIC_5_miss_lm$AIC)
mean(AIC_10_miss_lm$AIC)
sd(AIC_10_miss_lm$AIC)
mean(AIC_1_miss_gam$AIC)
sd(AIC_1_miss_gam$AIC)
mean(AIC_2_miss_gam$AIC)
sd(AIC_2_miss_gam$AIC)
mean(AIC_5_miss_gam$AIC)
sd(AIC_5_miss_gam$AIC)
mean(AIC_10_miss_gam$AIC)
sd(AIC_10_miss_gam$AIC)


#data mingling crosstable to list
library(tidyr)
final_df_1_missing_list <- pivot_longer(data = final_df_1_missing, cols = !c("year","iteration"),
                                        names_to = "taxon", values_to = "abundance")
final_df_1_missing_list<-as.data.frame(final_df_1_missing_list)
final_df_2_missing_list <- pivot_longer(data = final_df_2_missing, cols = !c("year","iteration"),
                                        names_to = "taxon", values_to = "abundance")
final_df_2_missing_list<-as.data.frame(final_df_2_missing_list)
final_df_5_missing_list <- pivot_longer(data = final_df_5_missing, cols = !c("year","iteration"),
                                        names_to = "taxon", values_to = "abundance")
final_df_5_missing_list<-as.data.frame(final_df_5_missing_list)
final_df_10_missing_list <- pivot_longer(data = final_df_10_missing, cols = !c("year","iteration"),
                                         names_to = "taxon", values_to = "abundance")
final_df_10_missing_list<-as.data.frame(final_df_10_missing_list)

#Mice loop#####
library(mice)

#1 gap filled
mice_1_filled <- data.frame()
iteration= unique(final_df_1_missing_list$iteration)
n <- iteration[1]
for(n in iteration){
  subset_df<- final_df_1_missing_list[final_df_1_missing_list$iteration == n,]
  
  # Extract the unique species from the subseted dataframe
  species <- unique(subset_df$taxon)
  for(j in 1:length(species)){
    # Subset the data frame to include only the current species
    subset_species_df <- subset_df %>% filter(taxon == species[j])
    # check if the species has less than 2  values
    if(sum(!is.na(subset_species_df$abundance[subset_species_df$taxon == species[j]])) < 2) {
      # if true, skip the loop for this species
      next
    }
    new_results <- data.frame()
    new_results <- complete(mice(subset_species_df, m = 5, maxit = 50, method = "pmm", seed = 123)) 
    mice_1_filled <- rbind(mice_1_filled, new_results)
    #write.table(Paride,paste0(iteration[i],"_", species[j],".xlsx"))
    cat("Paride is a turd and Fran will replace him because Tonda hates Paride",iteration[i])
  }
  
}

#2 gaps filled
mice_2_filled <- data.frame()
iteration= unique(final_df_2_missing_list$iteration)
n <- iteration[1]
for(n in iteration){
  subset_df<- final_df_2_missing_list[final_df_2_missing_list$iteration == n,]
  
  # Extract the unique species from the subseted dataframe
  species <- unique(subset_df$taxon)
  for(j in 1:length(species)){
    # Subset the data frame to include only the current species
    subset_species_df <- subset_df %>% filter(taxon == species[j])
    # check if the species has less than 2  values
    if(sum(!is.na(subset_species_df$abundance[subset_species_df$taxon == species[j]])) < 2) {
      # if true, skip the loop for this species
      next
    }
    new_results <- data.frame()
    new_results <- complete(mice(subset_species_df, m = 5, maxit = 50, method = "pmm", seed = 123)) 
    mice_2_filled <- rbind(mice_2_filled, new_results)
    #write.table(Paride,paste0(iteration[i],"_", species[j],".xlsx"))
    cat("Paride is a turd and Fran will replace him because Tonda hates Paride",iteration[i])
  }
  
}

#5 gaps filled
mice_5_filled <- data.frame()
iteration= unique(final_df_5_missing_list$iteration)
n <- iteration[1]
for(n in iteration){
  subset_df<- final_df_5_missing_list[final_df_5_missing_list$iteration == n,]
  
  # Extract the unique species from the subseted dataframe
  species <- unique(subset_df$taxon)
  for(j in 1:length(species)){
    # Subset the data frame to include only the current species
    subset_species_df <- subset_df %>% filter(taxon == species[j])
    # check if the species has less than 2  values
    if(sum(!is.na(subset_species_df$abundance[subset_species_df$taxon == species[j]])) < 2) {
      # if true, skip the loop for this species
      next
    }
    new_results <- data.frame()
    new_results <- complete(mice(subset_species_df, m = 5, maxit = 50, method = "pmm", seed = 123)) 
    mice_5_filled <- rbind(mice_5_filled, new_results)
    #write.table(Paride,paste0(iteration[i],"_", species[j],".xlsx"))
    cat("Paride is a turd and Fran will replace him because Tonda hates Paride",iteration[i])
  }
  
}

#10 gaps filled
mice_10_filled <- data.frame()
iteration= unique(final_df_10_missing_list$iteration)
n <- iteration[1]
for(n in iteration){
  subset_df<- final_df_10_missing_list[final_df_10_missing_list$iteration == n,]
  
  # Extract the unique species from the subseted dataframe
  species <- unique(subset_df$taxon)
  for(j in 1:length(species)){
    # Subset the data frame to include only the current species
    subset_species_df <- subset_df %>% filter(taxon == species[j])
    # check if the species has less than 2  values
    if(sum(!is.na(subset_species_df$abundance[subset_species_df$taxon == species[j]])) < 2) {
      # if true, skip the loop for this species
      next
    }
    new_results <- data.frame()
    new_results <- complete(mice(subset_species_df, m = 5, maxit = 50, method = "pmm", seed = 123)) 
    mice_10_filled <- rbind(mice_10_filled, new_results)
    #write.table(Paride,paste0(iteration[i],"_", species[j],".xlsx"))
    cat("Paride is a turd and Fran will replace him because Tonda hates Paride",iteration[i])
  }
  
}

#clean NAs
str(mice_1_filled)
mice_1_filled[is.na(mice_1_filled)] = 0 
str(mice_2_filled)
mice_2_filled[is.na(mice_2_filled)] = 0 
str(mice_5_filled)
mice_5_filled[is.na(mice_5_filled)] = 0 
str(mice_10_filled)
mice_10_filled[is.na(mice_10_filled)] = 0 



#back to wide format
str(mice_1_filled)
mice_1_filled$iteration <- as.character(mice_1_filled$iteration)
test1<- pivot_wider(data= mice_1_filled , id_cols =  c("year","iteration"),
                    names_from = "taxon", values_from = "abundance")
test12 <- test1 %>% select(3:142)
test12$Richness <- rowSums(test12 > 0)
test13 <- test1 %>% select(3:142)
test13$Total_Abundance <- apply(test13, 1, function(x) sum(x[which(x > 0)]))
final_df_1_filled<-cbind(test1$year,test1$iteration,test12$Richness, test13$Total_Abundance)
colnames(final_df_1_filled)<-c("year","iteration","richness","tot_abundance")
final_df_1_filled<-as.data.frame(final_df_1_filled)

str(mice_2_filled)
mice_2_filled$iteration <- as.character(mice_2_filled$iteration)
test2<- pivot_wider(data= mice_2_filled , id_cols =  c("year","iteration"),
                    names_from = "taxon", values_from = "abundance")
test22 <- test2 %>% select(3:142)
test22$Richness <- rowSums(test22 > 0)
test23 <- test2 %>% select(3:142)
test23$Total_Abundance <- apply(test23, 1, function(x) sum(x[which(x > 0)]))
final_df_2_filled<-cbind(test2$year,test2$iteration,test22$Richness, test23$Total_Abundance)
colnames(final_df_2_filled)<-c("year","iteration","richness","tot_abundance")
final_df_2_filled<-as.data.frame(final_df_2_filled)

str(mice_5_filled)
mice_5_filled$iteration <- as.character(mice_5_filled$iteration)
test5<- pivot_wider(data= mice_5_filled , id_cols =  c("year","iteration"),
                    names_from = "taxon", values_from = "abundance")
test52 <- test5 %>% select(3:142)
test52$Richness <- rowSums(test52 > 0)
test53 <- test5 %>% select(3:142)
test53$Total_Abundance <- apply(test53, 1, function(x) sum(x[which(x > 0)]))
final_df_5_filled<-cbind(test5$year,test5$iteration,test52$Richness, test53$Total_Abundance)
colnames(final_df_5_filled)<-c("year","iteration","richness","tot_abundance")
final_df_5_filled<-as.data.frame(final_df_5_filled)

str(mice_10_filled)
mice_10_filled$iteration <- as.character(mice_10_filled$iteration)
test10<- pivot_wider(data= mice_10_filled , id_cols =  c("year","iteration"),
                     names_from = "taxon", values_from = "abundance")
test102 <- test10 %>% select(3:142)
test102$Richness <- rowSums(test102 > 0)
test103 <- test10 %>% select(3:142)
test103$Total_Abundance <- apply(test103, 1, function(x) sum(x[which(x > 0)]))
final_df_10_filled<-cbind(test10$year,test10$iteration,test102$Richness, test103$Total_Abundance)
colnames(final_df_10_filled)<-c("year","iteration","richness","tot_abundance")
final_df_10_filled<-as.data.frame(final_df_10_filled)

#calculate annual species richness +SD
str(final_df_1_filled)
final_df_1_filled <- filter(final_df_1_filled, !is.na(richness))
final_df_1_filled$year<-as.numeric(final_df_1_filled$year)
final_df_1_filled$iteration<-as.numeric(final_df_1_filled$iteration)
final_df_1_filled$richness<-as.numeric(final_df_1_filled$richness)

final_df_2_filled <- filter(final_df_2_filled, !is.na(richness))
final_df_2_filled$year<-as.numeric(final_df_2_filled$year)
final_df_2_filled$iteration<-as.numeric(final_df_2_filled$iteration)
final_df_2_filled$richness<-as.numeric(final_df_2_filled$richness)

final_df_5_filled <- filter(final_df_5_filled, !is.na(richness))
final_df_5_filled$year<-as.numeric(final_df_5_filled$year)
final_df_5_filled$iteration<-as.numeric(final_df_5_filled$iteration)
final_df_5_filled$richness<-as.numeric(final_df_5_filled$richness)

final_df_10_filled <- filter(final_df_10_filled, !is.na(richness))
final_df_10_filled$year<-as.numeric(final_df_10_filled$year)
final_df_10_filled$iteration<-as.numeric(final_df_10_filled$iteration)
final_df_10_filled$richness<-as.numeric(final_df_10_filled$richness)

Richness_filled<-final_df_1_filled %>%
  dplyr::group_by(year) %>%  
  dplyr::summarize(mean_richness=mean(richness),
                   SD_richness=sd(richness))
Richness_filled
colnames(Richness_filled)<-c("year","mean_richness","sd_rich")
Richness_filled<-as.data.frame(Richness_filled)
plot10miss<-ggplot(Richness_filled, aes(x=year, y=mean_richness))+
  geom_point(color = "blue", size = 2)+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits = c(0, 50), breaks=seq(0, 50,5), labels = seq(0, 50,5))+
  geom_smooth()+
  theme_classic2()+
  theme_cleveland()
plot10miss

gam_filled = gam(mean_richness ~ s(year), 
                 data = Richness_filled, method="REML")
fit<-as.data.frame(gam_filled$fitted.values)
colnames(fit)[1]<-"annual_fit"
years <- seq(1991, 2016)
years <- data.frame(year = years)
fit1_filled<-cbind(years,fit)
#replotting after filling

#FILLED MODELS

#loop to fit every iteration individually and then averaging the fits
#1 filled
#MK trend
MK_1_fill <- as.data.frame(matrix(nrow=0,ncol=11))
columns=c("Corrected Zc","new P-value","N/N*","Original Z","old P.value","Tau","Sen s slope","old.variance","new.variance","S statistic","n")
colnames(MK_1_fill)<-columns

iteration<-unique(final_df_1_filled$iteration)
for(n in iteration){
  subset_df<- final_df_1_filled[final_df_1_filled$iteration == n,]
  MK<-My.mmkh(subset_df$richness)
  MK_1_fill<-rbind(MK_1_fill,MK)
  names(MK_1_fill) <- c("Corrected Zc","new P-value","N/N*","Original Z","old P.value","Tau","Sen s slope","old.variance","new.variance","S statistic","n")
  cat(".")
}
#gam
AIC_1_filled_gam <- as.data.frame(matrix(nrow=0,ncol=1))
columns="AIC"
colnames(AIC_1_filled_gam)<-columns

DEV_1_filled_gam <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(DEV_1_filled_gam)<-columns

RMSE_1_filled_gam <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(RMSE_1_filled_gam)<-columns



loop_1_filled <- as.data.frame(matrix(nrow=0,ncol=1)) 
columns="fit"
colnames(loop_1_missing)<-columns
iteration<-unique(final_df_1_filled$iteration)
for(n in iteration){
  subset_df<- final_df_1_filled[final_df_1_filled$iteration == n,]
  gam_1_filled<-gam(richness~ s(year), data =subset_df)
  fit<-gam_1_filled$fitted.values
  AIC<-AIC(gam_1_filled)
  AIC<-as.data.frame(AIC)
  colnames(AIC)[1]<-"AIC"
  AIC_1_filled_gam<-rbind(AIC_1_filled_gam,AIC)
  residuals <- residuals(gam_1_filled)
  mean_squared_error <- mean(residuals^2)
  rmse <- sqrt(mean_squared_error)
  RMSE_1_filled_gam<-rbind(RMSE_1_filled_gam,rmse)
  DEV<-((gam_1_filled$null.deviance - gam_1_filled$deviance)/gam_1_filled$null.deviance) * 100
  DEV<-as.data.frame(DEV)
  colnames(DEV)[1]<-"DEV"
  DEV_1_filled_gam<-rbind(DEV_1_filled_gam,DEV)
  fit<-as.data.frame(fit)
  colnames(fit)[1]<-"fit"
  loop_1_filled<-rbind(loop_1_filled,fit)
  #fit1miss_lm<-cbind(fit10miss$year,fit2)
  cat(".")
}
fit1_filled_gam<-cbind(final_df_1_filled$year,final_df_1_filled$iteration,loop_1_filled)
head(fit1_filled_gam)
colnames(fit1_filled_gam)<-c("year","iteration","fit")
head(fit1_filled_gam)

mean(DEV_1_filled_gam$DEV)
sd(DEV_1_filled_gam$DEV)

#lm
AIC_1_filled_lm <- as.data.frame(matrix(nrow=0,ncol=1))
columns="AIC"
colnames(AIC_1_filled_lm)<-columns

DEV_1_filled_lm <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(DEV_1_filled_lm)<-columns

RMSE_1_filled_lm <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(RMSE_1_filled_lm)<-columns



loop_1_filled <- as.data.frame(matrix(nrow=0,ncol=1)) 
columns="fit"
colnames(loop_1_missing)<-columns
iteration<-unique(final_df_1_filled$iteration)
for(n in iteration){
  subset_df<- final_df_1_filled[final_df_1_filled$iteration == n,]
  lm_1_filled<-lm(richness~ year, data =subset_df)
  fit<-lm_1_filled$fitted.values
  AIC<-AIC(lm_1_filled)
  AIC<-as.data.frame(AIC)
  colnames(AIC)[1]<-"AIC"
  AIC_1_filled_lm<-rbind(AIC_1_filled_lm,AIC)
  residuals <- residuals(lm_1_filled)
  mean_squared_error <- mean(residuals^2)
  rmse <- sqrt(mean_squared_error)
  RMSE_1_filled_lm<-rbind(RMSE_1_filled_lm,rmse)
  DEV<-deviance(lm_1_filled)/lm_1_filled$df.residual
  DEV<-as.data.frame(DEV)
  colnames(DEV)[1]<-"DEV"
  DEV_1_filled_lm<-rbind(DEV_1_filled_lm,DEV)
  fit<-as.data.frame(fit)
  colnames(fit)[1]<-"fit"
  loop_1_filled<-rbind(loop_1_filled,fit)
  #fit1miss_lm<-cbind(fit10miss$year,fit2)
  cat(".")
}
fit1_filled_lm<-cbind(final_df_1_filled$year,final_df_1_filled$iteration,loop_1_filled)
head(fit1_filled_lm)
colnames(fit1_filled_lm)<-c("year","iteration","fit")
head(fit1_filled_lm)

mean(DEV_1_filled_lm$DEV)
sd(DEV_1_filled_lm$DEV)

Annual<-fit1miss_lm %>%
  group_by(year) %>%  
  summarize(c(mean = mean(fit),SD=sd(fit)))
Annual
table(fit1miss_lm$fit)
#using fitted values with gaps
#lm
plot1fill_lm<-ggplot()+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits = c(20, 50), breaks=seq(20, 50,5), labels = seq(20, 50,5))+
  geom_smooth(data=fit1_filled_lm, aes(x=year, y=fit, group=iteration,colour =iteration),method = "loess", se=FALSE,size=0.25, alpha=0.5)+
  geom_smooth(data=original_full_fit_lm,aes(x=year, y=annual_fit),method = "loess",colour="red")+
  #geom_point(data=fit1_filled_lm,aes(x=year, y=fit, group=iteration),color = "blue", size = 3)+
  #geom_point(data=original_full_fit_lm,aes(x=year, y=annual_fit),color = "red", size = 3)+
  theme_classic2()+
  theme_cleveland()
plot1fill_lm

#gam
plot1fill_gam<-ggplot()+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits =c(25, 50), breaks=seq(25, 50,5), labels = seq(25, 50,5))+
  geom_smooth(data=fit1_filled_gam, aes(x=year, y=fit, group=iteration,colour =iteration),method = "gam", se=FALSE,size=0.25, alpha=0.5)+
  geom_smooth(data=original_full_fit,aes(x=year, y=annual_fit),method = "gam",se=FALSE,colour="red")+
  theme_classic2()+
  theme_cleveland()
plot1fill_gam

#raw data with loess and gaps
plot1fill_raw<-ggplot()+
  geom_smooth(data=final_df_1_gap, aes(x=year, y=richness, group=iteration),method = "loess",se=FALSE, size=0.25, alpha=0.5)+
  geom_point(data=final_df_1_gap,aes(x=year, y=richness, group=iteration),color = "blue", size = 3)+
  geom_smooth(data=df_no_gap,aes(x=year, y=richness),method = "loess", colour="red", se=FALSE)+
  geom_point(data=df_no_gap,aes(x=year, y=richness),color = "red", size = 1)+
  theme_classic2()+
  theme_cleveland()
plot1fill_raw

#2 filled
#MK trend
MK_2_fill <- as.data.frame(matrix(nrow=0,ncol=11))
columns=c("Corrected Zc","new P-value","N/N*","Original Z","old P.value","Tau","Sen s slope","old.variance","new.variance","S statistic","n")
colnames(MK_2_fill)<-columns

iteration<-unique(final_df_2_filled$iteration)
for(n in iteration){
  subset_df<- final_df_2_filled[final_df_2_filled$iteration == n,]
  MK<-My.mmkh(subset_df$richness)
  MK_2_fill<-rbind(MK_2_fill,MK)
  names(MK_2_fill) <- c("Corrected Zc","new P-value","N/N*","Original Z","old P.value","Tau","Sen s slope","old.variance","new.variance","S statistic","n")
  cat(".")
}
#gam
AIC_2_filled_gam <- as.data.frame(matrix(nrow=0,ncol=1))
columns="AIC"
colnames(AIC_2_filled_gam)<-columns

DEV_2_filled_gam <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(DEV_2_filled_gam)<-columns

RMSE_2_filled_gam <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(RMSE_2_filled_gam)<-columns



loop_2_filled <- as.data.frame(matrix(nrow=0,ncol=1)) 
columns="fit"
colnames(loop_2_filled)<-columns
iteration<-unique(final_df_2_filled$iteration)
for(n in iteration){
  subset_df<- final_df_2_filled[final_df_2_filled$iteration == n,]
  gam_2_filled<-gam(richness~ s(year), data =subset_df)
  fit<-gam_2_filled$fitted.values
  AIC<-AIC(gam_2_filled)
  AIC<-as.data.frame(AIC)
  colnames(AIC)[1]<-"AIC"
  AIC_2_filled_gam<-rbind(AIC_2_filled_gam,AIC)
  residuals <- residuals(gam_2_filled)
  mean_squared_error <- mean(residuals^2)
  rmse <- sqrt(mean_squared_error)
  RMSE_2_filled_gam<-rbind(RMSE_2_filled_gam,rmse)
  DEV<-((gam_2_filled$null.deviance - gam_2_filled$deviance)/gam_2_filled$null.deviance) * 100
  DEV<-as.data.frame(DEV)
  colnames(DEV)[1]<-"DEV"
  DEV_2_filled_gam<-rbind(DEV_2_filled_gam,DEV)
  fit<-as.data.frame(fit)
  colnames(fit)[1]<-"fit"
  loop_2_filled<-rbind(loop_2_filled,fit)
  #fit1miss_lm<-cbind(fit10miss$year,fit2)
  cat(".")
}

fit2_filled_gam<-cbind(final_df_2_filled$year,final_df_2_filled$iteration,loop_2_filled)
head(fit2_filled_gam)
colnames(fit2_filled_gam)<-c("year","iteration","fit")
head(fit2_filled_gam)
mean(DEV_2_filled_gam$DEV)
sd(DEV_2_filled_gam$DEV)

#lm
AIC_2_filled_lm <- as.data.frame(matrix(nrow=0,ncol=1))
columns="AIC"
colnames(AIC_2_filled_lm)<-columns

DEV_2_filled_lm <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(DEV_2_filled_lm)<-columns

RMSE_2_filled_lm <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(RMSE_2_filled_lm)<-columns



loop_2_filled <- as.data.frame(matrix(nrow=0,ncol=1)) 
columns="fit"
colnames(loop_2_filled)<-columns
iteration<-unique(final_df_2_filled$iteration)
for(n in iteration){
  subset_df<- final_df_2_filled[final_df_2_filled$iteration == n,]
  lm_2_filled<-lm(richness~ year, data =subset_df)
  fit<-lm_2_filled$fitted.values
  AIC<-AIC(lm_2_filled)
  AIC<-as.data.frame(AIC)
  colnames(AIC)[1]<-"AIC"
  AIC_2_filled_lm<-rbind(AIC_2_filled_lm,AIC)
  residuals <- residuals(lm_2_filled)
  mean_squared_error <- mean(residuals^2)
  rmse <- sqrt(mean_squared_error)
  RMSE_2_filled_lm<-rbind(RMSE_2_filled_lm,rmse)
  DEV<-deviance(lm_2_filled)/lm_2_filled$df.residual
  DEV<-as.data.frame(DEV)
  colnames(DEV)[1]<-"DEV"
  DEV_2_filled_lm<-rbind(DEV_2_filled_lm,DEV)
  fit<-as.data.frame(fit)
  colnames(fit)[1]<-"fit"
  loop_2_filled<-rbind(loop_2_filled,fit)
  #fit1miss_lm<-cbind(fit10miss$year,fit2)
  cat(".")
}
fit2_filled_lm<-cbind(final_df_2_filled$year,final_df_2_filled$iteration,loop_2_filled)
head(fit2_filled_lm)
colnames(fit2_filled_lm)<-c("year","iteration","fit")
head(fit2_filled_lm)
mean(DEV_2_filled_lm$DEV)
sd(DEV_2_filled_lm$DEV)

Annual<-fit2miss_lm %>%
  group_by(year) %>%  
  summarize(c(mean = mean(fit),SD=sd(fit)))
Annual
table(fit2miss_lm$fit)

#using fitted values with gaps
#lm
plot2fill_lm<-ggplot()+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits = c(20, 50), breaks=seq(20, 50,5), labels = seq(20, 50,5))+
  geom_smooth(data=fit2_filled_lm, aes(x=year, y=fit, group=iteration,colour =iteration),method = "loess", se=FALSE,size=0.25, alpha=0.5)+
  geom_smooth(data=original_full_fit_lm,aes(x=year, y=annual_fit),method = "loess",colour="red")+
  #geom_point(data=fit2_filled_lm,aes(x=year, y=fit, group=iteration),color = "blue", size = 3)+
  #geom_point(data=original_full_fit_lm,aes(x=year, y=annual_fit),color = "red", size = 3)+
  theme_classic2()+
  theme_cleveland()
plot2fill_lm

#gam
plot2fill_gam<-ggplot()+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits =c(25, 50), breaks=seq(25, 50,5), labels = seq(25, 50,5))+
  geom_smooth(data=fit2_filled_gam, aes(x=year, y=fit, group=iteration,colour =iteration),method = "gam", se=FALSE,size=0.25, alpha=0.5)+
  geom_smooth(data=original_full_fit,aes(x=year, y=annual_fit),method = "gam",se=FALSE,colour="red")+
  theme_classic2()+
  theme_cleveland()
plot2fill_gam

#raw data with loess and gaps
plot2fill_raw<-ggplot()+
  geom_smooth(data=final_df_2_gap, aes(x=year, y=richness, group=iteration),method = "loess",se=FALSE, size=0.25, alpha=0.5)+
  geom_point(data=final_df_2_gap,aes(x=year, y=richness, group=iteration),color = "blue", size = 3)+
  geom_smooth(data=df_no_gap,aes(x=year, y=richness),method = "loess", colour="red", se=FALSE)+
  geom_point(data=df_no_gap,aes(x=year, y=richness),color = "red", size = 1)+
  theme_classic2()+
  theme_cleveland()
plot2fill_raw

#5 filled
#MK trend
MK_5_fill <- as.data.frame(matrix(nrow=0,ncol=11))
columns=c("Corrected Zc","new P-value","N/N*","Original Z","old P.value","Tau","Sen s slope","old.variance","new.variance","S statistic","n")
colnames(MK_5_fill)<-columns

iteration<-unique(final_df_5_filled$iteration)
for(n in iteration){
  subset_df<- final_df_5_filled[final_df_5_filled$iteration == n,]
  MK<-My.mmkh(subset_df$richness)
  MK_5_fill<-rbind(MK_5_fill,MK)
  names(MK_5_fill) <- c("Corrected Zc","new P-value","N/N*","Original Z","old P.value","Tau","Sen s slope","old.variance","new.variance","S statistic","n")
  cat(".")
}

#gam
AIC_5_filled_gam <- as.data.frame(matrix(nrow=0,ncol=1))
columns="AIC"
colnames(AIC_5_filled_gam)<-columns

DEV_5_filled_gam <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(DEV_5_filled_gam)<-columns

RMSE_5_filled_gam <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(RMSE_5_filled_gam)<-columns



loop_5_filled <- as.data.frame(matrix(nrow=0,ncol=1)) 
columns="fit"
colnames(loop_5_filled)<-columns
iteration<-unique(final_df_5_filled$iteration)
for(n in iteration){
  subset_df<- final_df_5_filled[final_df_5_filled$iteration == n,]
  gam_5_filled<-gam(richness~ s(year), data =subset_df)
  fit<-gam_5_filled$fitted.values
  AIC<-AIC(gam_5_filled)
  AIC<-as.data.frame(AIC)
  colnames(AIC)[1]<-"AIC"
  AIC_5_filled_gam<-rbind(AIC_5_filled_gam,AIC)
  residuals <- residuals(gam_5_filled)
  mean_squared_error <- mean(residuals^2)
  rmse <- sqrt(mean_squared_error)
  RMSE_5_filled_gam<-rbind(RMSE_5_filled_gam,rmse)
  DEV<-((gam_5_filled$null.deviance - gam_5_filled$deviance)/gam_5_filled$null.deviance) * 100
  DEV<-as.data.frame(DEV)
  colnames(DEV)[1]<-"DEV"
  DEV_5_filled_gam<-rbind(DEV_5_filled_gam,DEV)
  fit<-as.data.frame(fit)
  colnames(fit)[1]<-"fit"
  loop_5_filled<-rbind(loop_5_filled,fit)
  #fit1miss_lm<-cbind(fit10miss$year,fit5)
  cat(".")
}
fit5_filled_gam<-cbind(final_df_5_filled$year,final_df_5_filled$iteration,loop_5_filled)
head(fit5_filled_gam)
colnames(fit5_filled_gam)<-c("year","iteration","fit")
head(fit5_filled_gam)
mean(DEV_5_filled_gam$DEV)
sd(DEV_5_filled_gam$DEV)

#lm
AIC_5_filled_lm <- as.data.frame(matrix(nrow=0,ncol=1))
columns="AIC"
colnames(AIC_5_filled_lm)<-columns

DEV_5_filled_lm <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(DEV_5_filled_lm)<-columns

RMSE_5_filled_lm <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(RMSE_5_filled_lm)<-columns



loop_5_filled <- as.data.frame(matrix(nrow=0,ncol=1)) 
columns="fit"
colnames(loop_5_filled)<-columns
iteration<-unique(final_df_5_filled$iteration)
for(n in iteration){
  subset_df<- final_df_5_filled[final_df_5_filled$iteration == n,]
  lm_5_filled<-lm(richness~ year, data =subset_df)
  fit<-lm_5_filled$fitted.values
  AIC<-AIC(lm_5_filled)
  AIC<-as.data.frame(AIC)
  colnames(AIC)[1]<-"AIC"
  AIC_5_filled_lm<-rbind(AIC_5_filled_lm,AIC)
  residuals <- residuals(lm_5_filled)
  mean_squared_error <- mean(residuals^2)
  rmse <- sqrt(mean_squared_error)
  RMSE_5_filled_lm<-rbind(RMSE_5_filled_lm,rmse)
  DEV<-deviance(lm_5_filled)/lm_5_filled$df.residual
  DEV<-as.data.frame(DEV)
  colnames(DEV)[1]<-"DEV"
  DEV_5_filled_lm<-rbind(DEV_5_filled_lm,DEV)
  fit<-as.data.frame(fit)
  colnames(fit)[1]<-"fit"
  loop_5_filled<-rbind(loop_5_filled,fit)
  #fit1miss_lm<-cbind(fit10miss$year,fit5)
  cat(".")
}
fit5_filled_lm<-cbind(final_df_5_filled$year,final_df_5_filled$iteration,loop_5_filled)
head(fit5_filled_lm)
colnames(fit5_filled_lm)<-c("year","iteration","fit")
head(fit5_filled_lm)
mean(DEV_5_filled_lm$DEV)
sd(DEV_5_filled_lm$DEV)

Annual<-fit5miss_lm %>%
  group_by(year) %>%  
  summarize(c(mean = mean(fit),SD=sd(fit)))
Annual
table(fit5miss_lm$fit)
#using fitted values with gaps
#lm
plot5fill_lm<-ggplot()+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits = c(20, 50), breaks=seq(20, 50,5), labels = seq(20, 50,5))+
  geom_smooth(data=fit5_filled_lm, aes(x=year, y=fit, group=iteration,colour =iteration),method = "loess", se=FALSE,size=0.25, alpha=0.5)+
  geom_smooth(data=original_full_fit_lm,aes(x=year, y=annual_fit),method = "loess",colour="red")+
  #geom_point(data=fit5_filled_lm,aes(x=year, y=fit, group=iteration),color = "blue", size = 3)+
  #geom_point(data=original_full_fit_lm,aes(x=year, y=annual_fit),color = "red", size = 3)+
  theme_classic2()+
  theme_cleveland()
plot5fill_lm

#gam
plot5fill_gam<-ggplot()+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits =c(25, 50), breaks=seq(25, 50,5), labels = seq(25, 50,5))+
  geom_smooth(data=fit5_filled_gam, aes(x=year, y=fit, group=iteration,colour =iteration),method = "gam", se=FALSE,size=0.25, alpha=0.5)+
  geom_smooth(data=original_full_fit,aes(x=year, y=annual_fit),method = "gam",se=FALSE,colour="red")+
  theme_classic2()+
  theme_cleveland()
plot5fill_gam

#raw data with loess and gaps
plot5fill_raw<-ggplot()+
  geom_smooth(data=final_df_5_gap, aes(x=year, y=richness, group=iteration),method = "loess",se=FALSE, size=0.25, alpha=0.5)+
  geom_point(data=final_df_5_gap,aes(x=year, y=richness, group=iteration),color = "blue", size = 3)+
  geom_smooth(data=df_no_gap,aes(x=year, y=richness),method = "loess", colour="red", se=FALSE)+
  geom_point(data=df_no_gap,aes(x=year, y=richness),color = "red", size = 1)+
  theme_classic2()+
  theme_cleveland()
plot5fill_raw

#10 filled
#MK trend
MK_10_fill <- as.data.frame(matrix(nrow=0,ncol=11))
columns=c("Corrected Zc","new P-value","N/N*","Original Z","old P.value","Tau","Sen s slope","old.variance","new.variance","S statistic","n")
colnames(MK_10_fill)<-columns

iteration<-unique(final_df_10_filled$iteration)
for(n in iteration){
  subset_df<- final_df_10_filled[final_df_10_filled$iteration == n,]
  MK<-My.mmkh(subset_df$richness)
  MK_10_fill<-rbind(MK_10_fill,MK)
  names(MK_10_fill) <- c("Corrected Zc","new P-value","N/N*","Original Z","old P.value","Tau","Sen s slope","old.variance","new.variance","S statistic","n")
  cat(".")
}

#gam
AIC_10_filled_gam <- as.data.frame(matrix(nrow=0,ncol=1))
columns="AIC"
colnames(AIC_10_filled_gam)<-columns

DEV_10_filled_gam <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(DEV_10_filled_gam)<-columns

RMSE_10_filled_gam <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(RMSE_10_filled_gam)<-columns



loop_10_filled <- as.data.frame(matrix(nrow=0,ncol=1)) 
columns="fit"
colnames(loop_10_filled)<-columns
iteration<-unique(final_df_10_filled$iteration)
for(n in iteration){
  subset_df<- final_df_10_filled[final_df_10_filled$iteration == n,]
  gam_10_filled<-gam(richness~ s(year), data =subset_df)
  fit<-gam_10_filled$fitted.values
  AIC<-AIC(gam_10_filled)
  AIC<-as.data.frame(AIC)
  colnames(AIC)[1]<-"AIC"
  AIC_10_filled_gam<-rbind(AIC_10_filled_gam,AIC)
  residuals <- residuals(gam_10_filled)
  mean_squared_error <- mean(residuals^2)
  rmse <- sqrt(mean_squared_error)
  RMSE_10_filled_gam<-rbind(RMSE_10_filled_gam,rmse)
  DEV<-((gam_10_filled$null.deviance - gam_10_filled$deviance)/gam_10_filled$null.deviance) * 100
  DEV<-as.data.frame(DEV)
  colnames(DEV)[1]<-"DEV"
  DEV_10_filled_gam<-rbind(DEV_10_filled_gam,DEV)
  fit<-as.data.frame(fit)
  colnames(fit)[1]<-"fit"
  loop_10_filled<-rbind(loop_10_filled,fit)
  #fit1miss_lm<-cbind(fit10miss$year,fit10)
  cat(".")
}
fit10_filled_gam<-cbind(final_df_10_filled$year,final_df_10_filled$iteration,loop_10_filled)
head(fit10_filled_gam)
colnames(fit10_filled_gam)<-c("year","iteration","fit")
head(fit10_filled_gam)
mean(DEV_10_filled_gam$DEV)
sd(DEV_10_filled_gam$DEV)

#lm
AIC_10_filled_lm <- as.data.frame(matrix(nrow=0,ncol=1))
columns="AIC"
colnames(AIC_10_filled_lm)<-columns

DEV_10_filled_lm <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(DEV_10_filled_lm)<-columns

RMSE_10_filled_lm <- as.data.frame(matrix(nrow=0,ncol=1))
columns="DEV"
colnames(RMSE_10_filled_lm)<-columns



loop_10_filled <- as.data.frame(matrix(nrow=0,ncol=1)) 
columns="fit"
colnames(loop_10_filled)<-columns
iteration<-unique(final_df_10_filled$iteration)
for(n in iteration){
  subset_df<- final_df_10_filled[final_df_10_filled$iteration == n,]
  lm_10_filled<-lm(richness~ year, data =subset_df)
  fit<-lm_10_filled$fitted.values
  AIC<-AIC(lm_10_filled)
  AIC<-as.data.frame(AIC)
  colnames(AIC)[1]<-"AIC"
  AIC_10_filled_lm<-rbind(AIC_10_filled_lm,AIC)
  residuals <- residuals(lm_10_filled)
  mean_squared_error <- mean(residuals^2)
  rmse <- sqrt(mean_squared_error)
  RMSE_10_filled_lm<-rbind(RMSE_10_filled_lm,rmse)
  DEV<-deviance(lm_10_filled)/lm_10_filled$df.residual
  DEV<-as.data.frame(DEV)
  colnames(DEV)[1]<-"DEV"
  DEV_10_filled_lm<-rbind(DEV_10_filled_lm,DEV)
  fit<-as.data.frame(fit)
  colnames(fit)[1]<-"fit"
  loop_10_filled<-rbind(loop_10_filled,fit)
  #fit1miss_lm<-cbind(fit10miss$year,fit10)
  cat(".")
}
fit10_filled_lm<-cbind(final_df_10_filled$year,final_df_10_filled$iteration,loop_10_filled)
head(fit10_filled_lm)
colnames(fit10_filled_lm)<-c("year","iteration","fit")
head(fit10_filled_lm)
mean(DEV_10_filled_lm$DEV)
sd(DEV_10_filled_lm$DEV)

Annual<-fit10miss_lm %>%
  group_by(year) %>%  
  summarize(c(mean = mean(fit),SD=sd(fit)))
Annual
table(fit10miss_lm$fit)
#using fitted values with gaps
#lm
plot10fill_lm<-ggplot()+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits = c(20, 50), breaks=seq(20, 50,5), labels = seq(20, 50,5))+
  geom_smooth(data=fit10_filled_lm, aes(x=year, y=fit, group=iteration,colour =iteration),method = "loess", se=FALSE,size=0.25, alpha=0.5)+
  geom_smooth(data=original_full_fit_lm,aes(x=year, y=annual_fit),method = "loess",colour="red")+
  #geom_point(data=fit10_filled_lm,aes(x=year, y=fit, group=iteration),color = "blue", size = 3)+
  #geom_point(data=original_full_fit_lm,aes(x=year, y=annual_fit),color = "red", size = 3)+
  theme_classic2()+
  theme_cleveland()
plot10fill_lm

#gam
plot10fill_gam<-ggplot()+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits =c(25, 50), breaks=seq(25, 50,5), labels = seq(25, 50,5))+
  geom_smooth(data=fit10_filled_gam, aes(x=year, y=fit, group=iteration,colour =iteration),method = "gam", se=FALSE,size=0.25, alpha=0.5)+
  geom_smooth(data=original_full_fit,aes(x=year, y=annual_fit),method = "gam",se=FALSE,colour="red")+
  theme_classic2()+
  theme_cleveland()
plot10fill_gam

#raw data with loess and gaps
plot10fill_raw<-ggplot()+
  geom_smooth(data=final_df_10_gap, aes(x=year, y=richness, group=iteration),method = "loess",se=FALSE, size=0.25, alpha=0.5)+
  geom_point(data=final_df_10_gap,aes(x=year, y=richness, group=iteration),color = "blue", size = 3)+
  geom_smooth(data=df_no_gap,aes(x=year, y=richness),method = "loess", colour="red", se=FALSE)+
  geom_point(data=df_no_gap,aes(x=year, y=richness),color = "red", size = 1)+
  theme_classic2()+
  theme_cleveland()
plot10fill_raw


#overall plot

#lm
plot_fill_overall_LM<-ggplot()+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits =c(25, 50), breaks=seq(25, 50,5), labels = seq(25, 50,5))+
  geom_smooth(data=fit10_filled_lm, aes(x=year, y=fit, group=iteration,colour =iteration),colour ="red",se=FALSE, size=0.25, alpha=0.5)+
  geom_smooth(data=fit5_filled_lm,aes(x=year, y=fit, group=iteration,colour =iteration),colour ="orange",se=FALSE, size=0.25, alpha=0.5)+
  geom_smooth(data=fit2_filled_lm, aes(x=year, y=fit, group=iteration,colour =iteration),colour ="yellow",se=FALSE, size=0.25, alpha=0.5)+
  geom_smooth(data=fit1_filled_lm, aes(x=year, y=fit, group=iteration,colour =iteration),colour ="green",se=FALSE, size=0.25, alpha=0.5)+
  geom_smooth(data=original_full_fit_lm,aes(x=year, y=annual_fit),colour="black",size=1)+
  theme_classic2()+
  theme_cleveland()
plot_fill_overall_LM

#gam
plot_fill_overall_GAM<-ggplot()+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits = c(25, 50), breaks=seq(25, 50,5), labels = seq(25, 50,5))+
  geom_smooth(data=fit10_filled_gam, aes(x=year, y=fit, group=iteration),colour ="red",se=FALSE, size=0.25, alpha=0.5)+
  geom_smooth(data=fit5_filled_gam, aes(x=year, y=fit, group=iteration),colour ="orange",se=FALSE, size=0.25, alpha=0.5)+
  geom_smooth(data=fit2_filled_gam, aes(x=year, y=fit, group=iteration),colour ="yellow",se=FALSE, size=0.25, alpha=0.5)+
  geom_smooth(data=fit1_filled_gam, aes(x=year, y=fit, group=iteration),colour ="green",se=FALSE, size=0.25, alpha=0.5)+
  geom_line(data=original_full_fit,aes(x=year, y=annual_fit),colour="black",size=1)+
  theme_classic2()+
  theme_cleveland()
plot_fill_overall_GAM

#raw data with loess and gaps
plot_fill_overall_raw<-ggplot()+
  scale_x_continuous(limits = c(1990, 2020), breaks=seq(1990,2020,5), labels = seq(1990,2020,5))+
  scale_y_continuous(limits = c(25, 50), breaks=seq(25, 50,5), labels = seq(25, 50,5))+
  geom_smooth(data=final_df_10_gap, aes(x=year, y=richness, group=iteration), colour = "red",method = "loess",se=FALSE, size=0.25, alpha=0.5)+  
  geom_smooth(data=final_df_5_gap, aes(x=year, y=richness, group=iteration), colour = "orange",method = "loess",se=FALSE, size=0.25, alpha=0.5)+  
  geom_smooth(data=final_df_2_gap, aes(x=year, y=richness, group=iteration), colour ="yellow",method = "loess",se=FALSE, size=0.25, alpha=0.5)+  
  geom_smooth(data=final_df_1_gap, aes(x=year, y=richness, group=iteration), colour="green",method = "loess",se=FALSE, size=0.25, alpha=0.5)+
  geom_smooth(data=df_no_gap,aes(x=year, y=richness),method = "loess", colour="black", se=FALSE)+
  geom_point(data=df_no_gap,aes(x=year, y=richness),color = "black", size = 1)+
  theme_classic2()+
  theme_cleveland()
plot_fill_overall_raw

library(cowplot)
plot3<-plot_grid(plot1fill_lm, plot1fill_gam,
                 plot2fill_lm, plot2fill_gam,
                 plot5fill_lm, plot5fill_gam,
                 plot10fill_lm, plot10fill_gam,
                 #labels = "AUTO", 
                 ncol = 2)
plot3
svg("plot_ALL_filled_NL_sep.svg")
plot3
dev.off()

plot4<-plot_grid(plot_fill_overall_LM,
                 plot_fill_overall_GAM,
                 plot_fill_overall_raw,
                 #labels = "AUTO", 
                 ncol = 1)
plot4
svg("plot_ALL_filled_NL.svg")
plot4
dev.off()


#Exporting Avg + SD of AICs                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
mean(AIC_1_filled_lm$AIC)
sd(AIC_1_filled_lm$AIC)
mean(AIC_2_filled_lm$AIC)
sd(AIC_2_filled_lm$AIC)
mean(AIC_5_filled_lm$AIC)
sd(AIC_5_filled_lm$AIC)
mean(AIC_10_filled_lm$AIC)
sd(AIC_10_filled_lm$AIC)
mean(AIC_1_filled_gam$AIC)
sd(AIC_1_filled_gam$AIC)
mean(AIC_2_filled_gam$AIC)
sd(AIC_2_filled_gam$AIC)
mean(AIC_5_filled_gam$AIC)
sd(AIC_5_filled_gam$AIC)
mean(AIC_10_filled_gam$AIC)
sd(AIC_10_filled_gam$AIC)

#averaging the richness values before AIC
#lm
Annual<-as.data.frame(final_df_1_gap %>%
                        group_by(year) %>%  
                        summarize(mean = mean(richness)))
avg1_miss<-lm(mean~year, data=Annual)
AIC(avg1_miss)
Annual<-as.data.frame(final_df_2_gap %>%
                        group_by(year) %>%  
                        summarize(mean = mean(richness)))
avg2_miss<-lm(mean~year, data=Annual)
AIC(avg2_miss)
Annual<-as.data.frame(final_df_5_gap %>%
                        group_by(year) %>%  
                        summarize(mean = mean(richness)))
avg5_miss<-lm(mean~year, data=Annual)
AIC(avg5_miss)
Annual<-as.data.frame(final_df_10_gap %>%
                        group_by(year) %>%  
                        summarize(mean = mean(richness)))
avg10_miss<-lm(mean~year, data=Annual)
AIC(avg10_miss)

AICs_missing_lm_NL<-rbind(AIC(avg1_miss),AIC(avg2_miss),AIC(avg5_miss),AIC(avg10_miss))
colnames(AICs_missing_lm_NL)<-c("AIC_lm_missing")
rownames(AICs_missing_lm_NL)<-c("1","2","5","10")

Annual<-as.data.frame(final_df_1_filled %>%
                        group_by(year) %>%  
                        summarize(mean = mean(richness)))
avg1_fill<-lm(mean~year, data=Annual)
AIC(avg1_fill)
Annual<-as.data.frame(final_df_2_filled %>%
                        group_by(year) %>%  
                        summarize(mean = mean(richness)))
avg2_fill<-lm(mean~year, data=Annual)
AIC(avg2_fill)
Annual<-as.data.frame(final_df_5_filled %>%
                        group_by(year) %>%  
                        summarize(mean = mean(richness)))
avg5_fill<-lm(mean~year, data=Annual)
AIC(avg5_fill)
Annual<-as.data.frame(final_df_10_filled %>%
                        group_by(year) %>%  
                        summarize(mean = mean(richness)))
avg10_fill<-lm(mean~year, data=Annual)
AIC(avg10_fill)

AICs_fill_lm_NL<-rbind(AIC(avg1_fill),AIC(avg2_fill),AIC(avg5_fill),AIC(avg10_fill))
colnames(AICs_fill_lm_NL)<-c("AIC_lm_filled")
rownames(AICs_fill_lm_NL)<-c("1","2","5","10")
#gam
Annual<-as.data.frame(final_df_1_gap %>%
                        group_by(year) %>%  
                        summarize(mean = mean(richness)))
avg1_miss<-gam(mean~s(year), data=Annual)
AIC(avg1_miss)
Annual<-as.data.frame(final_df_2_gap %>%
                        group_by(year) %>%  
                        summarize(mean = mean(richness)))
avg2_miss<-gam(mean~s(year), data=Annual)
AIC(avg2_miss)
Annual<-as.data.frame(final_df_5_gap %>%
                        group_by(year) %>%  
                        summarize(mean = mean(richness)))
avg5_miss<-gam(mean~s(year), data=Annual)
AIC(avg5_miss)
Annual<-as.data.frame(final_df_10_gap %>%
                        group_by(year) %>%  
                        summarize(mean = mean(richness)))
avg10_miss<-gam(mean~s(year), data=Annual)
AIC(avg10_miss)

AICs_missing_gam_NL<-rbind(AIC(avg1_miss),AIC(avg2_miss),AIC(avg5_miss),AIC(avg10_miss))
colnames(AICs_missing_gam_NL)<-c("AIC_gam_miss")
rownames(AICs_missing_gam_NL)<-c("1","2","5","10")

Annual<-as.data.frame(final_df_1_filled %>%
                        group_by(year) %>%  
                        summarize(mean = mean(richness)))
avg1_fill<-gam(mean~s(year), data=Annual)
AIC(avg1_fill)
Annual<-as.data.frame(final_df_2_filled %>%
                        group_by(year) %>%  
                        summarize(mean = mean(richness)))
avg2_fill<-gam(mean~s(year), data=Annual)
AIC(avg2_fill)
Annual<-as.data.frame(final_df_5_filled %>%
                        group_by(year) %>%  
                        summarize(mean = mean(richness)))
avg5_fill<-gam(mean~s(year), data=Annual)
AIC(avg5_fill)
Annual<-as.data.frame(final_df_10_filled %>%
                        group_by(year) %>%  
                        summarize(mean = mean(richness)))
avg10_fill<-gam(mean~s(year), data=Annual)
AIC(avg10_fill)

AICs_fill_gam_NL<-rbind(AIC(avg1_fill),AIC(avg2_fill),AIC(avg5_fill),AIC(avg10_fill))
colnames(AICs_fill_gam_NL)<-c("AIC_gam_fill")
rownames(AICs_fill_gam_NL)<-c("1","2","5","10")

AIC_averaged_NL<-cbind(AICs_missing_lm_NL,AICs_fill_lm_NL,AICs_missing_gam_NL,AICs_fill_gam_NL)

#MK stuff
colnames(MK_1_miss)[10]<-"S_stat"
colnames(MK_2_miss)[10]<-"S_stat"
colnames(MK_5_miss)[10]<-"S_stat"
colnames(MK_10_miss)[10]<-"S_stat"
colnames(MK_1_fill)[10]<-"S_stat"
colnames(MK_2_fill)[10]<-"S_stat"
colnames(MK_5_fill)[10]<-"S_stat"
colnames(MK_10_fill)[10]<-"S_stat"



#Meta regression plot
orig<-My.mmkh(Richness_0missing$mean_richness)
orig<-as.data.frame(t(orig))
colnames(orig)[10]<-"S_stat"
MK_orig <- rma(S_stat, old.variance, data = orig)


MK_trend_1_miss <- rma(S_stat, old.variance, data = MK_1_miss)
MK_trend_2_miss <- rma(S_stat, old.variance, data = MK_2_miss)
MK_trend_5_miss <- rma(S_stat, old.variance, data = MK_5_miss)
MK_trend_10_miss <- rma(S_stat, old.variance, data = MK_10_miss)
MK_trend_1_fill <- rma(S_stat, old.variance, data = MK_1_fill)
MK_trend_2_fill <- rma(S_stat, old.variance, data = MK_2_fill)
MK_trend_5_fill <- rma(S_stat, old.variance, data = MK_5_fill)
MK_trend_10_fill <- rma(S_stat, old.variance, data = MK_10_fill)

t0<-cbind(MK_orig$b,
          MK_orig$se,MK_orig$ci.lb,MK_orig$ci.ub)
t1<-cbind(MK_trend_1_miss$b,
          MK_trend_1_miss$se,MK_trend_1_miss$ci.lb,MK_trend_1_miss$ci.ub)
t2<-cbind(MK_trend_2_miss$b,
          MK_trend_2_miss$se,MK_trend_2_miss$ci.lb,MK_trend_2_miss$ci.ub)
t3<-cbind(MK_trend_5_miss$b,
          MK_trend_5_miss$se,MK_trend_5_miss$ci.lb,MK_trend_5_miss$ci.ub)
t4<-cbind(MK_trend_10_miss$b,
          MK_trend_10_miss$se,MK_trend_10_miss$ci.lb,MK_trend_10_miss$ci.ub)
q1<-cbind(MK_trend_1_fill$b,
          MK_trend_1_fill$se,MK_trend_1_fill$ci.lb,MK_trend_1_fill$ci.ub)
q2<-cbind(MK_trend_2_fill$b,
          MK_trend_2_fill$se,MK_trend_2_fill$ci.lb,MK_trend_2_fill$ci.ub)
q3<-cbind(MK_trend_5_fill$b,
          MK_trend_5_fill$se,MK_trend_5_fill$ci.lb,MK_trend_5_fill$ci.ub)
q4<-cbind(MK_trend_10_fill$b,
          MK_trend_10_fill$se,MK_trend_10_fill$ci.lb,MK_trend_10_fill$ci.ub)

pl_miss_NL<-rbind(t0,t1,t2,t3,t4)
rownames(pl_miss_NL)<-c("0","1","2","5","10")
pl_miss_NL <- cbind(rownames(pl_miss_NL), data.frame(pl_miss_NL, row.names=NULL))
colnames(pl_miss_NL)<-c("rep","estimate","se","ci.lb","ci.ub")
pl_miss_NL$rep<-as.numeric(pl_miss_NL$rep)
pl_fill_NL<-rbind(t0,q1,q2,q3,q4)
rownames(pl_fill_NL)<-c("0","1","2","5","10")
pl_fill_NL <- cbind(rownames(pl_fill_NL), data.frame(pl_fill_NL, row.names=NULL))
colnames(pl_fill_NL)<-c("rep","estimate","se","ci.lb","ci.ub")
pl_fill_NL$rep<-as.numeric(pl_fill_NL$rep)

MK_plot_miss<-ggplot(pl_miss_NL,aes(x=rep,y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=ci.lb,ymax=ci.ub,width=0.5))+
  theme_bw()+
  theme_cleveland()
MK_plot_miss

MK_plot_fill<-ggplot(pl_fill_NL,aes(x=rep,y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=ci.lb,ymax=ci.ub,width=0.5))+
  theme_bw()+
  theme_cleveland()
MK_plot_fill

plot_MK_trend_NL<-plot_grid(MK_plot_miss,MK_plot_fill,ncol=1)
plot_MK_trend_NL
svg("plot_MK_trends_both_NL.svg")
plot_MK_trend_NL
dev.off()

library(sjPlot)
p<-plot_model(MK_trend_1_miss, type = "est")
pMK_miss_1<-p +theme_sjplot()
p<-plot_model(MK_trend_2_miss, type = "est")
pMK_miss_2<-p +theme_sjplot()
p<-plot_model(MK_trend_5_miss, type = "est")
pMK_miss_5<-p +theme_sjplot()
p<-plot_model(MK_trend_10_miss, type = "est")
pMK_miss_10<-p +theme_sjplot()
detach("package:sjPlot", unload=TRUE)
library(cowplot)
plot_MK1_NL<-plot_grid(pMK_miss_1,pMK_miss_2,
                       pMK_miss_5,pMK_miss_10, ncol=2)
plot_MK1_NL
svg("plot_MK_trends_miss_NL.svg")
plot_MK1_NL
dev.off()

library(sjPlot)
p<-plot_model(MK_trend_1_fill, type = "est")
pMK_fill_1<-p +theme_sjplot()
p<-plot_model(MK_trend_2_fill, type = "est")
pMK_fill_2<-p +theme_sjplot()
p<-plot_model(MK_trend_5_fill, type = "est")
pMK_fill_5<-p +theme_sjplot()
p<-plot_model(MK_trend_10_fill, type = "est")
pMK_fill_10<-p +theme_sjplot()
detach("package:sjPlot", unload=TRUE)
library(cowplot)
plot_MK1_NL<-plot_grid(pMK_fill_1,pMK_fill_2,
                       pMK_fill_5,pMK_fill_10, ncol=2)
plot_MK1_NL
svg("plot_MK_trends_filled_NL.svg")
plot_MK1_NL
dev.off()

#RMSE
colnames(RMSE_1_miss_gam)[1]<-"value"
colnames(RMSE_2_miss_gam)[1]<-"value"
colnames(RMSE_5_miss_gam)[1]<-"value"
colnames(RMSE_10_miss_gam)[1]<-"value"
colnames(RMSE_1_miss_lm)[1]<-"value"
colnames(RMSE_2_miss_lm)[1]<-"value"
colnames(RMSE_5_miss_lm)[1]<-"value"
colnames(RMSE_10_miss_lm)[1]<-"value"

colnames(RMSE_1_filled_gam)[1]<-"value"
colnames(RMSE_2_filled_gam)[1]<-"value"
colnames(RMSE_5_filled_gam)[1]<-"value"
colnames(RMSE_10_filled_gam)[1]<-"value"
colnames(RMSE_1_filled_lm)[1]<-"value"
colnames(RMSE_2_filled_lm)[1]<-"value"
colnames(RMSE_5_filled_lm)[1]<-"value"
colnames(RMSE_10_filled_lm)[1]<-"value"

mean_1_miss_gam<-mean(RMSE_1_miss_gam$value)
se_1_miss_gam<-sd(RMSE_1_miss_gam$value)
mean_2_miss_gam<-mean(RMSE_2_miss_gam$value)
se_2_miss_gam<-sd(RMSE_2_miss_gam$value)
mean_5_miss_gam<-mean(RMSE_5_miss_gam$value)
se_5_miss_gam<-sd(RMSE_5_miss_gam$value)
mean_10_miss_gam<-mean(RMSE_10_miss_gam$value)
se_10_miss_gam<-sd(RMSE_10_miss_gam$value)
mean_1_miss_lm<-mean(RMSE_1_miss_lm$value)
se_1_miss_lm<-sd(RMSE_1_miss_lm$value)
mean_2_miss_lm<-mean(RMSE_2_miss_lm$value)
se_2_miss_lm<-sd(RMSE_2_miss_lm$value)
mean_5_miss_lm<-mean(RMSE_5_miss_lm$value)
se_5_miss_lm<-sd(RMSE_5_miss_lm$value)
mean_10_miss_lm<-mean(RMSE_10_miss_lm$value)
se_10_miss_lm<-sd(RMSE_10_miss_lm$value)

mean_1_fill_gam<-mean(RMSE_1_filled_gam$value)
se_1_fill_gam<-sd(RMSE_1_filled_gam$value)
mean_2_fill_gam<-mean(RMSE_2_filled_gam$value)
se_2_fill_gam<-sd(RMSE_2_filled_gam$value)
mean_5_fill_gam<-mean(RMSE_5_filled_gam$value)
se_5_fill_gam<-sd(RMSE_5_filled_gam$value)
mean_10_fill_gam<-mean(RMSE_10_filled_gam$value)
se_10_fill_gam<-sd(RMSE_10_filled_gam$value)
mean_1_fill_lm<-mean(RMSE_1_filled_lm$value)
se_1_fill_lm<-sd(RMSE_1_filled_lm$value)
mean_2_fill_lm<-mean(RMSE_2_filled_lm$value)
se_2_fill_lm<-sd(RMSE_2_filled_lm$value)
mean_5_fill_lm<-mean(RMSE_5_filled_lm$value)
se_5_fill_lm<-sd(RMSE_5_filled_lm$value)
mean_10_fill_lm<-mean(RMSE_10_filled_lm$value)
se_10_fill_lm<-sd(RMSE_10_filled_lm$value)


t0<-rbind(rmse_orig_lm, mean_1_miss_lm,mean_2_miss_lm,mean_5_miss_lm,mean_10_miss_lm)
t1<-rbind(rmse_orig_lm, se_1_miss_lm,se_2_miss_lm,se_5_miss_lm,se_10_miss_lm)
RMSE_lm<-cbind(t0,t1)
rownames(RMSE_lm)<-c("0","1","2","5","10")
RMSE_lm_miss <- cbind(rownames(RMSE_lm), data.frame(RMSE_lm, row.names=NULL))
colnames(RMSE_lm_miss)<-c("rep","mean","se")

t0<-rbind(rmse_orig_gam, mean_1_miss_gam,mean_2_miss_gam,mean_5_miss_gam,mean_10_miss_gam)
t1<-rbind(rmse_orig_gam, se_1_miss_gam,se_2_miss_gam,se_5_miss_gam,se_10_miss_gam)
RMSE_gam<-cbind(t0,t1)
rownames(RMSE_gam)<-c("0","1","2","5","10")
RMSE_gam_miss <- cbind(rownames(RMSE_gam), data.frame(RMSE_gam, row.names=NULL))
colnames(RMSE_gam_miss)<-c("rep","mean","se")

t0<-rbind(rmse_orig_lm, mean_1_fill_lm,mean_2_fill_lm,mean_5_fill_lm,mean_10_fill_lm)
t1<-rbind(rmse_orig_lm, se_1_fill_lm,se_2_fill_lm,se_5_fill_lm,se_10_fill_lm)
RMSE_lm<-cbind(t0,t1)
rownames(RMSE_lm)<-c("0","1","2","5","10")
RMSE_lm_fill <- cbind(rownames(RMSE_lm), data.frame(RMSE_lm, row.names=NULL))
colnames(RMSE_lm_fill)<-c("rep","mean","se")

t0<-rbind(rmse_orig_gam, mean_1_fill_gam,mean_2_fill_gam,mean_5_fill_gam,mean_10_fill_gam)
t1<-rbind(rmse_orig_gam, se_1_fill_gam,se_2_fill_gam,se_5_fill_gam,se_10_fill_gam)
RMSE_gam<-cbind(t0,t1)
rownames(RMSE_gam)<-c("0","1","2","5","10")
RMSE_gam_fill <- cbind(rownames(RMSE_gam), data.frame(RMSE_gam, row.names=NULL))
colnames(RMSE_gam_fill)<-c("rep","mean","se")

RMSE_plot_gam_miss<-ggplot(RMSE_gam_miss,aes(x=rep,y=mean))+
  geom_point()+
  geom_errorbar(aes(ymin=mean+se,ymax=mean-se,width=0.5))+
  theme_bw()+
  theme_cleveland()
RMSE_plot_gam_miss

RMSE_plot_gam_fill<-ggplot(RMSE_gam_fill,aes(x=rep,y=mean))+
  geom_point()+
  geom_errorbar(aes(ymin=mean+se,ymax=mean-se,width=0.5))+
  theme_bw()+
  theme_cleveland()
RMSE_plot_gam_fill

RMSE_plot_lm_miss<-ggplot(RMSE_lm_miss,aes(x=rep,y=mean))+
  geom_point()+
  geom_errorbar(aes(ymin=mean+se,ymax=mean-se,width=0.5))+
  theme_bw()+
  theme_cleveland()
RMSE_plot_lm_miss

RMSE_plot_lm_fill<-ggplot(RMSE_lm_fill,aes(x=rep,y=mean))+
  geom_point()+
  geom_errorbar(aes(ymin=mean+se,ymax=mean-se,width=0.5))+
  theme_bw()+
  theme_cleveland()
RMSE_plot_lm_fill

plot_RMSE_trend_NL<-plot_grid(RMSE_plot_gam_miss,RMSE_plot_gam_fill,
                              RMSE_plot_lm_miss,RMSE_plot_lm_fill,ncol=2)
plot_RMSE_trend_NL
svg("plot_RMSE_trends_both_NL.svg")
plot_RMSE_trend_NL
dev.off()
