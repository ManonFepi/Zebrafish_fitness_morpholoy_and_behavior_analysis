##########Script for "Environmentally relevant exposure to PFAS during zebrafish development impacts behavior, fitness, and the transcriptome" ##########

####Call the librairies####
library(tibble)
library(DESeq2)
library("writexl")
library(ggplot2)
library(openxlsx)
library(forcats)
library(datasets)
library(RColorBrewer)
library(FSA)
library(dplyr)
library(ggpubr)
library(nortest)
library(lme4)
library(arm)
library(stringr)
library("survival")
library("survminer")
library(factoextra)
library("multcompView")
library(ggsignif)
library(tidyverse)
library(plyr)
library(reshape2)
library(grid)
library(gridExtra)
library(ape)
library(scales)
library(data.table) # different syntax to merge data
install.packages("readxl")
library("readxl")
library(car)
library(pheatmap)
library(dendextend)
library(VennDiagram)
#-------------------------
####Set working directory####
setwd("~/Papers/EpiTrans F0/R analysis")
#-------------------------
#####Sup.Fig.3-FET####
#Import data
surv<-read.xlsx("FET kaplan.xlsx")

#A. Survival curve
survival <- survfit(Surv(Timepoint, Statut) ~ Condition, data = surv, type ="kaplan-meier")
ggsurvplot(survival, conf.int = TRUE, risk.table = FALSE, ylim=c(0.45,1), pval = TRUE,pval.coord = c(20,0.7), data = surv,xlab = "Hours post fertilization",ylab = "Good health probability",legend=c(0.2, 0.2),palette = c("cornflowerblue","darkgreen","green","#FFCC33","#FFFF00"))

#B.Forest plot
CoxMod<-coxph(Surv(Timepoint,Statut)~Condition,data=surv)
ggforest(CoxMod)
#-------------------------
####Sup.Fig.4.A Spawning success#####
#Import data
spawn=read.xlsx("Spawning success_29072022.xlsx")
spawn$Condition <- fct_relevel(spawn$Condition, c("Control", "PFOS_LC", "PFOS_HC", "PFBS_LC", "PFBS_HC"))

Control<-subset(spawn,Condition=="Control")
PFOS_LC<-subset(spawn,Condition=="PFOS_LC")
PFOS_HC<-subset(spawn,Condition=="PFOS_HC")
PFBS_LC<-subset(spawn,Condition=="PFBS_LC")
PFBS_HC<-subset(spawn,Condition=="PFBS_HC")

shapiro.test(Control$Percentage_good_repro)
shapiro.test(PFOS_LC$Percentage_good_repro)
shapiro.test(PFOS_HC$Percentage_good_repro)
shapiro.test(PFBS_LC$Percentage_good_repro)
shapiro.test(PFBS_HC$Percentage_good_repro)

t.test(Control$Percentage_good_repro, PFOS_LC$Percentage_good_repro)
t.test(Control$Percentage_good_repro, PFOS_HC$Percentage_good_repro) 
t.test(Control$Percentage_good_repro, PFBS_LC$Percentage_good_repro) 
t.test(Control$Percentage_good_repro, PFBS_HC$Percentage_good_repro) 

my_comparisons <- list( c("Control", "PFOS_LC"), c("Control", "PFOS_HC"), c("Control", "PFBS_LC"), c("Control", "PFBS_HC"))

Sup_Fig_4_A <- ggplot(spawn, aes(x = Condition, y = Percentage_good_repro, fill=Condition)) +
  geom_boxplot(width=0.3) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  theme_light() + 
  scale_fill_manual(values = c("cornflowerblue","#FFFF00","#FFCC33","green","darkgreen")) + 
  #ggtitle("Boxplot of the total number of eggs produced by reproduction (n=30)") + 
  ylab("Spawning sucess (%)") + #  stat_compare_means(comparisons = my_comparisons, method = "t.test") +
  theme(axis.text.x = element_text(size=18,angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  geom_pwc(method = "t.test", label = "{p.format}", hide.ns = F, tip.length = 0, p.adjust.method = "none", ref.group = "Control", label.size = 6)
Sup_Fig_4_A
#-------------------------
####Fig. 1.A Eggs production####
#Import data
eggs <- read.xlsx("Eggs production.xlsx")
#Eggs_total
eggs$Condition <- fct_relevel(eggs$Condition, c("Control", "PFOS_LC", "PFOS_HC", "PFBS_LC", "PFBS_HC"))
head(eggs)

hist(eggs$Eggs_total,prob=T,main="Eggs_total")
points(density(eggs$Eggs_total),type="l",col="red", lwd=2)
shapiro.test(eggs$Eggs_total)

#Subset to see difference between each condition according to sex
Control<-subset(eggs,Condition=="Control")
PFOS_LC<-subset(eggs,Condition=="PFOS_LC")
PFOS_HC<-subset(eggs,Condition=="PFOS_HC")
PFBS_LC<-subset(eggs,Condition=="PFBS_LC")
PFBS_HC<-subset(eggs,Condition=="PFBS_HC")

shapiro.test(Control$Eggs_total)
shapiro.test(PFOS_LC$Eggs_total)
shapiro.test(PFOS_HC$Eggs_total)
shapiro.test(PFBS_LC$Eggs_total)
shapiro.test(PFBS_HC$Eggs_total)

wilcox.test(Control$Eggs_total, PFOS_LC$Eggs_total)
wilcox.test(Control$Eggs_total, PFOS_HC$Eggs_total) 
wilcox.test(Control$Eggs_total, PFBS_LC$Eggs_total) 
wilcox.test(Control$Eggs_total, PFBS_HC$Eggs_total) 

Fig_1_A <- ggplot(eggs, aes(x = Condition, y = Eggs_total, fill=Condition)) +
  geom_boxplot(width=0.3) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  theme_light() + 
  scale_fill_manual(values = c("cornflowerblue","#FFFF00","#FFCC33","green","darkgreen")) + 
  #ggtitle("Boxplot of the total number of eggs produced by reproduction (n=30)") + 
  ylab("Number of eggs") +
  theme(axis.text.x = element_text(size=18,angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
geom_pwc(method = "wilcox.test", label = "{p.format} {p.signif}", hide.ns = TRUE, tip.length = 0, p.adjust.method = "none", ref.group = "Control", label.size = 6)
Fig_1_A
#-------------------------
####Sup.Fig.4.B Eggs quality####
#Eggs_good
ggdensity(eggs$Eggs_good, 
          main = "Density plot of good eggs",
          xlab = "Number of good eggs")
shapiro.test(eggs$Eggs_good)

wilcox.test(Control$Eggs_good, PFOS_LC$Eggs_good)
wilcox.test(Control$Eggs_good, PFOS_HC$Eggs_good) 
wilcox.test(Control$Eggs_good, PFBS_LC$Eggs_good) 
wilcox.test(Control$Eggs_good, PFBS_HC$Eggs_good) 

Nb_Good <- ggplot(eggs, aes(x = Condition, y = Eggs_good, fill=Condition)) +
  geom_boxplot(width=0.3) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  theme_light() + 
  scale_fill_manual(values = c("cornflowerblue","#FFFF00","#FFCC33","green","darkgreen")) + 
  #ggtitle("Boxplot of the total number of eggs produced by reproduction (n=30)") + 
  ylab("Number of eggs") + 
  theme(axis.text.x = element_text(size=18,angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  geom_pwc(method = "wilcox.test", label = "p = {p.format} {p.signif}", hide.ns = TRUE, tip.length = 0, p.adjust.method = "none", ref.group = "Control", label.size = 6)
Nb_Good

#Percentage_Eggs_good
ggdensity(eggs$Percentage_good, 
          main = "Density plot of good eggs",
          xlab = "Number of good eggs")
shapiro.test(eggs$Percentage_good)

wilcox.test(Control$Percentage_good, PFOS_LC$Percentage_good)
wilcox.test(Control$Percentage_good, PFOS_HC$Percentage_good) 
wilcox.test(Control$Percentage_good, PFBS_LC$Percentage_good) 
wilcox.test(Control$Percentage_good, PFBS_HC$Percentage_good) 

Sup_Fig_4_B <- ggplot(eggs, aes(x = Condition, y = Percentage_good, fill=Condition)) +
  geom_boxplot(width=0.3) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  theme_light() + 
  scale_fill_manual(values = c("cornflowerblue","#FFFF00","#FFCC33","green","darkgreen")) + 
  #ggtitle("Boxplot of the total number of eggs produced by reproduction (n=30)") + 
  ylab("Percentage of good eggs") + 
  theme(axis.text.x = element_text(size=18,angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  geom_pwc(method = "wilcox.test", label = "{p.format}", hide.ns = F, tip.length = 0, p.adjust.method = "none", ref.group = "Control", label.size = 6)
Sup_Fig_4_B
#-------------------------
####Sup.Fig.4.C BCF global####
rm(list = ls())
BodyCondition<-read.table("F0_onlyPairs_BodyConditionOrgans.txt",h=T)
BodyCondition$Condition <- fct_relevel(BodyCondition$Condition, c("Control", "PFOS_LC", "PFOS_HC", "PFBS_LC", "PFBS_HC"))

#Subset to see difference between each condition according to sex
Control<-subset(BodyCondition,Condition=="Control")
PFOS_LC<-subset(BodyCondition,Condition=="PFOS_LC")
PFOS_HC<-subset(BodyCondition,Condition=="PFOS_HC")
PFBS_LC<-subset(BodyCondition,Condition=="PFBS_LC")
PFBS_HC<-subset(BodyCondition,Condition=="PFBS_HC")

#Body Condition

hist(Control$Condition_Factor_K)
sd(Control$Condition_Factor_K)
shapiro.test(Control$Condition_Factor_K)
ggdensity(Control$Condition_Factor_K, 
          main = "Density plot of Total_distance moved (cm)",
          xlab = "Total distance moved cm")

hist(PFOS_HC$Condition_Factor_K,prob=T,main="Condition_Factor_K")
points(density(PFOS_HC$Condition_Factor_K),type="l",col="red", lwd=2)

shapiro.test(Control$Condition_Factor_K)
shapiro.test(PFOS_LC$Condition_Factor_K)
shapiro.test(PFOS_HC$Condition_Factor_K)
shapiro.test(PFBS_LC$Condition_Factor_K)
shapiro.test(PFBS_HC$Condition_Factor_K)

t.test(Control$Condition_Factor_K, PFOS_LC$Condition_Factor_K)
t.test(Control$Condition_Factor_K, PFOS_HC$Condition_Factor_K) 
t.test(Control$Condition_Factor_K, PFBS_LC$Condition_Factor_K) 
t.test(Control$Condition_Factor_K, PFBS_HC$Condition_Factor_K) 

my_comparisons <- list( c("Control", "PFOS_HC"), c("Control", "PFBS_HC"))

Sup_Fig_4_C <- ggplot(BodyCondition, aes(x = Condition, y = Condition_Factor_K, fill=Condition)) +
  geom_boxplot(width=0.3) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  theme_light() + 
  scale_fill_manual(values = c("cornflowerblue","#FFFF00","#FFCC33","green","darkgreen")) + 
  #ggtitle("Boxplot of the total number of eggs produced by reproduction (n=30)") + 
  ylab("BCF") + geom_hline(yintercept = 1, linetype = 2) +
  theme(axis.text.x = element_text(size=18,angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  geom_pwc(method = "t.test", label = "{p.format}", hide.ns = F, tip.length = 0, p.adjust.method = "none", ref.group = "Control", label.size = 6)
Sup_Fig_4_C
#-------------------------
####Fig. 1.B BCF males####
###Per conditions in males
males<-subset(BodyCondition,Sex=="Male")
summary(males)

males$Condition <- fct_relevel(males$Condition, c("Control", "PFOS_LC", "PFOS_HC", "PFBS_LC", "PFBS_HC"))

#subset males for pairwise test
malesC<-subset(males,Condition=="Control") #To select only the females
malesPFOSLC<-subset(males,Condition=="PFOS_LC") #To select only the females
malesPFOSHC<-subset(males,Condition=="PFOS_HC") #To select only the females
malesPFBSLC<-subset(males,Condition=="PFBS_LC") #To select only the females
malesPFBSHC<-subset(males,Condition=="PFBS_HC") #To select only the females

hist(males$Condition_Factor_K)
sd(males$Condition_Factor_K)
summary(males)
shapiro.test(rnorm(1000, mean = 1.0042  , sd = 0.1366105))

hist(PFOS_HC$Condition_Factor_K,prob=T,main="Distance moved")
points(density(PFOS_HC$Condition_Factor_K),type="l",col="red", lwd=2)

#pairwise test for significances
t.test(malesC$Condition_Factor_K, malesPFOSLC$Condition_Factor_K) 
t.test(malesC$Condition_Factor_K, malesPFOSHC$Condition_Factor_K) 
t.test(malesC$Condition_Factor_K, malesPFBSLC$Condition_Factor_K) 
t.test(malesC$Condition_Factor_K, malesPFBSHC$Condition_Factor_K) 

compare_means(Condition_Factor_K ~ Condition,  data = males, method = "t.test")

#Condition_Factor_K Control PFOS HC p = 0.0219
#Condition_Factor_K Control PFBS LC p = 0.0170

my_comparisons <- list( c("Control", "PFOS_HC"), c("Control", "PFBS_LC"))
Fig_2_B <- ggplot(males, aes(x = Condition, y = Condition_Factor_K, fill=Condition)) +
  geom_boxplot(width=0.3) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  theme_light() + 
  scale_fill_manual(values = c("cornflowerblue","#FFFF00","#FFCC33","green","darkgreen")) + 
  #ggtitle("Boxplot of the total number of eggs produced by reproduction (n=30)") + 
  ylab("BCF in males") + geom_hline(yintercept = 1, linetype = 2) +
  theme(axis.text.x = element_text(size=18,angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  geom_pwc(method = "t.test", label = "{p.format} {p.signif}", hide.ns = TRUE, tip.length = 0, p.adjust.method = "none", ref.group = "Control", label.size = 6)
Fig_2_B
#-------------------------
####Fig. 1.C BCF males vs. females####
#By sex
# Now make a interaction between two factors  
# on x axis   
BodyCondition$Factor1Factor2 <- interaction(BodyCondition$Condition,BodyCondition$Sex)  

my_comparisons <- list( c("Control.Female", "Control.Male"), c("PFOS_LC.Female", "PFOS_LC.Male"), c("PFOS_HC.Female", "PFOS_HC.Male"), c("PFBS_LC.Female", "PFBS_LC.Male"), c("PFBS_HC.Female","PFBS_HC.Male"))

# now Plot Boxplot with fill color according 
# to factor1 and factor2 
BodyCondition$Factor1Factor2 <- fct_relevel(BodyCondition$Factor1Factor2, c("Control.Female", "Control.Male","PFOS_LC.Female","PFOS_LC.Male","PFOS_HC.Female","PFOS_HC.Male","PFBS_LC.Female","PFBS_LC.Male","PFBS_HC.Female","PFBS_HC.Male"))
Fig_2_C<-ggplot(aes(y = Condition_Factor_K , x = Factor1Factor2), data = BodyCondition) +  
  geom_boxplot(aes(fill=Factor1Factor2))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  theme_light() +  stat_compare_means(comparisons = my_comparisons, label = "p.format", label.y = 1.60,size=5,method = "t.test",hide.ns = FALSE)+
  scale_fill_manual(values = c("cornflowerblue","cornflowerblue","#FFFF00","#FFFF00","#FFCC33","#FFCC33","green","green","darkgreen","darkgreen"))+
  ylab("BCF") + geom_hline(yintercept = 1, linetype = 2) + xlab("Condition & Sex")+
  theme(axis.text.x = element_text(size=18,angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))
Fig_2_C
#-------------------------
####Sup.Fig.4.D & E gonad somatic index females and males####
#Gonad somatic Index
###Per conditions in females
females<-subset(BodyCondition,Sex=="Female")
boxplot(females$Gonade_somatic)$out
Q <- quantile(females$Gonade_somatic, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(females$Gonade_somatic)
Gonade_somatic_without_outliers_females<- subset(females, females$Gonade_somatic > (Q[1] - 1.5*iqr) & females$Gonade_somatic < (Q[2]+1.5*iqr))
boxplot(Gonade_somatic_without_outliers_females$Gonade_somatic)

hist(Gonade_somatic_without_outliers_females$Gonade_somatic)
sd(Gonade_somatic_without_outliers_females$Gonade_somatic)
summary(Gonade_somatic_without_outliers_females)
shapiro.test(rnorm(1000, mean = 14.986  , sd = 3.397044))

hist(Gonade_somatic_without_outliers_females$Gonade_somatic,prob=T,main="Gonade somatic",ylim = c(0,0.15))
points(density(Gonade_somatic_without_outliers_females$Gonade_somatic),type="l",col="red", lwd=2)

femalesC<-subset(BodyCondition,Condition=="Control")
femalesPFOSLC<-subset(BodyCondition,Condition=="PFOS_LC")
femalesPFOSHC<-subset(BodyCondition,Condition=="PFOS_HC")
femalesPFBSLC<-subset(BodyCondition,Condition=="PFBS_LC")
femalesPFBSHC<-subset(BodyCondition,Condition=="PFBS_HC")

#pairwise test for significances
t.test(femalesC$Gonade_somatic, femalesPFOSLC$Gonade_somatic) 
t.test(femalesC$Gonade_somatic, femalesPFOSHC$Gonade_somatic) 
t.test(femalesC$Gonade_somatic, femalesPFBSLC$Gonade_somatic) 
t.test(femalesC$Gonade_somatic, femalesPFBSHC$Gonade_somatic) 

compare_means(Gonade_somatic ~ Condition,  data = Gonade_somatic_without_outliers_females, method = "t.test")

Gonade_somatic_females <- ggplot(Gonade_somatic_without_outliers_females, aes(x = Condition, y = Gonade_somatic, fill=Condition)) +
  geom_boxplot(width=0.3) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  theme_light() + 
  scale_fill_manual(values = c("cornflowerblue","#FFFF00","#FFCC33","green","darkgreen")) + 
  #ggtitle("F0: Gonade_somatic in females") + 
  ylab("Gonade somatic in females") + 
  theme(axis.text.x = element_text(size=18,angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  geom_pwc(method = "t.test", label = "{p.format}", hide.ns = F, tip.length = 0, p.adjust.method = "none", ref.group = "Control", label.size = 6)
Gonade_somatic_females

###Per conditions in males
boxplot(males$Gonade_somatic)$out
Q <- quantile(males$Gonade_somatic, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(males$Gonade_somatic)
Gonade_somatic_without_outliers<- subset(males, males$Gonade_somatic > (Q[1] - 1.5*iqr) & males$Gonade_somatic < (Q[2]+1.5*iqr))
boxplot(Gonade_somatic_without_outliers$Gonade_somatic)

hist(Gonade_somatic_without_outliers$Gonade_somatic)
sd(Gonade_somatic_without_outliers$Gonade_somatic)
summary(Gonade_somatic_without_outliers)
shapiro.test(rnorm(1000, mean = 1.1098  , sd = 0.3870416))

hist(Gonade_somatic_without_outliers$Gonade_somatic,prob=T,main="Gonade_somatic",ylim=c(0,1.5))
points(density(Gonade_somatic_without_outliers$Gonade_somatic),type="l",col="red", lwd=2)

#pairwise test for significances
t.test(malesC$Gonade_somatic, malesPFOSLC$Gonade_somatic) 
t.test(malesC$Gonade_somatic, malesPFOSHC$Gonade_somatic) 
t.test(malesC$Gonade_somatic, malesPFBSLC$Gonade_somatic) 
t.test(malesC$Gonade_somatic, malesPFBSHC$Gonade_somatic) 

compare_means(Gonade_somatic ~ Condition,  data = Gonade_somatic_without_outliers, method = "t.test")

Gonade_somatic_males <- ggplot(Gonade_somatic_without_outliers, aes(x = Condition, y = Gonade_somatic, fill=Condition)) +
  geom_boxplot(width=0.3) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  theme_light() + 
  scale_fill_manual(values = c("cornflowerblue","#FFFF00","#FFCC33","green","darkgreen")) + 
  #ggtitle("F0: Gonade_somatic in males") + 
  ylab("Gonade somatic in males") + 
  theme(axis.text.x = element_text(size=18,angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  geom_pwc(method = "t.test", label = "{p.format}", hide.ns = F, tip.length = 0, p.adjust.method = "none", ref.group = "Control", label.size = 6)
Gonade_somatic_males
#-------------------------
####Fig. 1. D Liver males####
#Liver somatic Index
###Per conditions in males
hist(males$Liver_somatic)
sd(males$Liver_somatic,na.rm = "TRUE")
Liver_somatic_males_wt_NA<-na.omit(males$Liver_somatic)
summary(males)
shapiro.test(rnorm(1000, mean = 0.5313  , sd = 0.2103225))

hist(males$Liver_somatic,prob=T,main="Liver_somatic")
points(density(Liver_somatic_males_wt_NA),type="l",col="red", lwd=2)

#pairwise test for significances
t.test(malesC$Liver_somatic, malesPFOSLC$Liver_somatic) 
t.test(malesC$Liver_somatic, malesPFOSHC$Liver_somatic) 
t.test(malesC$Liver_somatic, malesPFBSLC$Liver_somatic) 
t.test(malesC$Liver_somatic, malesPFBSHC$Liver_somatic) 

compare_means(Liver_somatic ~ Condition,  data = males, method = "t.test")

my_comparisons <- list(c("Control", "PFBS_HC"))

Fig_1_D <- ggplot(males, aes(x = Condition, y = Liver_somatic, fill=Condition)) +
  geom_boxplot(width=0.3) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  theme_light() + 
  scale_fill_manual(values = c("cornflowerblue","#FFFF00","#FFCC33","green","darkgreen")) + 
  #ggtitle("F0: Total Liver_somatic in males") + 
  ylab("Liver somatic in males") + # stat_compare_means(comparisons = my_comparisons, method = "t.test")+
  theme(axis.text.x = element_text(size=18,angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  geom_pwc(method = "t.test", label = "{p.format} {p.signif}", hide.ns = TRUE, tip.length = 0, p.adjust.method = "none", ref.group = "Control", label.size = 6)
Fig_1_D
#-------------------------
####Sup.Fig.4.F####
#Brain somatic Index
###Per conditions in females
hist(females$Brain_somatic)
sd(females$Brain_somatic,na.rm="TRUE")
summary(females)
shapiro.test(rnorm(1000, mean = 0.8265, sd = 0.3163376))

hist(females$Brain_somatic,prob=T,main="Brain_somatic")
Brain_somatic_females_wt_NA<-na.omit(females$Brain_somatic)
points(density(Brain_somatic_females_wt_NA),type="l",col="red", lwd=2)

#pairwise test for significances
t.test(femalesC$Brain_somatic, femalesPFOSLC$Brain_somatic) 
t.test(femalesC$Brain_somatic, femalesPFOSHC$Brain_somatic) 
t.test(femalesC$Brain_somatic, femalesPFBSLC$Brain_somatic) 
t.test(femalesC$Brain_somatic, femalesPFBSHC$Brain_somatic) 

compare_means(Brain_somatic ~ Condition,  data = females, method = "t.test")

my_comparisons <- list( c("Control", "PFOS_HC"), c("Control", "PFBS_HC"))

Fig_4_F <- ggplot(females, aes(x = Condition, y = Brain_somatic, fill=Condition)) +
  geom_boxplot(width=0.3) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  theme_light() + 
  scale_fill_manual(values = c("cornflowerblue","#FFFF00","#FFCC33","green","darkgreen")) + 
  #ggtitle("F0: Brain_somatic in females") + 
  ylab("Brain somatic in females") + 
  theme(axis.text.x = element_text(size=18,angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  geom_pwc(method = "t.test", label = "{p.format}", hide.ns = F, tip.length = 0, p.adjust.method = "none", ref.group = "Control", label.size = 6)
Fig_4_F
#-------------------------
####Fig. 1. E####
#By sex
# Now make a interaction between two factors  
# on x axis   
my_comparisons <- list( c("Control.Female", "Control.Male"), c("PFOS_LC.Female", "PFOS_LC.Male"), c("PFOS_HC.Female", "PFOS_HC.Male"), c("PFBS_LC.Female", "PFBS_LC.Male"), c("PFBS_HC.Female","PFBS_HC.Male"))

# now Plot Boxplot with fill color according 
# to factor1 and factor2 
BodyCondition$Factor1Factor2 <- fct_relevel(BodyCondition$Factor1Factor2, c("Control.Female", "Control.Male","PFOS_LC.Female","PFOS_LC.Male","PFOS_HC.Female","PFOS_HC.Male","PFBS_LC.Female","PFBS_LC.Male","PFBS_HC.Female","PFBS_HC.Male"))
Fig_1_E<-ggplot(aes(y = Brain_somatic , x = Factor1Factor2), data = BodyCondition) +  
  geom_boxplot(aes(fill=Factor1Factor2))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  theme_light() + 
  scale_fill_manual(values = c("cornflowerblue","cornflowerblue","#FFFF00","#FFFF00","#FFCC33","#FFCC33","green","green","darkgreen","darkgreen"))+
  ylab("Brain somatic") + geom_hline(yintercept = 1, linetype = 2) + xlab("Condition & Sex")+
  theme(axis.text.x = element_text(size=18,angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.format", label.y = 1.60,size=5,method = "t.test",hide.ns = FALSE)
Fig_1_E
#-------------------------
####Sup Fig. 5. A Larvae behavioral test checking####
rm(list = ls())

#Set the working directory
B1<- read.table("StatsLight120ALL.txt",h=T)
summary(B1)

# Remove all rows with too high values due to instrument error 
B<-B1[!(B1$Well=="F8" & B1$Plate=="Trial_1"),]
B<-B[!(B$Well=="E1" & B$Plate=="Trial_1"),]
B<-B[!(B$Well=="A6" & B$Plate=="Trial_1"),]
B<-B[!(B$Well=="B8" & B$Plate=="Trial_1"),]
B<-B[!(B$Well=="C8" & B$Plate=="Trial_1"),]
B<-B[!(B$Well=="D10" & B$Plate=="Trial_1"),]
B<-B[!(B$Well=="A12" & B$Plate=="Trial_1"),]
B<-B[!(B$Well=="H2" & B$Plate=="Trial_1"),]

B<-B[!(B$Well=="F9" & B$Plate=="Trial_2"),]
B<-B[!(B$Well=="E1" & B$Plate=="Trial_2"),]
B<-B[!(B$Well=="D10" & B$Plate=="Trial_2"),]
B<-B[!(B$Well=="A9" & B$Plate=="Trial_2"),]
B<-B[!(B$Well=="E2" & B$Plate=="Trial_2"),]
B<-B[!(B$Well=="F1" & B$Plate=="Trial_2"),]
B<-B[!(B$Well=="A5" & B$Plate=="Trial_2"),]
B<-B[!(B$Well=="B6" & B$Plate=="Trial_2"),]
B<-B[!(B$Well=="C8" & B$Plate=="Trial_2"),]
B<-B[!(B$Well=="E10" & B$Plate=="Trial_2"),]
B<-B[!(B$Well=="C11" & B$Plate=="Trial_2"),]

B<-B[!(B$Well=="C11" & B$Plate=="Trial_3"),]
B<-B[!(B$Well=="H5" & B$Plate=="Trial_3"),]
B<-B[!(B$Well=="H10" & B$Plate=="Trial_3"),]
B<-B[!(B$Well=="D11" & B$Plate=="Trial_3"),]
B<-B[!(B$Well=="A12" & B$Plate=="Trial_3"),]
B<-B[!(B$Well=="E12" & B$Plate=="Trial_3"),]
B<-B[!(B$Well=="E6" & B$Plate=="Trial_3"),]
B<-B[!(B$Well=="D9" & B$Plate=="Trial_3"),]
B<-B[!(B$Well=="C10" & B$Plate=="Trial_3"),]
B<-B[!(B$Well=="B11" & B$Plate=="Trial_3"),]
B<-B[!(B$Well=="C12" & B$Plate=="Trial_3"),]

B<-B[!(B$Well=="E9" & B$Plate=="Trial_4"),]
B<-B[!(B$Well=="D6" & B$Plate=="Trial_4"),]
B<-B[!(B$Well=="G8" & B$Plate=="Trial_4"),]
B<-B[!(B$Well=="F11" & B$Plate=="Trial_4"),]
B<-B[!(B$Well=="C1" & B$Plate=="Trial_4"),]
B<-B[!(B$Well=="D1" & B$Plate=="Trial_4"),]
B<-B[!(B$Well=="B8" & B$Plate=="Trial_4"),]
B<-B[!(B$Well=="D12" & B$Plate=="Trial_4"),]

B<-B[!(B$Well=="F1" & B$Plate=="Trial_5"),]
B<-B[!(B$Well=="G12" & B$Plate=="Trial_5"),]
B<-B[!(B$Well=="C10" & B$Plate=="Trial_5"),]
B<-B[!(B$Well=="C6" & B$Plate=="Trial_5"),]
B<-B[!(B$Well=="H8" & B$Plate=="Trial_5"),]
B<-B[!(B$Well=="F8" & B$Plate=="Trial_5"),]
B<-B[!(B$Well=="C11" & B$Plate=="Trial_5"),]

B<-B[!(B$Well=="E12" & B$Plate=="Trial_6"),]
B<-B[!(B$Well=="G6" & B$Plate=="Trial_6"),]
B<-B[!(B$Well=="B12" & B$Plate=="Trial_6"),]
B<-B[!(B$Well=="C10" & B$Plate=="Trial_6"),]
B<-B[!(B$Well=="D9" & B$Plate=="Trial_6"),]
B<-B[!(B$Well=="C11" & B$Plate=="Trial_6"),]
B<-B[!(B$Well=="C2" & B$Plate=="Trial_6"),]
B<-B[!(B$Well=="H2" & B$Plate=="Trial_6"),]
B<-B[!(B$Well=="A6" & B$Plate=="Trial_6"),]
B<-B[!(B$Well=="E7" & B$Plate=="Trial_6"),]

B<-B[!(B$Well=="F7" & B$Plate=="Trial_7"),]
B<-B[!(B$Well=="A6" & B$Plate=="Trial_7"),]

B<-B[!(B$Well=="B2" & B$Plate=="Trial_8"),]
B<-B[!(B$Well=="C1" & B$Plate=="Trial_8"),]
B<-B[!(B$Well=="E6" & B$Plate=="Trial_8"),]
B<-B[!(B$Well=="F11" & B$Plate=="Trial_8"),]
B<-B[!(B$Well=="G11" & B$Plate=="Trial_8"),]
B<-B[!(B$Well=="F1" & B$Plate=="Trial_8"),]
B<-B[!(B$Well=="G1" & B$Plate=="Trial_8"),]

B<-B[!(B$Well=="H10" & B$Plate=="Trial_9"),]
B<-B[!(B$Well=="D12" & B$Plate=="Trial_9"),]
B<-B[!(B$Well=="E6" & B$Plate=="Trial_9"),]
B<-B[!(B$Well=="G10" & B$Plate=="Trial_9"),]
B<-B[!(B$Well=="F1" & B$Plate=="Trial_9"),]
B<-B[!(B$Well=="A2" & B$Plate=="Trial_9"),]
B<-B[!(B$Well=="B5" & B$Plate=="Trial_9"),]
B<-B[!(B$Well=="F5" & B$Plate=="Trial_9"),]
B<-B[!(B$Well=="C8" & B$Plate=="Trial_9"),]
B<-B[!(B$Well=="D8" & B$Plate=="Trial_9"),]
B<-B[!(B$Well=="G9" & B$Plate=="Trial_9"),]
B<-B[!(B$Well=="D11" & B$Plate=="Trial_9"),]
B<-B[!(B$Well=="F12" & B$Plate=="Trial_9"),]

B<-B[!(B$Well=="F1" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="G7" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="C10" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="A7" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="F8" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="A11" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="H1" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="F2" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="G2" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="A6" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="C7" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="E7" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="A8" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="D8" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="B9" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="C11" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="D11" & B$Plate=="Trial_10"),]

B<-B[!(B$Well=="A6" & B$Plate=="Trial_11"),]
B<-B[!(B$Well=="H9" & B$Plate=="Trial_11"),]
B<-B[!(B$Well=="A8" & B$Plate=="Trial_11"),]
B<-B[!(B$Well=="H6" & B$Plate=="Trial_11"),]
B<-B[!(B$Well=="D7" & B$Plate=="Trial_11"),]
B<-B[!(B$Well=="E12" & B$Plate=="Trial_11"),]

B<-B[!(B$Well=="F10" & B$Plate=="Trial_12"),]
B<-B[!(B$Well=="F7" & B$Plate=="Trial_12"),]
B<-B[!(B$Well=="B1" & B$Plate=="Trial_12"),]
B<-B[!(B$Well=="C1" & B$Plate=="Trial_12"),]
B<-B[!(B$Well=="C2" & B$Plate=="Trial_12"),]
B<-B[!(B$Well=="E6" & B$Plate=="Trial_12"),]
B<-B[!(B$Well=="H10" & B$Plate=="Trial_12"),]
B<-B[!(B$Well=="C12" & B$Plate=="Trial_12"),]

B<-B[!(B$Well=="B6" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="A2" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="A1" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="B2" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="A5" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="C5" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="E5" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="G5" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="C7" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="D7" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="B8" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="C8" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="A12" & B$Plate=="Trial_13"),]

B<-B[!(B$Well=="H9" & B$Plate=="Trial_14"),]
B<-B[!(B$Well=="B9" & B$Plate=="Trial_14"),]
B<-B[!(B$Well=="B5" & B$Plate=="Trial_14"),]
B<-B[!(B$Well=="C5" & B$Plate=="Trial_14"),]
B<-B[!(B$Well=="D7" & B$Plate=="Trial_14"),]
B<-B[!(B$Well=="D9" & B$Plate=="Trial_14"),]
B<-B[!(B$Well=="G9" & B$Plate=="Trial_14"),]
B<-B[!(B$Well=="B10" & B$Plate=="Trial_14"),]
B<-B[!(B$Well=="C11" & B$Plate=="Trial_14"),]
B<-B[!(B$Well=="E11" & B$Plate=="Trial_14"),]
B<-B[!(B$Well=="F12" & B$Plate=="Trial_14"),]

B<-B[!(B$Well=="B11" & B$Plate=="Trial_15"),]
B<-B[!(B$Well=="G1" & B$Plate=="Trial_15"),]
B<-B[!(B$Well=="D8" & B$Plate=="Trial_15"),]
B<-B[!(B$Well=="D5" & B$Plate=="Trial_15"),]
B<-B[!(B$Well=="H6" & B$Plate=="Trial_15"),]
B<-B[!(B$Well=="A10" & B$Plate=="Trial_15"),]
B<-B[!(B$Well=="D10" & B$Plate=="Trial_15"),]
B<-B[!(B$Well=="D11" & B$Plate=="Trial_15"),]
B<-B[!(B$Well=="G11" & B$Plate=="Trial_15"),]
B<-B[!(B$Well=="G12" & B$Plate=="Trial_15"),]

B<-B[!(B$Well=="E6" & B$Plate=="Trial_16"),]
B<-B[!(B$Well=="B10" & B$Plate=="Trial_16"),]
B<-B[!(B$Well=="F10" & B$Plate=="Trial_16"),]
B<-B[!(B$Well=="E1" & B$Plate=="Trial_16"),]
B<-B[!(B$Well=="E5" & B$Plate=="Trial_16"),]
B<-B[!(B$Well=="C7" & B$Plate=="Trial_16"),]
B<-B[!(B$Well=="E8" & B$Plate=="Trial_16"),]
B<-B[!(B$Well=="D10" & B$Plate=="Trial_16"),]
B<-B[!(B$Well=="E11" & B$Plate=="Trial_16"),]
B<-B[!(B$Well=="D12" & B$Plate=="Trial_16"),]

B<-B[!(B$Well=="G9" & B$Plate=="Trial_17"),]
B<-B[!(B$Well=="F10" & B$Plate=="Trial_17"),]
B<-B[!(B$Well=="H11" & B$Plate=="Trial_17"),]
B<-B[!(B$Well=="E12" & B$Plate=="Trial_17"),]
B<-B[!(B$Well=="D7" & B$Plate=="Trial_17"),]
B<-B[!(B$Well=="E11" & B$Plate=="Trial_17"),]
B<-B[!(B$Well=="F11" & B$Plate=="Trial_17"),]
B<-B[!(B$Well=="F6" & B$Plate=="Trial_17"),]
B<-B[!(B$Well=="F7" & B$Plate=="Trial_17"),]
B<-B[!(B$Well=="D8" & B$Plate=="Trial_17"),]
B<-B[!(B$Well=="E8" & B$Plate=="Trial_17"),]
B<-B[!(B$Well=="D9" & B$Plate=="Trial_17"),]

B<-B[!(B$Well=="F8" & B$Plate=="Trial_18"),]
B<-B[!(B$Well=="C1" & B$Plate=="Trial_18"),]
B<-B[!(B$Well=="E10" & B$Plate=="Trial_18"),]
B<-B[!(B$Well=="H9" & B$Plate=="Trial_18"),]
B<-B[!(B$Well=="B9" & B$Plate=="Trial_18"),]
B<-B[!(B$Well=="D12" & B$Plate=="Trial_18"),]
B<-B[!(B$Well=="H11" & B$Plate=="Trial_18"),]

B<-B[!(B$Well=="E7" & B$Plate=="Trial_19"),]
B<-B[!(B$Well=="C11" & B$Plate=="Trial_19"),]
B<-B[!(B$Well=="A11" & B$Plate=="Trial_19"),]
B<-B[!(B$Well=="A12" & B$Plate=="Trial_19"),]
B<-B[!(B$Well=="D6" & B$Plate=="Trial_19"),]
B<-B[!(B$Well=="F12" & B$Plate=="Trial_19"),]
B<-B[!(B$Well=="H12" & B$Plate=="Trial_19"),]
B<-B[!(B$Well=="F7" & B$Plate=="Trial_19"),]
B<-B[!(B$Well=="G8" & B$Plate=="Trial_19"),]
B<-B[!(B$Well=="B7" & B$Plate=="Trial_19"),]
B<-B[!(B$Well=="D1" & B$Plate=="Trial_19"),]
B<-B[!(B$Well=="E12" & B$Plate=="Trial_19"),]

B<-B[!(B$Well=="C12" & B$Plate=="Trial_20"),]
B<-B[!(B$Well=="G8" & B$Plate=="Trial_20"),]
B<-B[!(B$Well=="G7" & B$Plate=="Trial_20"),]
B<-B[!(B$Well=="H9" & B$Plate=="Trial_20"),]
B<-B[!(B$Well=="G10" & B$Plate=="Trial_20"),]
B<-B[!(B$Well=="E10" & B$Plate=="Trial_20"),]
B<-B[!(B$Well=="D6" & B$Plate=="Trial_20"),]
B<-B[!(B$Well=="B9" & B$Plate=="Trial_20"),]
B<-B[!(B$Well=="F9" & B$Plate=="Trial_20"),]
B<-B[!(B$Well=="A8" & B$Plate=="Trial_20"),]
B<-B[!(B$Well=="B2" & B$Plate=="Trial_20"),]
B<-B[!(B$Well=="A12" & B$Plate=="Trial_20"),]

B1<-B
summary(B1) # check if all values are numeric and if range of values is ok

LON1<-B1[B1$Condition=="LON_1",]
LOFF<-B1[B1$Condition=="LOFF",]
LON2<-B1[B1$Condition=="LON_2",]

B1$Condition <- fct_relevel(B1$Condition, c("LON_1", "LOFF", "LON_2"))
B1$Expo <- fct_relevel(B1$Expo, c("Control", "PFOS_LC", "PFOS_HC", "PFBS_LC", "PFBS_HC"))

levels(B1$Condition)
levels(B1$Expo)

Sup_Fig_5_A <- ggplot(B1, aes(x = Expo, y = Dmcptot, fill=Condition)) + scale_y_continuous(limits=c(0,180))+
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)+
  theme_light() + 
  scale_fill_manual(values = c("yellow", "darkgrey","yellow")) + 
  xlab("Condition") + ylab("Total distance moved (cm)") + 
  theme(axis.text.x = element_text(size=18,angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))
Sup_Fig_5_A + guides(fill=guide_legend(title="Light exposure"))
#ggsave("bp.png", width = 16, height = 9, bp)

LOFF$Expo <- factor(LOFF$Expo,levels = c("Control", "PFOS_LC", "PFOS_HC", "PFBS_LC","PFBS_HC"))

# Tapping test
B2<- read.table("StatsTap120ALLN.txt",h=T)
summary (B2)

# remove all values, where embryos were not hatched, coagulated, or deformed
B<-B2[!(B2$Well=="F8" & B2$Plate=="Trial_1"),]
B<-B[!(B$Well=="E1" & B$Plate=="Trial_1"),]
B<-B[!(B$Well=="A6" & B$Plate=="Trial_1"),]
B<-B[!(B$Well=="B8" & B$Plate=="Trial_1"),]
B<-B[!(B$Well=="C8" & B$Plate=="Trial_1"),]
B<-B[!(B$Well=="D10" & B$Plate=="Trial_1"),]
B<-B[!(B$Well=="A12" & B$Plate=="Trial_1"),]

B<-B[!(B$Well=="F9" & B$Plate=="Trial_2"),]
B<-B[!(B$Well=="E1" & B$Plate=="Trial_2"),]
B<-B[!(B$Well=="D10" & B$Plate=="Trial_2"),]
B<-B[!(B$Well=="A9" & B$Plate=="Trial_2"),]
B<-B[!(B$Well=="E2" & B$Plate=="Trial_2"),]
B<-B[!(B$Well=="F1" & B$Plate=="Trial_2"),]
B<-B[!(B$Well=="A5" & B$Plate=="Trial_2"),]
B<-B[!(B$Well=="B6" & B$Plate=="Trial_2"),]
B<-B[!(B$Well=="C8" & B$Plate=="Trial_2"),]
B<-B[!(B$Well=="E10" & B$Plate=="Trial_2"),]
B<-B[!(B$Well=="C11" & B$Plate=="Trial_2"),]

B<-B[!(B$Well=="C11" & B$Plate=="Trial_3"),]
B<-B[!(B$Well=="H5" & B$Plate=="Trial_3"),]
B<-B[!(B$Well=="H10" & B$Plate=="Trial_3"),]
B<-B[!(B$Well=="D11" & B$Plate=="Trial_3"),]
B<-B[!(B$Well=="A12" & B$Plate=="Trial_3"),]
B<-B[!(B$Well=="E12" & B$Plate=="Trial_3"),]
B<-B[!(B$Well=="E6" & B$Plate=="Trial_3"),]
B<-B[!(B$Well=="D9" & B$Plate=="Trial_3"),]
B<-B[!(B$Well=="C10" & B$Plate=="Trial_3"),]
B<-B[!(B$Well=="B11" & B$Plate=="Trial_3"),]
B<-B[!(B$Well=="C12" & B$Plate=="Trial_3"),]

B<-B[!(B$Well=="E9" & B$Plate=="Trial_4"),]
B<-B[!(B$Well=="D6" & B$Plate=="Trial_4"),]
B<-B[!(B$Well=="G8" & B$Plate=="Trial_4"),]
B<-B[!(B$Well=="F11" & B$Plate=="Trial_4"),]
B<-B[!(B$Well=="E9" & B$Plate=="Trial_4"),]
B<-B[!(B$Well=="C1" & B$Plate=="Trial_4"),]
B<-B[!(B$Well=="D1" & B$Plate=="Trial_4"),]
B<-B[!(B$Well=="B8" & B$Plate=="Trial_4"),]
B<-B[!(B$Well=="D12" & B$Plate=="Trial_4"),]
B<-B[!(B$Well=="E6" & B$Plate=="Trial_4"),]
B<-B[!(B$Well=="D11" & B$Plate=="Trial_4"),]
B<-B[!(B$Well=="B12" & B$Plate=="Trial_4"),]

B<-B[!(B$Well=="F1" & B$Plate=="Trial_5"),]
B<-B[!(B$Well=="G12" & B$Plate=="Trial_5"),]
B<-B[!(B$Well=="C10" & B$Plate=="Trial_5"),]
B<-B[!(B$Well=="C6" & B$Plate=="Trial_5"),]
B<-B[!(B$Well=="H8" & B$Plate=="Trial_5"),]
B<-B[!(B$Well=="F8" & B$Plate=="Trial_5"),]
B<-B[!(B$Well=="C11" & B$Plate=="Trial_5"),]
B<-B[!(B$Well=="C2" & B$Plate=="Trial_5"),]
B<-B[!(B$Well=="G6" & B$Plate=="Trial_5"),]
B<-B[!(B$Well=="G9" & B$Plate=="Trial_5"),]


B<-B[!(B$Well=="E12" & B$Plate=="Trial_6"),]
B<-B[!(B$Well=="G6" & B$Plate=="Trial_6"),]
B<-B[!(B$Well=="B12" & B$Plate=="Trial_6"),]
B<-B[!(B$Well=="C10" & B$Plate=="Trial_6"),]
B<-B[!(B$Well=="D9" & B$Plate=="Trial_6"),]
B<-B[!(B$Well=="C11" & B$Plate=="Trial_6"),]
B<-B[!(B$Well=="C2" & B$Plate=="Trial_6"),]
B<-B[!(B$Well=="H2" & B$Plate=="Trial_6"),]
B<-B[!(B$Well=="A6" & B$Plate=="Trial_6"),]
B<-B[!(B$Well=="E7" & B$Plate=="Trial_6"),]
B<-B[!(B$Well=="A2" & B$Plate=="Trial_6"),]

B<-B[!(B$Well=="F7" & B$Plate=="Trial_7"),]
B<-B[!(B$Well=="A6" & B$Plate=="Trial_7"),]

B<-B[!(B$Well=="B2" & B$Plate=="Trial_8"),]
B<-B[!(B$Well=="C1" & B$Plate=="Trial_8"),]
B<-B[!(B$Well=="E6" & B$Plate=="Trial_8"),]
B<-B[!(B$Well=="F11" & B$Plate=="Trial_8"),]
B<-B[!(B$Well=="G11" & B$Plate=="Trial_8"),]
B<-B[!(B$Well=="F1" & B$Plate=="Trial_8"),]
B<-B[!(B$Well=="G1" & B$Plate=="Trial_8"),]

B<-B[!(B$Well=="H10" & B$Plate=="Trial_9"),]
B<-B[!(B$Well=="D12" & B$Plate=="Trial_9"),]
B<-B[!(B$Well=="E6" & B$Plate=="Trial_9"),]
B<-B[!(B$Well=="G10" & B$Plate=="Trial_9"),]
B<-B[!(B$Well=="F1" & B$Plate=="Trial_9"),]
B<-B[!(B$Well=="A2" & B$Plate=="Trial_9"),]
B<-B[!(B$Well=="B5" & B$Plate=="Trial_9"),]
B<-B[!(B$Well=="F5" & B$Plate=="Trial_9"),]
B<-B[!(B$Well=="C8" & B$Plate=="Trial_9"),]
B<-B[!(B$Well=="D8" & B$Plate=="Trial_9"),]
B<-B[!(B$Well=="G9" & B$Plate=="Trial_9"),]
B<-B[!(B$Well=="D11" & B$Plate=="Trial_9"),]
B<-B[!(B$Well=="F12" & B$Plate=="Trial_9"),]

B<-B[!(B$Well=="F1" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="G7" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="C10" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="A7" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="F8" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="A11" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="H1" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="F2" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="G2" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="A6" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="C7" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="E7" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="A8" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="D8" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="B9" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="C11" & B$Plate=="Trial_10"),]
B<-B[!(B$Well=="D11" & B$Plate=="Trial_10"),]

B<-B[!(B$Well=="A6" & B$Plate=="Trial_11"),]
B<-B[!(B$Well=="H9" & B$Plate=="Trial_11"),]
B<-B[!(B$Well=="A8" & B$Plate=="Trial_11"),]
B<-B[!(B$Well=="H6" & B$Plate=="Trial_11"),]
B<-B[!(B$Well=="D7" & B$Plate=="Trial_11"),]
B<-B[!(B$Well=="E12" & B$Plate=="Trial_11"),]
B<-B[!(B$Well=="G12" & B$Plate=="Trial_11"),]

B<-B[!(B$Well=="F10" & B$Plate=="Trial_12"),]
B<-B[!(B$Well=="F7" & B$Plate=="Trial_12"),]
B<-B[!(B$Well=="B1" & B$Plate=="Trial_12"),]
B<-B[!(B$Well=="C1" & B$Plate=="Trial_12"),]
B<-B[!(B$Well=="C2" & B$Plate=="Trial_12"),]
B<-B[!(B$Well=="E6" & B$Plate=="Trial_12"),]
B<-B[!(B$Well=="H10" & B$Plate=="Trial_12"),]
B<-B[!(B$Well=="C12" & B$Plate=="Trial_12"),]

B<-B[!(B$Well=="B6" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="A2" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="G2" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="A1" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="B2" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="A5" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="C5" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="E5" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="G5" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="C7" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="D7" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="B8" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="C8" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="A12" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="H7" & B$Plate=="Trial_13"),]
B<-B[!(B$Well=="G1" & B$Plate=="Trial_13"),]

B<-B[!(B$Well=="H9" & B$Plate=="Trial_14"),]
B<-B[!(B$Well=="B9" & B$Plate=="Trial_14"),]
B<-B[!(B$Well=="B5" & B$Plate=="Trial_14"),]
B<-B[!(B$Well=="C5" & B$Plate=="Trial_14"),]
B<-B[!(B$Well=="D7" & B$Plate=="Trial_14"),]
B<-B[!(B$Well=="D9" & B$Plate=="Trial_14"),]
B<-B[!(B$Well=="G9" & B$Plate=="Trial_14"),]
B<-B[!(B$Well=="B10" & B$Plate=="Trial_14"),]
B<-B[!(B$Well=="C11" & B$Plate=="Trial_14"),]
B<-B[!(B$Well=="E11" & B$Plate=="Trial_14"),]
B<-B[!(B$Well=="F12" & B$Plate=="Trial_14"),]
B<-B[!(B$Well=="A5" & B$Plate=="Trial_14"),]
B<-B[!(B$Well=="A12" & B$Plate=="Trial_14"),]


B<-B[!(B$Well=="B11" & B$Plate=="Trial_15"),]
B<-B[!(B$Well=="G1" & B$Plate=="Trial_15"),]
B<-B[!(B$Well=="D8" & B$Plate=="Trial_15"),]
B<-B[!(B$Well=="D5" & B$Plate=="Trial_15"),]
B<-B[!(B$Well=="H6" & B$Plate=="Trial_15"),]
B<-B[!(B$Well=="A10" & B$Plate=="Trial_15"),]
B<-B[!(B$Well=="D10" & B$Plate=="Trial_15"),]
B<-B[!(B$Well=="D11" & B$Plate=="Trial_15"),]
B<-B[!(B$Well=="G11" & B$Plate=="Trial_15"),]
B<-B[!(B$Well=="G12" & B$Plate=="Trial_15"),]

B<-B[!(B$Well=="E6" & B$Plate=="Trial_16"),]
B<-B[!(B$Well=="B10" & B$Plate=="Trial_16"),]
B<-B[!(B$Well=="F10" & B$Plate=="Trial_16"),]
B<-B[!(B$Well=="E1" & B$Plate=="Trial_16"),]
B<-B[!(B$Well=="E5" & B$Plate=="Trial_16"),]
B<-B[!(B$Well=="C7" & B$Plate=="Trial_16"),]
B<-B[!(B$Well=="E8" & B$Plate=="Trial_16"),]
B<-B[!(B$Well=="D10" & B$Plate=="Trial_16"),]
B<-B[!(B$Well=="E11" & B$Plate=="Trial_16"),]
B<-B[!(B$Well=="D12" & B$Plate=="Trial_16"),]

B<-B[!(B$Well=="G9" & B$Plate=="Trial_17"),]
B<-B[!(B$Well=="F10" & B$Plate=="Trial_17"),]
B<-B[!(B$Well=="H11" & B$Plate=="Trial_17"),]
B<-B[!(B$Well=="E12" & B$Plate=="Trial_17"),]
B<-B[!(B$Well=="D7" & B$Plate=="Trial_17"),]
B<-B[!(B$Well=="E11" & B$Plate=="Trial_17"),]
B<-B[!(B$Well=="F11" & B$Plate=="Trial_17"),]
B<-B[!(B$Well=="F6" & B$Plate=="Trial_17"),]
B<-B[!(B$Well=="F7" & B$Plate=="Trial_17"),]
B<-B[!(B$Well=="D8" & B$Plate=="Trial_17"),]
B<-B[!(B$Well=="E8" & B$Plate=="Trial_17"),]
B<-B[!(B$Well=="D9" & B$Plate=="Trial_17"),]

B<-B[!(B$Well=="F8" & B$Plate=="Trial_18"),]
B<-B[!(B$Well=="C1" & B$Plate=="Trial_18"),]
B<-B[!(B$Well=="E10" & B$Plate=="Trial_18"),]
B<-B[!(B$Well=="H9" & B$Plate=="Trial_18"),]
B<-B[!(B$Well=="B9" & B$Plate=="Trial_18"),]
B<-B[!(B$Well=="D12" & B$Plate=="Trial_18"),]
B<-B[!(B$Well=="H11" & B$Plate=="Trial_18"),]
B<-B[!(B$Well=="H12" & B$Plate=="Trial_18"),]

B<-B[!(B$Well=="E7" & B$Plate=="Trial_19"),]
B<-B[!(B$Well=="C11" & B$Plate=="Trial_19"),]
B<-B[!(B$Well=="A11" & B$Plate=="Trial_19"),]
B<-B[!(B$Well=="A12" & B$Plate=="Trial_19"),]
B<-B[!(B$Well=="D6" & B$Plate=="Trial_19"),]
B<-B[!(B$Well=="F12" & B$Plate=="Trial_19"),]
B<-B[!(B$Well=="H12" & B$Plate=="Trial_19"),]
B<-B[!(B$Well=="F7" & B$Plate=="Trial_19"),]
B<-B[!(B$Well=="G8" & B$Plate=="Trial_19"),]
B<-B[!(B$Well=="B7" & B$Plate=="Trial_19"),]
B<-B[!(B$Well=="D1" & B$Plate=="Trial_19"),]
B<-B[!(B$Well=="E12" & B$Plate=="Trial_19"),]
B<-B[!(B$Well=="G1" & B$Plate=="Trial_19"),]
B<-B[!(B$Well=="A1" & B$Plate=="Trial_19"),]

B<-B[!(B$Well=="C12" & B$Plate=="Trial_20"),]
B<-B[!(B$Well=="G8" & B$Plate=="Trial_20"),]
B<-B[!(B$Well=="G7" & B$Plate=="Trial_20"),]
B<-B[!(B$Well=="H9" & B$Plate=="Trial_20"),]
B<-B[!(B$Well=="G10" & B$Plate=="Trial_20"),]
B<-B[!(B$Well=="E10" & B$Plate=="Trial_20"),]
B<-B[!(B$Well=="D6" & B$Plate=="Trial_20"),]
B<-B[!(B$Well=="B9" & B$Plate=="Trial_20"),]
B<-B[!(B$Well=="F9" & B$Plate=="Trial_20"),]
B<-B[!(B$Well=="A8" & B$Plate=="Trial_20"),]
B<-B[!(B$Well=="B2" & B$Plate=="Trial_20"),]
B<-B[!(B$Well=="A12" & B$Plate=="Trial_20"),]

B<-B[!(B$Well=="B11" & B$Plate=="Trial_21"),]
B2<-B[!(B$Well=="H9" & B$Plate=="Trial_21"),]

summary(B2)

# subset off and on
AT1<-subset(B2, Condition== "AT1")
AT5<-subset(B2, Condition== "AT5")
AT10<-subset(B2, Condition== "AT10")
BT<-subset(B2, Condition== "BT")

B2$Condition <- fct_relevel(B2$Condition, c("BT","AT1","AT5","AT10"))
B2$Expo <- fct_relevel(B2$Expo, c("Control","PFOS_LC","PFOS_HC","PFBS_LC", "PFBS_HC"))

levels(B2$Condition)

Sup_Fig_5_B <- ggplot(B2, aes(x=Expo, y=Dmcptot, fill=Condition)) + scale_y_continuous(limits=c(0,1.5))+ # oder colour for outer lines
  geom_boxplot() + scale_fill_manual(name="Periods",breaks = c("BT", "AT1", "AT5","AT10"), 
                                     values=c("white", "magenta", "red", "darkred")) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)+
  # theme(legend.position = "none")+
  theme_light() + 
  xlab("Condition") + ylab("Distance moved (cm)") + 
  theme(axis.text.x = element_text(size=18,angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))
Sup_Fig_5_B
#-------------------------
####Fig. 2. up####
#Light test
#subset Concentrations for LOFF
LOFFC<-subset(LOFF,Expo=="Control")
LOFFPFOSLC<-subset(LOFF,Expo=="PFOS_LC")
LOFFPFOSHC<-subset(LOFF,Expo=="PFOS_HC")
LOFFPFBSLC<-subset(LOFF,Expo=="PFBS_LC")
LOFFPFBSHC<-subset(LOFF,Expo=="PFBS_HC")

#check distribution with histogram
hist(LOFF$Dmcptot,prob=T,main="Distance moved",ylim=c(0,0.025))
points(density(LOFF$Dmcptot),type="l",col="red", lwd=2)

shapiro.test(LOFFC$Dmcptot)
shapiro.test(LOFFPFOSLC$Dmcptot)
shapiro.test(LOFFPFOSHC$Dmcptot)
shapiro.test(LOFFPFBSLC$Dmcptot)
shapiro.test(LOFFPFBSHC$Dmcptot)

compare_means(Dmcptot ~ Expo,  data = LOFF, method = "wilcox.test")

# Create boxplot to graphically show your results
LOFF$Expo <- fct_relevel(LOFF$Expo, c("Control", "PFOS_LC", "PFOS_HC", "PFBS_LC", "PFBS_HC"))

Fig_2_A_up <- ggplot(LOFF, aes(x=Expo, y=Dmcptot, fill=Expo)) + # oder colour for outer lines
  geom_boxplot() +
  scale_fill_manual(values = c("cornflowerblue","#FFFF00","#FFCC33","green","darkgreen")) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)+
  # theme(legend.position = "none")+
  theme_light()+
  xlab("Condition")+
  ylab("Total distance moved (cm)")+
  theme(axis.text.x = element_text(size=18,angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  geom_pwc(method = "wilcox.test", label = "{p.format} {p.signif}", hide.ns = TRUE, tip.length = 0, p.adjust.method = "none", ref.group = "Control", label.size = 6)
Fig_2_A_up + guides(fill=guide_legend(title="Condition"))

#Tapping test
###Distribution
# histogram to graphically view how the data is distributed

hist(AT1$Dmcptot,prob=T,main="Distance moved")
points(density(AT1$Dmcptot),type="l",col="red", lwd=2)

hist(BT$Dmcptot,prob=T,main="Distance moved", ylim = c(0,8))
points(density(BT$Dmcptot),type="l",col="red", lwd=2)

hist(AT5$Dmcptot,prob=T,main="Distance moved")
points(density(AT5$Dmcptot),type="l",col="red", lwd=2)

hist(AT10$Dmcptot,prob=T,main="Distance moved")
points(density(AT10$Dmcptot),type="l",col="red", lwd=2)

#AT5
#subset data

AT5_Control<-AT5[AT5$Expo=="Control",]
AT5_PFOS_LC<-AT5[AT5$Expo=="PFOS_LC",]
AT5_PFOS_HC<-AT5[AT5$Expo=="PFOS_HC",]
AT5_PFBS_LC<-AT5[AT5$Expo=="PFBS_LC",]
AT5_PFBS_HC<-AT5[AT5$Expo=="PFBS_HC",]

#distribution

plot(density(AT5_Control$Dmcptot))
plot(density(AT5_PFOS_LC$Dmcptot))
plot(density(AT5_PFOS_HC$Dmcptot))
plot(density(AT5_PFBS_LC$Dmcptot))
plot(density(AT5_PFBS_HC$Dmcptot))

#check distribution with Shapiro Test
shapiro.test(AT5_Control$Dmcptot)
shapiro.test(AT5_PFOS_LC$Dmcptot)
shapiro.test(AT5_PFOS_HC$Dmcptot)
shapiro.test(AT5_PFBS_LC$Dmcptot)
shapiro.test(AT5_PFBS_HC$Dmcptot)

#all not normally distributed! --> Wilcox Test

#pairwise comparison
wilcox.test(AT5_Control$Dmcptot, AT5_PFOS_LC$Dmcptot)
wilcox.test(AT5_Control$Dmcptot, AT5_PFOS_HC$Dmcptot) 
wilcox.test(AT5_Control$Dmcptot, AT5_PFBS_LC$Dmcptot) 
wilcox.test(AT5_Control$Dmcptot, AT5_PFBS_HC$Dmcptot) 

compare_means(Dmcptot ~ Expo,  data = AT5, method = "wilcox.test")

AT5$Expo <- fct_relevel(AT5$Expo, c("Control", "PFOS_LC", "PFOS_HC", "PFBS_LC", "PFBS_HC"))

Fig_2_B_up<- ggplot(AT5, aes(x=Expo, y=Dmcptot, fill=Expo)) + # oder colour for outer lines
  geom_boxplot() +
  scale_fill_manual(values = c("cornflowerblue","#FFFF00","#FFCC33","green","darkgreen")) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2)+
  # theme(legend.position = "none")+
  theme_light()+
  xlab("Condition")+
  ylab("Total distance moved (cm)")+
  theme(axis.text.x = element_text(size=18,angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  geom_pwc(method = "wilcox.test", label = "{p.format} {p.signif}", hide.ns = TRUE, tip.length = 0, p.adjust.method = "none", ref.group = "Control", label.size = 6)
Fig_2_B_up+ guides(fill=guide_legend(title="Condition"))
#-------------------------
####Sup.Fig.5 B & C####
#Light test
#subset Concentrations for LON1
LON1C<-subset(LON1,Expo=="Control")
LON1PFOSLC<-subset(LON1,Expo=="PFOS_LC")
LON1PFOSHC<-subset(LON1,Expo=="PFOS_HC")
LON1PFBSLC<-subset(LON1,Expo=="PFBS_LC")
LON1PFBSHC<-subset(LON1,Expo=="PFBS_HC")

#check distribution with histogram
hist(LON1$Dmcptot,prob=T,main="Distance moved",ylim=c(0,0.05))
points(density(LON1$Dmcptot),type="l",col="red", lwd=2)

shapiro.test(LON1C$Dmcptot)
shapiro.test(LON1PFOSLC$Dmcptot)
shapiro.test(LON1PFOSHC$Dmcptot)
shapiro.test(LON1PFBSLC$Dmcptot)
shapiro.test(LON1PFBSHC$Dmcptot)

compare_means(Dmcptot ~ Expo,  data = LON1, method = "wilcox.test")

Sup_Fig_5_B <- ggplot(LON1, aes(x=Expo, y=Dmcptot, fill=Expo)) + # oder colour for outer lines
  geom_boxplot() +
  scale_fill_manual(values = c("cornflowerblue","#FFFF00","#FFCC33","green","darkgreen")) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)+
  # theme(legend.position = "none")+
  theme_light()+
  xlab("Condition")+
  ylab("Distance moved (cm)")+
  theme(axis.text.x = element_text(size=18,angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  labs(title = "Acclimation phase (LON1)") + 
  geom_pwc(method = "wilcox.test", hide.ns = FALSE, tip.length = 0, p.adjust.method = "none", ref.group = "Control", label.size = 6)
Sup_Fig_5_B + guides(fill=guide_legend(title="Condition"))

#subset Concentrations for LON2
LON2C<-subset(LON2,Expo=="Control")
LON2PFOSLC<-subset(LON2,Expo=="PFOS_LC")
LON2PFOSHC<-subset(LON2,Expo=="PFOS_HC")
LON2PFBSLC<-subset(LON2,Expo=="PFBS_LC")
LON2PFBSHC<-subset(LON2,Expo=="PFBS_HC")

#check distribution with histogram
hist(LON2$Dmcptot,prob=T,main="Distance moved",ylim=c(0,0.055))
points(density(LON2$Dmcptot),type="l",col="red", lwd=2)

shapiro.test(LON2C$Dmcptot)
shapiro.test(LON2PFOSLC$Dmcptot)
shapiro.test(LON2PFOSHC$Dmcptot)
shapiro.test(LON2PFBSLC$Dmcptot)
shapiro.test(LON2PFBSHC$Dmcptot)

compare_means(Dmcptot ~ Expo,  data = LON2, method = "wilcox.test")

Sup_Fig_5_C <- ggplot(LON2, aes(x=Expo, y=Dmcptot, fill=Expo)) + # oder colour for outer lines
  geom_boxplot() +
  scale_fill_manual(values = c("cornflowerblue","#FFFF00","#FFCC33","green","darkgreen")) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)+
  # theme(legend.position = "none")+
  theme_light()+
  xlab("Condition")+
  ylab("Distance moved (cm)")+
  labs(title = "Resting phase (LON2)") + 
  theme(axis.text.x = element_text(size=18,angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  geom_pwc(method = "wilcox.test", hide.ns = FALSE, tip.length = 0, p.adjust.method = "none", ref.group = "Control", label.size = 6)
Sup_Fig_5_C + guides(fill=guide_legend(title="Condition"))

#Tapping test
#BT
#subset data

BT_Control<-BT[BT$Expo=="Control",]
BT_PFOS_LC<-BT[BT$Expo=="PFOS_LC",]
BT_PFOS_HC<-BT[BT$Expo=="PFOS_HC",]
BT_PFBS_LC<-BT[BT$Expo=="PFBS_LC",]
BT_PFBS_HC<-BT[BT$Expo=="PFBS_HC",]

#check distribution with Shapiro Test
shapiro.test(BT_Control$Dmcptot)
shapiro.test(BT_PFOS_LC$Dmcptot)
shapiro.test(BT_PFOS_HC$Dmcptot)
shapiro.test(BT_PFBS_LC$Dmcptot)
shapiro.test(BT_PFBS_HC$Dmcptot)

#all not normally distributed! --> Wilcox Test

#pairwise comparison
wilcox.test(BT_Control$Dmcptot, BT_PFOS_LC$Dmcptot)
wilcox.test(BT_Control$Dmcptot, BT_PFOS_HC$Dmcptot) 
wilcox.test(BT_Control$Dmcptot, BT_PFBS_LC$Dmcptot) 
wilcox.test(BT_Control$Dmcptot, BT_PFBS_HC$Dmcptot) 

#no significance

BT$Expo <- fct_relevel(BT$Expo, c("Control", "PFOS_LC", "PFOS_HC", "PFBS_LC", "PFBS_HC"))

Sup_Fig_5_D <- ggplot(BT, aes(x=Expo, y=Dmcptot, fill=Expo)) + # oder colour for outer lines
  geom_boxplot() +
  scale_fill_manual(values = c("cornflowerblue","#FFFF00","#FFCC33","green","darkgreen")) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2)+
  # theme(legend.position = "none")+
  theme_light()+
  xlab("Condition")+
  ylab("Distance moved (cm)")+
  labs(title = "Before tapping (BT)") + 
  theme(axis.text.x = element_text(size=18,angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  geom_pwc(method = "wilcox.test", hide.ns = FALSE, tip.length = 0, p.adjust.method = "none", ref.group = "Control", label.size = 6)
Sup_Fig_5_D + guides(fill=guide_legend(title="Condition"))

#AT1
#subset data
AT1_Control<-AT1[AT1$Expo=="Control",]
AT1_PFOS_LC<-AT1[AT1$Expo=="PFOS_LC",]
AT1_PFOS_HC<-AT1[AT1$Expo=="PFOS_HC",]
AT1_PFBS_LC<-AT1[AT1$Expo=="PFBS_LC",]
AT1_PFBS_HC<-AT1[AT1$Expo=="PFBS_HC",]

plot(density(AT1_Control$Dmcptot))
plot(density(AT1_PFOS_LC$Dmcptot))
plot(density(AT1_PFOS_HC$Dmcptot))
plot(density(AT1_PFBS_LC$Dmcptot))
plot(density(AT1_PFBS_HC$Dmcptot))

shapiro.test(AT1_Control$Dmcptot)
shapiro.test(AT1_PFOS_LC$Dmcptot)
shapiro.test(AT1_PFOS_HC$Dmcptot)
shapiro.test(AT1_PFBS_LC$Dmcptot)
shapiro.test(AT1_PFBS_HC$Dmcptot)

#pairwise comparison
wilcox.test(AT1_Control$Dmcptot, AT1_PFOS_LC$Dmcptot)
wilcox.test(AT1_Control$Dmcptot, AT1_PFOS_HC$Dmcptot) 
wilcox.test(AT1_Control$Dmcptot, AT1_PFBS_LC$Dmcptot) 
wilcox.test(AT1_Control$Dmcptot, AT1_PFBS_HC$Dmcptot) 

compare_means(Dmcptot ~ Expo,  data = AT1, method = "wilcox.test")

#no significance

AT1$Expo <- fct_relevel(AT1$Expo, c("Control", "PFOS_LC", "PFOS_HC", "PFBS_LC", "PFBS_HC"))

Sup_Fig_5_E <- ggplot(AT1, aes(x=Expo, y=Dmcptot, fill=Expo)) + # oder colour for outer lines
  geom_boxplot() +
  scale_fill_manual(values = c("cornflowerblue","#FFFF00","#FFCC33","green","darkgreen")) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)+
  # theme(legend.position = "none")+
  theme_light()+
  xlab("Condition")+
  ylab("Distance moved (cm)")+
  labs(title = "After tapping 1 (AT1)") + 
  theme(axis.text.x = element_text(size=18,angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  geom_pwc(method = "wilcox.test", hide.ns = FALSE, tip.length = 0, p.adjust.method = "none", ref.group = "Control", label.size = 6)
Sup_Fig_5_E + guides(fill=guide_legend(title="Condition"))

#AT10
#subset data
AT10_Control<-AT10[AT10$Expo=="Control",]
AT10_PFOS_LC<-AT10[AT10$Expo=="PFOS_LC",]
AT10_PFOS_HC<-AT10[AT10$Expo=="PFOS_HC",]
AT10_PFBS_LC<-AT10[AT10$Expo=="PFBS_LC",]
AT10_PFBS_HC<-AT10[AT10$Expo=="PFBS_HC",]

shapiro.test(AT10_Control$Dmcptot)
shapiro.test(AT10_PFOS_LC$Dmcptot)
shapiro.test(AT10_PFOS_HC$Dmcptot)
shapiro.test(AT10_PFBS_LC$Dmcptot)
shapiro.test(AT10_PFBS_HC$Dmcptot)
#all not normally distributed! --> Wilcox Test

#pairwise comparison
wilcox.test(AT10_Control$Dmcptot, AT10_PFOS_LC$Dmcptot)
wilcox.test(AT10_Control$Dmcptot, AT10_PFOS_HC$Dmcptot) 
wilcox.test(AT10_Control$Dmcptot, AT10_PFBS_LC$Dmcptot) 
wilcox.test(AT10_Control$Dmcptot, AT10_PFBS_HC$Dmcptot) 

compare_means(Dmcptot ~ Expo,  data = AT10, method = "wilcox.test")

AT10$Expo <- fct_relevel(AT10$Expo, c("Control", "PFOS_LC", "PFOS_HC", "PFBS_LC", "PFBS_HC"))

Sup_Fig_5_F <- ggplot(AT10, aes(x=Expo, y=Dmcptot, fill=Expo)) + # oder colour for outer lines
  geom_boxplot() +
  scale_fill_manual(values = c("cornflowerblue","#FFFF00","#FFCC33","green","darkgreen")) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2)+
  # theme(legend.position = "none")+
  theme_light()+
  xlab("Condition")+
  ylab("Distance moved (cm)")+
  labs(title = "After tapping 1 (AT10)") + 
  theme(axis.text.x = element_text(size=18,angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  geom_pwc(method = "wilcox.test", hide.ns = FALSE, tip.length = 0, p.adjust.method = "none", ref.group = "Control", label.size = 6)
Sup_Fig_5_F + guides(fill=guide_legend(title="Condition"))
#-------------------------
####Fig. 2. C####
setwd("~/Papers/EpiTrans F0/R analysis") #the path to where you have your files 
behav<-read.table("Statistics-Jonas_F0_AdultBehavior(Improved).txt",h=T)
summary(behav) #to get basic information about your data
head(behav) #to see the head of your dataset

#put conditions in the right order
behav$Condition <- fct_relevel(behav$Condition, c("Control", "PFOS_LC", "PFOS_HC", "PFBS_LC", "PFBS_HC"))
behav$Condition_sex <- fct_relevel(behav$Condition_sex, c("Control_female", "Control_male", "PFOS_LC_female","PFOS_LC_male", "PFOS_HC_female", "PFOS_HC_male", "PFBS_LC_female", "PFBS_LC_male", "PFBS_HC_female", "PFBS_HC_male"))

#Subset to see difference between each condition according to sex
Control<-subset(behav,Condition=="Control")
PFOS_LC<-subset(behav,Condition=="PFOS_LC")
PFOS_HC<-subset(behav,Condition=="PFOS_HC")
PFBS_LC<-subset(behav,Condition=="PFBS_LC")
PFBS_HC<-subset(behav,Condition=="PFBS_HC")

Fig_2_C <- ggplot(behav, aes(x=Condition, y=Latency_to_first_entry_in_top_zone, fill=Condition)) + # oder colour for outer lines
  geom_boxplot() +
  scale_fill_manual(values = c("cornflowerblue","#FFFF00","#FFCC33","green","darkgreen"), labels = c("Control", "PFOS LC", "PFOS HC", "PFBS LC", "PFBS HC")) + 
  labs(title = "F0 - NTT: Latency to top zone", y = "Latency to top zone (s)", x = "Condition") + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  # theme(legend.position = "none")+
  theme_light()+
  theme(axis.text.x = element_text(size=18,angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  geom_pwc(method = "wilcox.test", label = "{p.format} {p.signif}", hide.ns = TRUE, tip.length = 0, p.adjust.method = "none", ref.group = "Control", label.size = 6)
Fig_2_C

# Wilcox test as alternative for non-normally distributed data

wilcox.test(Control$Latency_to_first_entry_in_top_zone, PFOS_LC$Latency_to_first_entry_in_top_zone)
wilcox.test(Control$Latency_to_first_entry_in_top_zone, PFOS_HC$Latency_to_first_entry_in_top_zone) 
wilcox.test(Control$Latency_to_first_entry_in_top_zone, PFBS_LC$Latency_to_first_entry_in_top_zone) 
wilcox.test(Control$Latency_to_first_entry_in_top_zone, PFBS_HC$Latency_to_first_entry_in_top_zone) 

#data:  Control and PFOS_HC
#W = 295.5, p-value = 0.02278
#-------------------------
####Fig. 2. D####
###Per conditions in females
females<-subset(behav,Sex=="Female") #To select only the females
summary(females)
females$Condition <- fct_relevel(females$Condition, c("Control", "PFOS_LC", "PFOS_HC", "PFBS_LC", "PFBS_HC"))

###Per conditions in males
males<-subset(behav,Sex=="Male") #To select only the males
summary(males)
males$Condition <- fct_relevel(males$Condition, c("Control", "PFOS_LC", "PFOS_HC", "PFBS_LC", "PFBS_HC"))

####Entries to top zone
#Boxplot entries to top zone per condition and to sex
Fig_2_D <- ggplot(behav, aes(x=Condition, y=Entries_to_top_zone, fill=Condition)) + # oder colour for outer lines
  geom_boxplot() +
  facet_grid(cols = vars(Sex))+
  scale_fill_manual(values = c("cornflowerblue","#FFFF00","#FFCC33","green","darkgreen"), labels = c("Control", "PFOS LC", "PFOS HC", "PFBS LC", "PFBS HC")) + 
  labs(title = "F0 - NTT: Entries to top zone", y = "Entries to top zone", x = "Condition") + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  # theme(legend.position = "none")+
  theme_light()+
  theme(axis.text.x = element_text(size=18,angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  geom_pwc(method = "wilcox.test", label = "{p.format} {p.signif}", hide.ns = TRUE, tip.length = 0, p.adjust.method = "none", y.position = 40, ref.group = "Control", label.size = 6)+
  theme(strip.text.x = element_text(size = 16))
Fig_2_D

#ggsave("Entries to top zone.png", width = 16, height = 9, p)
#-------------------------
####Fig. 2. E####
####Total duration in bottom zone
#Boxplot Total duration in bottom zone within condition 
Fig_2_E <- ggplot(behav, aes(x = Sex, y = Cumulative_duration_in_bottom_zone, fill=Sex)) +
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  facet_grid(cols = vars(Condition), labeller = labeller(cols = vars(Condition))) + 
  scale_fill_manual(values = c("grey", "red"), labels = c("Female", "Male")) + 
  labs(title = "F0 - NTT: Total duration in bottom zone", y = "Total duration in bottom zone (s)", x = "Condition") + 
  geom_pwc(method = "wilcox.test", label = "{p.format} {p.signif}", hide.ns = TRUE, tip.length = 0, p.adjust.method = "none", label.size = 6) +
  theme_light()+
  theme(axis.text.x = element_text(size=18,angle = 60, vjust = 1, hjust=1), 
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(color = "black", size = 20, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 20, angle = 90, face = "bold"))+
  theme(strip.text.x = element_text(size = 16))
Fig_2_E
#ggsave("test8.png", width = 16, height = 9, p)
#-------------------------
#####Epigenetic####
####Big heatmap####
epi_PFOS<-read.table("Epi_heatmap_PFOS_&_PFBS.txt",header=T, sep="\t")
colnames(epi_PFOS)<-c("Epigenetic function", "gene","2hpf PFOS LC", "2hpf PFOS HC", "5dpf PFOS LC", "5dpf PFOS HC", "20dpf PFOS LC", "20dpf PFOS HC",
                      "2hpf PFBS LC", "2hpf PFBS HC", "5dpf PFBS LC", "5dpf PFBS HC", "20dpf PFBS LC", "20dpf PFBS HC")
rownames(epi_PFOS) <- epi_PFOS$gene
geneExp_matrix <- as.matrix(epi_PFOS[3:14])

my_hclust_gene <- hclust(dist(geneExp_matrix), method = "complete")
as.dendrogram(my_hclust_gene) %>%
  plot(horiz = TRUE)
my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 1)
my_gene_col <- data.frame(my_gene_col)
my_gene_col$`Epigenetic function` <- epi_PFOS$`Epigenetic function`

annotation <- data.frame(Var1 = factor(1:10 %% 2 == 0, labels = c("Exp1", "Exp2")))
rownames(annotation) <- colnames(test) # check out the row names of annotation
pheatmap(test, annotation = annotation)

# Plot the heatmap
bk1 <- c(seq(-4,-0.1,by=0.1),-0.001)
bk2 <- c(0.001,seq(0.1,3,by=0.1))
bk <- c(bk1,bk2)  #combine the break limits for purpose of graphing

my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(bk1)-1.5),
                "white", "white",
                c(colorRampPalette(colors = c("tomato1","darkred"))(n = length(bk2)-1)))

pheatmap(geneExp_matrix, color = my_palette, breaks = bk, 
         cluster_rows = T, cluster_cols = T, margin = c(5,5),annotation_row = my_gene_col,fontsize_row = 2)

Epi_heatmap<-pheatmap(geneExp_matrix,color=my_palette,annotation_row = my_gene_col,fontsize_row = 4)
####Work on 5dpf####
Data_5dpf<-as.data.frame(epi_PFOS[,1])
row.names(Data_5dpf)<-rownames(epi_PFOS)
Data_5dpf<-cbind(Data_5dpf, epi_PFOS[,5])
Data_5dpf<-cbind(Data_5dpf, epi_PFOS[,6])
Data_5dpf<-cbind(Data_5dpf, epi_PFOS[,11])
Data_5dpf<-cbind(Data_5dpf, epi_PFOS[,12])
colnames(Data_5dpf)<-c("Epigenetic function","5dpf PFOS LC","5dpf PFOS HC", "5dpf PFBS LC","5dpf PFBS HC")

epi_chrom_remod<-subset(Data_5dpf,Data_5dpf$`Epigenetic function`=="chromatin remodeler")
geneExp_matrix <- as.matrix(epi_chrom_remod[,2:5])
pheatmap(geneExp_matrix, color = my_palette, breaks = bk, 
         cluster_rows = T, cluster_cols = T, margin = c(5,5),fontsize_row = 9)

epi_RNA_modification<-subset(epi_PFOS,epi_PFOS$`Epigenetic function`=="RNA modification")
geneExp_matrix <- as.matrix(epi_RNA_modification[3:14])
Epi_heatmap<-pheatmap(geneExp_matrix,color=colorRampPalette(c("navy", "white", "red"))(50),fontsize_row = 30)

epi_histone_eraser<-subset(epi_PFOS,epi_PFOS$`Epigenetic function`=="histone eraser")
geneExp_matrix <- as.matrix(epi_histone_eraser[3:14])
Epi_heatmap<-pheatmap(geneExp_matrix,color=colorRampPalette(c("navy", "white", "red"))(50),fontsize_row = 15)

#####Venn diagram####
twohpf_PFOS_LC_genes<-subset(epi_PFOS,epi_PFOS$`2hpf PFOS LC`!=0.000000)
twohpf_PFOS_LC_genes<-row.names(twohpf_PFOS_LC_genes)
twohpf_PFOS_HC_genes<-subset(epi_PFOS,epi_PFOS$`2hpf PFOS HC`!=0.000000)
twohpf_PFOS_HC_genes<-row.names(twohpf_PFOS_HC_genes)
fivedpf_PFOS_LC_genes<-subset(epi_PFOS,epi_PFOS$`5dpf PFOS LC`!=0.000000)
fivedpf_PFOS_LC_genes<-row.names(fivedpf_PFOS_LC_genes)
fivedpf_PFOS_HC_genes<-subset(epi_PFOS,epi_PFOS$`5dpf PFOS HC`!=0.000000)
fivedpf_PFOS_HC_genes<-row.names(fivedpf_PFOS_HC_genes)
twentydpf_PFOS_LC_genes<-subset(epi_PFOS,epi_PFOS$`20dpf PFOS LC`!=0.000000)
twentydpf_PFOS_LC_genes<-row.names(twentydpf_PFOS_LC_genes)
twentydpf_PFOS_HC_genes<-subset(epi_PFOS,epi_PFOS$`20dpf PFOS HC`!=0.000000)
twentydpf_PFOS_HC_genes<-row.names(twentydpf_PFOS_HC_genes)
twohpf_PFBS_LC_genes<-subset(epi_PFOS,epi_PFOS$`2hpf PFBS LC`!=0.000000)
twohpf_PFBS_LC_genes<-row.names(twohpf_PFBS_LC_genes)
twohpf_PFBS_HC_genes<-subset(epi_PFOS,epi_PFOS$`2hpf PFBS HC`!=0.000000)
twohpf_PFBS_HC_genes<-row.names(twohpf_PFBS_HC_genes)
fivedpf_PFBS_LC_genes<-subset(epi_PFOS,epi_PFOS$`5dpf PFBS LC`!=0.000000)
fivedpf_PFBS_LC_genes<-row.names(fivedpf_PFBS_LC_genes)
fivedpf_PFBS_HC_genes<-subset(epi_PFOS,epi_PFOS$`5dpf PFBS HC`!=0.000000)
fivedpf_PFBS_HC_genes<-row.names(fivedpf_PFBS_HC_genes)
twentydpf_PFBS_LC_genes<-subset(epi_PFOS,epi_PFOS$`20dpf PFBS LC`!=0.000000)
twentydpf_PFBS_LC_genes<-row.names(twentydpf_PFBS_LC_genes)
twentydpf_PFBS_HC_genes<-subset(epi_PFOS,epi_PFOS$`20dpf PFBS HC`!=0.000000)
twentydpf_PFBS_HC_genes<-row.names(twentydpf_PFBS_HC_genes)

PFOS_genes<-unique(c(twohpf_PFOS_LC_genes,twohpf_PFOS_HC_genes,fivedpf_PFOS_LC_genes,fivedpf_PFOS_HC_genes,twentydpf_PFOS_LC_genes,twentydpf_PFOS_HC_genes))
PFBS_genes<-unique(c(twohpf_PFBS_LC_genes,twohpf_PFBS_HC_genes,fivedpf_PFBS_LC_genes,fivedpf_PFBS_HC_genes,twentydpf_PFBS_LC_genes,twentydpf_PFBS_HC_genes))

PFOS_2hpf<-unique(c(twohpf_PFOS_LC_genes,twohpf_PFOS_HC_genes))
PFOS_5dpf<-unique(c(fivedpf_PFOS_LC_genes,fivedpf_PFOS_HC_genes))
PFOS_20dpf<-unique(c(twentydpf_PFOS_LC_genes,twentydpf_PFOS_HC_genes))

PFBS_2hpf<-unique(c(twohpf_PFBS_LC_genes,twohpf_PFBS_HC_genes))
PFBS_5dpf<-unique(c(fivedpf_PFBS_LC_genes,fivedpf_PFBS_HC_genes))
PFBS_20dpf<-unique(c(twentydpf_PFBS_LC_genes,twentydpf_PFBS_HC_genes))

#PFOS vs PFBS
venn.diagram(
  x = list(PFOS_genes, PFBS_genes),
  category.names = c("" , ""),
  filename = 'Epi_genes_PFOS_&_PFBS.png',
  output=TRUE,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  lwd = 2,
  lty = 'blank',
  fill = c("orange", "green"),
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer"
)

#PFOS all
x = list(twohpf_PFOS_LC_genes,fivedpf_PFOS_LC_genes,fivedpf_PFOS_HC_genes,twentydpf_PFOS_LC_genes,twentydpf_PFOS_HC_genes)
ggVennDiagram(x,category.names = c("2hpf PFOS LC" , "5dpf PFOS LC", "5dpf PFOS HC", "20dpf PFOS LC", "20dpf PFOS HC"),label_size = 6) + scale_fill_gradient(low="blue",high = "red")

#PFBS all
y = list(twohpf_PFBS_HC_genes,fivedpf_PFBS_LC_genes,fivedpf_PFBS_HC_genes,twentydpf_PFBS_LC_genes,twentydpf_PFBS_HC_genes)
ggVennDiagram(y,category.names = c("2hpf PFBS HC", "5dpf PFBS LC", "5dpf PFBS HC", "20dpf PFBS LC", "20dpf PFBS HC"),label_size = 6) + scale_fill_gradient(low="blue",high = "red")

#PFOS time
venn.diagram(
  x = list(PFOS_2hpf, PFOS_5dpf,PFOS_20dpf),
  category.names = c("" , "", ""),
  filename = 'Epi_genes_PFOS.png',
  output=TRUE,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  lwd = 2,
  lty = 'blank',
  fill = c("yellow", "orange","red"),
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer"
)

#PFBS time
venn.diagram(
  x = list(PFBS_2hpf, PFBS_5dpf,PFBS_20dpf),
  category.names = c("" , "", ""),
  filename = 'Epi_genes_PFBS.png',
  output=TRUE,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  lwd = 2,
  lty = 'blank',
  fill = c("lightgreen", "green","darkgreen"),
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer"
)
