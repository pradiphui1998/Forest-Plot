library(foreign)
data<- read.spss("data.sav",to.data.frame=TRUE)
data1=as.data.frame(data)
View(data1)
attach(data1)

#### For hazard ratio
library(survival)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(condSURV)

## For PFS
## Molecular Subgroup
data1$molecular_sugroup <- factor(data1$molecular_sugroup, levels = c("WNT","SHH","group3","group4"))

coxph(Surv(truncated_pfs, truncated_pfs_status) ~molecular_sugroup, data =data1) %>% 
  tbl_regression(exp = TRUE) 


### Age group
data1$age_group<- factor(data1$age_group, levels = c("Less than 3","Greater than 3"))

coxph(Surv(truncated_pfs, truncated_pfs_status) ~age_group, data =data1) %>% 
  tbl_regression(exp = TRUE) 


### Gender
data1$Gender<- factor(data1$Gender, levels = c("Male","Female"))

coxph(Surv(truncated_pfs, truncated_pfs_status) ~Gender, data =data1) %>% 
  tbl_regression(exp = TRUE) 


### Rec-EOR
data1$rec_EOR<- factor(data1$rec_EOR, levels = c("STR","GTR+NTR"))

coxph(Surv(truncated_pfs, truncated_pfs_status) ~rec_EOR, data =data1) %>% 
  tbl_regression(exp = TRUE) 


### Met-prep
data1$mets_presentation<- factor(data1$mets_presentation, levels = c("NO","YES"))

coxph(Surv(truncated_pfs, truncated_pfs_status) ~mets_presentation, data =data1) %>% 
  tbl_regression(exp = TRUE) 



#---------------------------------------------------------------------------------------
library(forestplot)
library(dplyr)

data <- structure(list(mean  = c(NA,2.31,2.81,2.89,0.48,0.87,0.90,2.55,NA),
                       lower = c(NA,1.16,1.33,1.45,0.26,0.56,0.59,1.72,NA),
                       upper = c(NA,4.58,5.94,5.78,0.89,1.35,1.38,3.78,NA)),
                  .Names = c("mean", "lower", "upper"),
                  row.names = c(NA, -9L),
                  class = "data.frame")
data


tabletext <- cbind(c("", "Molecular Subgroup","","", "Age", "Gender","Extended Resections","Metastasis at Presentation",NA),
                   c("","SHH","group3","group4","Greater than 3","Female","GTR+NTR","YES",NA),
                  
                   c("Hazard Ratio",2.31,2.81,2.89,0.48,0.87,0.90,2.55,NA),
                   c("95% CI","(1.16,4.58)","(1.33,5.94)","(1.45,5.78)","(0.26,0.89)","(0.56,1.35)","(0.59,1.38)","(1.72,3.78)",NA)
                   #c("P-Value",0.017,0.007,0.003,0.020,0.500,0.600,"<0.001","")
                   )
tabletext

data %>%
  forestplot(labeltext = tabletext,
             #is.summary = c(rep(TRUE, 1), rep(FALSE, 7)),## Make test bold to highlight summary rows
             #clip = c(0.1,5),## Lower and upper limit of clipping confidence intervals to arrows
             lty.ci=1,
             hrzl_lines=list("2"=gpar(lty=1,lwd=2),
                             "9"=gpar(lwd=2,colums=1:2,col="#000044")),
             xlog = TRUE,
            xlab="Hazard Ratio",
            graph.pos=4,
            align="l",
            graphwidth = "auto",
             title="Forest Plot for Progression Free Survival",
             boxsize=0.1,
             col = fpColors(box = "red",
                            line = "darkblue"))



#----------------------------------------------------------------------------------------
  
  
  
  
  
  
  
  
  
  
  




#---------------------------------------------------------------------------------------
  ## For OS
  ## Molecular Subgroup
  data1$molecular_sugroup <- factor(data1$molecular_sugroup, levels = c("WNT","SHH","group3","group4"))

coxph(Surv(truncated_OS,truncates_osstatus) ~molecular_sugroup, data =data1) %>% 
  tbl_regression(exp = TRUE) 


### Age group
data1$age_group<- factor(data1$age_group, levels = c("Less than 3","Greater than 3"))

coxph(Surv(truncated_OS,truncates_osstatus) ~age_group, data =data1) %>% 
  tbl_regression(exp = TRUE) 


### Gender
data1$Gender<- factor(data1$Gender, levels = c("Male","Female"))

coxph(Surv(truncated_OS,truncates_osstatus) ~Gender, data =data1) %>% 
  tbl_regression(exp = TRUE) 


### Rec-EOR
data1$rec_EOR<- factor(data1$rec_EOR, levels = c("STR","GTR+NTR"))

coxph(Surv(truncated_OS,truncates_osstatus) ~rec_EOR, data =data1) %>% 
  tbl_regression(exp = TRUE) 


### Met-prep
data1$mets_presentation<- factor(data1$mets_presentation, levels = c("NO","YES"))

coxph(Surv(truncated_OS,truncates_osstatus) ~mets_presentation, data =data1) %>% 
  tbl_regression(exp = TRUE) 






#---------------------------------------------------------------------------------------
library(forestplot)
library(dplyr)

data <- structure(list(mean  = c(NA,3.55,5.17,3.23,1.11,1.25,0.61,2.58,NA),
                       lower = c(NA,1.25,1.74,1.11,0.40,0.74,0.36,1.58,NA),
                       upper = c(NA,10.1,15.4,9.47,3.04,2.12,1.01,4.21,NA)),
                  .Names = c("mean", "lower", "upper"),
                  row.names = c(NA, -9L),
                  class = "data.frame")
data


tabletext <- cbind(c("", "Molecular Subgroup","","", "Age", "Gender","Extended Resections","Metastasis at Presentation",NA),
                   c("","SHH","group3","group4","Greater than 3","Female","GTR+NTR","YES",NA),
                   
                   c("Hazard Ratio",3.55,5.17,3.23,1.11,1.25,0.61,2.58,NA),
                   c("95% CI","(1.25,10.1)","(1.74,15.4)","(1.11,9.47)","(0.40,3.04)","(0.74,2.12)","(0.36,1.01)","(1.58,4.21)",NA)
                   #c("P-Value",0.017,0.007,0.003,0.020,0.500,0.600,"<0.001","")
)
tabletext

data %>%
  forestplot(labeltext = tabletext,
             #is.summary = c(rep(TRUE, 1), rep(FALSE, 7)),## Make test bold to highlight summary rows
             #clip = c(0.1, 3),## Lower and upper limit of clipping confidence intervals to arrows
             lty.ci=1,
             hrzl_lines=list("2"=gpar(lty=1,lwd=2),
                             "9"=gpar(lwd=2,colums=1:2,col="#000044")),
             xlog = TRUE,
             title="Forest Plot for Overall Survival",
             boxsize=0.1,
             xlab="Hazard Ratio",
             graphwidth = "auto",
             graph.pos= 4,
             align="l",
             new_page = TRUE,
             col = fpColors(box = "red",
                            line = "darkblue"))




