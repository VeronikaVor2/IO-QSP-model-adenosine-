####
#### Script by to make model-based simulations for AZD4635
#### Author:  Veronika Voronova

## clean workspace
rm(list = ls())
## step 1: upload libraries
library(tidyverse)
library(mlxR)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(mlxR); initMlxR("C:/ProgramData/Lixoft/MonolixSuite2019R1")
library(psych)
## function to convert mlxR simulations from short to long format
proc_out <- function(out, trt,ind){
  out_res <- as.data.frame(matrix(ncol = 8,nrow = 0))
  for (i in 1:(length(out) - 2)) {
    out_i <- out[[i]] %>% mutate(trt = trt,var = names(out)[i],value = .[,3]) %>% select(-3) %>% merge(.,ind,by = "id")
    colnames(out_res) <- colnames(out_i)    
    out_res <- rbind(out_res,as.data.frame(out_i))
  }
  return(out_res)
}
# define working directory - have to be changed to the file location directory
setwd("...") # please, select working directory
## read the model
model <- "A2A_clean_tv_sim.txt"
## read the parameters
param_0 <- read.csv("populationParameters.txt",sep = ",")
param_pop <- param_0$value
names(param_pop) <- param_0$parameter
## set additional parameters for alternative combinations to default values
param_pop <- c(param_pop,"cd73" = 0, "cytost" = 0, "cytotox" = 0, "vacc" = 0,"act" = 0,'Ag_fl' = 1)
## model-based simulations of TME (figure 2) -------------------------------------------------------------------

## define treatment groups
adm0   <- list(time = seq(from = 7,to = 22,by = 3.5), amount = 0,target = "Ac1")       ## placebo
adm1   <- list(time = seq(from = 7,to = 22,by = 3.5), amount = 0.125,target = "Ac1")  ## PD-L1
adm2   <- list(time = seq(from = 7,to = 22,by = 0.5), amount = 1.25,target = "Ad2")    ## AZD4635
## define covariates
ind <- list(name = c("Model","Study","sR","sL"))
mod <- data.frame(Study = c('CIV226','CIV227', 'CIV151', 'CIV258'), Model = c('CT26','MC38','MCA205','MCA205'))
##  30 animals per group for the simulations
param_ind <- expand.grid(id = 1:30,Study = mod$Study) %>% 
  merge(.,mod[,c("Study","Model")],by = "Study") %>% mutate(id = 1:nrow(.))
## define outputs
out <- list(list(name = c( "Tum","Ag","CTL","Ado_suppr","PDL1free","ISC","TKR","Ag_norm"), time = seq(0, 30, 0.1)))
### run simulations
res0  <- simulx(model = model,parameter = list(param_pop,param_ind),output = list(out),treatment = adm0)
res1  <- simulx(model = model,parameter = list(param_pop,param_ind),output = list(out),treatment = list(adm0,adm2)) 
res2  <- simulx(model = model,parameter = list(param_pop,param_ind),output = list(out),treatment = list(adm0,adm1)) 
res3  <- simulx(model = model,parameter = list(param_pop,param_ind),output = list(out),treatment = list(adm1,adm2)) 
## convert simulations to the single dataset

res_t <- do.call("rbind",list(proc_out(res0,"Vehicle",param_ind),
                            proc_out(res1,"AZD",param_ind),proc_out(res2,"PDL1",param_ind),proc_out(res3,"combo",param_ind)))

ressum <- res_t %>% group_by(time,trt,var,Study,Model) %>% summarize(ci05 = quantile(value,0.05),
                                                                     ci50 = quantile(value,0.5),
                                                                     ci95 = quantile(value,0.95)) %>% ungroup() %>% 
  mutate(trt = fct_relevel(trt, "Vehicle","AZD","PDL1","combo"), Study = fct_relevel(Study, "CIV226","CIV227","CIV258","CIV151"))

#### visualisation

labs_1 <- as_labeller(c(`CTL`="dTeff, cells/uL", `ISC`="ISC",`Ado_suppr`="Adenosine",`Ado_exp`="Adenosine, uM", 
                      `CD8_tot`="Total CD8, cells uL",
                      `Ag_norm`="Ag, 1/day",`PDL1free`="Free PD-L1",`PDL1`="total PD-(L)1" , `PRfunc`="IAR, -",`Tum`="Tumor vol, uL",
                      `Vehicle`="Vehicle",PD1Ab="PD-L1 mAb",`AZD4635`="AZD4635", `BID50`="PD-L1 Ab + \n 50 mg/kg AZD4635",
                      `CIV226`="CT26",`CIV227`="MC38",`CIV258`="MCA205: 1",`CIV151`="MCA205: 2",
                      `cd73`="Vmax ado",`vacc`="sL",`cytost`="sR"))

alt_col <- c('black','green3','blue','red')
ribbon_data <- ressum %>% group_by(var) %>% summarize(maxv = max(ci95))

f2a <- ggplot(ressum[ressum$var=="Tum",])+facet_grid(var~Study,scales="free",labeller=labs_1)+
  geom_rect(data = ribbon_data %>% filter(var == 'Tum'),aes(xmin = 7, xmax = 22,ymin = 0, ymax = maxv),alpha = 0.2) +
  geom_ribbon(aes(x=time, ymin=ci05,ymax=ci95,group=trt,fill=trt),alpha=0.2)+
  geom_line(aes(x=time, y=ci50,group=trt,color=trt,linetype = trt))+theme_bw(base_size=12)+
  xlab("Time, days")+ylab("")+
  scale_color_manual(values = alt_col) + scale_fill_manual(values = alt_col) +ggtitle("C.Tumor size dynamics")+
  scale_linetype_manual(values = c('solid','dashed','solid','dashed')) +
  theme(plot.title = element_text(size=16),legend.position = "none")

f2b <- ggplot(ressum[ressum$var %in% c("Ado_suppr","PDL1free","ISC"),])+facet_grid(var~Study,scales="free",labeller=labs_1)+
  geom_rect(data = ribbon_data %>% filter(var %in% c("Ado_suppr","PDL1free","ISC")),aes(xmin = 7, xmax = 22,ymin = 0, ymax = maxv),alpha = 0.2) +
  geom_line(aes(x=time, y=ci50,group=trt,color=trt,linetype = trt))+theme_bw(base_size=12)+
  geom_ribbon(aes(x=time, ymin=ci05,ymax=ci95,group=trt,fill=trt),alpha=0.2)+
  xlab("Time, days")+ylab("")+
  scale_color_manual(values = alt_col) + scale_fill_manual(values = alt_col)+ggtitle("B.Dynamics of immunosupressive components")+
  scale_linetype_manual(values = c('solid','dashed','solid','dashed')) +
  theme(plot.title = element_text(size=16),legend.position = "none")

f2c <- ggplot(ressum[ressum$var %in% c("CTL"),])+facet_grid(var~Study,scales="free",labeller=labs_1)+
  geom_rect(data = ribbon_data %>% filter(var == 'CTL'),aes(xmin = 7, xmax = 22,ymin = 0, ymax = maxv),alpha = 0.2) +
  geom_line(aes(x=time, y=ci50,group=trt,color=trt))+theme_bw(base_size=12)+
  geom_ribbon(aes(x=time, ymin=ci05,ymax=ci95,group=trt,fill=trt),alpha=0.2)+
  xlab("Time, days")+ylab("")+
  scale_color_manual(values = alt_col) + scale_fill_manual(values = alt_col)+ggtitle("A.Lymphocytes dynamics")+
  scale_linetype_manual(values = c('solid','dashed','solid','dashed')) +
  theme(plot.title = element_text(size=16),legend.position = "none")

### differnces in the estimated parameters between the individuals
#param_ind_0<-read.csv("D:\\Personal_folders\\Veronika\\1. IO\\1. Work\\10. AZD4635\\7. Manuscript\\MODELS\\ADOMLX_cleanfit\\run_00_final_refit\\IndividualParameters\\estimatedIndividualParameters.txt")
param_ind_0<-read.csv("estimatedIndividualParameters.txt")
param_ind<-param_ind_0 %>% mutate(sR_SAEM=1/sR_SAEM,sL_SAEM=1/sL_SAEM) %>%
  select(Study,sR_SAEM,sL_SAEM, Vado_SAEM)%>%gather(key="variable",value="value",-Study)

param_ind$Study<-factor(param_ind$Study,levels=c("CIV226","CIV227","CIV258","CIV151"))

signif<-data.frame(pval=c(
  wilcox.test(sL_SAEM~Study,data=param_ind_0[param_ind_0$Study%in%c("CIV258","CIV151"),])$p.value,
  wilcox.test(sL_SAEM~Study,data=param_ind_0[param_ind_0$Study%in%c("CIV226","CIV227"),])$p.value,
  wilcox.test(sR_SAEM~Study,data=param_ind_0[param_ind_0$Study%in%c("CIV258","CIV151"),])$p.value,
  wilcox.test(sR_SAEM~Study,data=param_ind_0[param_ind_0$Study%in%c("CIV226","CIV227"),])$p.value),
  variable=c("sL_SAEM","sL_SAEM","sR_SAEM","sR_SAEM"),x=c(3.5,1.5,3.5,1.5),y=c(1,1,0.03,0.03))

f2d<-ggplot(param_ind,aes(y=value,x=Study))+
  facet_wrap(~variable,ncol=1,scales="free",labeller=as_labeller(c("sL_SAEM"="1/sL, day/uL\n(nTeff flux sensitivity)",
                                                                   "sR_SAEM"="1/sR, day/uL\n(ISC flux sensitivity)",
                                                                   "Vado_SAEM"="Vmaxado, umol/day\n(Adenosine flux)")))+
  geom_boxplot()+geom_jitter(height = 0,width=0.2,size=1,alpha=0.5)+
  theme_bw(base_size=12)+xlab("")+ylab("")+
  scale_x_discrete(labels=c(`CIV226`="CT26",`CIV227`="MC38",`CIV258`="MCA205 \n 1",`CIV151`="MCA205 \n 2"))+ggtitle("D. Estimated cell flux\nparameters")+ 
  theme(plot.title = element_text(size=16),plot.margin=unit(c(0.1,0.1,5,0.1),"cm"))+
  geom_text(data=signif,aes(x=x,y=y,label=paste("p-val=",signif(pval,2))))

### evaluation of statistical differences between the groups on the day 30 (figure s4) ------------------------

my_comparisons <- list( c("PDL1", "combo"), c("AZD", "combo")) 
lb <- as_labeller(c(`CIV258`="CIV258 \n MCA205: 2",`CIV226`="CIV226 \n CT26",`CIV151`="CIV151 \n MCA205",`CIV227`="CIV227 \n MC38"))
##

figs5_2 <- ggplot(res_t %>% filter(time == 30 & var == 'Tum') %>% mutate(trt = as.factor(trt), 
                                                                         trt = fct_relevel(trt, 'Vehicle','AZD','PDL1','combo')),
                  aes(x = trt, y = value)) + facet_wrap(~ Study, labeller = lb) +
  geom_boxplot(outlier.alpha = 0.001) + geom_jitter(width = 0.1, height = 0) + xlab('') + ylab(expression(Tumor~volume~mm^{3})) +
  scale_x_discrete(labels = c('Vehicle', 'AZD4635','PD-L1 mAb', 'Combo\n50 mg/kg')) +theme_bw(base_size = 16) +
  stat_compare_means(comparisons = my_comparisons, method = 't.test',alternative = "two.sided") + ylim(0,3500)

######validation using CFM data - figure s5 ##### ---------------------------------------------------

valid_model<-res_t%>% filter(time>=23&time<23.1&var%in%c( "Ag_norm","CTL")&Study=="CIV227")%>%
  group_by(var)%>%mutate(median_pbo=median(value[trt=="Vehicle"]))%>%ungroup()%>%
  mutate(fold=value/median_pbo,process=recode(var,`Ag_norm`="APC",`CTL`="CTL")) %>%
  filter(trt %in% c("Vehicle","PDL1","AZD","combo"))%>%
  mutate(Group=recode(trt,"Vehicle"="Vehicle",'PDL1'="PD1Ab","AZD"="AZD4635",'combo'="BID50"),Variable="model")%>%
  select(Model,Variable,Group,fold,process)

valid_data<-read.csv("FlowValidation_AZpaper.csv")
valid_data<-valid_data%>% group_by(Variable)%>%mutate(mean_pbo=median(Value[Group=="Vehicle"]))%>%
  mutate(fold=Value/mean_pbo)%>%filter(Variable!="CD8_perc") %>%
  mutate(process=recode(Variable,`CD86_DC`="APC",`MHCII_DC`="APC",`CD86_M`="APC",`MHCII_M`="APC",`CD8_PD1_perc`="CTL"))%>%
  select(Model,Variable,Group,fold,process)

validation<-rbind(as.data.frame(valid_model),as.data.frame(valid_data))
validation$Group<-factor(validation$Group,levels=c("Vehicle","PD1Ab","AZD4635","BID50"))

## visualise validation result

f2d<-ggplot(validation)+facet_wrap(~process,scales="free",labeller=as_labeller(c("CTL"="CTL activation","APC"="APC activation")))+
  geom_boxplot(aes(x=Group ,fill=Variable,y=fold))+theme_bw(base_size=12)+
  ylab("vehicle-normalized value")+xlab("")+
  scale_fill_manual(values=c("salmon","lightpink","salmon","red","orange1","lightblue","lightblue"),
                    labels=c("PD1+CD8 cells, % from live (data)","CD86 on DC (MFI) (data)","CD86 on M (MFI) (data)",
                             "MHC II on DC (MFI) (data)","MHC II on M (MFI) (data)",
                             "dTeff or APC (model)"))+
  scale_x_discrete(labels=c("Vehicle","PD-L1 Ab","AZD4635","Combo"))+ggtitle("Model validation")+ 
  theme(plot.title = element_text(size=16))

#### evaluation of between-animal variability, fig 3 --------------------------------------------
## select syngeneic model for evaluation
param_ind<- data.frame(id=1:100,Study="CIV258", Model="MCA205")
out <- list(list(name =c( "Tum","Ag","CTL","Ado_suppr","PDL1","ISC","TKR","Ag_norm","CD8_tot","sL_out","sR_out"), time = seq(0, 30, 0.1)))

### run simulations
res_258  <- simulx(model=model,parameter=list(param_pop,param_ind),output=list(out),treatment=list(adm1,adm2))
res_t258<-proc_out(res_258 ,"combo",param_ind)

## define responders and non-responders
resp<-res_t258%>% filter(var=="Tum")%>% filter(time==max(time)|(time>=7&time<7.1)) %>% 
  spread(key=time,value=value)
colnames(resp)[c(6,7)]<-c("start","end")
resp<-resp %>% mutate(resp=ifelse(start>end,"resp","prog")) %>% dplyr::select(id,resp)
res_t258_2<-merge(res_t258,resp,by="id")

## make simulations for placebo treatment
res_258_pbo  <- simulx(model=model,parameter=list(param_pop,param_ind),output=list(out),treatment=adm0) 
res_t258_pbo<-proc_out(res_258_pbo ,"Vehicle",param_ind) %>%mutate(resp="vehicle")
res_t258_3<-rbind(res_t258_2,res_t258_pbo)

## merge all data
ressum_t258<-res_t258_3%>%group_by(time,resp,var) %>%
  summarize(ci05=quantile(value,0.05),ci50=quantile(value,0.5),ci95=quantile(value,0.95))%>% ungroup()%>%
  mutate(resp=factor(resp,levels=c("vehicle","resp","prog")))%>%
  filter(var %in% c("Tum","ISC","CTL","Ado_suppr","PDL1","Ag_norm","CD8_tot"))

## make TME 'snapshots'
ressum_t258_bp<-res_t258_3 %>%filter(var %in% c("Tum","ISC","CTL","Ado_suppr","PDL1","CD8_tot","Ag_norm")&
                                       (time>=14&time<14.1)|(time>=7&time<7.1))%>%
  mutate(resp=factor(resp,levels=c("vehicle","prog","resp")))

param_resp<-res_t258_3 %>% dplyr::select(id,var,value,resp) %>%filter(var %in% c("sL_out","sR_out"))%>% 
  unique()%>% filter(resp!="vehicle") %>% mutate(resp=factor(resp,levels=c("resp","prog")))

## visualise the result
f3a<-ggplot(ressum_t258[ressum_t258$var==("Tum"),])+# facet_wrap(~var,scales="free",labeller=labs_1)+
  geom_line(aes(x=time, y=ci50,group=resp,color=resp))+theme_bw(base_size=12)+
  geom_ribbon(aes(x=time, ymin=ci05,ymax=ci95,group=resp,fill=resp),alpha=0.2)+
  xlab("Time, days")+ylab("Tumor volume, uL")+
  scale_colour_manual(values=c("blue","forestgreen","darkred"),name="Group",labels=c("vehicle","non-progressors","progressors"))+ 
  scale_fill_manual(values=c("blue","forestgreen","darkred"),name="Group",labels=c("vehicle","non-progressors","progressors"))+
  ggtitle("A.Tumor  dynamics")+
  theme(plot.title = element_text(size=16),legend.position = "right")

f3b<-ggplot(param_resp)+facet_wrap(~var,scales="free",label=as_labeller(c("sL_out"="nTeff flux","sR_out"="ISC flux")))+ 
  geom_density(aes(x=1/value,fill=resp),position="identity",alpha=0.5)+
  theme_bw(base_size=12)+ggtitle("B.Estimated fluxes")+ xlab("")+
  theme(plot.title = element_text(size=16),legend.position = "none")+
  scale_fill_manual(values=c("forestgreen","darkred"),name="Group",labels=c("non-progressors","progressors"))

f3c<-ggplot(ressum_t258[(ressum_t258$var%in% c("CTL","ISC","PDL1","CD8_tot")),])+
  facet_wrap(~var,scales="free",labeller=labs_1,ncol=4)+
  geom_line(aes(x=time, y=ci50,group=resp,color=resp))+theme_bw(base_size=12)+
  geom_ribbon(aes(x=time, ymin=ci05,ymax=ci95,group=resp,fill=resp),alpha=0.2)+
  xlab("Time, days")+ylab("value")+
  scale_colour_manual(values=c("blue","forestgreen","darkred"),name="Group",labels=c("vehicle","non-progressors","progressors"))+ 
  scale_fill_manual(values=c("blue","forestgreen","darkred"),name="Group",labels=c("vehicle","non-progressors","progressors"))+
  ggtitle("C.Biomarker dynamics")+
  theme(plot.title = element_text(size=16),legend.position = "none")


f3d<-ggplot(ressum_t258_bp[(ressum_t258_bp$var%in% c("CTL","ISC","PDL1","CD8_tot")),])+
  facet_wrap(~var,scales="free",labeller=labs_1,ncol=4)+
  geom_boxplot(aes(x=as.factor(time),y=value,fill=resp),alpha=0.4)+theme_bw(base_size=12)+
  scale_fill_manual(values=c("blue","darkred","forestgreen"),name="Group",labels=c("vehicle","progressors","non-progressors"))+
  scale_x_discrete(labels=c("baseline","post-treat"))+xlab("")+
  ggtitle("D.Biomarker levels")+
  theme(plot.title = element_text(size=16),legend.position = "none")

############ sensitivity analysis (figure 4A) --------------------------------------
## create a function for simulations processing
proc_out2<-function(out, trt){
  out_res<-as.data.frame (matrix(ncol=4,nrow=0))
  for (i in 1:(length(out)-2)){
    out_i<-out[[i]]%>%mutate(trt=trt,var=names(out)[i],val=.[,2]) %>% select(-2)
    colnames(out_res)<-colnames(out_i)    
    out_res<-rbind(out_res,as.data.frame(out_i))
  }
  return(out_res)
}
#param_0<-read.csv("run_00_final_refit_scf1_tkr\\populationParameters.txt",sep=",")
param_0<-read.csv("populationParameters.txt",sep=",")

param_pop<-param_0$value
names(param_pop)<-param_0$parameter
## set other treatment parameters to default
param_pop<-c(param_pop,"cd73"=0, "cytost"=0, "cytotox"=0, "vacc"=1, 'act'=0,'Ag_fl'=1)
## set variabilities to 0 to have population simulations
param_pop[c('omega_sL','omega_sR','b1')] <- c(1e-6,1e-6,1e-6)
### define parameters to test
par_var <- c('sL_pop','sR_pop','kLn_pop','Vado_pop','Kp_pop', 'beff_pop','r_pop','Ag_fl')#,'vacc'

adm0   <- list(time = seq(from=7,to=22,by=3.5), amount = 0,target="Ac1")       ## placebo
adm1   <- list(time = seq(from=7,to=22,by=3.5), amount = 2*6.66,target="Ac1")  ## PD-L1
adm2   <- list(time = seq(from=7,to=22,by=0.5), amount = 1.25,target="Ad2")    ## AZD4635
## define range for the parameter change 
scF <- 2

## make model simulations with changed parameters
sensitivity <- as.data.frame(matrix(ncol=9, nrow=0))
for (i in 1:length(par_var)){
  param_ind_i <- data.frame(id=1:2,param_i=c(param_pop[par_var[i]]/scF,param_pop[par_var[i]]*scF),Study="CIV227", Model="MC38")
  # param_ind_i <- data.frame(id=1:2,param_i=c(param_pop[par_var[i]]/scF,param_pop[par_var[i]]*scF),Study="CIV258", Model="MCA205")
  colnames(param_ind_i)[2]<- par_var[i]
  param_ind_i2 <-param_ind_i %>% rename(par = par_var[i])
  param_pop_i <- param_pop[setdiff(names(param_pop),par_var[i])]
  res_im  <- simulx(model=model,parameter=list(param_ind_i,param_pop_i),
                    output=list(list(name =c( "Tum","Ag","CTL","ISC","PRfunc","PDL1free"), time = seq(0, 30, 0.1))),
                    treatment=list(adm0,adm2))
  res_ic  <- simulx(model=model,parameter=list(param_ind_i,param_pop_i),
                    output=list(list(name =c( "Tum","Ag","CTL","ISC","PRfunc","PDL1free"), time = seq(0, 30, 0.1))),
                    treatment=list(adm1,adm2))
  
  pred_im_i<-rbind(proc_out(res_im,'mono',param_ind_i2),proc_out(res_ic,'combo',param_ind_i2)) %>% 
    mutate(param=par_var[i])
  colnames(sensitivity) <- colnames(pred_im_i)
  sensitivity <- rbind(sensitivity,pred_im_i)
}

## simulations for default parameter values
res_imd  <- simulx(model=model,parameter=list(param_pop,data.frame(id=1,Study="CIV227", Model="MC38")),#Study="CIV227", Model="MC38")),
                   output=list(list(name =c( "Tum","Ag","CTL","ISC","PRfunc","PDL1free"), time = seq(0, 30, 0.1))),
                   treatment=list(adm0,adm2)) 
res_icd  <- simulx(model=model,parameter=list(param_pop,data.frame(id=1,Study="CIV227", Model="MC38")),
                   output=list(list(name =c( "Tum","Ag","CTL","ISC","PRfunc","PDL1free"), time = seq(0, 30, 0.1))),
                   treatment=list(adm1,adm2))
## extract variables 
pred_imd<-rbind(proc_out2(res_imd,'mono'),proc_out2(res_icd,'combo'))
pred_sum <- pred_imd%>%filter(var=='Tum'&time>29.9)
pred_sum2 <- pred_imd%>%filter(var=='CTL') %>% group_by(trt) %>% summarize(PR_mean=mean(val))
pred_sum3 <- pred_imd%>%filter(var=='ISC') %>% group_by(trt) %>% summarize(PR_mean=mean(val))

### summarize the result
sens_sum <- sensitivity %>%filter(var=='Tum'&time>29.9) %>%   #group_by(id, param,var) %>% summarize()
  merge(.,pred_sum[,c('val','trt')],by='trt') %>% mutate(ch_bln=(value-val)/val*100) %>% 
  mutate(id=recode(id,'1'='-50%','2'='+50%'),trt=as.factor(recode(trt,'mono'='AZD4635 alone','combo'='AZD4635 + anti-PD-L1 mAb')),
         xlabs=as.factor(gsub('_pop',"",param)),xlabs=fct_reorder(xlabs, desc(ch_bln)), 
         trt=fct_relevel(trt, 'AZD4635 alone','AZD4635 + anti-PD-L1 mAb' )) 
## visualise the result
fs4a <- ggplot(sens_sum, aes(xlabs, ch_bln, fill=as.factor(id))) +facet_wrap(~trt)+theme_bw(base_size=12)+
  geom_bar(position="identity", stat="identity",alpha=0.5,color='black') + coord_flip()+
  scale_x_discrete(labels=c('kLn'='kLn','beff'='beff','r'='r','sR'='sR','sL'= 'sL','Ag_fl'= expression("Ag"["norm"]),'Vado'=expression("Vmax"["Ado"]),'Kp'='Kp'))+
  xlab('parameter')+ylab('Tumor volume change from the default simulation on day 30, %')+
  scale_fill_manual(values=c('red','blue'),name='parameter\nchange')

############ alternative treatment combinations (figure 4B) --------------------------------------
## define model parameters
out<-list(list(name =c( "Tum","CTL","ISC"), time = seq(0, 28, 0.1)))
param_0<-read.csv("run01_finalref\\populationParameters.txt",sep=",")
param_pop<-param_0$value
names(param_pop)<-param_0$parameter
param_pop<-c(param_pop,"cytotox"=0,'Ag_fl'=1) #,"cd73"=0, "cytost"=0,  "vacc"=1, 'act'=0

#param_pop<-c(param_pop, "cytotox"=0) #"cd73"=0, "cytost"=0, , 'act'=0, "vacc"=0
noa <- 10 # number of animals # 
## create a dataframe, specifying parameters for various treatment options
param_ind <- data.frame(id=1:(5*noa),Study="CIV227", Model="MC38",
                        cd73=rep(c(0,1,0,0,0),each=noa),act=rep(c(0,0,1,0,0),each=noa),
                        cytost=rep(c(0,0,0,1,0),each=noa),vacc=rep(c(0,0,0,0,1),each=noa))
adm0   <- list(time = seq(from=7,to=22,by=3.5), amount = 0,target="Ac1")       ## placebo
adm1   <- list(time = seq(from=7,to=22,by=3.5), amount = 2*6.66,target="Ac1")  ## PD-L1
adm2   <- list(time = seq(from=7,to=22,by=0.5), amount = 1.25,target="Ad2")    ## AZD4635

## run 100 model simulations with 10 animals per group
res_s<-data.frame(matrix(ncol=9,nrow=0))
for (i in 1:100){
  res_ip  <- simulx(model=model,parameter=list(param_pop,param_ind),
                    output=out,treatment=adm0)
  res_im  <- simulx(model=model,parameter=list(param_pop,param_ind),
                    output=out,treatment=list(adm0,adm2))
  res_ic  <- simulx(model=model,parameter=list(param_pop,param_ind),
                    output=out,treatment=list(adm1,adm2))
  
  pred_im_i<-do.call('rbind',list(proc_out(res_ip,'plac',param_ind[,1:4]),
                                  proc_out(res_im,'mono',param_ind[,1:4]),proc_out(res_ic,'combo',param_ind[,1:4]))) %>% 
    mutate(simn=i)
  colnames(res_s) <- colnames(pred_im_i)
  res_s <- rbind(res_s,pred_im_i)
}
## process the simulations
res_s <- res_s %>%select(-cd73) %>%  mutate(varpar=ifelse(id %in% 1:10,'no',ifelse(id %in% 11:20,'cd73',ifelse(id %in% 21:30,'act',
                                                                                                               (ifelse(id %in% 31:40,'cytost','vacc'))))))
res_s_sum1 <- res_s %>%filter(var=="Tum"&time>27.9) %>% group_by(time,trt,var,simn,varpar) %>% 
  mutate(TV_m=geometric.mean(value)) %>% group_by(simn,var,varpar) %>%
  mutate(TGI=(1-TV_m/TV_m[trt=='plac'])*100) %>% 
  group_by(trt,var,varpar) %>% summarize(mean_tg=median(TGI),ci05=quantile(TGI,0.05),ci95=quantile(TGI,0.95)) %>% ungroup() %>% 
  mutate(varpar=fct_relevel(varpar,'act','cytost','cd73','vacc','no'),trt=fct_relevel(trt,'no','mono','combo'))

res_s_sum1a <- res_s %>%filter(var=="Tum"&time>27.9) %>% group_by(time,trt,var,simn,varpar) %>% 
  filter(value<100) %>% group_by(trt,simn,varpar) %>% count() %>% mutate(rr=n/10*100) %>% ## define number of responders
  group_by(trt,varpar) %>% summarize(mean_rr=mean(rr),ci05=quantile(rr,0.05),ci95=quantile(rr,0.95)) %>% ungroup() %>% 
  mutate(varpar=fct_relevel(varpar,'act','cytost','cd73','vacc','no'),trt=fct_relevel(trt,'no','mono','combo'))
## visualise the result
f4b <- ggplot(res_s_sum1%>% filter(trt!='plac'))+#facet_wrap(~var)+# +
  geom_tile(aes(x=trt,y=varpar,fill=mean_tg),color="black")+theme_bw(base_size=16)+
  scale_fill_gradientn(colours=c('red3','yellow',"green","forestgreen"),limits=c(0,100),name="TGI, %")+ #,"darkblue"
  geom_text(aes(x = trt, y = varpar,
                label =  paste(round(mean_tg),"\n (",round(ci05),"-",round(ci95),")",sep="")),
            size=5,color="black")+
  scale_y_discrete(labels=c('no'='+ vehicle','vacc'='+vaccine','cd73'='+CD73 Ab','cytost'='+ ISC bl.','act'='+ACT'))+
  scale_x_discrete(labels=c('mono'='AZD4635 \n alone', 'combo'='AZD4635+ \n anti-PD-L1 mAb'))+xlab('')+ylab('')

