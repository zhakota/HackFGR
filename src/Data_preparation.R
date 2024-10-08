library(dplyr)
library(ggplot2)
library(ggsignif)
library(utils)
library(ggsci)
library(wesanderson)
library(meta)
library(strucchange)

#Source data acquisition
soure_data <- read.csv("data/raw/outleirs_to_na.csv") 


mean_source <- soure_data %>% group_by(gender, gestation) %>% summarize(Mean_wg = mean(BodyWeight), 
                                                                       sd_wg = sd(BodyWeight)) %>% 
  as.data.frame() %>% filter(gestation<41)


mean_source[mean_source$gender==1,1] <-'male'
mean_source[mean_source$gender==2,1] <-'female'

#Calculating mean body weight estimates for source data
source_overall <- mean_source %>% group_by(gestation) %>% 
  summarize(Mean_BodyWeight =mean(Mean_wg), SD_BodyWeight =mean(sd_wg ),
            Gender='overall',Method='Pathology', Study='Source') %>% as.data.frame() %>% unique()

colnames(source_overall)[1] <-'GA._w'

mean_source_renamed <- data.frame(GA._w=mean_source$gestation, Mean_BodyWeight =mean_source$Mean_wg,
                                  SD_BodyWeight=mean_source$sd_wg, Gender=mean_source$gender , Method='Pathology', Study='Source')


#Uploading data from publications
summary_data_raw <- read.csv('All_data_comparisions.csv', sep=";")
summary_data_raw$Mean_BodyWeight <- as.numeric(summary_data_raw$Mean_BodyWeight)

summary_data_raw <- summary_data_raw[!(summary_data_raw$Study=='Salomon' & summary_data_raw$Gender=='overall'),]


#Findind studies without mean estimates
unique(summary_data_raw[,c(4,6)]) #Kiserud, Kramer

#Calculating mean values for elected studies
Kiserud_overall <- summary_data_raw[summary_data_raw$Study=='Kiserud',] %>% group_by(GA._w) %>% 
  summarize(Mean_BodyWeight =mean(Mean_BodyWeight ), SD_BodyWeight =mean(SD_BodyWeight ),
            Gender='overall',Method=Method, Study=Study) %>% as.data.frame() %>% unique()
Kiramer_overall <- summary_data_raw[summary_data_raw$Study=='Kramer',] %>% group_by(GA._w) %>% 
  summarize(Mean_BodyWeight =mean(Mean_BodyWeight ), SD_BodyWeight =mean(SD_BodyWeight ),
            Gender='overall',Method=Method, Study=Study) %>% as.data.frame() %>% unique()

Intergrowth_overall <- summary_data_raw[summary_data_raw$Study=='Intergrowth21' & summary_data_raw$Method=='Neonatology',] %>% group_by(GA._w) %>% 
  summarize(Mean_BodyWeight =mean(Mean_BodyWeight ), SD_BodyWeight =mean(SD_BodyWeight ),
            Gender='overall',Method=Method, Study=Study) %>% as.data.frame() %>% unique()

Salomon_overall <- summary_data_raw[summary_data_raw$Study=='Salomon',] %>% group_by(GA._w) %>% 
  summarize(Mean_BodyWeight =mean(Mean_BodyWeight ), SD_BodyWeight =mean(SD_BodyWeight ),
            Gender='overall',Method=Method, Study=Study) %>% as.data.frame() %>% unique()

summary_data_raw <- rbind(summary_data_raw,Kiserud_overall, Kiramer_overall, 
                          mean_source_renamed,source_overall, Intergrowth_overall, Salomon_overall)

summary_data_raw$GA._w <- factor(summary_data_raw$GA._w)
summary_data_raw$Method <-factor(summary_data_raw$Method)
summary_data_raw$Gender <-factor(summary_data_raw$Gender)


#Plotting mean body weight classified according to the used method
ggplot(summary_data_raw, aes(fill=Method, y=Mean_BodyWeight, x=GA._w)) +
  geom_boxplot()+
  facet_wrap(~Gender)+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=16),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 12, hjust = 1))+
  xlab('Gestation week')+
  ylab("Weight")+
  scale_fill_lancet()
#5//10


#Plotting mean body weight estimates for each study regarding gender and used methods
ggplot(summary_data_raw, aes( y=Mean_BodyWeight, x=GA._w, color=Gender)) + #shape=Gender,
  geom_point(size=2, alpha=0.7)+
  scale_color_lancet()+
  facet_wrap(~Study)+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=16),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 12, hjust = 1))+
  ylab('Weight')+
  xlab("Gestation week")
#9//13

#Plotting mean body weight estimates for each study regarding gender and used methods in one plot
ggplot(summary_data_raw, aes(shape=Method, y=Mean_BodyWeight, x=GA._w, color=Study )) +
  geom_point(size=2, alpha=0.7)+
  facet_wrap(~Gender)+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=16),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 12, hjust = 1))+
  ylab('Weight')+
  xlab("Gestation week")+
  scale_color_manual(values = wes_palette("Darjeeling2", 12, type = "continuous"))
#6//12




#Comparing differences between genders for each gestation week within a particular method using t-test
gender_pval_df = data.frame(
  GA._w=c(),Method=c(),Pval=c())

for (gestation in levels(summary_data_raw$GA._w)){
  for (method in levels(summary_data_raw$Method)){
    dat_male <- summary_data_raw[summary_data_raw$GA._==gestation & summary_data_raw$Method==method & summary_data_raw$Gender=='male',2]
    dat_female <- summary_data_raw[summary_data_raw$GA._==gestation & summary_data_raw$Method==method & summary_data_raw$Gender=='female',2]
    if (length(na.omit(dat_male))>1 & length(na.omit(dat_female))>1){
      pval <- t.test(dat_male,dat_female)$p.value} else {
        pval<-NA
      }
    gender_subdf<-data.frame(GA._w=gestation, Method=method,
                             Pval=pval)
    gender_pval_df<-rbind(gender_pval_df,gender_subdf)
  }
}

#Adjusting p-values
gender_pval_df$adj_pval_BH <- p.adjust(gender_pval_df$Pval, method='BH')
gender_pval_df$adj_pval_FDR <- p.adjust(gender_pval_df$Pval, method='fdr')

#Plotting distributions of p-values using pair-wise t-test when comparing mean body weight between genders
ggplot(na.omit(gender_pval_df), aes(y=GA._w, x=adj_pval_FDR, color=Method)) +
  geom_point(size=2)+
  geom_vline(xintercept = 0.05, size=1.1, color='darkgrey')+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=16),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 12, hjust = 1))+
  ylab('Gestation week')+
  xlab("P-value")+
  scale_color_lancet()
#5//9

#Performing pair-wise comparisons of mean body weight for pairs of analysis methods using t-test for a certain gender or mean estimates
summary_pval_df = data.frame(
  GA._w=c(),Method1=c(),Method2=c(),Pval=c(), Gender=c())

compare_combinations <- t(combn(levels(summary_data_raw$Method),2))


for (gestation in levels(summary_data_raw$GA._w)){
  for (comare_ind in 1:nrow(compare_combinations)){
    for (gender in levels(summary_data_raw$Gender)){
    method1 <- compare_combinations[comare_ind,1]
    method2 <- compare_combinations[comare_ind,2]
    dat1 <- summary_data_raw[summary_data_raw$GA._==gestation & summary_data_raw$Method==method1 & summary_data_raw$Gender==gender ,2]
    dat2 <- summary_data_raw[summary_data_raw$GA._==gestation & summary_data_raw$Method==method2 & summary_data_raw$Gender==gender,2]
    if (length(na.omit(dat1))>1 & length(na.omit(dat2))>1){
    pval <- t.test(dat1,dat2)$p.value} else {
      pval<-NA
    }
    comp_subdf<-data.frame(GA._w=gestation, Method1=method1, Method2=method2,
                           Pval=pval, Gender=gender)
    summary_pval_df<-rbind(summary_pval_df,comp_subdf)
    }
  }
}


#Adjusting p-values
summary_pval_df$adj_pval_BH <- p.adjust(summary_pval_df$Pval, method='BH')
summary_pval_df$adj_pval_FDR <- p.adjust(summary_pval_df$Pval, method='fdr')

summary_pval_df$pair <- paste(summary_pval_df$Method1,summary_pval_df$Method2, sep='_')


#Plotting distribution of p-values using t-test for comparing pairs of used methods per gestation weak
ggplot(summary_pval_df, aes( y=GA._w, x=adj_pval_FDR, color=pair)) +
  geom_point(size=2)+
  geom_vline(xintercept = 0.05, size=1.1, color='darkgrey')+
  facet_wrap(~Gender)+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=16),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 12, hjust = 1))+
  ylab('Gestation week')+
  xlab("P-value")+
  scale_color_lancet()
#5//10


conf_interval <- function(mean_value, SD, n){
  standard_error <- SD / sqrt(n)
  alpha = 0.3
  degrees_of_freedom = n
  t_score = qt(p=alpha/2, df=degrees_of_freedom,lower.tail=F)
  margin_error <- t_score * standard_error
  
  lower_bound <- mean_value - margin_error
  upper_bound <- mean_value + margin_error
  return(data.frame(Low=lower_bound, Upp=upper_bound))
  
}

#Preparing data for meta analysis
for_metaanalysis <- summary_data_raw %>% group_by(GA._w,Gender, Method) %>% 
  summarize(Mean_BodyWeight=mean(Mean_BodyWeight,na.rm=TRUE),SD_BodyWeight = mean(SD_BodyWeight,na.rm=TRUE), n_sample=n()) %>% 
  as.data.frame() %>% na.omit()

for_metaanalysis <- cbind(for_metaanalysis, conf_interval(for_metaanalysis$Mean_BodyWeight, 
                                                          for_metaanalysis$SD_BodyWeight,
                                                          for_metaanalysis$n_sample))



#Testing methods for certain gestation weeks
test_overall <- for_metaanalysis[for_metaanalysis$Gender=='overall',]

m38 <- metamean(n_sample , Mean_BodyWeight , SD_BodyWeight, 
               studlab=Method, data=test_overall[test_overall$GA._w=='33',])
m31 <- metamean(n_sample , Mean_BodyWeight , SD_BodyWeight, 
                studlab=Method, data=test_overall[test_overall$GA._w=='31',])

forest(m31)
forest(m38)
#5//11


#Preparing the dataset with meta analysis results when comparing methods within a certain gestation week 
meta_summary_df <- data.frame(
  GA._w=c(),Pval=c(), Gender=c())

for (gestation in levels(summary_data_raw$GA._w)){
  for (gender in levels(summary_data_raw$Gender)){
    extr_df <- for_metaanalysis[for_metaanalysis$GA._==gestation & for_metaanalysis$Gender==gender,]
    meta_res <- metamean(n_sample , Mean_BodyWeight , SD_BodyWeight, 
                         studlab=Method, data=extr_df)
    meta_subdf<-data.frame(GA._w=gestation,
                           Pval=meta_res$pval.Q, Gender=gender)
    meta_summary_df <-rbind(meta_summary_df,meta_subdf)
  }
}

#Adjusting p-values
meta_summary_df$P_adj <- p.adjust(meta_summary_df$Pval, method='fdr')

#Plotting the distribution of p-values when comparing methods using meta analysis
ggplot(meta_summary_df, aes( y=GA._w, x=P_adj, shape=Gender)) +
  geom_point(size=3, color = '#937028')+
  geom_vline(xintercept = 0.05, size=1.1, color='darkgrey')+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=16),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 12, hjust = 1))+
  ylab('Gestation week')+
  xlab("P-value")
#5//7

#Building forest plot for a certain gestation week and gender according to raw p-values
m33 <- metamean(n_sample , Mean_BodyWeight , SD_BodyWeight, 
                studlab=Method, data=for_metaanalysis[for_metaanalysis$GA._w=='33' & 
                                                        for_metaanalysis$Gender=='overall',])
m30 <- metamean(n_sample , Mean_BodyWeight , SD_BodyWeight, 
                studlab=Method, data=for_metaanalysis[for_metaanalysis$GA._w=='30' & 
                                                        for_metaanalysis$Gender=='overall',])
forest(m33)
forest(m30)
#5//11

#Mean weights per different methods
ggplot(for_metaanalysis, aes( y=Mean_BodyWeight, x=GA._w, 
                              fill=Method, ymin=Low, ymax=Upp)) +
  geom_point(size=3, alpha=0.7,shape=21)+
  facet_wrap(~Gender)+ 
  geom_errorbar(width=.2)+
  geom_line(data = for_metaanalysis, aes(color=Method, y=Mean_BodyWeight, x=GA._w, group=Method))+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=16),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 12, hjust = 1))+
  ylab('Weight')+
  xlab("Gestation week")+
  scale_color_lancet()+scale_fill_lancet()

#5//10


#Time series regression Chow test for overall values
method_pairs <-   t(combn(levels(factor(for_metaanalysis$Method)),2))
Chow_overall_df <- data.frame(Method1=c(), Method2=c(), Gender=c(),F_val=c(),Pval=c())

for (gender in levels(for_metaanalysis$Gender)){
  for (Method_index in 1:nrow(method_pairs)){
    method1 <- method_pairs[Method_index,2]
    method2 <- method_pairs[Method_index,1]
    weigt_vec1 <- for_metaanalysis[for_metaanalysis$Gender==gender &
                                     for_metaanalysis$Method==method1 ,4]
    weigt_vec2 <- for_metaanalysis[for_metaanalysis$Gender==gender &
                                     for_metaanalysis$Method==method2 ,4]
    
    if (length(weigt_vec1)>0 & length(weigt_vec2)>0){
      min_length <- min(c(length(weigt_vec1),length(weigt_vec2)))
      weigt_vec1 <- weigt_vec1[1:min_length]
      weigt_vec2 <- weigt_vec2[1:min_length]
      diff_vec <- c(diff(weigt_vec1), diff(weigt_vec2))
      diff_frame <- data.frame(x= 1:length(diff_vec), dif_vals = diff_vec)
      sc_test <- sctest(diff_frame$dif_vals ~ diff_frame$x, type="Chow", point=min_length)
      sc_df <- data.frame(Method1=method1, Method2=method2,
                          Gender=gender, F_val = sc_test$statistic, Pval = sc_test$p.value)
      Chow_overall_df <- rbind(Chow_overall_df, sc_df)
    }
  }
}

#Add data for the same method
for (gender in levels(for_metaanalysis$Gender)){
  for(method in levels(factor(for_metaanalysis$Method))){
    
    sc_df <- data.frame(Method1=method, Method2=method,
                        Gender=gender, F_val = 0, Pval = 1)
    Chow_overall_df <- rbind(Chow_overall_df, sc_df)
  }
}


#Draw heatmaps based on F-values
ggplot(data = Chow_overall_df, aes(Method1, Method2, fill = F_val))+
  geom_tile(color = 'black')+
  scale_fill_gradient(low = "#DCDCDC", high = "#A60B0B", limit = c(0,7.5), space = "Lab", 
                      name="F value")+
  geom_text(aes(Method1, Method2, label = round(Pval,3)), color = "black", size = 4.5)+
  theme_bw() + xlab('Method') + ylab('Method') +
  theme( axis.text.y = element_text(color='black', 
                                    size=10),
         axis.title.y=element_text(color="black", 
                                   size=14),
         panel.background = element_rect(fill = "#DCDCDC"),
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_line(size = 0.1, linetype = 'dashed',
                                         colour = "black"),
         #legend.position = "none",
         axis.text.x = element_text(color='black', 
                                    angle = 60, vjust = 1, 
                                    size = 10, hjust = 1),
         axis.title.x = element_text( color="black", 
                                      size=14))+
  coord_fixed()+
  facet_wrap(~Gender)


#Meta-analysis per individual study

for_metamean_no_groupping <- summary_data_raw %>%
  group_by(GA._w,Study,Method, Gender) %>% 
  summarize(Mean_BodyWeight=mean(Mean_BodyWeight,na.rm=TRUE),SD_BodyWeight = mean(SD_BodyWeight,na.rm=TRUE), 
            n_sample=n(), Method = Method, Gender=Gender) %>% 
  as.data.frame() %>% na.omit()

m24_no <- metamean(n_sample , Mean_BodyWeight , SD_BodyWeight, 
                   studlab=Study, control=list(maxiter=1000),
                   data=for_metamean_no_groupping[for_metamean_no_groupping$GA._w=='24' & 
                                                    for_metamean_no_groupping$Gender=='overall',])
m27_no <- metamean(n_sample , Mean_BodyWeight , SD_BodyWeight, 
                   studlab=Study, control=list(maxiter=1000),
                   data=for_metamean_no_groupping[for_metamean_no_groupping$GA._w=='27' & 
                                                    for_metamean_no_groupping$Gender=='overall',])
m30_no <- metamean(n_sample , Mean_BodyWeight , SD_BodyWeight, 
                   studlab=Study, control=list(maxiter=1000),
                   data=for_metamean_no_groupping[for_metamean_no_groupping$GA._w=='30' & 
                                                    for_metamean_no_groupping$Gender=='overall',])
m33_no <- metamean(n_sample , Mean_BodyWeight , SD_BodyWeight, 
                   studlab=Study, control=list(maxiter=1000),
                   data=for_metamean_no_groupping[for_metamean_no_groupping$GA._w=='33' & 
                                                    for_metamean_no_groupping$Gender=='overall',])
forest(m24_no)
forest(m27_no)
forest(m30_no)
forest(m33_no)
#5//11


#Compare means for each individual study
for_metamean <- cbind(for_metamean_no_groupping, conf_interval(for_metamean_no_groupping$Mean_BodyWeight, 
                                                               for_metamean_no_groupping$SD_BodyWeight,
                                                               for_metamean_no_groupping$n_sample))

metameans_means <- for_metamean %>% group_by(Gender,GA._w ) %>% summarize(Mean_BodyWeight=mean(Mean_BodyWeight),
                                                n_sample=n(),
                                                SD_BodyWeight=mean(SD_BodyWeight)) %>% 
  as.data.frame()

ggplot(data=for_metamean[for_metamean$GA._w=='33',], 
       aes(y=Study, x=Mean_BodyWeight, xmin=Low, xmax=Upp, fill=Method)) +
  geom_point(shape=21, size=5) + 
  geom_errorbarh(height=.2)+
  geom_vline(data = metameans_means[metameans_means$GA._w=='33',], aes(xintercept = Mean_BodyWeight))+
  facet_wrap(~Gender)+
  scale_fill_lancet()+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=16),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 12, hjust = 1))+
  ylab('Study')+
  xlab("Body weight")

#5//10


#Plot KS test p-values when comparing genders from source data
pvals_source <- read.csv("data/raw/auxology_comparison_by_gender.tsv", sep ='\t') 

ggplot(pvals_source, aes( y=Week, x=pvals_source$KS.p.val_b2g, color=Organ)) +
  geom_point(size=2)+
  geom_vline(xintercept = 0.05, size=1.1, color='darkgrey')+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=16),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 12, hjust = 1))+
  ylab('Gestation week')+
  xlab("P-value")+
  scale_color_manual(values = wes_palette("Darjeeling2", 11, type = "continuous"))
