
just_n_save<-function(file_name){
  ggsave(paste(file_name,'.png',sep=''),width = 10,height = 10,dpi=600)
}

plot_n_save<-function(dsets,dset_names,dir_suffix='',vari,legend_title,lower_q=NULL,upper_q=NULL,digits=3){
  
  pal<-(rcartocolor::carto_pal(7,'Emrld'))
  cols_pal<- colorRampPalette(
    colors = pal)
  point_cols<- cols_pal(20)
  
  dummy_frame<-data.frame(x=1:20,y=1:20,point_cols,point_vals=rep(NA,20))
  
  dir_name<-paste('./Figures/simplex_',vari,'_plots',dir_suffix,sep='')
  
  
  min_val<-min(all_vals)-(1e-10)
  max_val<-max(all_vals)+(1e-10)
  smallest<-NULL
  largest <-NULL
  l_out <- 21 ##The number of boundaries for our color pallete
  if(!is.null(lower_q)){
    smallest<-min_val
    min_val<-quantile(all_vals,lower_q)
    l_out<-l_out-1
  }
  if(!is.null(upper_q)){
    largest<-max_val
    max_val<-quantile(all_vals,upper_q)
    l_out<-l_out-1
  }
  
  
  val_breaks<-c(smallest,seq(min_val,max_val,length.out=l_out),largest)
  dummy_frame$point_vals<-val_breaks[-1]
  dummy_frame$point_vals[nrow(dummy_frame)]<-dummy_frame$point_vals[1]

  legend_plot<- ggplot(data=dummy_frame,mapping=aes(x=x,y=y,fill=point_vals),limits=c(min_val,max_val))+
    geom_point()+
    scale_fill_gradientn(colours=point_cols,
                         breaks=round(seq(min_val,max_val,length.out=5),digits = digits),
                         limits=round(c(min_val,max_val),digits=digits))+
    theme(
      legend.title=element_text(size=25),
      legend.key.size = unit(0.1,'npc'),
      legend.text = element_text(size=20))
  legend_plot$labels$fill<- legend_title
  legend_plot<-ggdraw(get_legend(legend_plot))
  

  dir.create(dir_name)
  
  print(legend_plot)
  ggsave(paste(dir_name,'/',vari,'_legend.png',sep=''),width = 10,height = 10,dpi=600)
   
  for(i in 1:length(dsets)){
    dset<-dsets[[i]]
    dset_name<-dset_names[i]
    
    dset_cols <-point_cols[as.numeric(cut(unlist(dset[,vari]),breaks=val_breaks))]
    dset_plot <-make.ggtern.points(dset,dset_cols)
    
    print(dset_plot)
    ggsave(paste(dir_name,'/',dset_name,'_',vari,'.png',sep=''),width = 10,height = 10,dpi=600)
  }
 }





##############################
##Line plots with changing s##
##############################


########################
## Generate LTT plots ##
########################

## NOTE! THESE PLOTS ARE NO LONGER USED
## The function for expected ltt is not correct


# #function for N(T)
# ltt <- function(time,r,s,n0=2){
#   numer <- r-s
#   denom1<- (n0*s+r-s)/n0
#   denom2<-exp(-(r-s)*time)
#   denom<-denom1*denom2-s
#   lineages<-numer/denom
#   return(lineages)
# }
# 
# times<-seq(0,1,by=0.0002)
# divs <-seq(-0.8,0.8,by=0.2)
# ltt_vals<-matrix(nrow = length(times)*length(divs),ncol = 3,data = NA)
# for(i in 1:length(divs)){
#   s<-divs[i]
#   ltt_vals[((i-1)*length(times)+1):((i)*length(times)),1]<-times
#   ltt_vals[((i-1)*length(times)+1):((i)*length(times)),2]<-ltt(times,12,s)
#   ltt_vals[((i-1)*length(times)+1):((i)*length(times)),3]<-rep(s,length(times))
# }
# ltt_vals<-as.data.frame(ltt_vals)
# colnames(ltt_vals)<-c('time','lineages','s')
# ltt_vals<-ltt_vals[ltt_vals$lineages>=0,]
# ltt_vals$s<-(ltt_vals$s)
# 
# ltt.plot<-ggplot(data=ltt_vals,mapping=aes(x=time,y=lineages,fill=factor(s,levels = (unique(s)))))+
#   geom_line() + xlim(0,1) +ylim(0,150) +
#   geom_line(aes(color=factor(s,levels = rev(unique(s)))),size=4)+
#   scale_color_manual(values=rev(line_col)) +
#   theme(
#    #legend.position = 'none',
#     legend.title=element_text(size=15),
#     legend.key.size = unit(0.1,'npc'),
#     legend.text = element_text(size=15),
#     axis.title = element_text(size=15))+
#   xlab("Time") +
#   ylab("Number of Lineages N(t)")+
#   guides(color=guide_legend(title= expression("s = ("*nu['+']*- nu['-']*')')))
# 
# # 
# ltt.plot
# leg <- ggdraw(get_legend(ltt.plot))
# leg
# ggsave('div_n_leg.png',width = 10,height = 10,dpi=600)


##################################################
##############  Load Core Datasets ###############
##################################################

source('Summarize_data.R')

##########################################################
##########################################################

library(rcartocolor)
library(cowplot)
library(gridGraphics)
pal<-(rcartocolor::carto_pal(7,'Emrld'))
cols_pal<- colorRampPalette(
  colors = pal)
point_cols<- cols_pal(20)

dummy_frame<-data.frame(x=1:20,y=1:20,point_cols,point_vals=rep(NA,20))



###########################
### Make summary tables ###
###########################

##Add identifier for each dataset
c_dat$ID<-rep('Complete',nrow(c_dat))
r_dat$ID<-rep('Extant-Only',nrow(c_dat))
c_gd_dat$ID<-rep('Distance Dependence',nrow(c_dat))
r_i_dat$ID<-rep('Incomplete Sampling',nrow(c_dat))

all_dat<-rbind(c_dat,r_dat,c_gd_dat,r_i_dat)
all_dat<-all_dat[,-c(1,2,16,17,18)] ##Get rid of unnessecary columns

all_sum <- all_dat %>%
  group_by(ID) %>%
  summarise(across(everything(),list(mean=mean,var=var,min=min,max=max)))

div_sum <- all_sum[,c(1,18:33)]

xtable(class_sum<-all_sum[,c(1,6:13,2:5,14:17)])


###################
### Level Plots ###
###################

ggplot(simp_frame[(!simp_frame$recon) & (simp_frame$rets<=10),],mapping=aes(x=rets,level,group=rets))+
  geom_violin(adjust=2)+scale_y_continuous(breaks=1:10)+
  scale_x_continuous(breaks=1:10)+
  labs(x='Number of Reticulations',y="Level")+
  theme(text =  element_text(size=20))
just_n_save('level_per_ret10')

ave_levl <- simp_frame[(!simp_frame$recon),] %>%
  group_by(rets) %>%
  summarize(lvl = mean(level))

ggplot(ave_levl,mapping=aes(x=rets,lvl))+
  geom_point()+
  labs(x='Number of Reticulations',y=" Average Level")+
  theme(text =  element_text(size=20))
just_n_save('level_per_ret_all')


##############################################################
##############################################################
######                                                  ######
###### Genetic Distance and Incomplete Sampling Section ######
######                                                  ######
##############################################################
##############################################################




################################
### Load incomplete sampling ###
################################


gathered_if<- if_dat %>% 
  select(frac,ave_rets,rets0,ratio,tb,tc,fus,nor) %>%
  gather(key=class,value=prop,-c(frac,ave_rets,ratio,rets0))


incom_by_rets<-incom_frame %>% 
  filter(rets<=10) %>%
  group_by(frac,rets) %>%
  summarize(tb = sum(tree_based)/n(),
            fus= sum(fu_stable)/n(),
            tc = sum(tree_child)/n(),
            nor= sum(normal)/n())


incom_dir<-'incom_plots'
dir.create(incom_dir)



ggplot(data=gathered_if, mapping = aes(x=frac,y=prop,color=factor(class,levels = c('tb','fus','tc','nor'))))+
  geom_point(size=5) +
  scale_color_discrete(name='Class',labels=c('Tree-Based','FU-Stable','Tree-Child','Normal'))+
  xlab(expression('Sampling Fraction '*rho))+
  ylab('Proportion in class')+
  theme(
    legend.title = element_text(size=30), #change legend title font size
    legend.text = element_text(size=20),
    axis.title=element_text(size=30,face="bold"))
ggsave(paste(incom_dir,'/incon_props.png',sep=''),width = 10,height = 10,dpi=600)



ggplot(data=gathered_if, mapping = aes(x=frac,y=rets0))+
  geom_line(size=4) +
  xlab(expression('Sampling Fraction '*rho))+
  ylab('Proportion without Reticulations')+
  theme(
    legend.title = element_text(size=30), #change legend title font size
    legend.text = element_text(size=20),
    axis.title=element_text(size=30,face="bold"))
ggsave(paste(incom_dir,'/incon_ret0.png',sep=''),width = 10,height = 10,dpi=600)

ggplot(data=gathered_if, mapping = aes(x=frac,y=ave_rets))+
  geom_line(size=4) +
  xlab(expression('Sampling Fraction '*rho))+
  ylab('Average Number of Reticulations')+
  theme(
    legend.title = element_text(size=30), #change legend title font size
    legend.text = element_text(size=20),
    axis.title=element_text(size=30,face="bold"))
ggsave(paste(incom_dir,'/incon_rets.png',sep=''),width = 10,height = 10,dpi=600)

ggplot(data=gathered_if, mapping = aes(x=frac,y=ratio))+
  geom_line(size=4) +
  xlab(expression('Sampling Fraction '*rho))+
  ylab('Average Reticulation Density')+
  theme(
    legend.title = element_text(size=30), #change legend title font size
    legend.text = element_text(size=20),
    axis.title=element_text(size=30,face="bold"))
ggsave(paste(incom_dir,'/incon_ratio.png',sep=''),width = 10,height = 10,dpi=600)


ggplot(data=incom_by_rets,mapping=aes(x=rets,y=tb,color=factor(frac,levels = sort(unique(frac)))))+
  geom_line(size=2)+
  scale_color_viridis_d(name='Sampling Fraction')+
  ylab('Proportion in Class')+
  xlab('Number of Reticulations')+
  theme(
    legend.title = element_text(size=30), #change legend title font size
    legend.text = element_text(size=20),
    axis.title=element_text(size=30,face="bold"))
ggsave(paste(incom_dir,'/incon_ret_tb.png',sep=''),width = 10,height = 10,dpi=600)

ggplot(data=incom_by_rets,mapping=aes(x=rets,y=fus,color=factor(frac,levels = sort(unique(frac)))))+
  geom_line(size=2)+
  scale_color_viridis_d(name='Sampling Fraction')+
  ylab('Proportion in Class')+
  xlab('Number of Reticulations')+
  theme(
    legend.title = element_text(size=30), #change legend title font size
    legend.text = element_text(size=20),
    axis.title=element_text(size=30,face="bold"))
ggsave(paste(incom_dir,'/incon_ret_fus.png',sep=''),width = 10,height = 10,dpi=600)

ggplot(data=incom_by_rets,mapping=aes(x=rets,y=tc,color=factor(frac,levels = sort(unique(frac)))))+
  geom_line(size=2)+
  scale_color_viridis_d(name='Sampling Fraction')+
  ylab('Proportion in Class')+
  xlab('Number of Reticulations')+
  theme(
    legend.title = element_text(size=30), #change legend title font size
    legend.text = element_text(size=20),
    axis.title=element_text(size=30,face="bold"))
ggsave(paste(incom_dir,'/incon_ret_tc.png',sep=''),width = 10,height = 10,dpi=600)

ggplot(data=incom_by_rets,mapping=aes(x=rets,y=nor,color=factor(frac,levels = sort(unique(frac)))))+
  geom_line(size=2)+
  scale_color_viridis_d(name='Sampling Fraction')+
  ylab('Proportion in Class')+
  xlab('Number of Reticulations')+
  theme(
    legend.title = element_text(size=30), #change legend title font size
    legend.text = element_text(size=20),
    axis.title=element_text(size=30,face="bold"))
ggsave(paste(incom_dir,'/incon_ret_nor.png',sep=''),width = 10,height = 10,dpi=600)


if_prop<- (if_dat[,c('frac','gen_prop','degen_prop','neu_prop')])


temp_dat<-data.frame(xx=1/3,yy=1/3,zz=1/3)
ggtern(if_prop, mapping=aes(x=gen_prop,y=degen_prop,z=neu_prop))+
  geom_point(mapping=ggtern::aes(color=frac),size=6)+
  geom_point(temp_dat,mapping=aes(x=xx,y=yy,z=zz),size=8,show.legend = F,shape='diamond')+
  theme_custom(
    col.T = '#DDCC77',
    col.L = '#CC6677',
    col.R = '#AA4499',
  )+
  theme_nomask() +
  theme(tern.axis.arrow = element_line(size = 4),
        tern.panel.grid.minor.T = element_line(color = "gray80"),
        tern.panel.grid.major.T = element_line(color = "gray80"),
        tern.panel.grid.minor.L = element_line(color = "gray80"),
        tern.panel.grid.major.L = element_line(color = "gray80"),
        tern.panel.grid.minor.R = element_line(color = "gray80"),
        tern.panel.grid.major.R = element_line(color = "gray80"),
        tern.panel.background = element_rect(fill = "gray92"),
        axis.title = element_text(size=25,face="bold"),
        tern.axis.text=element_text(size=13.5),
        tern.axis.text.T=element_text(hjust=1,angle=-30),
        tern.axis.text.R=element_text(hjust=1.5,angle=-90),
        tern.axis.text.L=element_text(hjust=0.15,angle=30),
        tern.axis.arrow.text=element_text(size=25,face="bold"),
        legend.title=element_text(size=25),
        
        legend.text = element_text(size=20)
  )+
  theme_arrowlarge()+
  labs(x = 'm-type', xarrow='Proportion',
       z = 'n-type', zarrow='Proportion',
       y = 'y-type', yarrow='Proportion',
       color='Sampling \n Fraction'
  ) + 
  scale_L_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  scale_R_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  scale_T_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  theme_latex(TRUE)+
  theme_rotate(210)

just_n_save('sampling_frac_props')

#############################
### Load Genetic Distance ###
#############################

gathered_gd<- gd_dat %>% 
  select(stre,ave_rets,rets0,ratio,tb,tc,fus,nor) %>%
  gather(key=class,value=prop,-c(stre,ave_rets,ratio,rets0))



gd_by_rets<-gd_frame %>% 
  filter(rets<=10) %>%
  group_by(rets,stre) %>%
  summarize(tb = sum(tree_based)/n(),
            fus= sum(fu_stable)/n(),
            tc = sum(tree_child)/n(),
            nor= sum(normal)/n())


gd_dir<-'gd_plots'
dir.create(gd_dir)






ggplot(data=gathered_gd, mapping = aes(x=stre,y=prop,color=factor(class,levels = c('tb','fus','tc','nor'))))+
  geom_point(size=5) +
  scale_color_discrete(name='Class',labels=c('Tree-Based','FU-Stable','Tree-Child','Normal'))+
  xlab(expression('Strength of Distance Dependence '*delta))+
  ylab('Proportion in class')+
  theme(
    legend.title = element_text(size=30), #change legend title font size
    legend.text = element_text(size=20),
    axis.title=element_text(size=30,face="bold"))
ggsave(paste(gd_dir,'/gd_props.png',sep=''),width = 10,height = 10,dpi=600)



ggplot(data=gathered_gd, mapping = aes(x=stre,y=rets0))+
  geom_line(size=4) +
  xlab(expression('Distance Dependence Strength '*delta))+
  ylab('Proportion without Reticulations')+
  theme(
    legend.title = element_text(size=30), #change legend title font size
    legend.text = element_text(size=20),
    axis.title=element_text(size=30,face="bold"))
ggsave(paste(gd_dir,'/gd_ret0.png',sep=''),width = 10,height = 10,dpi=600)

ggplot(data=gathered_gd, mapping = aes(x=stre,y=ave_rets))+
  geom_line(size=4) +
  xlab(expression('Distance Dependence Strength '*delta))+
  ylab('Average Number of Reticulations')+
  theme(
    legend.title = element_text(size=30), #change legend title font size
    legend.text = element_text(size=20),
    axis.title=element_text(size=30,face="bold"))
ggsave(paste(gd_dir,'/gd_rets.png',sep=''),width = 10,height = 10,dpi=600)

ggplot(data=gathered_gd, mapping = aes(x=stre,y=ratio))+
  geom_line(size=4) +
  xlab(expression('Distance Dependence Strength '*delta))+
  ylab('Average Reticulation Density')+
  theme(
    legend.title = element_text(size=30), #change legend title font size
    legend.text = element_text(size=20),
    axis.title=element_text(size=30,face="bold"))
ggsave(paste(gd_dir,'/gd_ratio.png',sep=''),width = 10,height = 10,dpi=600)


ggplot(data=gd_by_rets,mapping=aes(x=rets,y=tb,color=factor(stre,levels = sort(unique(stre)))))+
  geom_line(size=2)+
  scale_color_viridis_d(name='Distance Dependence')+
  ylab('Proportion in Class')+
  xlab('Number of Reticulations')+
  theme(
    legend.title = element_text(size=30), #change legend title font size
    legend.text = element_text(size=20),
    axis.title=element_text(size=30,face="bold"))
ggsave(paste(gd_dir,'/gd_ret_tb.png',sep=''),width = 10,height = 10,dpi=600)

ggplot(data=gd_by_rets,mapping=aes(x=rets,y=fus,color=factor(stre,levels = sort(unique(stre)))))+
  geom_line(size=2)+
  scale_color_viridis_d(name='Distance Dependence')+
  ylab('Proportion in Class')+
  xlab('Number of Reticulations')+
  theme(
    legend.title = element_text(size=30), #change legend title font size
    legend.text = element_text(size=20),
    axis.title=element_text(size=30,face="bold"))
ggsave(paste(gd_dir,'/gd_ret_fus.png',sep=''),width = 10,height = 10,dpi=600)

ggplot(data=gd_by_rets,mapping=aes(x=rets,y=tc,color=factor(stre,levels = sort(unique(stre)))))+
  geom_line(size=2)+
  scale_color_viridis_d(name='Distance Dependence')+
  ylab('Proportion in Class')+
  xlab('Number of Reticulations')+
  theme(
    legend.title = element_text(size=30), #change legend title font size
    legend.text = element_text(size=20),
    axis.title=element_text(size=30,face="bold"))
ggsave(paste(gd_dir,'/gd_ret_tc.png',sep=''),width = 10,height = 10,dpi=600)

ggplot(data=gd_by_rets,mapping=aes(x=rets,y=nor,color=factor(stre,levels = sort(unique(stre)))))+
  geom_line(size=2)+
  scale_color_viridis_d(name='Distance Dependence')+
  ylab('Proportion in Class')+
  xlab('Number of Reticulations')+
  theme(
    legend.title = element_text(size=30), #change legend title font size
    legend.text = element_text(size=20),
    axis.title=element_text(size=30,face="bold"))
ggsave(paste(gd_dir,'/gd_ret_nor.png',sep=''),width = 10,height = 10,dpi=600)

