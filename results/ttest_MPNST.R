#install.packages('dplyr')
#install.packages('MASS')
if(!require(dplyr))
  install.packages('dplyr')
if(!require(MASS))
  install.packages("MASS")
if(!require(stringr))
  install.packages("stringr")
if(!require(ggplot2))
  install.packages('ggplot2')
if(!require(purrr))
  install.packages("purrr")
library(MXM)
library(stringr)

syn=MXM::synapseLogin()
df = MXM::querySynapseTable("syn22279826")


#independent t test
community = df %>% filter(net2_type == 'community')

##added this to help plotting
community$category = rep('Tumor',nrow(community))
community$category[grep('xenograft',community$net1)]<-'Xenograft'

ttest_comm = function(x){
  data = community %>% filter(net2 == x)
  tumor = data %>% filter(str_detect(net1,'tumor'))
  xenograft = data %>% filter(str_detect(net1,'xenograft'))
  model = t.test(tumor$distance, xenograft$distance, paired = FALSE)
  return (model)
}


community_list = unique(community$net2)
#models = rep(0,45)
#for (val in community_list){
#  print(paste("Community number", val))
#  model = ttest_comm(val)
#  print(model)
#}

pvals<-community_list%>%
  purrr::map(~ttest_comm(.)$p.value)%>%as.numeric()%>%unlist()

#now adding larger data frame to plot
res.df<-data.frame(Community=unlist(community_list),PValue=pvals)
community <- community%>%rename(Community=net2)%>%left_join(res.df)

ggplot(subset(community,PValue<0.3),aes(x=Community,y=distance,fill=category))+geom_boxplot()

