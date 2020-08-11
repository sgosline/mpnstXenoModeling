#install.packages('dplyr')
#install.packages('MASS')
if(!require(dplyr))
  install.packages('dplyr')
if(!require(MASS))
  install.packages("MASS")
if(!require(stringr))
  install.packages("stringr")

library(MXM)
library(stringr)

syn=MXM::synapseLogin()
df = MXM::querySynapseTable("syn22279826")


#independent t test
community = df %>% filter(net2_type == 'community')


ttest_comm = function(x){
  data = community %>% filter(net2 == x)
  tumor = data %>% filter(str_detect(net1,'tumor'))
  xenograft = data %>% filter(str_detect(net1,'xenograft'))
  model = t.test(tumor$distance, xenograft$distance, paired = FALSE)
  return (model)
}


community_list = unique(community$net2)
#models = rep(0,45)
for (val in community_list){
  print(paste("Community number", val))
  model = ttest_comm(val)
  print(model)
}

