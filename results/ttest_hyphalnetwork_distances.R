#install.packages('dplyr')
#install.packages('MASS')
if(!require(dplyr))
  install.packages('dplyr')
if(!require(MASS))
  install.packages("MASS")

library(MXM)


syn=MXM::synapseLogin()
data = MXM::querySynapseTable("syn22279826")


#independent t test
community = MPNST %>% filter(net2_type == 'community')
tumor = community %>% filter(net1 == 'WU505_tumor')
xe = community %>% filter(net1 == 'WU545_xenograft')
model = t.test(tumor$distance, xe$distance, paired = FALSE)
model


tumor2 = community %>% filter(net1 == 'JHU002_tumor')
xe2 = community %>% filter(net1 == 'JHU002_xenograft')
model2 = t.test(tumor2$distance, xe2$distance, paired = FALSE)
model2


tumor3 = community %>% filter(net1 == 'JHU023_tumor')
xe3 = community %>% filter(net1 == 'JHU023_xenograft')
model3 = t.test(tumor3$distance, xe3$distance, paired = FALSE)
model3


tumor4 = community %>% filter(net1 == 'JHU031_tumor')
xe4 = community %>% filter(net1 == 'JHU031_xenograft')
model4 = t.test(tumor4$distance, xe4$distance, paired = FALSE)
model4
