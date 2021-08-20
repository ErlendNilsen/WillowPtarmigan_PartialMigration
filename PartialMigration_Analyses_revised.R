
##########################################################################################################
#########################################################################################################
### Loading some libraries;
library(tidyverse)
library(lme4)
library(glmmTMB)
library(DHARMa)
library(here)
library(rptR)
library(AICcmodavg)
library(ggpubr)
library(ggalluvial)

##########################################################################################################
##########################################################################################################
#### Revisions R1;

# Test only first year after capture for prediction 1 and 3
# Added elevation data for nests
# Changed colours ++ on figures.


#########################################################################################################
#Loading files - and mutating some variables;

ndft<- as_tibble(read.csv(file= here("data/ndft.csv"), header = TRUE, sep=",")) #Partial migration data

ndft <- ndft %>% mutate(dist2=dist/1000,
                        Id=as.factor(RingNR),
                        Weigth2=scale(Weigth, center=TRUE),
                        dist3=if_else(dist2==0, 0.01, dist2),
                        age2=if_else(Age=="ad", "Adult", "Juvenile"))


##### Creating data set  including birds only first year after capture, as suggested by referee on V1

ndft2 <- ndft %>% group_by(RingNR) %>%
  slice(which.min(Year))


#################################################################################################
#################################################################################################
## Reading the nest data file;

ns <- read.csv(file= here("NestSuccess_rev.csv"), header = TRUE, sep=",") #Nesting Success data

ns <- ns %>% mutate(Id=as.factor(RingNR),
                    Weigth2=scale(Weigth, center=TRUE))


ns2 <- ns %>% group_by(RingNR) %>%
  slice(which.min(Year))

## Assessing difference in altitude between migrants and residents;

M1 <- lm(alt_dem10~as.factor(BOL)-1, data=ns)
M1b <- lm(alt_dem10~as.factor(BOL), data=ns)

M2 <- lm(alt_dem10~1, data=ns)

AICc(c(M1, M2))



########################################################################################################
########################################################################################################
### Testing prediction 1 a (P(migration))

GLMM1 <- glmmTMB(BOL ~ Age + (1 | Id ), family = "binomial", data = ndft)
GLMM2 <- glmmTMB(BOL ~ Weigth2 + (1 | Id ), family = "binomial", data = ndft)
GLMM3 <- glmmTMB(BOL ~ Weigth2 + Age + (1 | Id ), family = "binomial", data = ndft)
GLMM4 <- glmmTMB(BOL ~ Weigth2 * Age + (1 | Id ), family = "binomial", data = ndft)
GLMM5 <- glmmTMB(BOL ~ 1 + (1 | Id ), family = "binomial", data = ndft)

cand_1 <- list(GLMM1, GLMM2, GLMM3, GLMM4, GLMM5)
## Table 4
AICcmodavg::aictab(cand_1,modnames=c("Age", "Weigth", "Weigth+Age", "Weigth*Age", "Null"))


################################################################################################
#### Diagnostic;

simulationOutputGLMM4 <- simulateResiduals(fittedModel = GLMM4, n=1000)
testUniformity(simulationOutput = simulationOutputGLMM4)

#################################################################################################
##### Including birds only first year after capture, as suggested by referee on V1
#####

GLM1 <- glm(BOL ~ Age, family = "binomial", data = ndft2)
GLM2 <- glm(BOL ~ Weigth2, family = "binomial", data = ndft2)
GLM3 <- glm(BOL ~ Weigth2+ Age, family = "binomial", data = ndft2)
GLM4 <- glm(BOL ~ Weigth2 * Age, family = "binomial", data = ndft2)
GLM5 <- glm(BOL ~ 1, family = "binomial", data = ndft2)

cand_1 <- list(GLM1, GLM2, GLM3, GLM4, GLM5)
## Table 4
AICcmodavg::aictab(cand_1,modnames=c("Age", "Weigth", "Weigth+Age", "Weigth*Age", "Null"))


#################################################################################################
#################################################################################################
#### Prediction 1b - Distance moved rather than P(migration);

LMM1 <- glmmTMB(log(dist3) ~ Age + (1 | Id ), data = ndft, REML =FALSE)
LMM2 <- glmmTMB(log(dist3) ~ Weigth2 + (1 | Id ), data = ndft, REML =FALSE)
LMM3 <- glmmTMB(log(dist3) ~ Weigth2 + Age + (1 | Id ), data = ndft, REML =FALSE)
LMM4 <- glmmTMB(log(dist3) ~ Weigth2 * Age + (1 | Id ), data = ndft, REML =FALSE)
LMM5 <- glmmTMB(log(dist3) ~ 1 + (1 | Id ), data = ndft, REML =FALSE)

cand_1 <- list(LMM1, LMM2, LMM3, LMM4, LMM5)
## Table 5
AICcmodavg::aictab(cand_1, modnames=c("Age", "Weigth", "Weigth+Age", "Weigth*Age", "Null"))

################################################################################################
##### Diagnostics

simulationOutputLMM4 <- simulateResiduals(fittedModel = LMM4, n=1000)
testUniformity(simulationOutput = simulationOutputLMM4)

################################################################################################
###### Including only first year movement

LM1 <- lm(log(dist3) ~ Age, data = ndft2)
LM2 <- lm(log(dist3) ~ Weigth2, data = ndft2)
LM3 <- lm(log(dist3) ~ Weigth2 + Age, data = ndft2)
LM4 <- lm(log(dist3) ~ Weigth2 * Age, data = ndft2)
LM5 <- lm(log(dist3) ~ 1, data = ndft2)

cand_1 <- list(LM1, LM2, LM3, LM4, LM5)
## Table 5
AICcmodavg::aictab(cand_1, modnames=c("Age", "Weigth", "Weigth+Age", "Weigth*Age", "Null"))


#################################################################################################
#################################################################################################
## Prediction 2: Repeatability in movement (migration) behaviour

temp <- ndft %>% group_by (RingNR) %>%
  summarize(antall=n()) %>%
  filter(antall>1) %>%
  left_join(., ndft)


Agreement_rep <- rptGaussian(log(dist3) ~  1 + (1 | Id ), grname = "Id", data = temp)
summary(Agreement_rep)
Adjusted_rep <- rptGaussian(log(dist3) ~  Age + (1 | Id ), grname = "Id", data = temp)
summary(Adjusted_rep)

##################################################################################
### Using binial response - not reported due to model convergence failure
Agreement_repB <- rpt(BOL ~  1 + (1 | Id ), grname = "Id", data = temp, datatype="Binary")
summary(Agreement_repB)
Adjusted_repB <- rpt(BOL ~  Age + (1 | Id ), grname = "Id", data = temp, datatype="Binary")
summary(Adjusted_repB)





#################################################################################################
### Testing prediction 3a; about clutch size

tt1 <- glmmTMB(N_egg ~ 1 + (1 | Id ), data = ns, family = compois)
tt2 <- glmmTMB(N_egg ~ Age + (1 | Id ), data = ns, family = compois)
tt3 <- glmmTMB(N_egg ~ Weigth2 + (1 | Id ), data = ns, family = compois)
tt4 <- glmmTMB(N_egg ~ BOL + (1 | Id ), data = ns, family = compois)
tt5 <- glmmTMB(N_egg ~ Age + Weigth2 + (1 | Id ), data = ns, family = compois)
tt6 <- glmmTMB(N_egg ~ Age + BOL + (1 | Id ), data = ns, family = compois)
tt7 <- glmmTMB(N_egg ~ Weigth2 + BOL +(1 | Id ), data = ns, family = compois)
tt8 <- glmmTMB(N_egg ~ Weigth2 + BOL + Age + (1 | Id ), data = ns, family = compois)

cand_2 <- list(tt1,tt2,tt3,tt4,tt5, tt6, tt7, tt8)
AICcmodavg::aictab(cand_2, modnames=c("Null", "Age", "Weigth", "BOL", "Age+Weigth",
                                      "Age+BOL", "Weigth+BOL", "Weigth+BOL+Age"))

################################################################################################
##### Diagnostics

simulationOutputBS1 <- simulateResiduals(fittedModel = tt2, n=1000)
testUniformity(simulationOutput = simulationOutputBS1)


################################################################################################
#### Only first year of data;

lm1 <- glm(N_egg ~ 1, data = ns2, family = poisson)
lm2 <- glm(N_egg ~ Age, data = ns2, family = poisson)
lm3 <- glm(N_egg ~ Weigth2, data = ns2, family = poisson)
lm4 <- glm(N_egg ~ BOL, data = ns2, family = poisson)
lm5 <- glm(N_egg ~ Age + Weigth2, data = ns2, family = poisson)
lm6 <- glm(N_egg ~ Age + BOL, data = ns2, family = poisson)
lm7 <- glm(N_egg ~ Weigth2 + BOL, data = ns2, family = poisson)
lm8 <- glm(N_egg ~ Weigth2 + BOL + Age, data = ns2, family = poisson)

cand_2 <- list(lm1,lm2,lm3,lm4,lm5, lm6, lm7, lm8)
AICcmodavg::aictab(cand_2, modnames=c("Null", "Age", "Weigth", "BOL", "Age+Weigth",
                                      "Age+BOL", "Weigth+BOL", "Weigth+BOL+Age"))




#################################################################################################
### Testing prediction 3b;

NS1 <- glmmTMB(Nest_state ~ 1 + (1 | Id ), data = ns, family = "binomial")
NS2 <- glmmTMB(Nest_state ~ Age + (1 | Id ), data = ns, family = "binomial")
NS3 <- glmmTMB(Nest_state ~ Weigth2 + (1 | Id ), data = ns, family = "binomial")
NS4 <- glmmTMB(Nest_state ~ BOL + (1 | Id ), data = ns, family = "binomial")
NS5 <- glmmTMB(Nest_state ~ Age + Weigth2 + (1 | Id ), data = ns, family = "binomial")
NS6 <- glmmTMB(Nest_state ~ Age + BOL + (1 | Id ), data = ns, family = "binomial")
NS7 <- glmmTMB(Nest_state ~ Weigth2 + BOL + (1 | Id ), data = ns, family = "binomial")
NS8 <- glmmTMB(Nest_state ~ Weigth2 + BOL + Age + (1 | Id ), data = ns, family = "binomial")

cand_2 <- list(NS1, NS2,NS3, NS4, NS5, NS6, NS7, NS8)
# Table 7
AICcmodavg::aictab(cand_2,  modnames=c("Null", "Age", "Weigth", "BOL", "Age+Weigth",
                                       "Age+BOL", "Weigth+BOL", "Weigth+BOL+Age"))

################################################################################################
##### Diagnostics

simulationOutputNS1 <- simulateResiduals(fittedModel = NS1, n=1000)
testUniformity(simulationOutput = simulationOutputNS1)

################################################################################################
##### Using only first year of data

N1 <- glm(Nest_state ~ 1, data = ns2, family = "binomial")
N2 <- glm(Nest_state ~ Age, data = ns2, family = "binomial")
N3 <- glm(Nest_state ~ Weigth2, data = ns2, family = "binomial")
N4 <- glm(Nest_state ~ BOL, data = ns2, family = "binomial")
N5 <- glm(Nest_state ~ Age + Weigth2, data = ns2, family = "binomial")
N6 <- glm(Nest_state ~ Age + BOL, data = ns, family = "binomial")
N7 <- glm(Nest_state ~ Weigth2 + BOL, data = ns2, family = "binomial")
N8 <- glm(Nest_state ~ Weigth2 + BOL + Age, data = ns2, family = "binomial")

cand_2 <- list(N1, N2,N3, N4, N5, N6, N7, N8)
# Table 7
AICcmodavg::aictab(cand_2,  modnames=c("Null", "Age", "Weigth", "BOL", "Age+Weigth",
                                       "Age+BOL", "Weigth+BOL", "Weigth+BOL+Age"))




#################################################################################################
#################################################################################################
#################################################################################################
#Plotting FIGURE 3 in manuscript;

ndft <- ndft %>% mutate(
  dist3=dist/1000, BOL3 = if_else(BOL==1, "Migratory", "Resident"))

# A) Hist - movement distances
gp1 <- ggplot(ndft, aes(x=dist3, color = BOL3, fill = BOL3)) +
  geom_histogram(bins = 20, position="identity", alpha = 0.6, show.legend = TRUE) +
  xlab("Distance migrated (km)") + ylab("Frequency")+
  scale_fill_manual(values = c("#fdae6b", "#756bb1"))+
  scale_colour_manual(values = c("black", "black")) +
  theme(legend.position = c(0.7, 0.7),
        legend.title = element_blank(), text=element_text(size=17))

# B) Hist - movement distances : study sites
gp2 <- ggplot(ndft, aes(x=dist3, fill = BOL3, color = BOL3))+
  facet_wrap(vars(Area))+
  geom_histogram(position = "identity", bins = 20,alpha = 0.6, show.legend = F) +
  scale_fill_manual(values = c("#fdae6b", "#756bb1"))+
  scale_colour_manual(values = c("black", "black"))+
  xlab("Distance migrated (km)") + ylab("Frequency") +
  theme(legend.position = c(0.9, 0.9),
        legend.title = element_blank(), text=element_text(size=17))

#c) Weigth vs distance: juveniles
gp3 <- ggplot(filter(ndft, age2 == "Juvenile"), aes(x=Weigth, y=dist))+
  geom_point(position = "identity", aes(color = BOL3), alpha = 0.6, size=5)+
  scale_colour_manual(values = c("#fdae6b", "#756bb1"))+
  ylab("Migr. dist (m)")+ xlab("Individual weight (g)")+
  geom_hline(yintercept=1276, color = "black", size = 1)+
  geom_vline(xintercept=590, color = "black", size = 1, alpha =.5, linetype = "dashed")+
  annotate("text", x=545, y = 0, label = "Residents", size = 5, color = "black") +
  annotate("text", x=545, y = 29000, label = "Migrants", size = 5, color = "black")+
  xlim(515,670)+
  annotate("label", x=630, y = 29000, label = "Juvenile", size = 5, color = "black") +
  theme(legend.position = "none", text=element_text(size=17))

#D) Weigth vs distance: adults
gp4 <-ggplot(filter(ndft, age2 == "Adult"), aes(x=Weigth, y=dist))+
  geom_point(position = "identity", aes(col = BOL3), alpha = 0.6, size=5)+
  scale_colour_manual(values = c("#fdae6b", "#756bb1"))+
  ylab("Migr. dist (m)")+
  xlab("Individual weight (g)")+
  geom_hline(yintercept=1276, color = "black", size = 1)+
  geom_vline(xintercept=600, color = "black", size = 1, alpha =.5, linetype = "dashed")+
  annotate("text", x=545, y = -1500, label = "Residents", size = 5, color = "black") +
  annotate("text", x=545, y = 45000, label = "Migrants", size = 5, color = "black")+
  annotate("label", x=630, y = 45000, label = "Adult", size = 5, color = "black")+
  theme(legend.position = "none", text=element_text(size=17))

ggarrange(gp1, gp2, gp3, gp4, labels = c("A", "B", "C", "D"))
ggsave("Fig3.jpg", plot=last_plot())


##################################################################################################
##################################################################################################
##### figure 4 in manscript;

Fig4 <- ggplot(ndft2, aes(y=BOL, x=Weigth)) + #Erlend sitt plot
  facet_wrap(vars(age2)) +
  ylab("P(migratory)") +
  xlab("Weight") +
  geom_smooth(method="glm", method.args = list(family = "binomial"),
              color="darkorange", alpha=.5) +
  theme_bw() +
  theme(text=element_text(size=20))

Fig4
ggsave("Fig4.jpg", plot=last_plot())

#################################################################################################
#################################################################################################
###### Figure 5 in ms;

temp <- ndft %>% group_by (RingNR) %>%
  summarize(antall=n()) %>%
  filter(antall>1) %>%
  left_join(., ndft) %>%
  mutate(tBOL = ifelse(BOL == 0, 0, 1))

####
temp2_2 <-  temp %>% filter(antall==2) %>%
  group_by(RingNR) %>%
  summarise(x=sum(tBOL)) %>%
  mutate(sBOL = ifelse(x == 1,0,1)) %>%
  mutate(antall = ifelse(x==0,3,3))

temp3_2 <-  temp %>% filter(antall==3) %>%
  group_by(RingNR) %>%
  summarise(x=sum(tBOL)) %>%
  mutate(sBOL = ifelse(x == 2,0,1)) %>%
  mutate(antall = ifelse(x==0,4,4))

temp4_2 <-  temp %>% filter(antall>3) %>%
  group_by(RingNR) %>%
  summarise(x=sum(tBOL)) %>%
  mutate(sBOL = ifelse(x == 5,1,1)) %>%
  mutate(antall = ifelse(x==1,5,5))

temp2 <- bind_rows(temp2_2, temp3_2, temp4_2) %>%
  mutate(sBOL= if_else(sBOL==1, "Repeating", "Non-repeating"))


temp_2 <-Reduce(function(...) merge(..., all=TRUE), list(temp2_2, temp3_2, temp4_2)) %>%
  arrange(desc(antall))

temp_2$sBOL[temp_2$sBOL == 1] <- "Repeating n=10"
temp_2$sBOL[temp_2$sBOL == 0] <- "Non-repeating n=4"

Fig5 <- ggplot(temp2, aes(axis1 = desc(antall), axis2 = sBOL)) +
  geom_alluvium(aes(fill=sBOL), width = 1/15, show.legend = FALSE, color = "gray45", alpha = 0.6) +
  geom_stratum(width = 1/15, fill = "gray60") +
  scale_x_discrete(limits = c("Consecutive seasons", "Repeatability"), expand = c(.1, .05)) +
  scale_fill_manual(values = c("#fdae6b", "#756bb1"))+
  annotate("label", x=1, y = 7.5, label = "4", size = 8, color = "black") +
  annotate("label", x=1, y = 12, label = "5", size = 8, color = "black") +
  annotate("label", x=1, y = 2.5, label = "3", size = 8, color = "black") +
  annotate("label", x=1.91, y = 5, label = "Repeating", size = 8, color = "black") +
  annotate("label", x=1.89, y = 11.5, label = "Non-repeating", size = 8, color = "black") +
  theme_bw() +
  theme(axis.ticks.y=element_blank(),
        axis.text = element_text(size = 20, face = "bold", margin = margin(20)),
        axis.text.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        axis.text.x = element_text(margin = margin(t = -20)),
        axis.ticks = element_blank(),
        axis.title = element_blank())

Fig5
ggsave("Fig5.jpg", plot=last_plot())


