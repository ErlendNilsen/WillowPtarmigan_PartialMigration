
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

#### Revision R2;
# Present analyses based on winter-summer transition first year after capture in main text,


#########################################################################################################
#Loading files - and mutating some variables;

ndft<- as_tibble(read.csv(file= here("data/ndft.csv"), header = TRUE, sep=",")) #Partial migration data

ndft <- ndft %>% mutate(dist2=dist/1000,
                        Id=as.factor(RingNR),
                        Weigth2=scale(Weigth, center=TRUE),
                        dist3=if_else(dist2==0, 0.01, dist2),
                        age2=if_else(Age=="ad", "Adult", "Juvenile"))


##### Creating data set  only winter to summer-transition, and only first year after capture (for birds with more than one year of data)

ndft2 <- ndft %>% group_by(RingNR) %>%
         filter(Season2=="WS") %>%
  slice(which.min(Year))


#################################################################################################
#################################################################################################
## Reading the nest data file;

ns <- as_tibble(read.csv(file= here("data/NestSuccess_rev2.csv"), header = TRUE, sep=",")) #Nesting Success data

ns <- ns %>% mutate(Id=as.factor(RingNR),
                    Weigth2=scale(Weigth, center=TRUE))

## Only firt year after capture

ns2 <- ns %>% group_by(RingNR) %>%
  slice(which.min(Year))

## Assessing difference in altitude between migrants and residents;

M1 <- lm(alt_dem10~as.factor(BOL)-1, data=ns)
M1b <- lm(alt_dem10~as.factor(BOL), data=ns)

M2 <- lm(alt_dem10~1, data=ns)

AICc(c(M1, M2))

########################################################################################################
########################################################################################################
### summary statistics:

## N_birds in main analysis (only first year of data - only winter to summer transition)

N_birds_main <- ndft2 %>% group_by(Year) %>%
                summarize(Number=n())

N_birds_all <- ndft %>% group_by(Year) %>%
  summarize(Number=n())

N_nests_all <- ns %>% group_by(Year) %>%
          summarize(Number=n())

N_nests_main <- ns2 %>% group_by(Year) %>%
  summarize(Number=n())

### Table 2 - updated with only first year of data;

N_migrants_main <- ndft2 %>% group_by(Year, BOL) %>%
                  summarize(Number=n())

N_migrants_main2 <- ndft2 %>% group_by(BOL) %>%
  summarize(Number=n())

N_migrants_main3 <- ndft2 %>% group_by(Year) %>%
  summarize(Total=n())

N_migrants_main4 <- ndft2 %>% filter(BOL==0) %>% group_by(Year) %>%
  summarize(Residents=n())

N_migrants_main5 <- ndft2 %>% filter(BOL==1) %>% group_by(Year) %>%
  summarize(Migrants=n())

Table2 <- N_migrants_main4 %>% left_join(., N_migrants_main5) %>% left_join(., N_migrants_main3) %>%
          mutate(Perc_migratory=Migrants/Total)


### Table 3:

## Weigth:

Weigth_table <- ndft2 %>% group_by(Age) %>%
                mutate(vekt2=vekt*100) %>%
                summarize(Min=min(vekt2),
                          Mean=mean(vekt2),
                          Median=median(vekt2),
                          Max=max(vekt2),
                          N=n())

Movement_table <- ndft2 %>% group_by(Age) %>%
  summarize(Min=min(dist2),
            Mean=mean(dist2),
            Median=median(dist2),
            Max=max(dist2),
            N=n())


########################################################################################################
########################################################################################################
### Testing prediction 1 a (P(migration))

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


###############################################################################################
#### All data - reported in Appendix

GLMM1 <- glmmTMB(BOL ~ Age + (1 | Id ), family = "binomial", data = ndft)
GLMM2 <- glmmTMB(BOL ~ Weigth2 + (1 | Id ), family = "binomial", data = ndft)
GLMM3 <- glmmTMB(BOL ~ Weigth2 + Age + (1 | Id ), family = "binomial", data = ndft)
GLMM4 <- glmmTMB(BOL ~ Weigth2 * Age + (1 | Id ), family = "binomial", data = ndft)
GLMM5 <- glmmTMB(BOL ~ 1 + (1 | Id ), family = "binomial", data = ndft)

cand_1 <- list(GLMM1, GLMM2, GLMM3, GLMM4, GLMM5)
## Table A5
AICcmodavg::aictab(cand_1,modnames=c("Age", "Weigth", "Weigth+Age", "Weigth*Age", "Null"))


################################################################################################
#### Diagnostic;

simulationOutputGLMM4 <- simulateResiduals(fittedModel = GLMM4, n=1000)
testUniformity(simulationOutput = simulationOutputGLMM4)


#################################################################################################
#################################################################################################
#### Prediction 1b - Distance moved rather than P(migration);
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

###############################################################################################
#### Including all data

LMM1 <- glmmTMB(log(dist3) ~ Age + (1 | Id ), data = ndft, REML =FALSE)
LMM2 <- glmmTMB(log(dist3) ~ Weigth2 + (1 | Id ), data = ndft, REML =FALSE)
LMM3 <- glmmTMB(log(dist3) ~ Weigth2 + Age + (1 | Id ), data = ndft, REML =FALSE)
LMM4 <- glmmTMB(log(dist3) ~ Weigth2 * Age + (1 | Id ), data = ndft, REML =FALSE)
LMM5 <- glmmTMB(log(dist3) ~ 1 + (1 | Id ), data = ndft, REML =FALSE)

cand_1 <- list(LMM1, LMM2, LMM3, LMM4, LMM5)
## Table A6
AICcmodavg::aictab(cand_1, modnames=c("Age", "Weigth", "Weigth+Age", "Weigth*Age", "Null"))

################################################################################################
##### Diagnostics

simulationOutputLMM4 <- simulateResiduals(fittedModel = LMM4, n=1000)
testUniformity(simulationOutput = simulationOutputLMM4)


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
#################################################################################################
### Testing prediction 3a; about clutch size
################################################################################################
#### Only first year of data;

################################################################################################
### mpcmp::glm.cmp, due to strong underdispersion (which is common in data on clutch size)

#ns3 <- ns2 %>% filter(!is.na(N_egg))


COM_glm1 <- mpcmp::glm.cmp(N_egg ~ 1, data = ns2)
COM_glm2 <- mpcmp::glm.cmp(N_egg ~ Age, data = ns2)
COM_glm3 <- mpcmp::glm.cmp(N_egg ~ Weigth2, data = ns2)
COM_glm4 <- mpcmp::glm.cmp(N_egg ~ BOL, data = ns2)
COM_glm5 <- mpcmp::glm.cmp(N_egg ~ Age + Weigth2, data = ns2)
COM_glm6 <- mpcmp::glm.cmp(N_egg ~ Age + BOL, data = ns2)
COM_glm7 <- mpcmp::glm.cmp(N_egg ~ Weigth2 + BOL, data = ns2)
COM_glm8 <- mpcmp::glm.cmp(N_egg ~ Weigth2 + BOL + Age, data = ns2)

##Amnually copute AICc and associated stats, as AICcmodavg::aictab is not not implemented for
## mpcmp::glm.cmp

cand_2 <- tibble(Mod_names=c("Null", "Age", "Weigth", "BOL", "Age+Weigth",
                             "Age+BOL", "Weigth+BOL", "Weigth+BOL+Age"),
            AIC=c((summary(COM_glm1))$aic,(summary(COM_glm2))$aic,(summary(COM_glm3))$aic,(summary(COM_glm4))$aic,
            (summary(COM_glm5))$aic, (summary(COM_glm6))$aic, (summary(COM_glm7))$aic, (summary(COM_glm8))$aic),
            n=dim(ns3)[1], k=c(2, 3, 3, 3, 4, 4, 4, 5),
            Correction= ((2*k)*(k+1) / (n-k-1)),
            AICc=AIC+Correction,
            deltaAICc=AICc-min(AICc),
            Weigth1 = exp(-0.5 *deltaAICc),
            w_i=Weigth1/sum(Weigth1),
            Cum_wi=cumsum(w_i)) %>%
            select(-AIC, -Correction, -Weigth1, -n)

## Table 6
View(cand_2)

################################################################################################
#### All data

tt1 <- glmmTMB(N_egg ~ 1 + (1 | Id ), data = ns, family = compois)
tt2 <- glmmTMB(N_egg ~ Age + (1 | Id ), data = ns, family = compois)
tt3 <- glmmTMB(N_egg ~ Weigth2 + (1 | Id ), data = ns, family = compois)
tt4 <- glmmTMB(N_egg ~ BOL + (1 | Id ), data = ns, family = compois)
tt5 <- glmmTMB(N_egg ~ Age + Weigth2 + (1 | Id ), data = ns, family = compois)
tt6 <- glmmTMB(N_egg ~ Age + BOL + (1 | Id ), data = ns, family = compois)
tt7 <- glmmTMB(N_egg ~ Weigth2 + BOL +(1 | Id ), data = ns, family = compois)
tt8 <- glmmTMB(N_egg ~ Weigth2 + BOL + Age + (1 | Id ), data = ns, family = compois)

cand_2 <- list(tt1,tt2,tt3,tt4,tt5, tt6, tt7, tt8)

## Table A7
AICcmodavg::aictab(cand_2, modnames=c("Null", "Age", "Weigth", "BOL", "Age+Weigth",
                                      "Age+BOL", "Weigth+BOL", "Weigth+BOL+Age"))

################################################################################################
##### Diagnostics

simulationOutputBS1 <- simulateResiduals(fittedModel = tt2, n=1000)
testUniformity(simulationOutput = simulationOutputBS1)


#################################################################################################
### Testing prediction 3b;

################################################################################################
##### Using only first year of data

N1 <- glm(Nest_state ~ 1, data = ns2, family = "binomial")
N2 <- glm(Nest_state ~ Age, data = ns2, family = "binomial")
N3 <- glm(Nest_state ~ Weigth2, data = ns2, family = "binomial")
N4 <- glm(Nest_state ~ BOL, data = ns2, family = "binomial")
N5 <- glm(Nest_state ~ Age + Weigth2, data = ns2, family = "binomial")
N6 <- glm(Nest_state ~ Age + BOL, data = ns2, family = "binomial")
N7 <- glm(Nest_state ~ Weigth2 + BOL, data = ns2, family = "binomial")
N8 <- glm(Nest_state ~ Weigth2 + BOL + Age, data = ns2, family = "binomial")

cand_2 <- list(N1, N2,N3, N4, N5, N6, N7, N8)
# Table 7
AICcmodavg::aictab(cand_2,  modnames=c("Null", "Age", "Weigth", "BOL", "Age+Weigth",
                                       "Age+BOL", "Weigth+BOL", "Weigth+BOL+Age"))

###############################################################################################
### All data

NS1 <- glmmTMB(Nest_state ~ 1 + (1 | Id ), data = ns, family = "binomial")
NS2 <- glmmTMB(Nest_state ~ Age + (1 | Id ), data = ns, family = "binomial")
NS3 <- glmmTMB(Nest_state ~ Weigth2 + (1 | Id ), data = ns, family = "binomial")
NS4 <- glmmTMB(Nest_state ~ BOL + (1 | Id ), data = ns, family = "binomial")
NS5 <- glmmTMB(Nest_state ~ Age + Weigth2 + (1 | Id ), data = ns, family = "binomial")
NS6 <- glmmTMB(Nest_state ~ Age + BOL + (1 | Id ), data = ns, family = "binomial")
NS7 <- glmmTMB(Nest_state ~ Weigth2 + BOL + (1 | Id ), data = ns, family = "binomial")
NS8 <- glmmTMB(Nest_state ~ Weigth2 + BOL + Age + (1 | Id ), data = ns, family = "binomial")

cand_2 <- list(NS1, NS2,NS3, NS4, NS5, NS6, NS7, NS8)
# Table A8
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
### Revision: Only include winter - to - summer transition & only first year of data



ndft2 <- ndft2 %>% mutate(
  dist3=dist/1000, BOL3 = if_else(BOL==1, "Migratory", "Resident"))

# A) Hist - movement distances
gp1 <- ggplot(ndft2, aes(x=dist3, color = BOL3, fill = BOL3)) +
  geom_histogram(bins = 20, position="identity", alpha = 0.6, show.legend = TRUE) +
  xlab("Distance migrated (km)") + ylab("Frequency")+
  scale_fill_manual(values = c("#fdae6b", "#756bb1"))+
  scale_colour_manual(values = c("black", "black")) +
  theme(legend.position = c(0.7, 0.7),
        legend.title = element_blank(), text=element_text(size=17))

# B) Hist - movement distances : study sites
gp2 <- ggplot(ndft2, aes(x=dist3, fill = BOL3, color = BOL3))+
  facet_wrap(vars(Area))+
  geom_histogram(position = "identity", bins = 20,alpha = 0.6, show.legend = F) +
  scale_fill_manual(values = c("#fdae6b", "#756bb1"))+
  scale_colour_manual(values = c("black", "black"))+
  xlab("Distance migrated (km)") + ylab("Frequency") +
  theme(legend.position = c(0.9, 0.9),
        legend.title = element_blank(), text=element_text(size=17))

#c) Weigth vs distance: juveniles
gp3 <- ggplot(filter(ndft2, age2 == "Juvenile"), aes(x=Weigth, y=dist))+
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
gp4 <-ggplot(filter(ndft2, age2 == "Adult"), aes(x=Weigth, y=dist))+
  geom_point(position = "identity", aes(col = BOL3, ), alpha = 0.6, size=5)+
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
##### figure 4 in manuscript;

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


