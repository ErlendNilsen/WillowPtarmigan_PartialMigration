

### Loading libraries;
library(tidyverse)
library(lme4)
library(glmmTMB)
library(DHARMa)
library(here)
library(rptR)
library(AICcmodavg)


#########################################################################################################
#Loading files - and mutating some variables;

ndft<-read.csv(file= here("data/ndft.csv"), header = TRUE, sep=",") #Partial migration data

ndft <- ndft %>% mutate(dist2=dist/1000,
                        Id=as.factor(RingNR),
                        Weigth2=scale(Weigth, center=TRUE),
                        dist3=if_else(dist2==0, 0.01, dist2))


########################################################################################################
########################################################################################################
### Testing prediction 1 a (P(migration))

GLMM1 <- glmer(BOL ~ Age + (1 | Id ), family = "binomial", data = ndft, nAGQ = 25)
GLMM2 <- glmer(BOL ~ Weigth2 + (1 | Id ), family = "binomial", data = ndft, nAGQ = 25)
GLMM3 <- glmer(BOL ~ Weigth2 + Age + (1 | Id ), family = "binomial", data = ndft, nAGQ = 25)
GLMM4 <- glmer(BOL ~ Weigth2 * Age + (1 | Id ), family = "binomial", data = ndft, nAGQ = 25)
GLMM5 <- glmer(BOL ~ 1 + (1 | Id ), family = "binomial", data = ndft, nAGQ = 25)

cand_1 <- list(GLMM1, GLMM2, GLMM3, GLMM4, GLMM5)
AICcmodavg::aictab(cand_1,modnames=c("Age", "Weigth", "Weigth+Age", "Weigth*Age", "Null"))


################################################################################################
#### Diagnostic;

simulationOutputGLMM4 <- simulateResiduals(fittedModel = GLMM4, n=1000)
testUniformity(simulationOutput = simulationOutputGLMM4)


#################################################################################################
#################################################################################################
#### Prediction 1b - Distance moved rather than P(migration);


LMM1 <- lmer(log(dist3) ~ Age + (1 | Id ), data = ndft, REML =FALSE)
LMM2 <- lmer(log(dist3) ~ Weigth2 + (1 | Id ), data = ndft, REML =FALSE)
LMM3 <- lmer(log(dist3) ~ Weigth2 + Age + (1 | Id ), data = ndft, REML =FALSE)
LMM4 <- lmer(log(dist3) ~ Weigth2 * Age + (1 | Id ), data = ndft, REML =FALSE)
LMM5 <- lmer(log(dist3) ~ 1 + (1 | Id ), data = ndft, REML =FALSE)

cand_1 <- list(LMM1, LMM2, LMM3, LMM4, LMM5)
AICcmodavg::aictab(cand_1, modnames=c("Age", "Weigth", "Weigth+Age", "Weigth*Age", "Null"))

################################################################################################
##### Diagnostics

simulationOutputLMM4 <- simulateResiduals(fittedModel = LMM4, n=1000)
testUniformity(simulationOutput = simulationOutputLMM4)


#################################################################################################
#################################################################################################
## Prediction 2: Repeatability in movement (migration) behaviour

Agreement_rep <- rptGaussian(log(dist3) ~  1 + (1 | Id ), grname = "Id", data = ndft)
summary(Agreement_rep)
Adjusted_rep <- rptGaussian(log(dist3) ~  Age + (1 | Id ), grname = "Id", data = ndft)
summary(Adjusted_rep)

#################################################################################################
#################################################################################################
## Reading data to test predictions about nesting success;

## Reading the data file;

ns <- read.csv(file= here("NestSuccess.csv"), header = TRUE, sep=",") #Nesting Success data

ns <- ns %>% mutate(Id=as.factor(RingNR),
                    Weigth2=scale(Weigth, center=TRUE))


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
AICcmodavg::aictab(cand_2,  modnames=c("Null", "Age", "Weigth", "BOL", "Age+Weigth",
                                       "Age+BOL", "Weigth+BOL", "Weigth+BOL+Age"))

################################################################################################
##### Diagnostics

simulationOutputNS1 <- simulateResiduals(fittedModel = NS1, n=1000)
testUniformity(simulationOutput = simulationOutputNS1)

