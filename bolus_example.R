rm(list = ls())
library(glmmTMB)
library(dplyr)
library(xtable)
library(mvtnorm)
library(Matrix)
library(parallel)
library(foreach)
library(doParallel)
library(GGally)
library(lme4)
library(glm2)
library(faraway)
library(R2BayesX)
library(tidyr)


########################-------------Bolus Data-------------###############################
library(cold)
dat = bolus
setwd("C:/ANU/PhD/3rd")
dat$id = as.factor(dat$id)

glmmtmb_mod = glmmTMB(y  ~ 1 + time + (1 + time |id), data = dat,family = poisson)
glmmtmb_ranef = glmmtmb_mod %>% ranef %>% as.data.frame
glmmtmb_ranef
glmmtmb_mod
lme_mod = glmer(y  ~ 1 + time + (1 + time |id), data = dat,family = poisson)
lme_ranef = lme_mod %>% ranef %>% as.data.frame
lme_ranef
lme_mod

varcovar = VarCorr(glmmtmb_mod)$cond[[1]]
varcovar_lme = VarCorr(lme_mod)$id
X_1 = glmmtmb_mod$frame %>% filter(id=="1") %>% dplyr::select(time)
X_1 = rep(1,nrow(X_1)) %>% cbind(X_1) %>% as.matrix

#########---------- normal scale mixture quantiles ------------##########
S = matrix(nrow=10000,ncol=ncol(X_1))
for (i in 1:10000) {
  b_1 = rmvnorm(n=1,sigma = varcovar_lme) %>% t %>% as.vector
  eta_1 = X_1%*%lme_mod@beta + X_1%*%b_1
  W_1 = diag(nrow=nrow(X_1))
  diag(W_1) = poisson()$var(poisson()$linkinv(eta_1))
  S[i,] = rmvnorm(1,sigma = solve(crossprod(X_1,W_1)%*%X_1))
}
int_interval = quantile(S[,1],probs = c(0.025,0.975))
treat_interval = quantile(S[,2],probs = c(0.025,0.975))
lme_sds_int = lme_ranef %>% filter(term == '(Intercept)') %>% .$condsd
lme_sds_treat = lme_ranef %>% filter(term == 'time') %>% .$condsd
glmmtmb_sds_int = glmmtmb_ranef %>% filter(term == '(Intercept)') %>% .$condsd
glmmtmb_sds_treat = glmmtmb_ranef %>% filter(term == 'time') %>% .$condsd

patient = 1:65
all_groups = tibble(patient,lme_sds_int,glmmtmb_sds_int,lme_sds_treat,glmmtmb_sds_treat) %>%
  pivot_longer(c(lme_sds_int,glmmtmb_sds_int,lme_sds_treat,glmmtmb_sds_treat),
               names_to = "package",values_to = "estimate") %>%
  mutate(type = rep(c('Intercept','Intercept','Slope','Slope'),65)) %>%
  mutate(Software = rep(c('lme4','glmmTMB','lme4','glmmTMB'),65)) %>%
  mutate(intlength = 2*1.96*estimate) %>%
  mutate(unconditional = rep(c(int_interval[2]-int_interval[1],int_interval[2]-int_interval[1],
                               treat_interval[2]-treat_interval[1],treat_interval[2]-treat_interval[1]),65)) %>%
  mutate(paired = rep(1:130,each=2))

pdf('intlength.pdf', width=7, height=3.5)
all_groups %>%
  ggplot(aes(x = patient, y = intlength)) +
  facet_wrap(~ type, ncol = 2, scales = "free") +
  geom_point(aes(shape=Software, color = Software), size = 2) +
  labs(x = "Patient Number" , y = "Interval Length", col = "Package", shape = "Package") +
  geom_hline(aes(yintercept=unconditional)) +
  geom_line(aes(group = paired),linetype = 3) +
  scale_color_manual(name = "Package", labels = c("glmmTMB", "lme4"), values = c(2,4)) +
  scale_shape_manual(name = "Package", labels = c("glmmTMB", "lme4"), values = c(19,17)) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 65, by = 5))
dev.off()

plot(glmmtmb_ranef$condval,lme_ranef$condval,
     ylab = "lme4 Predicted Random Effect", xlab = "glmmTMB Predicted Random Effect")
abline(0,1,col="red")








