## Simulation on effect of deduplication, when using SNIP
## Allele frequency of 0.1
## Last updated: January 5, 2022


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

load(paste0(getwd(), "/Res_Sim_Effect_af01.RData"))

library(snipR)
library(tidyverse)
library(pROC)

nfam <- 10000
pDup.seq <- seq(0, 0.5, 0.05)
pFH.seq <- seq(0.5, 1, 0.1)

dat.fam <- dat.fam %>% mutate(BRCA = as.numeric(BRCA1 + BRCA2 > 0))

n.cv <- 100
res.cv <- vector("list", length(pDup.seq))
names(res.cv) <- paste0("pDup_", pDup.seq)
for(i in 1:length(pDup.seq)){
  res.cv[[i]] <- vector("list", length(pFH.seq))
  names(res.cv[[i]]) <- paste0("pFH_", pFH.seq)
}

## Performance metrics ##
set.seed(123)
cnames <- as.vector(outer(c("OE", "AUC", "rBS"),
                          c("cv", "val.ext", "train.ext"),
                          paste, sep = "_"))
start <- Sys.time()
for(k in 1:length(pDup.seq)){
  for(l in 1:length(pFH.seq)){
    print(c(k, l))
    dat.summ.added <- rbind(dat.summ, dat.dup.summ[[k]][[l]])
    
    dat.added <- rbind(dat.fam, dat.dup[[k]][[l]])
    
    dat.summ.added.dedup <- snip(dat.added, "RequestID", "isProband",
                                 c("isAffBC", "isAffOC", "CurAge", "AgeBC", "AgeOC", "isDead",
                                   "Gender", "BRCA1", "BRCA2"),
                                 priority = list(var = "ID", min = TRUE),
                                 seed = 99)
    id.dedup <- dat.summ.added.dedup$dupsReps$thresh_5
    dat.summ.added <- dat.summ.added %>% filter(RequestID %in% id.dedup)
    
    ## random cross-validation prediction error ##
    res.cv[[k]][[l]] <- setNames(data.frame(matrix(NA, n.cv, length(cnames))),
                                 cnames)
    
    for(i in 1:100){
      ind.cv <- sample(nrow(dat.summ.added), floor(nrow(dat.summ.added)/2))
      ind.boot.ext <- sample(nrow(dat.summ.ext), replace = TRUE)
      ind.boot.added <- sample(nrow(dat.summ.added), replace = TRUE)
      mod.cv <- glm(BRCA ~ proBC + proOC + relBC + relOC + proAge,
                    data = dat.summ.added[ind.cv, ], family = "binomial")
      mod.train.ext <- glm(BRCA ~ proBC + proOC + relBC + relOC + proAge,
                           data = dat.summ.ext, family = "binomial")
      pred.cv <- predict(mod.cv, dat.summ.added[-ind.cv, ], type = "response")
      pred.val.ext <- predict(mod.cv, dat.summ.ext[ind.boot.ext, ], type = "response")
      pred.train.ext <- predict(mod.train.ext, dat.summ.added[ind.boot.added, ], type = "response")
      
      #### Performance metrics ####
      
      ## training and validating on data with duplicates
      res.cv[[k]][[l]]$OE_cv[i] <- sum(dat.summ.added$BRCA[-ind.cv]) / sum(pred.cv)
      res.cv[[k]][[l]]$AUC_cv[i] <- auc(BRCA ~ pred.cv, data = dat.summ.added[-ind.cv, ])
      res.cv[[k]][[l]]$rBS_cv[i] <- sqrt(mean((dat.summ.added$BRCA[-ind.cv] - pred.cv)^2))
      
      ## training on data with duplicates, validating on external data set without duplicates ##
      res.cv[[k]][[l]]$OE_val.ext[i] <- sum(dat.summ.ext$BRCA[ind.boot.ext]) / sum(pred.val.ext)
      res.cv[[k]][[l]]$AUC_val.ext[i] <- auc(BRCA ~ pred.val.ext, data = dat.summ.ext[ind.boot.ext, ])
      res.cv[[k]][[l]]$rBS_val.ext[i] <- sqrt(mean((dat.summ.ext$BRCA[ind.boot.ext] - pred.val.ext)^2))
      
      ## training on external data without duplicates, validating on data with duplicates ##
      res.cv[[k]][[l]]$OE_train.ext[i] <- sum(dat.summ.added$BRCA[ind.boot.added]) / sum(pred.train.ext)
      res.cv[[k]][[l]]$AUC_train.ext[i] <- auc(BRCA ~ pred.train.ext, data = dat.summ.added[ind.boot.added, ])
      res.cv[[k]][[l]]$rBS_train.ext[i] <- sqrt(mean((dat.summ.added$BRCA[ind.boot.added] - pred.train.ext)^2))
    }
  }
}
print(difftime(Sys.time(), start, units = "secs"))

save(res.cv, file = paste0(getwd(), "/Res_Sim_Effect_af01_dedup.RData"))

getCI <- function(res, ndigits, comb = FALSE){
  if(comb == TRUE){
    paste0(format(round(mean(res), ndigits), nsmall = ndigits),
           ", (", format(round(sort(res)[3], ndigits), nsmall = ndigits),
           ", ", format(round(sort(res)[97], ndigits), nsmall = ndigits),
           ")")
  } else{
    c(format(round(mean(res), ndigits), nsmall = ndigits),
      format(round(sort(res)[3], ndigits), nsmall = ndigits),
      format(round(sort(res)[97], ndigits), nsmall = ndigits))
  }
}

ndigits <- 3
res.oe <- res.auc <- res.rbs <-
  lapply(setNames(vector("list", 3), c("cv", "val.ext", "train.ext")),
         function(x){
           z <- data.frame(matrix(NA, length(pDup.seq), length(pFH.seq)*3))
           rownames(z) <- paste0("pDup_", pDup.seq)
           colnames(z) <- apply(expand.grid(c("mean", "lo", "up"),
                                            paste0("pFH_", pFH.seq))[, c(2, 1)], 1,
                                function(x) paste(x[1], x[2], sep = "_"))
           return(z)
         })

for(i in 1:length(pDup.seq)){
  for(j in 1:length(pFH.seq)){
    ## training and validating on data with duplicates ##
    res.oe$cv[i, (3*(j-1)+1):(3*j)] <- getCI(res.cv[[i]][[j]]$OE_cv, ndigits)
    res.auc$cv[i, (3*(j-1)+1):(3*j)] <- getCI(res.cv[[i]][[j]]$AUC_cv, ndigits)
    res.rbs$cv[i, (3*(j-1)+1):(3*j)] <- getCI(res.cv[[i]][[j]]$rBS_cv, ndigits)
    
    ## training on data with duplicates, validating on external data set without duplicates ##
    res.oe$val.ext[i, (3*(j-1)+1):(3*j)] <- getCI(res.cv[[i]][[j]]$OE_val.ext, ndigits)
    res.auc$val.ext[i, (3*(j-1)+1):(3*j)] <- getCI(res.cv[[i]][[j]]$AUC_val.ext, ndigits)
    res.rbs$val.ext[i, (3*(j-1)+1):(3*j)] <- getCI(res.cv[[i]][[j]]$rBS_val.ext, ndigits)
    
    ## training on external data without duplicates, validating on data with duplicates ##
    res.oe$train.ext[i, (3*(j-1)+1):(3*j)] <- getCI(res.cv[[i]][[j]]$OE_train.ext, ndigits)
    res.auc$train.ext[i, (3*(j-1)+1):(3*j)] <- getCI(res.cv[[i]][[j]]$AUC_train.ext, ndigits)
    res.rbs$train.ext[i, (3*(j-1)+1):(3*j)] <- getCI(res.cv[[i]][[j]]$rBS_train.ext, ndigits)
  }
}

## Plots ##
oe.mean <- lapply(lapply(res.oe, function(x) x %>% select(ends_with("mean"))),
                  function(x){
                    pivot_longer(x, cols = starts_with("pFH")) %>%
                      transmute(value = as.numeric(value),
                                pDup = rep(pDup.seq, each = length(pFH.seq)),
                                pFH = rep(pFH.seq, length(pDup.seq)))
                  })
auc.mean <- lapply(lapply(res.auc, function(x) x %>% select(ends_with("mean"))),
                   function(x){
                     pivot_longer(x, cols = starts_with("pFH")) %>%
                       transmute(value = as.numeric(value),
                                 pDup = rep(pDup.seq, each = length(pFH.seq)),
                                 pFH = rep(pFH.seq, length(pDup.seq)))
                   })
rbs.mean <- lapply(lapply(res.rbs, function(x) x %>% select(ends_with("mean"))),
                   function(x){
                     pivot_longer(x, cols = starts_with("pFH")) %>%
                       transmute(value = as.numeric(value),
                                 pDup = rep(pDup.seq, each = length(pFH.seq)),
                                 pFH = rep(pFH.seq, length(pDup.seq)))
                   })

metrics.mean <- rbind(oe.mean$cv %>% mutate(metric = "O/E", type = "CV with duplicates"),
                      oe.mean$val.ext %>% mutate(metric = "O/E", type = "Train with duplicates, validate without duplicates"),
                      oe.mean$train.ext %>% mutate(metric = "O/E", type = "Train without duplicates, validate with duplicates"),
                      auc.mean$cv %>% mutate(metric = "AUC", type = "CV with duplicates"),
                      auc.mean$val.ext %>% mutate(metric = "AUC", type = "Train with duplicates, validate without duplicates"),
                      auc.mean$train.ext %>% mutate(metric = "AUC", type = "Train without duplicates, validate with duplicates"),
                      rbs.mean$cv %>% mutate(metric = "rBS", type = "CV with duplicates"),
                      rbs.mean$val.ext %>% mutate(metric = "rBS", type = "Train with duplicates, validate without duplicates"),
                      rbs.mean$train.ext %>% mutate(metric = "rBS", type = "Train without duplicates, validate with duplicates"))

ggplot(metrics.mean, aes(pDup, value)) +
  geom_line(aes(color = factor(pFH))) +
  geom_point(aes(color = factor(pFH))) +
  labs(x = "Proportion of families in original data with duplicates", y = "") +
  scale_color_discrete(name = "Among duplicates,\nproportion with \nextensive family \nhistory") +
  facet_grid(factor(metric, levels = c("O/E", "AUC", "rBS")) ~ type,
             scales = "free",
             switch = "y") +
  theme(axis.title.y = element_blank(), # remove the default y-axis title, "wt"
        strip.background.y = element_rect(fill = 'transparent'), # replace the strip backgrounds with transparent
        strip.placement = 'outside', # put the facet strips on the outside
        strip.text.y = element_text(angle = 180)) # rotate the y-axis text (optional)


metrics.mean_dedup <- metrics.mean

load(paste0(getwd(), "/Res_Sim_Effect_af01.RData"))

oe.mean <- lapply(lapply(res.oe, function(x) x %>% select(ends_with("mean"))),
                  function(x){
                    pivot_longer(x, cols = starts_with("pFH")) %>%
                      transmute(value = as.numeric(value),
                                pDup = rep(pDup.seq, each = length(pFH.seq)),
                                pFH = rep(pFH.seq, length(pDup.seq)))
                  })
auc.mean <- lapply(lapply(res.auc, function(x) x %>% select(ends_with("mean"))),
                   function(x){
                     pivot_longer(x, cols = starts_with("pFH")) %>%
                       transmute(value = as.numeric(value),
                                 pDup = rep(pDup.seq, each = length(pFH.seq)),
                                 pFH = rep(pFH.seq, length(pDup.seq)))
                   })
rbs.mean <- lapply(lapply(res.rbs, function(x) x %>% select(ends_with("mean"))),
                   function(x){
                     pivot_longer(x, cols = starts_with("pFH")) %>%
                       transmute(value = as.numeric(value),
                                 pDup = rep(pDup.seq, each = length(pFH.seq)),
                                 pFH = rep(pFH.seq, length(pDup.seq)))
                   })

metrics.mean <- rbind(oe.mean$cv %>% mutate(metric = "O/E", type = "CV with duplicates"),
                      oe.mean$val.ext %>% mutate(metric = "O/E", type = "Train with duplicates, validate without duplicates"),
                      oe.mean$train.ext %>% mutate(metric = "O/E", type = "Train without duplicates, validate with duplicates"),
                      auc.mean$cv %>% mutate(metric = "AUC", type = "CV with duplicates"),
                      auc.mean$val.ext %>% mutate(metric = "AUC", type = "Train with duplicates, validate without duplicates"),
                      auc.mean$train.ext %>% mutate(metric = "AUC", type = "Train without duplicates, validate with duplicates"),
                      rbs.mean$cv %>% mutate(metric = "rBS", type = "CV with duplicates"),
                      rbs.mean$val.ext %>% mutate(metric = "rBS", type = "Train with duplicates, validate without duplicates"),
                      rbs.mean$train.ext %>% mutate(metric = "rBS", type = "Train without duplicates, validate with duplicates"))


metrics.comb <- rbind(metrics.mean_dedup %>% mutate(Deduplicated = "Yes"),
                      metrics.mean %>% mutate(Deduplicated = "No"))

ggplot(metrics.comb %>% filter(pFH == 0.7), aes(pDup, value)) +
  geom_line(aes(color = factor(Deduplicated))) +
  geom_point(aes(color = factor(Deduplicated))) +
  geom_hline(data = metrics.comb %>% filter(pFH == 0.7, pDup == 0, Deduplicated == "No"),
             aes(yintercept = value), linetype = "dashed") +
  labs(x = "Proportion of families in original data with duplicates", y = "") +
  scale_color_discrete(name = "Deduplicated") +
  facet_grid(factor(metric, levels = c("O/E", "AUC", "rBS")) ~ type,
             scales = "free",
             switch = "y") +
  theme(axis.title.y = element_blank(), # remove the default y-axis title, "wt"
        strip.background.y = element_rect(fill = 'transparent'), # replace the strip backgrounds with transparent
        strip.placement = 'outside', # put the facet strips on the outside
        strip.text.y = element_text(angle = 180)) # rotate the y-axis text (optional)

