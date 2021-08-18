## Simulation on effect of deduplication
## Allele frequency of 0.01
## Last updated: August 18, 2021

library(tidyverse)
library(data.table)
library(pROC)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
sim.files <- list.files(paste0(getwd(), "/Generating Families Functions"))
for(i in 1:length(sim.files)){
  source(paste0(getwd(), "/Generating Families Functions/", sim.files[i]))
}

library(BayesMendel)
library(abind)

## loading BayesMendel colorectal and endometrial cancer penetrances and
## hazard for death from other causes

mutations <- c("BRCA1", "BRCA2")
cancers <- c("BC", "OC")
genos <- c("B00", "B10", "B01")


## getting baseline hazards -- assuming they are the MMRpro hazards
penF <- list(BC = penet.brca.net$fFX[, genos],
             OC = penet.brca.net$fFY[, genos])
penM <- list(BC = penet.brca.net$fMX[, genos],
             OC = penet.brca.net$fMY[, genos])

## using allele frequencies of 0.01
af <- setNames(rep(0.01, 2), mutations)
## baseline penetrance (MMRpro)
CP <- genCancerPen(mutations, cancers, penF, penM, maxK = length(mutations), age.last = 95)


## Generating original data, without duplicates ##
nfam <- 10000
start <- Sys.time()
for(i in 1:nfam){
  if(i %% 100 == 0) print(i)
  nSibsPatern <- sample(0:3, 2, replace = TRUE)
  nSibsMatern <- sample(0:3, 2, replace = TRUE)
  nSibs <- sample(0:3, 2, replace = TRUE)
  nGrandchild <- matrix(sample(0:3, 2 * (sum(nSibs) + 1), replace = TRUE), sum(nSibs) + 1, 2)
  fam <- sim.simFam(nSibsPatern, nSibsMatern, nSibs,
                    nGrandchild, af, CP, age.max = 94, includeGeno = TRUE)
  fam <- fam %>% mutate(FamID = i, RequestID = i, Duplicate = 0, nDuplicates = 0)
  if(i == 1){
    dat.fam <- fam
  } else{
    dat.fam <- rbind(dat.fam, fam)
  }
}

save(dat.fam, file = "SimFamEffect_Orig_001.RData")


## summary
dat.summ <- dat.fam %>%
  group_by(RequestID) %>%
  mutate(proBC = sum((isAffBC == 1) * (isProband == 1), na.rm = TRUE),
         proOC = sum((isAffOC == 1) * (isProband == 1), na.rm = TRUE),
         relBC = mean((isAffBC == 1) * (ID != 1), na.rm = TRUE),
         relOC = mean((isAffOC == 1) * (ID != 1), na.rm = TRUE),
         proAge = max(pmax(cbind(AgeBC * (isProband == 1), AgeOC * (isProband == 1), na.rm = TRUE))),
         AgeBC = sum((isAffBC == 1) * (isProband == 1) * AgeBC, na.rm = TRUE),
         AgeOC = sum((isAffOC == 1) * (isProband == 1) * AgeOC, na.rm = TRUE),
         relBCOC = mean((isAffBC >= 1 | isAffOC >= 1) * (ID != 1), na.rm = TRUE),
         BRCA = as.numeric(BRCA1 + BRCA2 > 0),
         FamSize = n()) %>%
  filter(isProband == 1) %>%
  select(FamID, RequestID, Duplicate, nDuplicates, Gender, proBC, proOC, relBC, relOC, proAge,
         AgeBC, AgeOC, relBCOC, BRCA1, BRCA2, BRCA, FamSize)

## cross-validation error on original data set, without duplicates ##
n.cv <- 100
res.cv.orig <- setNames(data.frame(matrix(NA, n.cv, 3)),
                        c("OE", "AUC", "rBS"))
set.seed(12345)
for(i in 1:100){
  ind.i <- sample(1:nrow(dat.summ), nrow(dat.summ)/2)
  mod.i <- glm(BRCA ~ proBC + proOC + relBC + relOC + proAge,
               data = dat.summ[ind.i, ], family = "binomial")
  pred.i <- predict(mod.i, dat.summ[-ind.i, ], type = "response")
  res.cv.orig$OE[i] <- sum(dat.summ$BRCA[-ind.i]) / sum(pred.i)
  res.cv.orig$AUC[i] <- auc(BRCA ~ pred.i, data = dat.summ[-ind.i, ])
  res.cv.orig$rBS[i] <- sqrt(mean((dat.summ$BRCA[-ind.i] - pred.i)^2))
}
summary(res.cv.orig)


### Adding duplicates ###

keyVar.bin <- c("Gender", "isAffBC", "isAffOC", "BRCA1", "BRCA2")
keyVar.cont <- c("AgeBC","AgeOC")

## Generating a database of duplicates ##
famid.db <- sample(unique(dat.fam$FamID), nfam, replace = TRUE)
set.seed(123)
start <- Sys.time()
for(i in 1:nfam){
  if(i %% 100 == 0) print(i)
  fam.dup <- simulateErrors(filter(dat.fam, FamID == famid.db[i]),
                            1.5, keyVar.bin, keyVar.cont, rep(1, length(c(keyVar.bin, keyVar.cont))))[[1]]
  if(i == 1){
    ## changing RequestID while keeping the same FamID
    fam.dup$RequestID <- max(dat.fam$RequestID) + 1
    dat.dup.db <- fam.dup
  } else{
    fam.dup$RequestID <- max(dat.dup.db$RequestID) + 1
    dat.dup.db <- rbind(dat.dup.db, fam.dup)
  }
}
print(difftime(Sys.time(), start, units = "secs"))

## summary of duplicate families database
dat.dup.db.summ <- dat.dup.db %>%
  mutate(across(starts_with(c("isAff", "Age", "BRCA", "Gender")), as.numeric)) %>%
  group_by(RequestID) %>%
  mutate(proBC = sum((isAffBC == 1) * (isProband == 1), na.rm = TRUE),
         proOC = sum((isAffOC == 1) * (isProband == 1), na.rm = TRUE),
         relBC = mean((isAffBC == 1) * (ID != 1), na.rm = TRUE),
         relOC = mean((isAffOC == 1) * (ID != 1), na.rm = TRUE),
         proAge = max(pmax(cbind(AgeBC * (isProband == 1), AgeOC * (isProband == 1), na.rm = TRUE))),
         AgeBC = sum((isAffBC == 1) * (isProband == 1) * AgeBC, na.rm = TRUE),
         AgeOC = sum((isAffOC == 1) * (isProband == 1) * AgeOC, na.rm = TRUE),
         relBCOC = mean((isAffBC >= 1 | isAffOC >= 1) * (ID != 1), na.rm = TRUE),
         BRCA = as.numeric(BRCA1 + BRCA2 > 0),
         FamSize = n()) %>%
  filter(isProband == 1) %>%
  select(FamID, RequestID, Duplicate, nDuplicates, Gender, proBC, proOC, relBC, relOC, proAge,
         AgeBC, AgeOC, relBCOC, BRCA1, BRCA2, BRCA, FamSize) %>%
  ungroup()

## RequestIDs based on family history
rid.fh <- dat.dup.db.summ %>%
  filter(relBCOC > quantile(dat.dup.db.summ$relBCOC, 0.9)) %>%
  select(RequestID) %>%
  unlist() %>%
  as.numeric()
rid.nfh <- setdiff(unique(dat.dup.db.summ$RequestID), rid.fh)
summary(dat.dup.db.summ %>% filter(RequestID %in% rid.fh) %>% select(relBCOC))
summary(dat.dup.db.summ %>% filter(RequestID %in% rid.nfh) %>% select(relBCOC))


#### Simulating effect of deduplication on random cross-validated performance
#### for logistic regression prediction model. Sensitivity analysis for the
#### proportion of duplicates in the data and proportion of duplicates who have
#### high family history
pDup.seq <- seq(0, 0.5, 0.05)
pFH.seq <- seq(0.5, 1, 0.1)
dat.dup <- dat.dup.summ <- res.cv <- vector("list", length(pDup.seq))
names(dat.dup) <- names(dat.dup.summ) <- names(res.cv) <- paste0("pDup_", pDup.seq)
for(i in 1:length(pDup.seq)){
  dat.dup[[i]] <- dat.dup.summ[[i]] <- res.cv[[i]] <- vector("list", length(pFH.seq))
  names(dat.dup[[i]]) <- names(dat.dup.summ[[i]]) <- names(res.cv[[i]]) <- paste0("pFH_", pFH.seq)
}

start <- Sys.time()
set.seed(12345)
for(k in 1:length(pDup.seq)){
  for(l in 1:length(pFH.seq)){
    print(c(k, l))
    
    nDup <- floor(pDup.seq[k] * nfam)
    nFH <- floor(pFH.seq[l] * nDup)
    
    if(nDup > 0){
      rid.dup <- c(sample(rid.fh, nFH, replace = TRUE),
                   sample(rid.nfh, nDup - nFH, replace = TRUE))
      rid.dup.unq <- unique(rid.dup)
      rid.dup.nunq <- rid.dup[duplicated(rid.dup)]
      ind.dup <- c()
      for(i in 1:length(rid.dup.unq)){
        ind.dup <- c(ind.dup, which(dat.dup.db$RequestID == rid.dup.unq[i]))
      }
      dat.dup.kl <- dat.dup.db[ind.dup, ]
      for(i in 1:length(rid.dup.nunq)){
        if(i == 1){
          dup.nunq <- dat.dup.db %>%
            filter(RequestID == rid.dup.nunq[1]) %>%
            mutate(RequestID = max(dat.dup.db$RequestID) + i)
        } else{
          dup.nunq <- rbind(dup.nunq, dat.dup.db %>%
                              filter(RequestID == rid.dup.nunq[i]) %>%
                              mutate(RequestID = max(dat.dup.db$RequestID) + i))
        }
      }
      
      dat.dup[[k]][[l]] <- rbind(dat.dup.kl, dup.nunq)
      
      dat.dup[[k]][[l]] <- dat.dup[[k]][[l]] %>%
        mutate(across(starts_with(c("isAff", "Age", "BRCA", "Gender")), as.numeric)) %>%
        mutate(BRCA = as.numeric(BRCA1 + BRCA2 > 0))
      
      ## Summary of duplicates ##
      dat.dup.summ[[k]][[l]] <- dat.dup[[k]][[l]] %>%
        group_by(RequestID) %>%
        mutate(proBC = sum((isAffBC == 1) * (isProband == 1)),
               proOC = sum((isAffOC == 1) * (isProband == 1)),
               relBC = mean((isAffBC == 1) * (ID != 1)),
               relOC = mean((isAffOC == 1) * (ID != 1)),
               proAge = max(pmax(cbind(AgeBC * (isProband == 1), AgeOC * (isProband == 1)))),
               ageBC = sum((isAffBC == 1) * (isProband == 1) * AgeBC),
               ageOC = sum((isAffOC == 1) * (isProband == 1) * AgeOC),
               relBCOC = mean((isAffBC == 1 | isAffOC == 1) * (ID != 1))) %>%
        filter(isProband == 1) %>%
        select(FamID, Duplicate, RequestID, Gender, proBC, proOC, relBC, relOC, proAge,
               ageBC, ageOC, relBCOC, BRCA1, BRCA2, BRCA)
    }
  }
}
print(difftime(Sys.time(), start, units = "secs"))

### Generating external validation data, separate from original data set ###
nfam.ext <- 10000
set.seed(99)
start <- Sys.time()
for(i in 1:nfam.ext){
  if(i %% 100 == 0) print(i)
  nSibsPatern <- sample(0:3, 2, replace = TRUE)
  nSibsMatern <- sample(0:3, 2, replace = TRUE)
  nSibs <- sample(0:3, 2, replace = TRUE)
  nGrandchild <- matrix(sample(0:3, 2 * (sum(nSibs) + 1), replace = TRUE), sum(nSibs) + 1, 2)
  fam <- sim.simFam(nSibsPatern, nSibsMatern, nSibs,
                    nGrandchild, af, CP, age.max = 94, includeGeno = TRUE)
  fam <- fam %>% mutate(FamID = i, RequestID = i, Duplicate = 0, nDuplicates = 0)
  if(i == 1){
    dat.fam.ext <- fam
  } else{
    dat.fam.ext <- rbind(dat.fam.ext, fam)
  }
}
print(difftime(Sys.time(), start, units = "secs"))

dat.summ.ext <- dat.fam.ext %>%
  group_by(RequestID) %>%
  mutate(proBC = sum((isAffBC == 1) * (isProband == 1), na.rm = TRUE),
         proOC = sum((isAffOC == 1) * (isProband == 1), na.rm = TRUE),
         relBC = mean((isAffBC == 1) * (ID != 1), na.rm = TRUE),
         relOC = mean((isAffOC == 1) * (ID != 1), na.rm = TRUE),
         proAge = max(pmax(cbind(AgeBC * (isProband == 1), AgeOC * (isProband == 1), na.rm = TRUE))),
         AgeBC = sum((isAffBC == 1) * (isProband == 1) * AgeBC, na.rm = TRUE),
         AgeOC = sum((isAffOC == 1) * (isProband == 1) * AgeOC, na.rm = TRUE),
         relBCOC = mean((isAffBC >= 1 | isAffOC >= 1) * (ID != 1), na.rm = TRUE),
         BRCA = as.numeric(BRCA1 + BRCA2 > 0),
         FamSize = n()) %>%
  filter(isProband == 1) %>%
  select(FamID, RequestID, Duplicate, nDuplicates, Gender, proBC, proOC, relBC, relOC, proAge,
         AgeBC, AgeOC, relBCOC, BRCA1, BRCA2, BRCA, FamSize)


## Performance metrics ##
set.seed(123)
cnames <- as.vector(outer(c("OE", "AUC", "rBS"),
                          c("cv", "val.ext", "train.ext",
                            "cv.ss", "val.ext.ss", "train.ext.ss"),
                          paste, sep = "_"))
start <- Sys.time()
for(k in 1:length(pDup.seq)){
  for(l in 1:length(pFH.seq)){
    print(c(k, l))
    dat.summ.added <- rbind(dat.summ, dat.dup.summ[[k]][[l]])
    
    ## removing families to keep sample size constant
    dat.summ.added.ss <- dat.summ.added %>%
      filter(RequestID %in% sample(unique(dat.summ.added$RequestID), nfam))
  
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
      
      ## keeping sample size constant regardless of pDup
      ind.cv.ss <- sample(1:nrow(dat.summ.added.ss), floor(nrow(dat.summ.added.ss)/2))
      ind.boot.added.ss <- sample(nrow(dat.summ.added.ss), replace = TRUE)
      mod.cv.ss <- glm(BRCA ~ proBC + proOC + relBC + relOC + proAge,
                       data = dat.summ.added.ss[ind.cv.ss, ], family = "binomial")
      pred.cv.ss <- predict(mod.cv.ss, dat.summ.added.ss[-ind.cv.ss, ], type = "response")
      pred.val.ext.ss <- predict(mod.cv.ss, dat.summ.ext[ind.boot.ext, ], type = "response")
      pred.train.ext.ss <- predict(mod.train.ext, dat.summ.added.ss[ind.boot.added.ss, ], type = "response")
      
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
      
      ### Constant sample size ###
      ## training and validating on data with duplicates
      res.cv[[k]][[l]]$OE_cv.ss[i] <- sum(dat.summ.added.ss$BRCA[-ind.cv.ss]) / sum(pred.cv.ss)
      res.cv[[k]][[l]]$AUC_cv.ss[i] <- auc(BRCA ~ pred.cv.ss, data = dat.summ.added.ss[-ind.cv.ss, ])
      res.cv[[k]][[l]]$rBS_cv.ss[i] <- sqrt(mean((dat.summ.added.ss$BRCA[-ind.cv.ss] - pred.cv.ss)^2))
      
      ## training on data with duplicates, validating on external data set without duplicates ##
      res.cv[[k]][[l]]$OE_val.ext.ss[i] <- sum(dat.summ.ext$BRCA[ind.boot.ext]) / sum(pred.val.ext.ss)
      res.cv[[k]][[l]]$AUC_val.ext.ss[i] <- auc(BRCA ~ pred.val.ext.ss, data = dat.summ.ext[ind.boot.ext, ])
      res.cv[[k]][[l]]$rBS_val.ext.ss[i] <- sqrt(mean((dat.summ.ext$BRCA[ind.boot.ext] - pred.val.ext.ss)^2))
      
      ## training on external data without duplicates, validating on data with duplicates ##
      res.cv[[k]][[l]]$OE_train.ext.ss[i] <- sum(dat.summ.added.ss$BRCA[ind.boot.added.ss]) / sum(pred.train.ext.ss)
      res.cv[[k]][[l]]$AUC_train.ext.ss[i] <- auc(BRCA ~ pred.train.ext.ss, data = dat.summ.added.ss[ind.boot.added.ss, ])
      res.cv[[k]][[l]]$rBS_train.ext.ss[i] <- sqrt(mean((dat.summ.added.ss$BRCA[ind.boot.added.ss] - pred.train.ext.ss)^2))
    }
  }
}
print(difftime(Sys.time(), start, units = "secs"))

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

## Mean performance metrics ##
ndigits <- 3
res.oe <- res.auc <- res.rbs <- res.oe.ss <- res.auc.ss <- res.rbs.ss <-
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
    
    ### Constant sample size ###
    ## training and validating on data with duplicates ##
    res.oe.ss$cv[i, (3*(j-1)+1):(3*j)] <- getCI(res.cv[[i]][[j]]$OE_cv.ss, ndigits)
    res.auc.ss$cv[i, (3*(j-1)+1):(3*j)] <- getCI(res.cv[[i]][[j]]$AUC_cv.ss, ndigits)
    res.rbs.ss$cv[i, (3*(j-1)+1):(3*j)] <- getCI(res.cv[[i]][[j]]$rBS_cv.ss, ndigits)
    
    ## training on data with duplicates, validating on external data set without duplicates ##
    res.oe.ss$val.ext[i, (3*(j-1)+1):(3*j)] <- getCI(res.cv[[i]][[j]]$OE_val.ext.ss, ndigits)
    res.auc.ss$val.ext[i, (3*(j-1)+1):(3*j)] <- getCI(res.cv[[i]][[j]]$AUC_val.ext.ss, ndigits)
    res.rbs.ss$val.ext[i, (3*(j-1)+1):(3*j)] <- getCI(res.cv[[i]][[j]]$rBS_val.ext.ss, ndigits)
    
    ## training on external data without duplicates, validating on data with duplicates ##
    res.oe.ss$train.ext[i, (3*(j-1)+1):(3*j)] <- getCI(res.cv[[i]][[j]]$OE_train.ext.ss, ndigits)
    res.auc.ss$train.ext[i, (3*(j-1)+1):(3*j)] <- getCI(res.cv[[i]][[j]]$AUC_train.ext.ss, ndigits)
    res.rbs.ss$train.ext[i, (3*(j-1)+1):(3*j)] <- getCI(res.cv[[i]][[j]]$rBS_train.ext.ss, ndigits)
  }
}

res.auc.wid <- lapply(res.auc, function(x) x %>%
                        mutate(across(everything(), as.numeric)) %>%
                        transmute(pFH_0.5_width = pFH_0.5_up - pFH_0.5_lo,
                                  pFH_0.6_width = pFH_0.6_up - pFH_0.6_lo,
                                  pFH_0.7_width = pFH_0.7_up - pFH_0.7_lo,
                                  pFH_0.8_width = pFH_0.8_up - pFH_0.8_lo,
                                  pFH_0.9_width = pFH_0.9_up - pFH_0.9_lo,
                                  pFH_1_width = pFH_1_up - pFH_1_lo))
res.auc.wid <- lapply(res.auc.wid, function(x){
  rownames(x) <- paste("pDup", pDup.seq, sep = "_")
  return(x)
})

res.oe.wid <- lapply(res.oe, function(x) x %>%
                       mutate(across(everything(), as.numeric)) %>%
                       transmute(pFH_0.5_width = pFH_0.5_up - pFH_0.5_lo,
                                 pFH_0.6_width = pFH_0.6_up - pFH_0.6_lo,
                                 pFH_0.7_width = pFH_0.7_up - pFH_0.7_lo,
                                 pFH_0.8_width = pFH_0.8_up - pFH_0.8_lo,
                                 pFH_0.9_width = pFH_0.9_up - pFH_0.9_lo,
                                 pFH_1_width = pFH_1_up - pFH_1_lo))
res.oe.wid <- lapply(res.oe.wid, function(x){
  rownames(x) <- paste("pDup", pDup.seq, sep = "_")
  return(x)
})

res.oe.ss.wid <- lapply(res.oe.ss, function(x) x %>%
                          mutate(across(everything(), as.numeric)) %>%
                          transmute(pFH_0.5_width = pFH_0.5_up - pFH_0.5_lo,
                                    pFH_0.6_width = pFH_0.6_up - pFH_0.6_lo,
                                    pFH_0.7_width = pFH_0.7_up - pFH_0.7_lo,
                                    pFH_0.8_width = pFH_0.8_up - pFH_0.8_lo,
                                    pFH_0.9_width = pFH_0.9_up - pFH_0.9_lo,
                                    pFH_1_width = pFH_1_up - pFH_1_lo))
res.oe.ss.wid <- lapply(res.oe.ss.wid, function(x){
  rownames(x) <- paste("pDup", pDup.seq, sep = "_")
  return(x)
})


res.auc.ss.wid <- lapply(res.auc.ss, function(x) x %>%
                           mutate(across(everything(), as.numeric)) %>%
                           transmute(pFH_0.5_width = pFH_0.5_up - pFH_0.5_lo,
                                     pFH_0.6_width = pFH_0.6_up - pFH_0.6_lo,
                                     pFH_0.7_width = pFH_0.7_up - pFH_0.7_lo,
                                     pFH_0.8_width = pFH_0.8_up - pFH_0.8_lo,
                                     pFH_0.9_width = pFH_0.9_up - pFH_0.9_lo,
                                     pFH_1_width = pFH_1_up - pFH_1_lo))
res.auc.ss.wid <- lapply(res.auc.ss.wid, function(x){
  rownames(x) <- paste("pDup", pDup.seq, sep = "_")
  return(x)
})

lapply(res.oe, function(x) x %>% select(ends_with("mean")))
lapply(res.auc, function(x) x %>% select(ends_with("mean")))
lapply(res.rbs, function(x) x %>% select(ends_with("mean")))




## Saving data ##
save(dat.fam, dat.summ, dat.dup.db, dat.dup.db.summ,
     dat.dup, dat.dup.summ, dat.fam.ext, dat.summ.ext,
     res.cv, res.cv.orig,
     res.oe, res.auc, res.rbs, res.oe.ss, res.auc.ss, res.rbs.ss,
     file = "Res_Sim_Effect_af001.RData")


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

