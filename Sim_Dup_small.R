## Simulation on smaller families
## Last updated: February 16, 2022

library(PanelPRO)
library(tidyverse)
library(PedUtils)
library(snipR)

genes <- c("BRCA1", "BRCA2", "MLH1", "MSH2", "MSH6", "CDKN2A")
cancers <- c("Breast", "Ovarian", "Colorectal", "Endometrial", "Pancreas", "Melanoma")



## Generating original data, without duplicates ##
nfam <- 5000
set.seed(123)
for(i in 1:nfam){
  if(i %% 100 == 0) print(i)
  
  ## sampling number of relatives
  nSibsPatern <- rpois(2, 0.5)
  nSibsMatern <- rpois(2, 0.5)
  nSibs <- rpois(2, 0.5)
  nGrandchild <- matrix(rpois(2 * (sum(nSibs) + 1), 0.5), sum(nSibs) + 1, 2)
  
  fam <- sim.runSimFam(nSibsPatern, nSibsMatern, nSibs, nGrandchild, PanelPRODatabase,
                       genes, cancers, includeGeno = TRUE)
  fam <- fam %>%
    mutate(FamID = i, RequestID = i, Duplicate = 0, nDuplicates = 0) %>%
    select(-riskmod, -interAge)
  fam$famSize <- nrow(fam)
  
  if(i == 1){
    dat.fam <- fam
  } else{
    dat.fam <- rbind(dat.fam, fam)
  }
}

## Randomly selecting some family members to have genetic testing ##
## 1% of people will have genetic testing results ##
for(i in 1:length(genes)){
  dat.fam[dat.fam[, genes[i]] == 0, genes[i]] <- 2
}
ind.gt <- sample(1:nrow(dat.fam), floor(0.01 * nrow(dat.fam)))
dat.fam[setdiff(1:nrow(dat.fam), ind.gt), genes] <- 0

save(dat.fam, file = "/Users/thuang/Dropbox (Partners HealthCare)/BayesMendel Dropbox folders/Deduplication Paper/Simulation/snipR-Simulations/snipR-Simulations/sim_dat_small.RData")


## Generating duplicates
source('/Users/thuang/Dropbox (Partners HealthCare)/BayesMendel Dropbox folders/Deduplication Paper/Simulation/snipR-Simulations/snipR-Simulations/Generating Families Functions/simulateFamily.R')

keyVar.bin <- setdiff(family %>% select(starts_with("isAff")) %>% names(), "isAffAny")
keyVar.cont <- setdiff(family %>% select(starts_with("Age")) %>% names(), "AgeAny")
keyVars <- c(keyVar.bin, keyVar.cont)
keyVars.female <- c("isAffOC", "isAffENDO")

## Generating a database of duplicates ##
## 200 families are duplicated once, 200 families are duplicated twice, ...,
## 200 families are duplicated 5 times
ndup <- floor(0.1 * nfam)
famid.dup.unique <- sample(unique(dat.fam$FamID), ndup)
dups.groups <- split(famid.dup.unique, ceiling((1:ndup) / (ndup / 5)))
famid.dup <- unlist(sapply(1:5, function(x) rep(dups.groups[[x]], each = x)))
set.seed(123)
for(i in 1:length(famid.dup)){
  fam.dup <- simulateErrors(filter(dat.fam, FamID == famid.dup[i]),
                            1.5, keyVar.bin, keyVar.cont, rep(1, length(c(keyVar.bin, keyVar.cont))))[[1]]
  if(i == 1){
    ## changing RequestID while keeping the same FamID
    fam.dup$RequestID <- max(dat.fam$RequestID) + 1
    dat.dup <- fam.dup
    # dat.dup.db <- rbind(dat.fam, fam.dup)
  } else{
    fam.dup$RequestID <- max(dat.dup$RequestID) + 1
    dat.dup <- rbind(dat.dup, fam.dup)
  }
}

dat.tot <- rbind(dat.fam, dat.dup)

for(i in 1:length(famid.dup)){
  dat.tot$Duplicate[dat.tot$FamID == dat.tot$FamID[dat.tot$RequestID == nfam + i][1]] <- 1
  dat.tot$nDuplicates[dat.tot$FamID == dat.tot$FamID[dat.tot$RequestID == nfam + i][1]] <-
    dat.tot$nDuplicates[dat.tot$FamID == dat.tot$FamID[dat.tot$RequestID == nfam + i][1]] + 1
}

save(dat.fam, dat.tot, file = "/Users/thuang/Dropbox (Partners HealthCare)/BayesMendel Dropbox folders/Deduplication Paper/Simulation/snipR-Simulations/snipR-Simulations/sim_dat_small.RData")


## Running SNIP
res.dup <- snip(pedigrees = dat.tot, requestID = "RequestID",
                isProband = "isProband",
                keyVars = keyVars, keyVars.female = keyVars.female,
                windowSN = 20,
                priority = list(var = "AgeBC", min = TRUE),
                seed = 99)

## assessing results
res.metrics <- summaryMetrics(res.dup, "RequestID", "FamID", 1:7)

## displaying results
library(xtable)
res.metrics[3:7, ] %>%
  select(-threshType) %>%
  xtable() %>%
  print(include.rownames = FALSE)

save(res.dup, res.metrics, file = "/Users/thuang/Dropbox (Partners HealthCare)/BayesMendel Dropbox folders/Deduplication Paper/Simulation/snipR-Simulations/snipR-Simulations/res_sim_small.RData")

