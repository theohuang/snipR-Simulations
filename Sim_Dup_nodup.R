## Simulation to assess deduplication performance for different parameter configurations ##
## No duplicates in data ##
## Last updated: January 11, 2022 ##


a1 <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

n.a1 <- 15

set.seed(12345)
library(snipR)
library(tidyverse)

# load("pedigrees.RData")
load(paste0(getwd(), "/Deduplication/Sim_Params/Sim_Dup_Dat.RData"))
load(paste0(getwd(), "/Deduplication/Sim_Params/Sim_Dup_Keyvars.RData"))


## Grid search parameters
num_neighbors <- seq(5, 20, by = 5) # Number of neighbors (window size)
sn_reps <- seq(1, 9, by = 4) # Number of repetitions for the sort keys
key_length <- c(3, 10, 20, 29) # Length of each sort key (removing 2 because males can't get the 2 female cancers, so those variables won't have any variation)
nreps <- 5 # Number of repetitions for each parameter configuration

thresh.i <- 1:7 # thresholds for intersection score
thresh.g <- seq(0.75, 0.95, 0.1) # percentiles for greedy match score
thresh.all <- thresh.i

cnames <- c("numNeighbors", "snReps", "keyLength", "threshold", "threshType", "iter", "runtime",
            "pairPrecision", "pairRecall", "pairF1",
            "clusterPrecision", "clusterRecall", "clusterF1", "GMD")


res.param <- setNames(data.frame(matrix(NA, nreps * length(num_neighbors) * length(sn_reps) *
                                          length(key_length), length(cnames))),
                      cnames)
res.obj <- vector("list", nreps * length(num_neighbors) * length(sn_reps) *
                    length(key_length))
ct <- 0
ct.list <- 0
for(rr in 1:nreps){
  for (ii in 1:length(num_neighbors)) {
    for(jj in 1:length(sn_reps)){
      for(kk in 1:length(key_length)){
        ct.list <- ct.list + 1
        for(ll in 1:length(thresh.all)){
          ct <- ct + 1
          if(thresh.all[ll] %in% thresh.i){
            res.param[ct, c("numNeighbors", "snReps", "keyLength", "threshold", "threshType", "iter")] <- 
              c(num_neighbors[ii], sn_reps[jj], key_length[kk], thresh.all[ll], "intersection", rr)
            if(ll == 1){
              names(res.obj)[ct.list] <- paste(num_neighbors[ii], sn_reps[jj], key_length[kk], rr, sep = "_")
            }
          } else{
            res.param[ct, c("numNeighbors", "snReps", "keyLength", "threshold", "threshType", "iter")] <- 
              c(num_neighbors[ii], sn_reps[jj], key_length[kk], thresh.all[ll], "greedy", rr)
            if(ll == 1){
              names(res.obj)[ct.list] <- paste(num_neighbors[ii], sn_reps[jj], key_length[kk], rr, sep = "_")
            }
          }
        }
      }
    }
  }
}
res.param <- res.param %>% mutate(across(c(numNeighbors, snReps, keyLength, threshold, iter),
                                         as.numeric))

res.param.nothresh <- res.param %>% filter(threshold == 1)

res.keyVars <- res.param.nothresh %>%
  select(numNeighbors, snReps, keyLength, iter) %>%
  mutate(filled = NA)
res.keyVars[keyVars] <- 0

## rows in res.dup to run
if(a1 == ceiling(nrow(res.keyVars) / n.a1)){
  rows <- (n.a1 * (a1 - 1) + 1):nrow(res.keyVars) 
} else{
  rows <- (n.a1 * (a1 - 1) + 1):(n.a1 * a1)
}

dat.nodup <- dat.tot %>%
  filter(RequestID %in% unique(dat.tot$FamID))

start <- Sys.time()
for(ii in rows){
  start_ii <- Sys.time()
  print(ii)
  nbs <- res.param.nothresh$numNeighbors[ii]
  reps <- res.param.nothresh$snReps[ii]
  kl <- res.param.nothresh$keyLength[ii]
  it <- res.param.nothresh$iter[ii]
  
  ## deduplicating
  res.dup <- snip(pedigrees = dat.nodup, requestID = "RequestID",
                 isProband = "isProband",
                 keyVars = keyVars, keyVars.female = keyVars.female,
                 keyWt = NULL, blockVar = NULL, repSN = reps, windowSN = nbs,
                 keyLength = rep(kl, reps), method = "both",
                 thresh = list(intersection = thresh.i, greedy = thresh.g),
                 priority = list(var = "AgeBC", min = TRUE),
                 seed = 99 * a1 + ii)
  
  ## assessing results
  res.metrics <- summaryMetrics(res.dup, "RequestID", "FamID",
                                thresh.i, thresh.g)
  
  ## storing the results
  res.param$runtime[which(res.param$numNeighbors == nbs &
                            res.param$snReps == reps &
                            res.param$keyLength == kl & 
                            res.param$iter == it)] <-
    as.numeric(difftime(Sys.time(), start_ii, units = "secs"))
  for(j in 1:length(thresh.all)){
    res.param[which(res.param$numNeighbors == nbs &
                      res.param$snReps == reps &
                      res.param$keyLength == kl & 
                      res.param$iter == it &
                      near(res.param$threshold, thresh.all[j])),
              c("pairPrecision", "pairRecall", "pairF1",
                "clusterPrecision", "clusterRecall", "clusterF1", "GMD")] <-
      res.metrics %>% filter(near(threshold, thresh.all[j])) %>% select(-threshold, -threshType)
  }
}
difftime(Sys.time(), start, units = "secs")

res.param <- res.param %>% filter(!is.na(runtime))
save(res.param,
     file = paste0(getwd(), "/Deduplication/Sim_Params/Res_Summary/No_Dup/res_sim_params_nodup_summ_", a1, ".RData"))