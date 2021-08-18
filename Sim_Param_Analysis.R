## Sim Param Analysis ##
## Last Updated: August 18, 2021 ##

library(snipR)
library(tidyverse)
library(ggplot2)
library(reshape2)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## loading results ##
path <- paste0(getwd(), "/Sim_Param_Results")
files <- list.files(path, full.names = TRUE)

for(i in 1:length(files[!file.info(files)$isdir])){
  load(paste0(path, "/res_sim_params_summ_", i, ".RData"))
  if(i == 1){
    res.param.tot <- res.param
  } else{
    res.param.tot <- rbind(res.param.tot, res.param)
  }
}
res.param <- res.param.tot

## collecting results ##
res.param.avg <- res.param %>%
  group_by(numNeighbors, snReps, keyLength, threshold, threshType) %>%
  summarize(runtime = mean(runtime),
            pairPrecision = mean(pairPrecision),
            pairRecall = mean(pairRecall),
            pairF1 = mean(pairF1),
            clusterPrecision = mean(clusterPrecision),
            clusterRecall = mean(clusterRecall),
            clusterF1 = mean(clusterF1),
            GMD = mean(GMD)) %>%
  ungroup()


res.param.avg %>%
  slice_max(pairF1, n = 10, with_ties = FALSE) %>%
  print(width = Inf)
res.param.avg %>%
  slice_max(clusterF1, n = 10, with_ties = FALSE) %>%
  print(width = Inf)
res.param.avg %>%
  slice_min(GMD, n = 10, with_ties = FALSE) %>%
  print(width = Inf)

## summary table ##
library(xtable)
print(xtable(res.param.avg %>%
               slice_min(GMD, n = 10, with_ties = FALSE) %>%
               select(-threshType) %>%
               as.data.frame(),
             digits = c(0, 0, 0, 0, 0, 2, 3, 3, 3, 3, 3, 3, 1)),
      include.rownames = FALSE)



res.long <- reshape2::melt(res.param.avg,
                           id.vars = c("numNeighbors", "snReps", "keyLength", "threshold"),
                           measure.vars = c("pairPrecision", "pairRecall", "pairF1",
                                            "clusterPrecision", "clusterRecall", "clusterF1",
                                            "GMD"),
                           variable.name = "metric")

## pairwise F1, cluster F1, and GMD in one plot
ggplot(res.long %>% filter(threshold %in% 4:7,
                           metric %in% c("pairF1", "clusterF1", "GMD")) %>%
         mutate(value = case_when(metric == "GMD" ~ value /
                                    max(res.long %>%
                                          filter(threshold %in% 4:7, metric == "GMD") %>%
                                          select(value), na.rm = TRUE),
                                  TRUE ~ value)),
       aes(numNeighbors, value)) +
  geom_line(aes(color = factor(snReps), linetype = factor(metric))) +
  geom_point(aes(color = factor(snReps), shape = factor(metric)), size = 1.7) +
  labs(y = "", x = "Sliding window size", linetype = "Metric", shape = "Metric") +
  scale_color_discrete(name = "Number of sort \nkey iterations") +
  scale_linetype_discrete(labels = c("Pairwise F1", "Cluster F1", "GMD")) +
  scale_shape_discrete(labels = c("Pairwise F1", "Cluster F1", "GMD")) +
  facet_grid(keyLength ~ threshold,
             labeller = labeller(keyLength = function(x) paste0("Key length: ", x),
                                 threshold = function(x) paste0("Threshold: ", x)))

ggplot(res.long %>% filter(threshold %in% 4:7, metric == "pairF1"), aes(numNeighbors, value)) +
  geom_line(aes(color = factor(snReps))) +
  geom_point(aes(color = factor(snReps)), size = 1) +
  labs(y = "Pairwise F1", x = "Sliding window size") +
  scale_color_discrete(name = "Number of sort \nkey iterations") +
  scale_linetype_discrete(name = "Metric") +
  facet_grid(keyLength ~ threshold,
             labeller = labeller(keyLength = function(x) paste0("Key length: ", x),
                                 threshold = function(x) paste0("Threshold: ", x)))
ggplot(res.long %>% filter(threshold %in% 4:7, metric == "clusterF1"), aes(numNeighbors, value)) +
  geom_line(aes(color = factor(snReps))) +
  geom_point(aes(color = factor(snReps)), size = 1) +
  labs(y = "Cluster F1", x = "Sliding window size") +
  scale_color_discrete(name = "Number of sort \nkey iterations") +
  scale_linetype_discrete(name = "Metric") +
  facet_grid(keyLength ~ threshold,
             labeller = labeller(keyLength = function(x) paste0("Key length: ", x),
                                 threshold = function(x) paste0("Threshold: ", x)))
ggplot(res.long %>% filter(threshold %in% 4:7, metric == "GMD"), aes(numNeighbors, value)) +
  geom_line(aes(color = factor(snReps))) +
  geom_point(aes(color = factor(snReps)), size = 1) +
  labs(y = "GMD", x = "Sliding window size") +
  scale_color_discrete(name = "Number of sort \nkey iterations") +
  scale_linetype_discrete(name = "Metric") +
  facet_grid(keyLength ~ threshold,
             labeller = labeller(keyLength = function(x) paste0("Key length: ", x),
                                 threshold = function(x) paste0("Threshold: ", x)))


## Plotting runtimes
ggplot(res.param.avg %>% filter(threshold %in% 1:7), aes(numNeighbors, runtime)) +
  geom_line(aes(color = factor(snReps), linetype = factor(keyLength))) +
  geom_point(aes(color = factor(snReps)), size = 1) +
  scale_color_discrete(name = "Number of sort \nkey iterations") +
  scale_linetype_discrete(name = "Key Length") +
  labs(x = "Sliding window size", y = "Runtime (secs)")

