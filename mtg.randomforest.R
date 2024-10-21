library(mikropml)
library(tidyverse)
library(phyloseq)
library(microbiome)
library(microViz)

###Making rf format using MTG feature table###
#Use MTG featuretable as input file
mtg <- read.csv("./yonsei.mtg.csv")
#select classification group
gc.mtg <- gc.mtg %>% filter(Group %in% c("HC", "GC"))

#select taxonomic level and aggregate certain taxonomic level
#The BLAST column is the BLAST result at the Species level of the MTG ASV.
gc.agg <- gc.mtg %>%
  group_by(BLAST, SampleID2) %>%
  summarise(Oral_Total = sum(Oral_count, na.rm = TRUE), Feces_Total = sum(Feces_count, na.rm = TRUE), .groups = 'drop')

# Assigning a Group to a sample
id.group <- mtg %>% select(SampleID2, Group)
gc.agg <- left_join(gc.agg, id.group, by ="SampleID2")

# remove duplication
gc.agg <- gc.agg %>% 
  group_by(SampleID2) %>%
  distinct(BLAST, .keep_all = TRUE)

gc.agg %>% colnames()

# Relative abundance of MTG transmission bacteria
gc.agg <- gc.agg %>% 
  mutate(Rel.abund.mouth = (Oral_Total/15000)*100,
         Rel.abund.feces = (Feces_Total/15000)*100)

# Pivot the data to wide format
gc.agg %>% colnames()
wide.agg <- gc.agg %>%
  select(SampleID2, Group, BLAST, Rel.abund.mouth) %>%
  spread(key = BLAST, value = Rel.abund.mouth)
write.csv(wide.agg, "ogc.rf.csv")

###Making rf format using Core feature table###
#Use core otu table extracted from phyloseq as input file
d <- readRDS("./Zhang.GC.CRC.RDS")

#Select classification group
gc <- d %>% ps_filter(Group %in% c("HC", "GC"))
oral.d <- gc %>% ps_filter(Type == "oral")
fecal.d <- gc %>% ps_filter(Type == "feces")

#Create an Core OTU table and merge it with ASV BLAST results
grm <- core(oral.d, detection = 2, prevalence = 10/100)
grf <- core(fecal.d, detection = 2, prevalence = 10/100)
d1 <- grm %>% otu_tibble()
write.csv(d1, "ogc.core.blast.csv")

#The file core.crcs.blast.xlsx is the BLAST results for the ASVs extracted from the Core OTU table.
d1 <- read.csv("./ogc.core.blast.csv")
d2 <- readxl::read_excel("./core.gc-crc.blast.xlsx", sheet = 3)
d3 <- left_join(d1, d2, by = "FeatureID")

#Aggregate species level
d3 <- d3 %>%
  group_by(BLAST) %>%
  summarise(across(-FeatureID, sum))
d3<- t(d3)
write.csv(d3, "ogc.core.blast.rf.csv")

###Merging internal and external set to make mikropml format ###
yonsei.rf <- read.csv("./")
yonsei.rf <- yonsei.rf %>% filter(Group %in% c("HC", "GC"))
valid1 <- read.csv("./")
valid1 <- valid1 %>% filter(Group %in% c("HC", "GC"))

yonsei.list <- yonsei.rf %>% colnames()
valid1.list <- valid1  %>% colnames()

#Identify common features between internal and external data
common <- intersect(yonsei.list, valid1.list)
#Extract only common features
rf <- yonsei.rf %>% select(all_of(common))
valid1 <- valid1 %>% select(all_of(common))

#Specify internal data set as t1 and external data set as v1
rf$set <- "t1"
valid1$set <- "v1"
#Merging internal and external set
total <- rbind(rf, valid1)
write.csv(total, "totalogc.rf.csv")

### Parallel processing ###
library(foreach)
library(future)
library(future.apply)
library(doFuture)
doFuture::registerDoFuture()
future::plan(future::multicore, workers = 40)

#read input file
df <- read.csv("./totalogc.rf.csv")
#assign internal & external data set
grps <- df$set
#remove set column
colnames(df)
df <- df[ , -86]
#preprocessing
pp <- preprocess_data(dataset = df, outcome_colname = "Group")
ppdf <- pp$dat_transformed

#case weight
set.seed(1004)
train_set_indices <- get_partition_indices(ppdf %>% pull(Group),training_frac = 0.70)

case_weights_dat <- ppdf %>%
  count(Group) %>%
  mutate(p = n / sum(n)) %>%
  select(Group, p) %>%
  right_join(ppdf, by = "Group") %>%
  select(Group, p) %>%
  mutate(
    row_num = row_number(),
    in_train = row_num %in% train_set_indices
  ) %>%
  filter(in_train)


head(case_weights_dat)
tail(case_weights_dat)
nrow(case_weights_dat) / nrow(ppdf)

#Customizing hyperparameters
get_hyperparams_list(ppdf, "rf")
new_hp <- list(
  alpha = 0,
  lambda = c(0.0001, 0.001, 0.01, 0.1, 1, 10)
)

new_hp <- list(
  maxdepth = c(1, 2, 4, 8, 16, 30)
)

new_hp <- list(
  mtry = c(8, 17, 34)
)
new_hp

###Evaluate modeling###
results_weighted <- future.apply::future_lapply(seq(100, 199), function(seed) {
  ml_result <- run_ml(ppdf,
                      "rf",
                      outcome_colname = "Group",
                      seed = seed,
                      training_frac = case_weights_dat %>% pull(row_num),
                      weights = case_weights_dat %>% pull(p),
                      groups = grps,
                      find_feature_importance = TRUE,
                      group_partitions = list(
                        train = c("t1"),
                        test = c("v1")
                      )
  )
}, future.seed = TRUE)

###hyperparameter data###
train_data <- do.call(rbind, lapply(results_weighted[1:100], function(x) data.frame(x$trained_model$results)))
write.csv(train_data, "hyper.ogc.core.csv")

###performance data###
performance_data <- do.call(rbind, lapply(results_weighted[1:100], function(x) data.frame(x$performance)))
write.csv(performance_data, "perfo.ogc.core.csv")

###feature_importance###
feature_importance <- do.call(rbind, lapply(results_weighted[1:100], function(x) data.frame(x$feature_importance)))
write.csv(feature_importance, "feature.ogc.csv")

#making ROC curve plot 
get_sensspec_seed <- function(seed) {
  ml_result <- run_ml(ppdf,
                      "rf",
                      outcome_colname = "Group",
                      seed = seed,
                      training_frac = case_weights_dat %>% pull(row_num),
                      weights = case_weights_dat %>% pull(p),
                      groups = grps,
                      group_partitions = list(
                        train = c("t1"),
                        test = c("v1")
                      )
  )
  sensspec <- calc_model_sensspec(
    ml_result$trained_model,
    ml_result$test_data,
    "Group"
  ) %>%
    dplyr::mutate(seed = seed)
  return(sensspec)
}
sensspec_dat <- purrr::map_dfr(seq(100, 199), get_sensspec_seed)

# calculate mean sensitivity over specificity

roc_dat <- calc_mean_roc(sensspec_dat)

roc_dat %>% plot_mean_roc(ribbon_fill = "#8A2BE2",
                          line_color = "black") 

# calculate mean precision over recall
prc_dat <- calc_mean_prc(sensspec_dat)

baseline_prec <- calc_baseline_precision(ppdf, "Group", "HC")

prc_dat %>%
  plot_mean_prc(baseline_precision = baseline_prec,
                ribbon_fill = "#8A2BE2",
                line_color = "black")