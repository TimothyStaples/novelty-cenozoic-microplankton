# ################################# ####
# Predicting Novel Communities      ####
# Author:    Timothy L Staples      ####
# Collaborators: John Pandolfi      #### 
#                Wolfgang Kiessling ####
# ################################# ####
# Global attributes & working directories ####

rm(list=ls())

# set working directory to folder where this file is located
setwd("/home/timothy/Dropbox/Tim/Post-doc/Research projects/novel_comms_neptune")

# I only have limited Dropbox space so I keep large files in a different folder.
# This is the path. If you want to keep everything together, make this the same
# as the working directory path.
large_file_path <- "/home/timothy/University files - offline/Post-doc/large data file storage/predicting_novel_comms"

# functions ####

# source functions from 'functions' sub-folder
sapply(list.files("./functions", pattern="\\.R", full.names=TRUE), 
       source)

package.loader(c("mgcv", "vegan", "lme4", "nlme", 
                 "DHARMa", "merTools", "shape",
                 "multcomp", "maptools", "sp", 
                 "divDyn", "plotrix", "raster",
                 "rgeos", "fun", "analogue",
                 "brms"))

# A little custom function to add dates to output files
date.wrap <- function(string, ext){
  paste0(string, " ", Sys.Date(), ext)
}

# ####
# IMPORT NEPTUNE MARINE PALEOECOLOGICAL DATA ####

nano <- read.csv(paste0(large_file_path,
                        "/NSB_1010/nannoplankton_2019-10-10_01-06-29.csv"), sep="\t", header=T)

colnames(nano) <- gsub("Resolved\\.", "", colnames(nano))

foram <- read.csv(paste0(large_file_path,
                         "/NSB_1010/foraminifera_2019-10-10_01-38-18.csv"), sep="\t", header=T)
colnames(foram) <- gsub("Resolved\\.", "", colnames(foram))

radio <- read.csv(paste0(large_file_path,
                         "/NSB_1010/radiolarians_2019-10-10_01-46-07.csv"), sep="\t", header=T)
colnames(radio) <- gsub("Resolved\\.", "", colnames(radio))

diatom <- read.csv(paste0(large_file_path,
                         "/NSB_1010/diatoms_2019-10-10_01-58-49.csv"), sep="\t", header=T)
colnames(diatom) <- gsub("Resolved\\.", "", colnames(diatom))

#### Analysis ####
# LONGHURST-LEVEL NOVELTY, TAXA AS RANDOM, LOW RICHNESS COMMS INCLUDED ####
#                     novel detection framework ####

longhurst.novel.list <- lapply(list(nano, foram, radio, diatom),
                               function(x){
                                 
                                 neptune.novelty(dataset = x,
                                                 bin.width = 0.1,
                                                 taxon.res = "Taxon.ID",
                                                 bin.cutoff = 10,
                                                 taxa.cutoff = 10,
                                                 novel.alpha = 0.05,
                                                 novel.metric = "jaccard",
                                                 ssmat.type="pa",
                                                 rich.cutoff=0,
                                                 plot=TRUE,
                                                 longhurst = TRUE)
                                 
                               })

names(longhurst.novel.list) <- c("nano", "foram", "radio", "diatom")

# remove first n time-points (edge effect)
longhurst.novel.cut <- lapply(longhurst.novel.list,
                              function(x){
                                x$novel <- cut.novel(x$novel, 5)
                                return(x)
                              })

# remove problematic samples
longhurst.novel.cut <- cut.complete(longhurst.novel.cut,
                                    cutoff = 0.05,
                                    plot.name = "longhurst_completeness_boxplot.pdf")

saveRDS(longhurst.novel.cut, "./outputs/processedNovelData.rds")

#                     probability of novelty ####

# Probability of shift occurrence
novel.prob.list <- novel.probability(longhurst.novel.cut)

figure1.plot.concept(example.ssmat = longhurst.novel.list[[1]]$ssmats[[length(longhurst.novel.list[[1]]$novel)-21]],
                    plot.name = "longhurst_complete_3col_concept",
                    example.xlims = rev(c(2,51)),
                    novel.lab.pos = 4,
                    ts.cut = 15)
close.screen(all.screens=TRUE)

figure1B.plot(prob.model.list = novel.prob.list,
              plot.name = "longhurst_1B")

#                     observed vs expected transition probabilities ####

# Comparing observed and expected probability of transitions
obs.exp.df <- estimate.observed.expected(prob.model.list = novel.prob.list,
                                         novel.list = longhurst.novel.cut)

obs.exp.df$obs.Estimate.back <- plogis(obs.exp.df$obs.Estimate)
obs.exp.df$obs.lower <- plogis(obs.exp.df$obs.Estimate - 1.96 * obs.exp.df$`obs.Std. Error`)
obs.exp.df$obs.upper <- plogis(obs.exp.df$obs.Estimate + 1.96 * obs.exp.df$`obs.Std. Error`)

write.csv(obs.exp.df, date.wrap(paste0("./outputs/",
                                       " obs vs exp transition"),
                                            ".csv"))

figure2.plot(trans.df = obs.exp.df, 
             plot.name = "longhurst_complete",
             ylims=log(c(0.35,15)))

#                     taxa turnover models ####

# Demographic models
demographic.df <- taxa.turnover.models(novel.list = longhurst.novel.cut,
                                       test.edge.effects=FALSE,
                                       ext.cutoff = 10,
                                       orig.cutoff = 10)

figure3.plot(turnover.models = demographic.df, 
             plot.name = "longhurst_complete",
             left.xlims = c(0,0.2),
             right.xlims = c(0,0.35),
             ylims = c(0,0.83))

close.screen(all.screens=TRUE)
dev.off()
figure3.longplot(turnover.models = demographic.df, 
                 plot.name = "longhurst_complete")

longhurst.analyses <- list(data = longhurst.novel.cut,
                           novel.prob = novel.prob.list,
                           novel.trans = obs.exp.df,
                           demo = demographic.df)

saveRDS(longhurst.analyses, paste0("./outputs/longhurst_analyses.rds"))

#                     taxa-specific obs vs exp probabilities ####

obs.exp.df.taxa <- estimate.observed.expected.fixed(prob.model.list = novel.prob.list,
                                                    novel.list = longhurst.novel.cut)

figure2.plot.fixed(trans.df = obs.exp.df.taxa, 
                   plot.name = "longhurst_fixed_taxa",
                   xlims=log(c(1.2e-4,1)), 
                   ylims=log(c(0.25,35)))

#                     taxa-specific taxa turnover ####

demographic.df <- taxa.turnover.models.fixed(novel.list = longhurst.novel.cut,
                                             test.edge.effects = FALSE,
                                             ext.cutoff = 10,
                                             orig.cutoff = 10)

figure3.plot.fixed(turnover.models = demographic.df,
                   plot.name = "longhurst_fixed_taxa",
                   left.xlim=c(0,0.35), 
                   left.ylim=c(0,0.9),
                   right.xlim=c(0,0.49),
                   right.ylim=c(0,0.99))

#### Analysis variants #### ####

# Site-LEVEL NOVELTY, TAXA AS RANDOM, LOW RICHNESS COMMS EXCLUDED ####

neptune.novelty.list <- lapply(list(nano, foram, radio, diatom),
                                function(x){
                                  
                                  neptune.novelty(dataset = x,
                                                  bin.width = 0.1,
                                                  taxon.res = "Taxon.ID",
                                                  bin.cutoff = 10,
                                                  taxa.cutoff = 10,
                                                  novel.alpha = 0.05,
                                                  novel.metric = "jaccard",
                                                  ssmat.type="pa",
                                                  rich.cutoff=0,
                                                  plot=FALSE,
                                                  longhurst=FALSE)
                                  
                                })
names(neptune.novelty.list) <- c("nano", "foram", "radio", "diatom")

# remove first n time-points (edge effect)
neptune.novelty.cut <- lapply(neptune.novelty.list,
                               function(x){
                                 x$novel <- cut.novel(x$novel, 5)
                                 return(x)
                               })

# remove problematic samples
neptune.novelty.cut <- cut.complete(neptune.novelty.cut,
                                     cutoff = 0.05,
                                     plot.name = "site_completeness_boxplot.pdf")

# Probability of shift occurrence
novel.prob.list <- novel.probability(neptune.novelty.cut)

figure1.plot.simple.3col(prob.model.list = novel.prob.list,
                         example.ssmat = neptune.novelty.cut[[4]]$ssmats[[1]],
                         example.xlims = rev(c(0,20)),
                         plot.name = "site")

# Comparing observed and expected probability of transitions
obs.exp.df <- estimate.observed.expected(prob.model.list = novel.prob.list,
                                         novel.list = neptune.novelty.cut)
figure2.plot(obs.exp.df, 
             plot.name = "site",
             ylims = log(c(0.3,20)))

# Demographic models
demographic.df <- taxa.turnover.models(novel.list = neptune.novelty.cut,
                                       test.edge.effects=FALSE,
                                       ext.cutoff = 10,
                                       orig.cutoff = 10)
figure3.plot(demographic.df, 
             "site",
             left.xlims = c(0,0.275),
             right.xlims = c(0,0.35),
             ylims = c(0,0.8))

saveRDS(list(data = neptune.novelty.cut,
             novel.prob = novel.prob.list,
             novel.trans = obs.exp.df,
             demo = demographic.df),
        paste0(large_file_path, "/site_analyses.rds"))

# LONGHURST-LEVEL NOVELTY, TAXA AS RANDOM, LOW RICHNESS COMMS EXCLUDED ####

neptune.low.rich.list <- lapply(list(nano, foram, radio, diatom),
                                function(x){
                                  
                                  neptune.novelty(dataset = x,
                                                  bin.width = 0.1,
                                                  taxon.res = "Taxon.ID",
                                                  bin.cutoff = 10,
                                                  taxa.cutoff = 10,
                                                  novel.alpha = 0.05,
                                                  novel.metric = "jaccard",
                                                  ssmat.type="pa",
                                                  rich.cutoff=4,
                                                  plot=FALSE,
                                                  longhurst=TRUE)
                                  
                                })
names(neptune.low.rich.list) <- c("nano", "foram", "radio", "diatom")

# remove first n time-points (edge effect)
neptune.low.rich.cut <- lapply(neptune.low.rich.list,
                            function(x){
                              x$novel <- cut.novel(x$novel, 5)
                              return(x)
                            })

# remove problematic samples
neptune.low.rich.cut <- cut.complete(neptune.low.rich.cut,
                                  cutoff = 0.05,
                                  plot.name = "low_richness_completeness_boxplot.pdf")

# Probability of shift occurrence
novel.prob.list <- novel.probability(neptune.low.rich.cut)

figure1.plot.simple.3col(prob.model.list = novel.prob.list,
                    example.ssmat = neptune.low.rich.cut[[4]]$ssmats[[1]],
                    example.xlims = rev(c(0,20)),
                    plot.name = "longhurst_low_richness")

# Comparing observed and expected probability of transitions
obs.exp.df <- estimate.observed.expected(prob.model.list = novel.prob.list,
                                              novel.list = neptune.low.rich.cut)
figure2.plot(obs.exp.df, 
             plot.name = "longhurst_low_richness",
             ylims = log(c(0.5,20)))

# Demographic models
demographic.df <- taxa.turnover.models(novel.list = neptune.low.rich.cut,
                                            test.edge.effects=FALSE,
                                            ext.cutoff = 10,
                                            orig.cutoff = 10)
figure3.plot(demographic.df, 
             "longhurst_low_richness",
             left.xlims = c(0,0.275),
             right.xlims = c(0,0.35),
             ylims = c(0,0.8))

saveRDS(list(data = neptune.low.rich.cut,
             novel.prob = novel.prob.list,
             novel.trans = obs.exp.df,
             demo = demographic.df),
        paste0(large_file_path, "/low_richness_analyses.rds"))

# LONGHURST NOVELTY, TAXA AS RANDOM, LOW RICHNESS INCLUDED, VARYING BIN WIDTHS ####

width.analyses <- lapply(seq(0.1,0.5,0.1), function(width){
  
  print(width)
  # Identify novel communities
  neptune.novel.list <- lapply(list(nano, foram, radio, diatom),
                               function(x){
                                 
                                 neptune.novelty(dataset = x,
                                                 bin.width = width,
                                                 taxon.res = "Taxon.ID",
                                                 bin.cutoff = 10,
                                                 taxa.cutoff = 10,
                                                 novel.alpha = 0.05,
                                                 novel.metric = "jaccard",
                                                 ssmat.type="pa",
                                                 rich.cutoff=0,
                                                 plot=FALSE,
                                                 longhurst=TRUE)
                                 
                               })
  names(neptune.novel.list) <- c("nano", "foram", "radio", "diatom")
  
  # remove first n time-points (edge effect)
  neptune.novel.cut <- lapply(neptune.novel.list,
                              function(x){
                                x$novel <- cut.novel(x$novel, 5)
                                return(x)
                              })
  
  neptune.novel.cut <- cut.complete(neptune.novel.cut,
                                    cutoff = 0.05,
                                    plot.name = paste0("bin_width (", width, ")_completeness_boxplot.pdf"))
  
  # Probability of shift occurrence
  novel.prob.list <- novel.probability(neptune.novel.cut)
 
  # Comparing observed and expected probability of transition
  obs.exp.df <- estimate.observed.expected(prob.model.list = novel.prob.list,
                                                novel.list = neptune.novel.cut)

  # Demographic models
  demographic.df <- taxa.turnover.models(novel.list = neptune.novel.cut,
                                         test.edge.effects=FALSE,
                                         ext.cutoff = 10,
                                         orig.cutoff = 10)
  
  
  return(list(data = neptune.novel.cut,
              novel.prob = novel.prob.list,
              novel.trans = obs.exp.df,
              demo = demographic.df))
  
  
})
names(width.analyses) <- seq(0.1,0.5,0.1)

lapply(1:length(width.analyses), function(n){
  print(n)
  temp.analyses <- width.analyses[[n]]
  temp.name <- names(width.analyses)[n]

  figure1.plot.simple.3col(prob.model.list = temp.analyses$novel.prob,
                      example.ssmat = temp.analyses$data[[4]]$ssmats[[1]],
                      plot.name = paste0("varying bin width: ", temp.name),
                      example.xlims = rev(c(0,20)))

  figure2.plot(trans.df = temp.analyses$novel.trans,
               plot.name = paste0("varying bin width: ", temp.name),
               ylims = log(c(0.45, 15)))
  
  figure3.plot(temp.analyses$demo, 
               paste0("varying bin width: ", temp.name), 
               left.xlims=c(0,c(0.4,0.4,0.4,0.4,0.4)[n]),
               right.xlims = c(0,c(0.4,0.4,0.4,0.5,0.5)[n]),
               ylims = c(0,0.8))
  
})

saveRDS(width.analyses,
        paste0(large_file_path, "/bin_width_analyses.rds"))

# LONGHURST NOVELTY, TAXA AS RANDOM, LOW RICHNESS INCLUDED, 5% PACMAN EXCLUSION ####

# Import NSB data with 5% Pacman exclusions
nano.pac <- read.csv(paste0(large_file_path,
                        "/NSB_1010/nannoplankton_pacman_2019-11-26_02-37-38.csv"), sep="\t", header=T)
colnames(nano.pac) <- gsub("Resolved\\.", "", colnames(nano.pac))

foram.pac <- read.csv(paste0(large_file_path,
                         "/NSB_1010/foraminifera_pacman_2019-11-26_05-12-53.csv"), sep="\t", header=T)
colnames(foram.pac) <- gsub("Resolved\\.", "", colnames(foram.pac))

radio.pac <- read.csv(paste0(large_file_path,
                         "/NSB_1010/radiolarians_pacman_2019-11-26_04-16-50.csv"), sep="\t", header=T)
colnames(radio.pac) <- gsub("Resolved\\.", "", colnames(radio.pac))

diatom.pac <- read.csv(paste0(large_file_path,
                          "/NSB_1010/diatom_pacman_2019-11-26_04-32-40.csv"), sep="\t", header=T)
colnames(diatom.pac) <- gsub("Resolved\\.", "", colnames(diatom.pac))

# Identify novel communities
neptune.pac.list <- lapply(list(nano.pac, foram.pac, radio.pac, diatom.pac),
                             function(x){
                               
                               neptune.novelty(dataset = x,
                                               bin.width = 0.1,
                                               taxon.res = "Taxon.ID",
                                               bin.cutoff = 10,
                                               taxa.cutoff = 10,
                                               novel.alpha = 0.05,
                                               novel.metric = "jaccard",
                                               ssmat.type="pa",
                                               rich.cutoff=0,
                                               plot=FALSE,
                                               longhurst=TRUE)
                               
                             })
names(neptune.pac.list) <- c("nano", "foram", "radio", "diatom")

# remove first n time-points (edge effect)
neptune.pac.cut <- lapply(neptune.pac.list,
                            function(x){
                              x$novel <- cut.novel(x$novel, 5)
                              return(x)
                            })

neptune.pac.cut <- cut.complete(neptune.pac.cut,
                                  cutoff = 0.05,
                                  plot.name = "pacman_completeness_boxplot.pdf")

# Probability of shift occurrence
novel.prob.list <- novel.probability(neptune.pac.cut)

figure1.plot.simple.3col(prob.model.list = novel.prob.list,
                    example.ssmat = neptune.pac.cut[[4]]$ssmats[[1]],
                    example.xlims = rev(c(2.5,30)),
                    plot.name = "pacman")

# Comparing observed and expected probability of transition
obs.exp.df <- estimate.observed.expected(prob.model.list = novel.prob.list,
                                         novel.list = neptune.pac.cut)

figure2.plot(trans.df = obs.exp.df, 
             plot.name = "pacman",
             ylims = log(c(0.45, 30)))

# Demographic models
demographic.df <- taxa.turnover.models(novel.list = neptune.pac.cut,
                                       test.edge.effects=FALSE,
                                       ext.cutoff = 10,
                                       orig.cutoff = 10)

figure3.plot(demographic.df, "pacman", 
             left.xlims=c(0, 0.25),
             right.xlims = c(0, 0.4),
             ylims = c(0, 1))

saveRDS(list(data = neptune.pac.cut,
             novel.prob = novel.prob.list,
             novel.trans = obs.exp.df,
             demo = demographic.df),
        paste0(large_file_path, "/pacman_analyses.rds"))

# ####
# SUPPLEMENTARY ANALYSES ####
#       Correlation between dissimilarity metrics ####

comb.df <- do.call("rbind", lapply(longhurst.novel.list, 
                                   function(x){do.call("rbind", x$novel)}))

cor.diss <- cor.test(comb.df$seq.dist, comb.df$raw.min.dist)
summary(cor.diss)

#       Probability based on time-series position ####

test.time.series.sensitivity(longhurst.novel.list,
                             longhurst.novel.cut)

#       Probability based on richness ####

test.richness.sensitivity(all.novel.list = longhurst.novel.cut)

#       Novel community probability with different thresholds ####

thresholds <- c(seq(0.01,0.1,0.005))
threshold.list <- lapply(thresholds, function(cutoff){
  
mp.novel.list <- lapply(list(nano, foram, radio, diatom),
                                  function(x){
                                    
                                    x.slim <- x[,c("Taxon.ID", 
                                                   "Site", 
                                                   "Longitude",
                                                   "Latitude",
                                                   "Age..Ma..Gradstein.et.al..2012")]
                                    
                                    neptune.novelty(dataset = x.slim,
                                                    bin.width = 0.1,
                                                    taxon.res = "Taxon.ID",
                                                    bin.cutoff = 10,
                                                    taxa.cutoff = 10,
                                                    novel.alpha = cutoff,
                                                    novel.metric = "jaccard",
                                                    ssmat.type="pa",
                                                    rich.cutoff=0,
                                                    plot=FALSE,
                                                    longhurst=TRUE)
                                    
                                    })

mp.novel.list <- lapply(mp.novel.list,
                            function(x){
                              x$novel <- cut.novel(x$novel, 5)
                              return(x)
                              })
names(mp.novel.list) <- c("nano", "foram", "radio", "diatom")
  
return(do.call("rbind", lapply(1:4, function(n){
    
    temp.novel <- mp.novel.list[[n]]$novel
    temp.novel <- do.call("rbind", temp.novel)
    temp.novel$taxa <- c("nano", "foram", "radio", "diatom")[n]
    return(temp.novel)
  })))
  
})
threshold.prob.list <- lapply(threshold.list, function(x){

  x$site <- paste0(x$taxa, ":", x$site)
  
  novel.ts.freq <- do.call("rbind", lapply(unique(x$site),
                                  function(site){
                                    
                                    temp <- x[x$site == site,]
                                    
                                    temp <- data.frame(site = site,
                                                       t(sapply(c("back", "instant", "cumul", "novel"),
                                                                function(y){sum(temp$cat == y)})))
                                    temp$taxa <- x$taxa[x$site == site][1]
                                    temp$lat <- x$data$lat[x$site == site][1]
                                    temp$long <- x$data$long[x$site == site][1]
                                    
                                    return(temp)
                                  }))
  novel.ts.freq <- as.data.frame(novel.ts.freq)
  
  # modelling the probability of each classification occurring.
  novel.ts.freq$non.novel <- rowSums(novel.ts.freq[,c("back", "cumul", "instant")])
  novel.ts.freq$non.back <- rowSums(novel.ts.freq[,c("novel", "cumul", "instant")])
  novel.ts.freq$non.instant <- rowSums(novel.ts.freq[,c("novel", "cumul", "back")])
  novel.ts.freq$non.cumul <- rowSums(novel.ts.freq[,c("novel", "back", "instant")])
  
  ## modeling the probability of our two GAM tests.
  novel.ts.freq$all.instant <- rowSums(novel.ts.freq[,c("novel", "instant")])
  novel.ts.freq$all.cumul <- rowSums(novel.ts.freq[,c("novel", "cumul")])
  novel.ts.freq$non.all.instant <- rowSums(novel.ts.freq[,c("cumul", "back")])
  novel.ts.freq$non.all.cumul <- rowSums(novel.ts.freq[,c("instant", "back")])
  
  novel.ts.freq$taxa <- as.factor(novel.ts.freq$taxa)
  
  fixed.prob.models <- fixed.taxa.prob.models.long(novel.freq.df = novel.ts.freq,
                                                   test.model=TRUE)
  
  random.prob.models <- random.taxa.prob.models(novel.freq.df = novel.ts.freq,
                                                test.model=TRUE)
  
  return(list(data = novel.ts.freq,
              fixed.prob.models = fixed.prob.models,
              random.prob.models = random.prob.models))
  
})

saveRDS(threshold.prob.list, "./outputs/alpha_threshold_prob_list.rds")

threshold.prob.list <- readRDS("./outputs/alpha_threshold_prob_list.rds")

probs = thresholds
varying.alpha.test(probs = thresholds,
                   threshold.prob.list = threshold.prob.list,
                   circle.1 = 0.01,
                   circle.2 = 0.05,
                   circle.3 = 0.1,
                   ylims = c(0.005,max(thresholds)+0.0035))
  
  
#       Simulation/Sensitivity test of novel detection method ####

# this test assumes all communities are 'background', and identifies how often
# novelty arises under these conditions

null.sim.results <- novel.sensitivity.test.null(observed.data = longhurst.analyses$data,
                                                probability.models = longhurst.analyses$novel.prob$random.prob.models,
                                                turnover.models = longhurst.analyses$demo,
                                                nsims = 999)

null.sim.results <- novel.sensitivity.test.controls(observed.data = longhurst.analyses$data,
                                                    probability.models = longhurst.analyses$novel.prob$random.prob.models,
                                                    turnover.models = longhurst.analyses$demo,
                                                    nsims = 999)

null.sim.results$cat <- as.factor(null.sim.results$cat)

sim.df <- do.call("rbind", lapply(unique(null.sim.results$sim),
                                function(site){
                                  print(site)
                                  temp <- null.sim.results[null.sim.results$sim==site,]
                                  
                                  temp.cut <- temp[-(1:5),]
                                  
                                  temp.cut <- data.frame(site = site,
                                                     t(sapply(c("back", "instant", "cumul", "novel"),
                                                              function(y){sum(temp.cut$cat == y)})))

                                  return(temp.cut)
                                }))

sim.df$non.novel <- rowSums(sim.df[,c("back", "cumul", "instant")])
sim.df$non.back <- rowSums(sim.df[,c("novel", "cumul", "instant")])
sim.df$non.instant <- rowSums(sim.df[,c("novel", "cumul", "back")])
sim.df$non.cumul <- rowSums(sim.df[,c("novel", "back", "instant")])

## modeling the probability of our two GAM tests.
sim.df$all.instant <- rowSums(sim.df[,c("novel", "instant")])
sim.df$all.cumul <- rowSums(sim.df[,c("novel", "cumul")])
sim.df$non.all.instant <- rowSums(sim.df[,c("cumul", "back")])
sim.df$non.all.cumul <- rowSums(sim.df[,c("instant", "back")])

sim.models <- lapply(1:4, function(n){
  
    success.var = c("back", "instant", "cumul", "novel")[n]
    failure.var = c("non.back", "non.instant", "non.cumul", "non.novel")[n]
    
    print(success.var)
    
    temp.df <- sim.df
    temp.df$success = temp.df[, success.var]
    temp.df$failure = temp.df[, failure.var]
    temp.df$trials <- temp.df$success + temp.df$failure
    
    taxa.prob.m <- glm(cbind(success, failure) ~ 1, data=temp.df, family=binomial)
    
    # group predictions
    pred.df <- as.data.frame(summary(taxa.prob.m)$coefficients)
    pred.df$taxa.rand <- summary(taxa.prob.m)$varcor$taxa[1,1]
    
    return(list(model=taxa.prob.m,
                pred.df = pred.df))
  
  
})

obs.df <- do.call("rbind", lapply(1:4, function(n){
  
  temp <- do.call('rbind', longhurst.analyses$data[[n]]$novel)
  temp$taxa <- c("nano", "foram", "radio", "diatom")[n]
  return(temp)
    
}))

t(sapply(novel.prob.list$random.prob.models, 
       function(x){
         temp.x <- summary(x$model)$coefficients
         cbind(plogis(temp.x[,1]),
               plogis(temp.x[,1] + 1.96 * temp.x[,2]),
               plogis(temp.x[,1] - 1.96 * temp.x[,2]))
       }))

t(sapply(sim.models, 
         function(x){
           temp.x <- summary(x$model)$coefficients
           cbind(plogis(temp.x[,1]),
                 plogis(temp.x[,1] + 1.96 * temp.x[,2]),
                 plogis(temp.x[,1] - 1.96 * temp.x[,2]))
         }))

write.csv(sim.results, "./outputs/novel simulation results (9999).csv")

#       Test of novel detection method with 'synthetic communities' ####

syn.test()

#       Novelty with different dissimilarity indices (Figure S4) ####

metric.vect <- c("bray", "jaccard", "gower", "kulczynski", "SQchord", "chord")
nano.metrics <- lapply(metric.vect,
                       function(metric){
                         print(metric)
                         lapply(list(nano, foram, radio, diatom),
                                function(x){
                                  
                                  neptune.novelty(dataset = x,
                                                  bin.width = 0.1,
                                                  taxon.res = "Taxon.ID",
                                                  bin.cutoff = 10,
                                                  taxa.cutoff = 10,
                                                  novel.alpha = 0.05,
                                                  novel.metric = metric,
                                                  ssmat.type="pa",
                                                  rich.cutoff=0,
                                                  plot=FALSE,
                                                  longhurst = TRUE)
                                  
                                })
                         
                       })

metric.prob <- lapply(nano.metrics, function(sub.df){
  novel.probability(sub.df)
})

pdf(date.wrap("./plots/metric comparison: combined dissimilarity pairs",".pdf"), 
    height=7, width=7.5, useDingbats = FALSE)

dissim.labels = c("Sorenson",
                  "Jaccard",
                  "Gower",
                  "Kulczynski",
                  "Chord^2",
                  "Chord")

novel.match <- lapply(c("instant", "cumul", "novel"),
                      function(cat){
                        print(cat)
                        # for each metrics
                        temp.novel <- do.call("cbind", lapply(nano.metrics, function(x){
                          
                          # for each taxa group
                          temp.novel <- do.call("rbind", lapply(x, function(x1){
                            do.call("rbind", x1$novel)
                          }))
                          
                          temp.novel$cat == cat
                          
                        }))
                        
                        col.comb <- expand.grid(1:ncol(temp.novel),
                                                1:ncol(temp.novel))
                        col.comb$metric1 <- metric.vect[col.comb[,1]]
                        col.comb$metric2 <- metric.vect[col.comb[,2]]
                        
                        col.comb$sum1 <- colSums(temp.novel[,col.comb[,1]] &
                                                   !temp.novel[,col.comb[,2]])
                        
                        col.comb$sum2 <- colSums(temp.novel[,col.comb[,2]] &
                                                   !temp.novel[,col.comb[,1]])
                        
                        col.comb$both.sum <- colSums(temp.novel[,col.comb[,1]] &
                                                       temp.novel[,col.comb[,2]])
                        
                        col.comb$neither.sum <- colSums(!temp.novel[,col.comb[,1]] &
                                                          !temp.novel[,col.comb[,2]] &
                                                          rowSums(temp.novel) > 0)
                        
                        return(col.comb)
                        
                      })

# Instantaneous novelty
temp.match <- novel.match[[1]]
temp.prob <- do.call("rbind", lapply(metric.prob, function(x){
  x$random.prob.models[[4]]$pred.df
}))

par(mar=c(0,0,0,0), mfrow=c(length(metric.vect), length(metric.vect)),
    oma=c(2,4,4,7), tcl = -0.25, las=1, mgp=c(3,0.5,0))

sapply(1:length(metric.vect), function(j){
  
  sapply(1:length(metric.vect), function(k){
    
    if(j==k){
      
      plot.new()
      text(x=0.5, y=0.6, 
           labels=paste0(sprintf("%.1f", plogis(temp.prob[j,"Estimate"])*100),
                         "%"),
           cex=2)
      text(x=0.5, y=0.3, 
           labels=paste0("(", 
                         sprintf("%.1f", plogis(temp.prob[j,"l-95% CI"])*100),
                         " - ",
                         sprintf("%.1f", plogis(temp.prob[j,"u-95% CI"])*100),
                         "%)"),
           cex=1.5)
      
      
      box()}
    
    if(j > k){
      
      sub.match <- temp.match[temp.match$Var1==j & temp.match$Var2 == k,]
      
      stats <- unlist(sub.match[1,c(7,5,6,8)])
      
      labels <- c("Both", "Row", "Col", "Neither")[stats > 0]
      cols <- colorRampPalette(c("red","white"))(4)[stats>0]
      stats<-stats[stats > 0]
      
      x<-barplot(cbind(stats), xlim=c(-0.25,1.75),
                 ylim=c(0, sum(stats)*1.1), yaxs="i", 
                 names.arg="", width=0.75,
                 axes=FALSE,
                 col=cols)
      text(x=x, y= cumsum(stats) - 0.5 * stats,
           labels=ifelse(stats / sum(stats) > 0.1,
                         sprintf("%.2f", round(stats / sum(stats), 2)), 
                         ""),
           font=2, cex=0.9)
      
      text(x=x+(0.75/2), y= cumsum(stats) - 0.5 * stats,
           labels=labels, 
           adj=0, pos=4, offset=0.2)
      box()
      
    }
    
    if(j < k){
      
      plot.new()
      
    }
    
    if(j==length(metric.vect)){
      
      mtext(side=1, line=0.1,
            text=dissim.labels[k],
            cex=0.8)
    }
    
    if(k==1){
      mtext(side=2, line=0.2,
            text=dissim.labels[j],
            cex=0.8, las=0)
      
    }
    
  })
  
})

# Cumulative novelty
temp.match <- novel.match[[2]]
temp.prob <- do.call("rbind", lapply(metric.prob, function(x){
  x$random.prob.models[[5]]$pred.df
}))

par(mar=c(0,0,0,0), mfrow=c(length(metric.vect), length(metric.vect)),
    oma=c(2,4,4,7), tcl = -0.25, las=1, mgp=c(3,0.5,0))

sapply(1:length(metric.vect), function(j){
  
  sapply(1:length(metric.vect), function(k){
    
    if(j==k){
      
      plot.new()
      text(x=0.5, y=0.6, 
           labels=paste0(sprintf("%.1f", plogis(temp.prob[j,"Estimate"])*100),
                         "%"),
           cex=2)
      text(x=0.5, y=0.3, 
           labels=paste0("(", 
                         sprintf("%.1f", plogis(temp.prob[j,"l-95% CI"])*100),
                         " - ",
                         sprintf("%.1f", plogis(temp.prob[j,"u-95% CI"])*100),
                         "%)"),
           cex=1.5)
      
      
      box()}
    
    if(j > k){
      
      sub.match <- temp.match[temp.match$Var1==j & temp.match$Var2 == k,]
      
      stats <- unlist(sub.match[1,c(7,5,6,8)])
      
      labels <- c("Both", "Row", "Col", "Neither")[stats > 0]
      cols <- colorRampPalette(c("blue","white"))(4)[stats>0]
      stats<-stats[stats > 0]
      
      x<-barplot(cbind(stats), xlim=c(-0.25,1.75),
                 ylim=c(0, sum(stats)*1.1), yaxs="i", 
                 names.arg="", width=0.75,
                 axes=FALSE,
                 col=cols)
      text(x=x, y= cumsum(stats) - 0.5 * stats,
           labels=ifelse(stats / sum(stats) > 0.1,
                         sprintf("%.2f", round(stats / sum(stats), 2)), 
                         ""),
           font=2, cex=0.9,
           col=c("white", "black", "black", "black"))
      
      text(x=x+(0.75/2), y= cumsum(stats) - 0.5 * stats,
           labels=labels, 
           adj=0, pos=4, offset=0.2)
      box()
      
    }
    
    if(j < k){
      
      plot.new()
      
    }
    
    if(j==length(metric.vect)){
      
      mtext(side=1, line=0.1,
            text=dissim.labels[k],
            cex=0.8)
    }
    
    if(k==1){
      mtext(side=2, line=0.2,
            text=dissim.labels[j],
            cex=0.8, las=0)
      
    }
    
  })
  
})

# Novel communities
# Cumulative novelty
temp.match <- novel.match[[3]]
temp.prob <- do.call("rbind", lapply(metric.prob, function(x){
  x$random.prob.models[[6]]$pred.df
}))

par(mar=c(0,0,0,0), mfrow=c(length(metric.vect), length(metric.vect)),
    oma=c(2,4,4,7), tcl = -0.25, las=1, mgp=c(3,0.5,0))

sapply(1:length(metric.vect), function(j){
  
  sapply(1:length(metric.vect), function(k){
    
    if(j==k){
      
      plot.new()
      text(x=0.5, y=0.6, 
           labels=paste0(sprintf("%.1f", plogis(temp.prob[j,"Estimate"])*100),
                         "%"),
           cex=2)
      text(x=0.5, y=0.3, 
           labels=paste0("(", 
                         sprintf("%.1f", plogis(temp.prob[j,"l-95% CI"])*100),
                         " - ",
                         sprintf("%.1f", plogis(temp.prob[j,"u-95% CI"])*100),
                         "%)"),
           cex=1.5)
      
      
      box()}
    
    if(j > k){
      
      sub.match <- temp.match[temp.match$Var1==j & temp.match$Var2 == k,]
      
      stats <- unlist(sub.match[1,c(7,5,6,8)])
      
      labels <- c("Both", "Row", "Col", "Neither")[stats > 0]
      cols <- colorRampPalette(c("orange","white"))(4)[stats>0]
      stats<-stats[stats > 0]
      
      x<-barplot(cbind(stats), xlim=c(-0.25,1.75),
                 ylim=c(0, sum(stats)*1.1), yaxs="i", 
                 names.arg="", width=0.75,
                 axes=FALSE,
                 col=cols)
      text(x=x, y= cumsum(stats) - 0.5 * stats,
           labels=ifelse(stats / sum(stats) > 0.1,
                         sprintf("%.2f", round(stats / sum(stats), 2)), 
                         ""),
           font=2, cex=0.9)
      
      text(x=x+(0.75/2), y= cumsum(stats) - 0.5 * stats,
           labels=labels, 
           adj=0, pos=4, offset=0.2)
      box()
      
    }
    
    if(j < k){
      
      plot.new()
      
    }
    
    if(j==length(metric.vect)){
      
      mtext(side=1, line=0.1,
            text=dissim.labels[k],
            cex=0.8)
    }
    
    if(k==1){
      mtext(side=2, line=0.2,
            text=dissim.labels[j],
            cex=0.8, las=0)
      
    }
    
  })
  
})
dev.off()

#       Novel climate correlation ####

# this function requires the raw data from Zachos et al 2001
# can be found here https://www.ncdc.noaa.gov/paleo-search/study/8674
# See function script to see where and how to name the data file.

novel.climate.correlation(novel.list = longhurst.novel.cut)

#       Sampling map ####

# this requires additional shape files etc and R packages, as well as 
# post-R processing as rgl() had issues displaying axis text correctly
# What I'm trying to say is... good luck!

comb.dat <- do.call("rbind", lapply(list(nano.novel$data, foram.novel$data, 
                                         diatom.novel$data, radio.novel$data), 
                                    function(x){
                                      
                                      x[,c("site", "lat", "long", "age")]
                                      
                                    }))
site.dat <- do.call("rbind", lapply(split(comb.dat, f=comb.dat$site), function(x){
  
  data.frame(site = x$site[1],
             lat = x$lat[1],
             long = x$long[1],
             start = max(x$age, na.rm=TRUE),
             end = min(x$age, na.rm=TRUE))
  
}))

write.csv(site.dat, "./outputs/site.data.csv")

site.dat <- read.csv("./outputs/site.data.csv")

library(rgl)
library(maptools)
library(raster)
library(plot3D)
library(rgdal)
library(rgeos)
library(plotrix)

#zoom<-par3d()$zoom
#userMatrix<-par3d()$userMatrix
#windowRect<-par3d()$windowRect
#write.csv(userMatrix, "./outputs/3dmat.csv")

library(rworldmap)
countries

outline <- readShapeSpatial("/home/timothy/Dropbox/Tim/Post-doc/Research projects/novel_comms/raw.datafiles/ne_110m_land.shp")
longhurst <- readShapeSpatial("/home/timothy/Dropbox/Tim/Post-doc/Research projects/novel_comms/raw.datafiles/longhurst_v4_2010/Longhurst_world_v4_2010.shp")
longhurstpoly <- disaggregate(longhurst)

#Polygons <- slot(coastsCoarse,"polygons")
open3d(zoom = 0.35532, userMatrix = as.matrix(read.csv("./outputs/3dmat.csv")[,-1]), 
       windowRect=c(2821,24,4680,1080))

plot3d(c(-250,250), c(-130,130), c(0,65), xaxs="i", yaxs="i", axes=FALSE, 
       xlab="", ylab="", zlab="")
aspect3d(1.65, 1, 0.1)
highlevel()  # To trigger display
polygon3d(x=c(-180,-180,180,180),
          y=c(-90,90,90,-90),
          z=rep(0, 4), fill=TRUE,
          col="white",lit=FALSE)

for (i in seq_along(outline@polygons)){
  print(i)
  #temp <- slot(slot(Polygons[[i]], "Polygons")[[1]], "coords")
  temp <- slot(slot(outline@polygons[[i]], "Polygons")[[1]], "coords")
  
  polygon3d(x=temp[,1], y=temp[,2], z=rep(0.01, dim(temp)[1]), fill=TRUE,
            col="grey75",lit=FALSE)
  #polygon3d(x=temp[,1], y=temp[,2], z=rep(0.01, dim(temp)[1]),
   #         lit=FALSE)
  
}

for (i in seq_along(longhurstpoly)){
  print(i)
  #temp <- slot(slot(Polygons[[i]], "Polygons")[[1]], "coords")
  temp <- slot(slot(longhurstpoly@polygons[[i]], "Polygons")[[1]], "coords")
  
  polygon3d(x=temp[,1], y=temp[,2], z=rep(0.01, dim(temp)[1]), fill=FALSE,
            lit=FALSE, lwd=3)
  
}


# Long axis
sapply(seq(-180,180,30), function(z){
  lines3d(x=rep(z, 2), y=c(-90,-93), z=c(0,0), lwd=2)
  # text3d(x=z, y=-98.5, z=0,
  #          texts=parse(text=paste0(z, "*degree")), adj=0.25, cex=1.5)
  
  # add in grey line across plot
  # if(z > -180 & z < 180){
  #   lines3d(x=rep(z, 2), y=c(-90,90), z=c(0,0), lwd=1, col="grey60")
  # }
})
lines3d(x=c(-180,180), y=c(-90,-90), z=c(0.02,0.02), lwd=1)
lines3d(x=c(-180,180), y=c(90,90), z=c(0.02,0.02), lwd=2)

# Lat axis
sapply(seq(-90,90,15), function(z){
  
  lines3d(x=c(-180, -184), y=rep(z, 2), z=c(0,0), lwd=2)
  lines3d(x=c(180, 184), y=rep(z, 2), z=c(0,0), lwd=2)
  # 
  # text3d(x=-195, y=z, z=0,
  #        texts=parse(text=paste0(z, "*degree")), adj=0.15, cex=1.5)
  # text3d(x=187.5, y=z, z=0,
  #        texts=parse(text=paste0(z, "*degree")), adj=0.15, cex=1.5)
  
  # # add in grey line across plot
  # if(z < 90 & z > -90){
  #   lines3d(x=c(-180,180), y=rep(z, 2), z=c(0,0), lwd=1, col="grey60")
  # }
})
lines3d(x=c(180,180), y=c(-90,90), z=c(0.02,0.02), lwd=2)
lines3d(x=c(-180,-180), y=c(-90,90), z=c(0.02,0.02), lwd=2)

# Cylinders
sapply(1:dim(site.dat)[1], function(n){
  
  print(n)
  
  temp.dat <- site.dat[n,]
  
  radius<-0.75
  scaling <- c(360/360, (360/180)/1.65)
  degs <- seq(0, pi*2, length.out=10)
  cyl.x = (radius * cos(degs)) * scaling[1]
  cyl.y = (radius * sin(degs)) * scaling[2]
  
  # generate several cylinders with different shading to represent 10MY
  tens <- c(seq(10,60,10))
  bottoms <- c(temp.dat$start, rev(tens[tens > temp.dat$end & tens < temp.dat$start]), temp.dat$end)
  z.place <- abs(bottoms - bottoms[1])
  
  cols <- data.frame(top=c(10,20,30,40,50,60,65),
                     bottom = c(0,10,20,30,40,50,60),
                     col = colorRampPalette(c(rgb(0.85,0.85,1), "cornflowerblue", "darkblue"),
                                            bias=3)(7),
                     stringsAsFactors = FALSE)
  
  cyl.coords <- do.call("rbind", lapply(2:length(bottoms), function(z){
    
    coords <- cbind(rep(temp.dat$long, 20),
                    rep(temp.dat$lat, 20),
                    seq(z.place[z-1], z.place[z], length.out=20))  
    
    cyls <- cylinder3d(center=coords, closed= -2, radius=1,
                       section=cbind(cyl.x, cyl.y))
    
    temp.col <- cols$col[cols$top == bottoms[z-1] |
                           cols$bottom == bottoms[z]]
    
    if(length(temp.col)==0){
      temp.col <- cols$col[cols$top > bottoms[z-1] &
                             cols$bottom < bottoms[z]]
    }
    
    shade3d(cyls, col = temp.col, lit=FALSE)
    #shade3d(addNormals(subdivision3d(cyls, depth = 2)), col = "green")
    return(coords)
  }))
  
  # plot some lines to highlight cylinders
  polygon3d(x=c(cyl.coords[1,1] + 1.15*cyl.x),
            y=c(cyl.coords[1,2] + 1.15*cyl.y),
            z=rep(max(cyl.coords[,3])+0.01, length(cyl.x)),
            border="black",lit=FALSE, lwd=3, fill=FALSE)
  
  polygon3d(x=c(cyl.coords[1,1] + 1.175*cyl.x),
            y=c(cyl.coords[1,2] + 1.175*cyl.y),
            z=rep(0.01, length(cyl.x)),
            border="black",lit=FALSE, lwd=3, fill=FALSE)
  
  lines3d(x=c(cyl.coords[1,1] + 1.4*cyl.x[1],
              cyl.coords[1,1] + 1.4*cyl.x[1]),
          y=c(cyl.coords[1,2] + 1.4*cyl.y[1],
              cyl.coords[1,2] + 1.4*cyl.y[1]),
          z=c(0, max(cyl.coords[,3])+0.025), lit=FALSE, lwd=3,
          border="black")
  
  lines3d(x=c(cyl.coords[1,1] + 1.3*cyl.x[6],
              cyl.coords[1,1] + 1.3*cyl.x[6]),
          y=c(cyl.coords[1,2] + 1.3*cyl.y[6],
              cyl.coords[1,2] + 1.3*cyl.y[6]),
          z=c(0, max(cyl.coords[,3])+0.02), lit=FALSE, lwd=3,
          border="black")
  
  plot(x=NULL, y=NULL, xlim=c(140, 143), ylim=c(0.5,3))
  text(x=cyl.coords[1,1] + cyl.x,
       y=cyl.coords[1,2] + cyl.y,
       labels=1:10)
  
  # segments3d(x=rep(cyl.coords[1,1] + cyl.x[1],2),
  #            y=rep(cyl.coords[1,2] + cyl.y[1],2),
  #            z=c(0,max(cyl.coords[,3])), lwd=6)
  # 
  # segments3d(x=rep(cyl.coords[1,1] + cyl.x[6],2),
  #            y=rep(cyl.coords[1,2] + cyl.y[6],2),
  #            z=c(0,max(cyl.coords[,3])), lwd=6)
  # 
  # segments3d(x=cyl.coords[1,1] + cyl.x,
  #         y=cyl.coords[1,2] + cyl.y,
  #         z=rep(max(cyl.coords[,3]), length(cyl.y)), lwd=6)
  
})

cols <- data.frame(top=c(10,20,30,40,50,60,65),
                   bottom = c(0,10,20,30,40,50,60),
                   col = colorRampPalette(c(rgb(0.85,0.85,1), "cornflowerblue", "darkblue"),
                                          bias=3)(7),
                   stringsAsFactors = FALSE)

legend3d("topright", 
         legend = rev(c("60 - 65",
                        "50 - 60",
                        "40 - 50",
                        "30 - 40",
                        "20 - 30",
                        "10 - 20",
                        "0 - 10")), fill=cols$col, 
         cex=2, inset=c(0.04), bty="n",
         title = "Millions of years\nbefore present")

rgl.postscript(date.wrap("./plots/persp3dd", ".pdf"),"pdf", drawText=FALSE) 

# rgl.postscript("./plots/persp3dd.eps","eps", drawText=FALSE) 
# rgl.postscript("./plots/persp3dd.svg","svg", drawText=FALSE) 
# rgl.postscript("./plots/persp3dd.ps","ps", drawText=FALSE) 
# rgl.postscript("./plots/persp3dd.tex","tex", drawText=FALSE) 
# rgl.postscript("./plots/persp3dd.pgf","pgf", drawText=FALSE) 

rgl.snapshot(filename = './plots/3dplot legend.png', fmt = 'png')



#       Test greater max K in novel framework GAMs ####

longhurst.novel.list.maxk <- lapply(list(nano, foram, radio, diatom),
                               function(x){
                                 
                                 neptune.novelty(dataset = x,
                                                 bin.width = 0.1,
                                                 taxon.res = "Taxon.ID",
                                                 bin.cutoff = 10,
                                                 taxa.cutoff = 10,
                                                 novel.alpha = 0.05,
                                                 novel.metric = "jaccard",
                                                 ssmat.type="pa",
                                                 rich.cutoff=0,
                                                 gam.max.k = 10,
                                                 plot=TRUE,
                                                 longhurst = TRUE)
                                 
                               })
names(longhurst.novel.list.maxk) <- c("nano", "foram", "radio", "diatom")

maxk.novel <- do.call("rbind", lapply(1:length(longhurst.novel.list.maxk), function(n){
temp <- do.call("rbind", longhurst.novel.list.maxk[[n]]$novel)
temp$taxa <- (1:4)[n]
return(temp)
}))

summary(maxk.novel$seq.edf[!duplicated(maxk.novel$site, ":",maxk.novel$taxa)])

reg.novel <- do.call("rbind", lapply(longhurst.novel.list, function(x){
  
  do.call("rbind", x$novel)
}))


#       Trends over time ####

time.trends(novel.list = longhurst.novel.cut,
            max.k = 20)

#       Species accumulation curves ####

# run for each taxa for each longhurst province,
# then model averages over time

sp.accum.df <- do.call("rbind", lapply(1:length(longhurst.novel.list), function(n1){
  
  x <- longhurst.novel.list[[n1]]
  taxa <- c("nano", "foram", "radio", "diatom")[n1]
  
  do.call("rbind", lapply(1:length(x$ssmats), function(n){
    
    ss <- x$ssmats[[n]]
    site <- x$sites[[n]]
  
    real.mat <- ss[order(rowSums(ss)),][nrow(ss):1,]
  
    test.mat <- rbind(rep(0, ncol(real.mat)), apply(real.mat[-nrow(real.mat),], 2, cumsum))
  
    accum <- data.frame(sp.accum = cumsum(rowSums(real.mat > test.mat)),
                        site = site,
                        taxa = taxa)
    accum$bin <- 1:nrow(accum)
    return(accum)
    
  }))
  
}))

# now model average, with site and taxa offets
sc.bin <- scale(sp.accum.df$bin)
sp.accum.df$sc.bin <- as.vector(sc.bin)
sp.accum.df$taxa <- as.factor(sp.accum.df$taxa)

sp.accum.m.base <- glm.nb(sp.accum ~ log(bin), data=sp.accum.df)
sp.accum.m <- glmer.nb(sp.accum ~ log(bin) + (1|taxa) + (1|site), data=sp.accum.df)
sp.accum.taxa <- glmer.nb(sp.accum ~ log(bin)  * taxa + (1|site), data=sp.accum.df)
summary(sp.accum.m)
summary(sp.accum.taxa)

pred.df <- data.frame(bin = 1:max(sp.accum.df$bin))
taxa.pred <- data.frame(bin = rep(pred.df$bin, 4),
                        taxa = rep(levels(sp.accum.df$taxa), each=nrow(pred.df)))
  

pred.df$sc.bin <- (pred.df$bin - attr(sc.bin, "scaled:center")) / attr(sc.bin, "scaled:scale")
pred.df <- cbind(pred.df, fit=predict(sp.accum.m, newdata=pred.df, re.form=NA, se.fit=TRUE))
taxa.pred <- cbind(taxa.pred, fit=predict(sp.accum.taxa, newdata=taxa.pred, re.form=NA, se.fit=TRUE))

pdf("./plots/species accumulation plot.pdf", height=3, width=4, useDingbats = FALSE)
par(mar=c(2,3,0.5,0.5), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL, y=NULL, xlim=c(1,515), ylim=c(7,400), log="y", axes=FALSE)

axis(side=1, mgp=c(3,0.1,0))
axis(side=2, at=c(seq(1,10,1), 
                  seq(10,100,10), 
                  seq(100,1000,100)),
                  labels=NA, tcl=-0.125)
axis(side=2, at=c(1,10,100,200,300))
mtext(side=1, line=1, text="Time bin")
mtext(side=2, line=2, las=0, text="Cumulative sampling bin richness")

lapply(split(sp.accum.df, f=paste0(sp.accum.df$site, sp.accum.df$taxa)), function(x){
  temp.col <- col2rgb(c("darkgreen", "blue", "orange", "red")[x$taxa[1]]) /255
  lines(x$sp.accum, col=rgb(temp.col[1], temp.col[2], temp.col[3], 0.35), lwd=0.5)
})
lapply(1:4, function(n){
  x <- split(taxa.pred, f=taxa.pred$taxa)[[n]]
  lines(exp(x$fit) ~ x$bin, lwd=2, col = c("red", "blue", "darkgreen", "orange")[n])
  text(y=exp(rev(x$fit)[1]), x=rev(x$bin)[1], pos=4,
       labels=c("D", "F", "N", "R")[n], font=2, col = c("red", "blue", "darkgreen", "orange")[n])
  
})
lines(exp(pred.df$fit) ~ pred.df$bin, lwd=4)
box()
dev.off()

# ####