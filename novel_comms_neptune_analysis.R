# ################################# ####
# Predicting Novel Communities      ####
# Author:    Timothy L Staples      ####
# Collaborators: John Pandolfi      #### 
#                Wolfgang Kiessling ####
# ################################# ####
# Global attributes & working directories ####

rm(list=ls())

# set working directory to folder where this file is located
setwd("/home/timothy/Dropbox/Tim/Post-doc/Research projects/Upload_projects/novelty-cenozoic-microplankton")

# functions ####

# source functions from 'functions' sub-folder
sapply(list.files("./functions", pattern="\\.R", full.names=TRUE), source)

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
#### Analysis ####
#                     import processed novelty data ####

longhurst.novel.cut <- readRDS("./outputs/processedNovelData.rds")

#                     probability of novelty ####

# Probability of novelty occurrence function
novel.prob.list <- novel.probability(longhurst.novel.cut)

figure1.plot.concept(example.ssmat = longhurst.novel.cut[[3]]$ssmats[[length(longhurst.novel.cut[[3]]$novel)-21]],
                    plot.name = "longhurst_complete_3col_concept",
                    example.xlims = rev(c(0,10)),
                    novel.lab.pos = 4,
                    ts.cut = 15)

figure1B.plot(prob.model.list = novel.prob.list,
              plot.name = "longhurst_1B")

#                     observed vs expected transition probabilities ####

# Comparing observed and expected probability of transitions
obs.exp.df <- estimate.observed.expected(prob.model.list = novel.prob.list,
                                         novel.list = longhurst.novel.cut)

obs.exp.df$obs.Estimate.back <- plogis(obs.exp.df$obs.Estimate)
obs.exp.df$obs.lower <- plogis(obs.exp.df$obs.Estimate - 1.96 * obs.exp.df$`obs.Std. Error`)
obs.exp.df$obs.upper <- plogis(obs.exp.df$obs.Estimate + 1.96 * obs.exp.df$`obs.Std. Error`)

write.csv(obs.exp.df, date.wrap(paste0("./outputs/", " obs vs exp transition"), ".csv"))

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

saveRDS(list(data = longhurst.novel.cut,
             novel.prob = novel.prob.list,
             novel.trans = obs.exp.df,
             demo = demographic.df),
        paste0("./outputs/longhurst_analyses1.rds"))


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
