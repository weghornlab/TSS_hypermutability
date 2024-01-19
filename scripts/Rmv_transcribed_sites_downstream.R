options(scipen=999)
library("dplyr")
library("data.table")
library("stringr")
args = commandArgs(trailingOnly=TRUE)
# args[1] = "chr1"
CHR = args[1]

set.seed(8)
     
### Dataset ###
##postTTS = fread(paste0("../results/BEDS/", CHR, "_transcripts_postTTS_FANTOM5.bed"), header = F)
# postTTS = fread(paste0("../results/BEDS/", CHR, "_transcripts_all_postTTS_FANTOM5.bed"), header = F)
## for postTTS rmv sites transcribed

# upstream = fread(paste0("../results/BEDS/", CHR, "_transcripts_upstream_FANTOM5.bed"), header = F)
## for upstream rmv sites transcribed

downstream = fread(paste0("../results/BEDS/", CHR, "_transcripts_downstream_FANTOM5.bed"), header = F)
## for downstream rmv sites > 50000 and sites transcribed by more than one transcript

mult_trans = fread("../results/sites_multiply_transcribed.bed", header = F)

# postTTS_2 = anti_join(postTTS, downstream, by = c("V2", "V3"))
# upstream_2 = anti_join(upstream, downstream, by = c("V2", "V3"))
downstream_2 = anti_join(downstream, mult_trans, by = c("V1", "V2", "V3"))


### downstream filtering ###
transcribed = nrow(downstream)
transcribed_once = nrow(downstream_2)
downstream_2 = subset(downstream_2, V7 <= 50)
transcribed_once_50kb = nrow(downstream_2)
# nrow(inner_join(postTTS_2, upstream_2, by = c("V2", "V3"))) ## there are many sites that are both and many upstreams and postTTS which are repeated (coming from different transcripts)

write.table(downstream_2, file = paste0("../results/BEDS/", CHR, "_transcripts_downstream_nontranscribed_FANTOM5.bed"), quote = F, row.names = F, col.names = F, sep = "\t", append = F)


# ### upstream filtering ###
# upstream_2 = upstream_2[order(upstream_2[,'V3'],upstream_2[,'V9']),] #take the sites closest to the nearest TSS
# upstream_3 = upstream_2[!duplicated(upstream_2$V3),] #there are sites that are at the same distance in different transcripts but I just take a random one

#   ## Sanity check chr22
#   # length((upstream_2$V3))
#   # length(unique(upstream_2$V3))
#   # length((upstream_3$V3))
#   # length(unique(upstream_3$V3))
#   # subset(upstream_2, V3 == "22002040")
#   # subset(upstream_3, V3 == "22002040")

#         ## very long and partial solution
#         # df = plyr::ldply(tapply(upstream_4$V5, upstream_4$V3, list), rbind)
#         # my_cols = colnames(df)[-1]
#         # df = df[my_cols]


# ### postTTS filtering ###
# postTTS_2 = postTTS_2[order(postTTS_2[,'V3'],postTTS_2[,'V9'],postTTS_2[,'V10']),] #take the sites closest to the nearest TTS, if the same then take the ones closest to the nearest TSS (min V10)
# postTTS_3 = postTTS_2[!duplicated(postTTS_2$V3),] #there are sites that are at the same distance in different transcripts but I just take a random one
  
#   ## Sanity check chr22
#   # length((postTTS_2$V3))
#   # length(unique(postTTS_2$V3))
#   # length((postTTS_3$V3))
#   # length(unique(postTTS_3$V3))

#   # length((downstream_2$V3))
#   # length(unique(downstream_2$V3))


# ### Merging files ###
# upstream_3$V10 = NA
# downstream_2$V10 = NA
# postTTS_3$V11 = "postTTS"
# upstream_3$V11 = "upstream"
# downstream_2$V11 = "downstream"

# # nrow(inner_join(postTTS_3, upstream_3, by = c("V2", "V3")))

# upstream_postTTS = rbind(upstream_3, postTTS_3)
# upstream_postTTS = upstream_postTTS[order(upstream_postTTS[,'V3'],upstream_postTTS[,'V9']),] #take the sites closest to the nearest TSS or TTS
# upstream_postTTS_2 = upstream_postTTS[!duplicated(upstream_postTTS$V3),] #there are sites that are at the same distance in different transcripts but I just take a random one

# upstream_postTTS_downstream = rbind(upstream_postTTS_2, downstream_2)
# upstream_postTTS_downstream = upstream_postTTS_downstream[order(upstream_postTTS_downstream[,'V3']),] 
#   # length((upstream_postTTS_downstream$V3))
#   # length(unique(upstream_postTTS_downstream$V3))


# ### Results to print ###
# posTTS_sites = nrow(postTTS)
# posTTS_2_sites = nrow(postTTS_3)
# posTTS_3_sites = nrow(subset(upstream_postTTS_downstream, V11 == "postTTS"))

# upstream_sites = nrow(upstream)
# upstream_2_sites = nrow(upstream_3)
# upstream_3_sites = nrow(subset(upstream_postTTS_downstream, V11 == "upstream"))

# results = c(CHR, transcribed, transcribed_once, transcribed_once_50kb, upstream_sites, upstream_2_sites, upstream_3_sites, posTTS_sites, posTTS_2_sites, posTTS_3_sites)

# write.table(upstream_postTTS_downstream, file = paste0("../results/BEDS/", CHR, "_transcripts_and_co_FANTOM5.bed"), quote = F, row.names = F, col.names = F, sep = "\t", append = F)
# write.table(results, file = paste0("../results/analyzed_sites_FANTOM5.bed"), quote = F, row.names = F, col.names = F, sep = "\t", append = T)

