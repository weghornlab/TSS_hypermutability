options(scipen=999)
library("dplyr")
library("data.table")
library("stringr")
library("tidyr")
args = commandArgs(trailingOnly=TRUE)
# args[1] = "chr22"
# args[2] = "CRG36_Exons"
CHR = args[1]

set.seed(8)
     
### Dataset ###
dt = fread(paste0("../results/BEDS/", CHR, "_transcripts_and_co_FANTOM5_", args[2], ".bed"), header = F)

dt = unique(dplyr::select(dt, V1, V11, V5, V7, V8, V9))
colnames(dt) <- c("chr", "position", "geneID", "1000", "100", "10")

dt_interrogated <- unique(gather(dt, bin_length, bin, `1000`:`10`, factor_key=TRUE))
dt_interrogated$bin = as.character(dt_interrogated$bin)
dt_interrogated$bin_length = as.numeric(as.character(dt_interrogated$bin_length))
head(dt_interrogated)
rm(list = "dt")

## downstream
downstream = fread(paste0("../results/BEDS/", CHR, "_transcripts_downstream_FANTOM5.bed"), header = F)
head(downstream)
## colnames(dt) <- c("chr", "start", "end", "transcript_type", "geneID", "strand", "bin1000", "bin100", "bin10")
## downstream_long =  unique(gather(downstream, bin_size, bin_number, bin:bin, factor_key=TRUE))
## head(downstream_long)

		dw_min = data.frame(tapply(downstream$V2, paste0(downstream$V5, "_",downstream$V7), min))
		dw_min$names <- rownames(dw_min)
		dw_min = dw_min %>%
		  separate(names, c("geneID", "bin"), "_")

		dw_max = data.frame(tapply(downstream$V2, paste0(downstream$V5, "_",downstream$V7), max))
		dw_max$names <- rownames(dw_max)
		dw_max = dw_max %>%
		  separate(names, c("geneID", "bin"), "_")

	dw = inner_join(dw_min, dw_max, by = c("geneID", "bin"))
	colnames(dw) <- c("start", "geneID", "bin", "end")
	dw$chr = CHR
	dw$position = "downstream_TSS"
	dw_1000 = dplyr::select(dw, chr, start, end, geneID, bin, position)

		dw_min = data.frame(tapply(downstream$V2, paste0(downstream$V5, "_",downstream$V8), min))
		dw_min$names <- rownames(dw_min)
		dw_min = dw_min %>%
		  separate(names, c("geneID", "bin"), "_")

		dw_max = data.frame(tapply(downstream$V2, paste0(downstream$V5, "_",downstream$V8), max))
		dw_max$names <- rownames(dw_max)
		dw_max = dw_max %>%
		  separate(names, c("geneID", "bin"), "_")

	dw = inner_join(dw_min, dw_max, by = c("geneID", "bin"))
	colnames(dw) <- c("start", "geneID", "bin", "end")
	dw$chr = CHR
	dw$position = "downstream_TSS"
	dw_100 = dplyr::select(dw, chr, start, end, geneID, bin, position)

	# 	dw_min = data.frame(tapply(downstream$V2, paste0(downstream$V5, "_",downstream$V9), min))
	# 	dw_min$names <- rownames(dw_min)
	# 	dw_min = dw_min %>%
	# 	  separate(names, c("geneID", "bin"), "_")

	# 	dw_max = data.frame(tapply(downstream$V2, paste0(downstream$V5, "_",downstream$V9), max))
	# 	dw_max$names <- rownames(dw_max)
	# 	dw_max = dw_max %>%
	# 	  separate(names, c("geneID", "bin"), "_")

	# dw = inner_join(dw_min, dw_max, by = c("geneID", "bin"))
	# colnames(dw) <- c("start", "geneID", "bin", "end")
	# dw$chr = CHR
	# dw$position = "downstream"
	# dw_10 = dplyr::select(dw, chr, start, end, geneID, bin, position)


## upstream TTS
rm(list = c("downstream", "dw", "dw_max", "dw_min"))
upstream = fread(paste0("../results/BEDS/", CHR, "_transcripts_upTTS_FANTOM5.bed"), header = F)
head(upstream)

		up_min = data.frame(tapply(upstream$V2, paste0(upstream$V5, "_",upstream$V7), min))
		up_min$names <- rownames(up_min)
		up_min = up_min %>%
		  separate(names, c("geneID", "bin"), "_")

		up_max = data.frame(tapply(upstream$V2, paste0(upstream$V5, "_",upstream$V7), max))
		up_max$names <- rownames(up_max)
		up_max = up_max %>%
		  separate(names, c("geneID", "bin"), "_")

	up = inner_join(up_min, up_max, by = c("geneID", "bin"))
	colnames(up) <- c("start", "geneID", "bin", "end")
	up$chr = CHR
	up$position = "upstream_TTS"
	upTTS_1000 = dplyr::select(up, chr, start, end, geneID, bin, position)

		up_min = data.frame(tapply(upstream$V2, paste0(upstream$V5, "_",upstream$V8), min))
		up_min$names <- rownames(up_min)
		up_min = up_min %>%
		  separate(names, c("geneID", "bin"), "_")

		up_max = data.frame(tapply(upstream$V2, paste0(upstream$V5, "_",upstream$V8), max))
		up_max$names <- rownames(up_max)
		up_max = up_max %>%
		  separate(names, c("geneID", "bin"), "_")

	up = inner_join(up_min, up_max, by = c("geneID", "bin"))
	colnames(up) <- c("start", "geneID", "bin", "end")
	up$chr = CHR
	up$position = "upstream_TTS"
	upTTS_100 = dplyr::select(up, chr, start, end, geneID, bin, position)

	# 	up_min = data.frame(tapply(upstream$V2, paste0(upstream$V5, "_",upstream$V9), min))
	# 	up_min$names <- rownames(up_min)
	# 	up_min = up_min %>%
	# 	  separate(names, c("geneID", "bin"), "_")

	# 	up_max = data.frame(tapply(upstream$V2, paste0(upstream$V5, "_",upstream$V9), max))
	# 	up_max$names <- rownames(up_max)
	# 	up_max = up_max %>%
	# 	  separate(names, c("geneID", "bin"), "_")

	# up = inner_join(up_min, up_max, by = c("geneID", "bin"))
	# colnames(up) <- c("start", "geneID", "bin", "end")
	# up$chr = CHR
	# up$position = "upstream"
	# up_10 = dplyr::select(up, chr, start, end, geneID, bin, position)


## upstream
rm(list = c("upstream", "up", "up_max", "up_min"))
upstream = fread(paste0("../results/BEDS/", CHR, "_transcripts_upstream_FANTOM5.bed"), header = F)
head(upstream)

		up_min = data.frame(tapply(upstream$V2, paste0(upstream$V5, "_",upstream$V7), min))
		up_min$names <- rownames(up_min)
		up_min = up_min %>%
		  separate(names, c("geneID", "bin"), "_")

		up_max = data.frame(tapply(upstream$V2, paste0(upstream$V5, "_",upstream$V7), max))
		up_max$names <- rownames(up_max)
		up_max = up_max %>%
		  separate(names, c("geneID", "bin"), "_")

	up = inner_join(up_min, up_max, by = c("geneID", "bin"))
	colnames(up) <- c("start", "geneID", "bin", "end")
	up$chr = CHR
	up$position = "upstream_TSS"
	up_1000 = dplyr::select(up, chr, start, end, geneID, bin, position)

		up_min = data.frame(tapply(upstream$V2, paste0(upstream$V5, "_",upstream$V8), min))
		up_min$names <- rownames(up_min)
		up_min = up_min %>%
		  separate(names, c("geneID", "bin"), "_")

		up_max = data.frame(tapply(upstream$V2, paste0(upstream$V5, "_",upstream$V8), max))
		up_max$names <- rownames(up_max)
		up_max = up_max %>%
		  separate(names, c("geneID", "bin"), "_")

	up = inner_join(up_min, up_max, by = c("geneID", "bin"))
	colnames(up) <- c("start", "geneID", "bin", "end")
	up$chr = CHR
	up$position = "upstream_TSS"
	up_100 = dplyr::select(up, chr, start, end, geneID, bin, position)

	# 	up_min = data.frame(tapply(upstream$V2, paste0(upstream$V5, "_",upstream$V9), min))
	# 	up_min$names <- rownames(up_min)
	# 	up_min = up_min %>%
	# 	  separate(names, c("geneID", "bin"), "_")

	# 	up_max = data.frame(tapply(upstream$V2, paste0(upstream$V5, "_",upstream$V9), max))
	# 	up_max$names <- rownames(up_max)
	# 	up_max = up_max %>%
	# 	  separate(names, c("geneID", "bin"), "_")

	# up = inner_join(up_min, up_max, by = c("geneID", "bin"))
	# colnames(up) <- c("start", "geneID", "bin", "end")
	# up$chr = CHR
	# up$position = "upstream"
	# up_10 = dplyr::select(up, chr, start, end, geneID, bin, position)


## postTTS
rm(list = c("upstream", "up", "up_max", "up_min"))
postTTS_file = fread(paste0("../results/BEDS/", CHR, "_transcripts_all_postTTS_FANTOM5.bed"), header = F)
head(postTTS_file)

		postTTS_min = data.frame(tapply(postTTS_file$V2, paste0(postTTS_file$V5, "_",postTTS_file$V7), min))
		postTTS_min$names <- rownames(postTTS_min)
		postTTS_min = postTTS_min %>%
		  separate(names, c("geneID", "bin"), "_")

		postTTS_max = data.frame(tapply(postTTS_file$V2, paste0(postTTS_file$V5, "_",postTTS_file$V7), max))
		postTTS_max$names <- rownames(postTTS_max)
		postTTS_max = postTTS_max %>%
		  separate(names, c("geneID", "bin"), "_")

	postTTS = inner_join(postTTS_min, postTTS_max, by = c("geneID", "bin"))
	colnames(postTTS) <- c("start", "geneID", "bin", "end")
	postTTS$chr = CHR
	postTTS$position = "downstream_TTS"
	postTTS_1000 = dplyr::select(postTTS, chr, start, end, geneID, bin, position)


		postTTS_min = data.frame(tapply(postTTS_file$V2, paste0(postTTS_file$V5, "_",postTTS_file$V8), min))
		postTTS_min$names <- rownames(postTTS_min)
		postTTS_min = postTTS_min %>%
		  separate(names, c("geneID", "bin"), "_")

		postTTS_max = data.frame(tapply(postTTS_file$V2, paste0(postTTS_file$V5, "_",postTTS_file$V8), max))
		postTTS_max$names <- rownames(postTTS_max)
		postTTS_max = postTTS_max %>%
		  separate(names, c("geneID", "bin"), "_")

	postTTS = inner_join(postTTS_min, postTTS_max, by = c("geneID", "bin"))
	colnames(postTTS) <- c("start", "geneID", "bin", "end")
	postTTS$chr = CHR
	postTTS$position = "downstream_TTS"
	postTTS_100 = dplyr::select(postTTS, chr, start, end, geneID, bin, position)


	# 	postTTS_min = data.frame(tapply(postTTS_file$V2, paste0(postTTS_file$V5, "_",postTTS_file$V9), min))
	# 	postTTS_min$names <- rownames(postTTS_min)
	# 	postTTS_min = postTTS_min %>%
	# 	  separate(names, c("geneID", "bin"), "_")

	# 	postTTS_max = data.frame(tapply(postTTS_file$V2, paste0(postTTS_file$V5, "_",postTTS_file$V9), max))
	# 	postTTS_max$names <- rownames(postTTS_max)
	# 	postTTS_max = postTTS_max %>%
	# 	  separate(names, c("geneID", "bin"), "_")

	# postTTS = inner_join(postTTS_min, postTTS_max, by = c("geneID", "bin"))
	# colnames(postTTS) <- c("start", "geneID", "bin", "end")
	# postTTS$chr = CHR
	# postTTS$position = "postTTS"
	# postTTS_10 = dplyr::select(postTTS, chr, start, end, geneID, bin, position)

rm(list = c("postTTS_file", "postTTS", "postTTS_max", "postTTS_min"))


### Merging ###
dw_1000$bin_length = 1000
dw_100$bin_length = 100
# dw_10$bin_length = 10
upTTS_1000$bin_length = 1000
upTTS_100$bin_length = 100
# up_10$bin_length = 10
up_1000$bin_length = 1000
up_100$bin_length = 100
# up_10$bin_length = 10
postTTS_1000$bin_length = 1000
postTTS_100$bin_length = 100
# postTTS_10$bin_length = 10

# dt_all = rbind(dw_1000, dw_100, dw_10, up_1000, up_100, up_10, postTTS_1000, postTTS_100, postTTS_10)
dt_all = rbind(dw_1000, dw_100, upTTS_1000, upTTS_100, up_1000, up_100, postTTS_1000, postTTS_100)

head(dt_all)
head(dt_interrogated)

dt_interrogated_coord = inner_join(dt_all, dt_interrogated, by = c("chr", "position", "geneID", "bin", "bin_length"))

## write.table(dt_interrogated_coord, file = paste0("../results/Bin_coordinates_FANTOM5.bed"), quote = F, row.names = F, col.names = F, sep = "\t", append = T)

dt_interrogated_coord1000 = subset(dt_interrogated_coord, bin_length == "1000" & start > 0)
write.table(dt_interrogated_coord1000, file = paste0("../results/Bin_coordinates_FANTOM5_1000_", args[2],"_last.bed"), quote = F, row.names = F, col.names = F, sep = "\t", append = T)

dt_interrogated_coord100 = subset(dt_interrogated_coord, bin_length == "100" & start > 0)
write.table(dt_interrogated_coord100, file = paste0("../results/Bin_coordinates_FANTOM5_100_", args[2],"_last.bed"), quote = F, row.names = F, col.names = F, sep = "\t", append = T)

# dt_interrogated_coord10 = subset(dt_interrogated_coord, bin_length == "10" & start > 0)
# write.table(dt_interrogated_coord10, file = paste0("../results/Bin_coordinates_FANTOM5_10_", args[2],".bed"), quote = F, row.names = F, col.names = F, sep = "\t", append = T)
