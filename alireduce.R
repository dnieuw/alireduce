#!/usr/bin/env Rscript
cat("Loading R libraries\n\n")

suppressMessages(require(lubridate))
suppressMessages(require(stringr))
suppressMessages(require(data.table))
suppressMessages(require(ape))
suppressMessages(require(Biostrings))
suppressMessages(require(optparse))

option_list <- list(
  make_option(c("-c", "--clusterby"),
              type="character",
              default="week", 
              help = "Timescale to cluster by (week/day/month) [default \"%default\"]"),
  make_option(c("-n", "--nt_diff"), 
              type="integer", 
              default=5,
              metavar="number",
              help="Number of differences to set as a threshold for sequence clustering [default %default]"),
  make_option(c("-s", "--start"), 
              type="character",
              default="2019-01-01",
              help="Start date to include sequences (YYYY-mm-dd) [default %default]"),
  make_option(c("-e", "--end"), 
              type="character", 
              default=today(),
              help="End date to include sequences (YYYY-mm-dd) [default %default]"),
  make_option(c("-f", "--forceall"), 
              type="integer", 
              default=0,
              metavar="number",
              help="Set to 1 to force all sequences without correct 
              date to be included at the end of the alignment [default %default]"),
  make_option(c("-r", "--randomsample"), 
              type="integer", 
              default=0,
              metavar="number",
              help="Set to anything more than 0 to to perform random 
              sampling to that number of sequences for clustering 
              (only use for huge alignments) [default %default]"),
  make_option(c("--usecountry"), 
              type="integer", 
              default=0,
              metavar="number",
              help="Set to 1 to perform random sampling while taking 
              into account country of origin. Make sure that your sequence 
              header is formatted as such: hCoV-19/Netherlands/dog-un-EMC-1/2020|EPI_ISL_722296|2020-10-20|Europe [default %default]"),
  make_option(c("-i", "--input"),
              type="character",
              help = "Input alignment"),
  make_option(c("-o", "--output"),
              type="character",
              help = "Output alignment")
)

opt <- parse_args(OptionParser(option_list=option_list))

clusterby <- opt$clusterby
nt_diff <- opt$nt_diff
date_start <- ymd(opt$start)
date_end <- ymd(opt$end)
force_all <- opt$forceall
randomsample <- opt$randomsample
usecountry <- ifelse(opt$usecountry==1,TRUE,FALSE)
infile <- opt$input
outfile <- opt$output

cat("Loading sequence alignment\n\n")
aln <- readDNAStringSet(infile, "fasta")

metadata <- data.table(seq_id=names(aln), sample_date=str_extract(names(aln),"[0-9]+-[0-9]+-[0-9]+"))
metadata[!is.na(dmy(sample_date, quiet = T)),sample_date:=as.character(dmy(sample_date, quiet = T))]
metadata[,sample_date:=ymd(sample_date, quiet = T)]

#If any records do not have a correct date:
if (sum(is.na(metadata$sample_date)>0)) {
  if (force_all>0){
    #Indicate the number of records without date
    cat("[", paste(sum(is.na(metadata$sample_date)),"] records do not have a YYYY-mm-dd or dd-mm-YYYY formatted date and are added without clustering\n\n"))
    unclustered <- metadata[is.na(sample_date)]
  } else {
    #Indicate the number of removed records
    cat("[", paste(sum(is.na(metadata$sample_date)),"] records are removed because they do not have a YYYY-mm-dd or dd-mm-YYYY formatted date in their seq ID\n\n"))
    unclustered <- NULL
  }
} else {
  unclustered <- NULL
}

#Remove records without date
metadata <- metadata[!is.na(sample_date)]
metadata <- metadata[sample_date>=date_start&sample_date<=date_end]

#Determine how to cluster
if (clusterby=="day") {
  metadata[,dategroup:=sample_date]
} else if (clusterby=="week") {
  metadata[,dategroup:=paste(lubridate::week(sample_date),lubridate::year(sample_date), sep="-")]
}else if (clusterby=="month") {
  metadata[,dategroup:=paste(lubridate::month(sample_date),lubridate::year(sample_date), sep="-")]
}

cat("Clustering [", nrow(metadata), "] sequences from [", as.character(date_start), "] to [", as.character(date_end), "] by [", clusterby, "] using a [", nt_diff, "] nucleotide difference threshold:\n\n")

if (randomsample > 0) {
  cat("Because randomsample is [",randomsample,"] we downsample each considered cluster to [",randomsample,"] sequences, only do this when you are clustering a huge alignement!\n\n")
}

if (usecountry){
  #Split by /
  splitname <- metadata[,str_split(seq_id,pattern = "/", simplify = T)]
  #Check if name is in format x/country/y/z
  if (ncol(splitname) != 4){
    cat("Too many columns after splitting on /, make sure all fasta headers are in format x/country/y/z. Results may be incorrect.\n\n")
  }
  metadata[,country:=str_split(seq_id,pattern = "/", simplify = T)[,2]]
  cat("Found countries: [",paste(metadata[,.N,country][,country], collapse=","),"] to take into account when downsampling\n\n")
}

reduce_daily_sequences <- function(daily, nt_threshold) {
  
  #Return if only 1 sequence was produced on that day
  if (nrow(daily)==1) {
    progress <<- progress + 1
    cat('\r',paste(progress,total, sep = "/"))
    flush.console() 
    return(daily)
  }
  
  #Transform sequence selection to DNAbin class for further analysis
  daily_subset <- as.DNAbin(as.matrix(aln[daily$seq_id,]))
  
  #extract variable sites from the sequence
  variable_sites <- which(sapply(1:ncol(daily_subset), function(x) !any(base.freq(daily_subset[,x])>0.99)))
  
  #Return if all sequences on that day have no variable sites
  if (length(variable_sites)==0) {
    progress <<- progress + 1
    cat('\r',paste(progress,total, sep = "/"))
    flush.console() 
    return(daily)
  }
  
  if (randomsample > 0){
    set.seed(1337)
    random_set <- sample(nrow(daily_subset), min(randomsample,nrow(daily_subset)))
    daily_subset <- daily_subset[random_set,]
    daily <- daily[random_set,]
  }
  
  #Calculate all vs all nucleotide distance based on variable sites
  variable_dist <- dist.dna(daily_subset[,variable_sites], model = "N", pairwise.deletion = T)
  
  #Cluster based on nucleotide difference
  seq_clusters <- cutree(hclust(variable_dist),h = nt_threshold)
  
  daily[,cluster:=seq_clusters]
  
  return_most_common <- function(seq_id) {
    seq_dist <- as.matrix(variable_dist)[seq_id,,drop=F]
    sum_seq_dist <- rowSums(seq_dist)
    most_common <- order(sum_seq_dist)[1]
    return(most_common)
  }
  
  #Find most common (least total mutation difference) sequence per cluster
  daily <- daily[,.SD[return_most_common(seq_id)],cluster]
  daily <- daily[,colnames(metadata), with=F]
  
  progress <<- progress + 1
  cat('\r',paste(progress,total, sep = "/"))
  flush.console() 
  
  return(daily)
}

if (usecountry){
  metadata[,metadata_cluster:=paste(country,dategroup, sep = "/")]
} else {
  metadata[,metadata_cluster:=dategroup]
}

uniq_clusters <- metadata[,.N,metadata_cluster][,metadata_cluster]

progress <- 0
total <- length(uniq_clusters)

result <- rbindlist(lapply(uniq_clusters, function(x) {
  subset <- metadata[metadata_cluster==x]
  reduce_daily_sequences(subset, nt_diff)
}))

cat("\n\nClustering resulted in [", nrow(result), "] sequences\n\n")

if (is.null(unclustered)) {
  filtered_seq_id <- result[,seq_id]
} else {
  cat("\n\nAdded [", nrow(unclustered), "] unclustered sequences\n\n")
  filtered_seq_id <- c(result[,seq_id], unclustered[,seq_id])
}

writeXStringSet(aln[filtered_seq_id],filepath = outfile)

cat("DONE!\n")
