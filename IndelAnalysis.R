library(ggplot2)
library(reshape2)

setwd("/Users/li/Documents/Research/Kharchenko_Lab/Miseq_data/20171006/analysis/eGFP")
## pattern = paste0("*",tail(unlist(strsplit(getwd(),"/")),1),"*",".csv")
pattern = tail(unlist(strsplit(getwd(),"/")),1)
csvFileNames = list.files(pattern = "*.csv")
##csvFileNames = list.files(pattern = "*.needleall.sam.csv")
sampleNames = c("WT_AGAT", "RR_AGAT", "RVR_AGAT", "WT_TTTA", "RR_TCCA", "RVR_TATT", "RR_ACCC", 
                "RVR_TACA", "RVR_AATA", "Control_Control")
csvFileList = list()
## Read each csv file into a list; store the sample name into a vector
for (i in 1:length(csvFileNames)){
  csvFileList[[i]] = read.csv(csvFileNames[i],stringsAsFactors = F)
}
names(csvFileList) = sampleNames

## Check out top 15 most frequent indels 
for (i in 1:length(csvFileList)){
  print(names(csvFileList)[i])
  print(head(sort(table(csvFileList[[i]][1]),decreasing=TRUE),15))
}

## Top N most frequent indel size/position 
indelTop = list(numeric(),numeric(),numeric())
indelTopFreq = list(numeric(),numeric(),numeric())
topN = 5

for (i in csvFileList){
  for (j in c(3:5)){
    temp = sort(table(i[,j]),decreasing=T)
    tempTop = names(temp)[1:topN]
    tempTopFreq = temp[tempTop]/sum(temp)
    indelTop[[j-2]] = c(indelTop[[j-2]], tempTop)
    indelTopFreq[[j-2]] = c(indelTopFreq[[j-2]], unname(tempTopFreq))
  }
} 

sampleNameDf = rep(names(csvFileList),rep(topN,length(csvFileList)))
indelTopDf = data.frame(sampleNameDf, indelTop[[1]], indelTopFreq[[1]], 
                        indelTop[[2]], indelTopFreq[[2]], 
                        indelTop[[3]], indelTopFreq[[3]])
colnames(indelTopDf) = c("sample", "indel_size", "indel_size_freq", 
                         "indel_position", "indel_position_freq", 
                         "indel_size_in_guide",  "indel_size_in_guide_freq")

## Max range of indel size/position for xlim of histogram
maxRange = c(0,0,0)
xlim=list()
for (i in csvFileList){
  for (j in c(3:5)){
    currentRange = max(i[,j]) - min(i[,j])
    if (currentRange > maxRange[j-2]){
      maxRange[j-2] = currentRange
      xlim[[j-2]] = c(min(i[,j]),max(i[,j]))
    }
  }
}

## Generate output path
outputPath = paste0(getwd(),"/","figures")
if (!file.exists(outputPath)){
  dir.create(outputPath)
}
if (pattern == "Cpf1"){
  title = "Guide Template"
}else{
  title = "eGFP Target"
}

## Indel frequency
metaFileNames = list.files(pattern = "*.meta")
metaFileList = list()
metaCombined = NULL
# Read each meta file into a list and then combine into a single data frame
for (i in 1:length(csvFileList)){
  metaFileList[[i]] = read.csv(metaFileNames[i],stringsAsFactors = F)
  metaCombined = rbind(metaCombined,metaFileList[[i]][1,])
}
colnames(metaCombined) = c("sumReads","sumMutations","sumDeletions","sumInsertions", "mutationFreq", "uniqueIndel")
metaCombined$sample = factor(names(csvFileList), level = names(csvFileList))
metaCombined$deletionFreq = metaCombined$sumDeletions/metaCombined$sumReads
metaCombined$insertionFreq = metaCombined$sumInsertions/metaCombined$sumReads
# Plot
delIns = metaCombined[,c("deletionFreq", "insertionFreq", "sample")]
colnames(delIns) = c("Deletion", "Insertion", "Sample")
delIns = melt(delIns)
if (grepl("Cpf1",pattern)){
  position = c(0.8,0.1)
}else{
  position = c(0.8,0.6)
}
filename = paste0(outputPath, "/indel_freq.tiff")
tiff(filename, width=1600, height=1200, res=300, units ="px", pointsize=8, compression="lzw")
ggplot(data=delIns,aes(x=Sample,y=value, fill=variable)) + geom_col() + coord_flip() + ggtitle(title) + ylim(0,0.4) +
        ylab ("Indel Frequency") + xlab("Constructs") + theme(plot.title = element_text(hjust = 0.5,face="bold",
        size=20), axis.text=element_text(size=10), axis.title=element_text(size=15,face="bold"),
        legend.position = position) + guides(fill=guide_legend(title="Mutation type"))
dev.off()

## num of unique indels 
filename = paste0(outputPath, "/uniqueIndel.tiff")
tiff(filename, width=1600, height=1200, res=300, units ="px", pointsize=8, compression="lzw")
ggplot(data=metaCombined,aes(x=sample,y=uniqueIndel)) + geom_col() + coord_flip() + ggtitle(title) + ylim(0,1000) +
        ylab ("Number of unique indels") + xlab("Constructs") + theme(plot.title = element_text(hjust = 0.5,
        face="bold", size=20), axis.text=element_text(size=10), axis.title=element_text(size=15,face="bold"))
dev.off()

## Concatenate all datasets
dataCombined = NULL
for (i in 1:length(csvFileList)){
  temp = csvFileList[[i]]
  temp$sample = names(csvFileList)[i]
  dataCombined = rbind(dataCombined,temp)
}
dataCombined$sample = factor(dataCombined$sample,level=names(csvFileList))

## Histogram for indel size
filename = paste0(outputPath,"/indel_size_dist.tiff")
tiff(filename, width=1600, height=1200, res=300, units ="px", pointsize=8, compression="lzw")
ggplot(data=dataCombined, aes(size)) + geom_histogram() + facet_wrap(~sample) + ggtitle(title) + 
        theme(plot.title = element_text(hjust = 0.5,face="bold", size=20), axis.text=element_text(size=6), 
              axis.title=element_text(size=15,face="bold")) + xlim(range(dataCombined$size)) + 
              xlab("Size of indels (bp)")
dev.off()

## Histogram for indel position
filename = paste0(outputPath,"/","indel_pos_dist.tiff")
tiff(filename, width=1600, height=1200, res=300, units ="px", pointsize=8, compression="lzw")
ggplot(data=dataCombined, aes(position)) + geom_histogram() + facet_wrap(~sample) + ggtitle(title) +
  theme(plot.title = element_text(hjust = 0.5,face="bold", size=20), axis.text=element_text(size=6), 
        axis.title=element_text(size=15,face="bold")) + xlim(range(dataCombined$position)) + 
        xlab("Position of indels relative to 3' of guide (bp)")
dev.off()

## Histogram for deletions in guide sequence
filename = paste0(outputPath,"/indel_size_guide_dist.tiff")
tiff(filename, width=1600, height=1200, res=300, units ="px", pointsize=8, compression="lzw")
ggplot(data=dataCombined, aes(delGuide)) + geom_histogram() + facet_wrap(~sample) + ggtitle(title) + 
  theme(plot.title = element_text(hjust = 0.5,face="bold", size=20), axis.text=element_text(size=6), 
        axis.title=element_text(size=15,face="bold")) + xlim(range(dataCombined$delGuide)) + 
        xlab("Deletions in guide sequence (bp)")
dev.off()

