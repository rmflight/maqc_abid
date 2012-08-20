# investigate what I should be using for unique identifiers in Affy and Agilent

affyDir <- "GSE5350_MAQC_AFX_123456_120CELs"
nRep <- 5
pDat <- data.frame(Sample=c(rep("A",nRep), rep("B", nRep)))
fileNames <- paste(rep("AFX_1",nRep*2), c(rep("_A", nRep), rep("_B", nRep)), seq(1, nRep), ".CEL", sep="" )
affyCel <- ReadAffy(filenames=fileNames[1:2], celfile.path=affyDir)
affyDat <- exprs(rma(affyCel))

affyMap <- read.table("GPL570-13270.txt", header=T, skip=16, sep="\t", quote="", strip.white=T, colClasses="character", stringsAsFactors=F)

sum(duplicated(affyMap$ID))
sum(rownames(affyDat) %in% affyMap$ID)
nrow(affyMap)

agilDir <- "GSE5350_MAQC_AG1_123_60TXTs"
nRep <- 5
fileNames <- paste(rep("AG1_1",nRep*2), c(rep("_A", nRep), rep("_B", nRep)), seq(1, nRep), ".txt", sep="" )
agilDat <- read.maimages(files=fileNames[1:2], path=agilDir, green.only=T, source="agilent")

agilMap <- read.table("GPL1708-20418.txt", header=T, skip=20, sep="\t", quote="", strip.white=T, colClasses="character", stringsAsFactors=F, comment.char="")

nrow(agilDat$E)

sum(duplicated(agilDat$genes$ProbeName)) # this is no good

agilOrg <- read.table(file.path(agilDir, fileNames[1]), header=T, skip=9, sep="\t", quote="", strip.white=T, stringsAsFactors=F, comment.char="")
max(agilOrg$FeatureNum)

agilOrg$FeatureNum <- as.character(agilOrg$FeatureNum)
agilDat$genes$FeatureNum <- as.character(agilDat$genes$Row * agilDat$genes$Col)


ranFeat <- sample(agilDat$genes$FeatureNum, 20)

agilDat$genes$ProbeName[(agilDat$genes$FeatureNum %in% ranFeat)]