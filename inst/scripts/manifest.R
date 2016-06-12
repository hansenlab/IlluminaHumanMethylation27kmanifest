library(minfi)

manifest <- "../../../files_27k/humanmethylation27_270596_v1-2.csv"
stopifnot(file.exists(manifest))


read.manifest.27k <- function(file) {
	temp <- read.csv(file)
	control.line <- which(temp[,1]=="[Controls]")
	assay.line   <- which(temp[,1]=="[Assay]")
	colNames <- temp[assay.line+1,]
	colClasses <- rep("character", length(colNames))

	manifest <- read.csv(file, header = TRUE,
	                           sep = ",", comment.char = "", quote = "", 
	                           skip = assay.line+1, colClasses = colClasses,
	                           nrows = control.line - (assay.line + 1), fill = TRUE)

	TypeI <- manifest[ c("Name", "AddressA_ID", "AddressB_ID", "Color_Channel",
	                         "Next_Base", "AlleleA_ProbeSeq", "AlleleB_ProbeSeq")]
	TypeI <- TypeI[TypeI$Name != "",]
	names(TypeI)[c(2,3,4,5,6,7)] <- c("AddressA", "AddressB", "Color", "NextBase",
	                                  "ProbeSeqA", "ProbeSeqB")
	TypeI <- as(TypeI, "DataFrame")
	TypeI$ProbeSeqA <- DNAStringSet(TypeI$ProbeSeqA)
	TypeI$ProbeSeqB <- DNAStringSet(TypeI$ProbeSeqB)
	TypeI$NextBase <- DNAStringSet(TypeI$NextBase)
	TypeI$nCpG <- as.integer(oligonucleotideFrequency(TypeI$ProbeSeqB,
	                                                  width = 2)[, "CG"] - 1)
	TypeI$nCpG[TypeI$nCpG < 0] <- 0L

	TypeII <- manifest[manifest$Infinium_Design_Type == "II",
	                    c("Name", "AddressA_ID", "AlleleA_ProbeSeq")]
	names(TypeII)[c(2,3)] <- c("AddressA", "ProbeSeqA")
	TypeII <- as(TypeII, "DataFrame")
	TypeII$ProbeSeqA <- BStringSet(TypeII$ProbeSeqA)
	TypeII$nCpG <- as.integer(letterFrequency(TypeII$ProbeSeq, letters = "R"))

	controls <- read.table(file, skip = control.line+1,
	                           sep = ",", comment.char = "", quote = "",
	                           colClasses = c(rep("character", 4)))

	TypeControl <- controls[, 1:4]
	names(TypeControl) <- c("Address", "Type", "Color", "ExtendedType")
	TypeControl <- as(TypeControl, "DataFrame")

	snps <- TypeControl[TypeControl$Type == "Genotyping",]
	TypeControl <- TypeControl[TypeControl$Type != "Genotyping",]
	rsname <- sub("_[AB]", "", snps$ExtendedType)
	snps.sp <- split(snps, rsname)
	snps.sp <- lapply(names(snps.sp), function(rs) {
	    snp <- snps.sp[[rs]]
	    DataFrame(Name = rs,
	              AddressA = snp[grep("_A", snp$ExtendedType), "Address"],
	              AddressB = snp[grep("_B", snp$ExtendedType), "Address"],
	              Color = "Unknown")
	})
	TypeSnpI <- do.call(rbind, snps.sp)
	TypeSnpII <- TypeSnpI[0,]
	     

	list(manifestList = list(TypeI = TypeI, TypeII = TypeII, TypeControl = TypeControl,
	TypeSnpI = TypeSnpI, TypeSnpII = TypeSnpII),
	manifest = manifest, controls = controls)  
}




maniTmp <- read.manifest.27k(manifest)

## Manifest package
maniList <- maniTmp$manifestList
IlluminaHumanMethylation27kmanifest <- IlluminaMethylationManifest(TypeI = maniList$TypeI,
                                                                    TypeII = maniList$TypeII,
                                                                    TypeControl = maniList$TypeControl,
                                                                    TypeSnpI = maniList$TypeSnpI,
                                                                    TypeSnpII = maniList$TypeSnpII,
                                                                    annotation = "IlluminaHumanMethylation27k")
stopifnot(validObject(IlluminaHumanMethylation27kmanifest))
save(IlluminaHumanMethylation27kmanifest, compress = "xz",
     file = "../../data/IlluminaHumanMethylation27kmanifest.rda")






