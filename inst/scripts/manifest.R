library(minfi)

manifest <- "../../../IlluminaHumanMethylation27k_files/data/HumanMethylation27_270596_v.1.2.bpm"
stopifnot(file.exists(manifest))
maniTmp <- minfi:::read.manifest.27k(manifest)

## Manifest package
maniList <- maniTmp$manifestList

# Adding colors to Type I SNP probes as discovered experimentally by Tim Triche:
maniList$TypeSnpI[match(
        c("rs1019916", "rs10457834", "rs1416770", "rs1941955", "rs2125573", 
        "rs2235751", "rs2521373", "rs264581", "rs2804694", "rs2959823", 
        "rs5931272", "rs6546473", "rs739259", "rs798149", "rs845016", 
        "rs866884"), maniList$TypeSnpI$Name), "Color"] <- 
        c("Red", "Red", "Red", "Red", "Red", "Grn", "Grn", "Red", "Red", 
		    "Red", "Red", "Grn", "Red", "Grn", "Grn", "Red")

IlluminaHumanMethylation27kmanifest <- IlluminaMethylationManifest(TypeI = maniList$TypeI,
                                                                    TypeII = maniList$TypeII,
                                                                    TypeControl = maniList$TypeControl,
                                                                    TypeSnpI = maniList$TypeSnpI,
                                                                    TypeSnpII = maniList$TypeSnpII,
                                                                    annotation = "IlluminaHumanMethylation27k")
stopifnot(validObject(IlluminaHumanMethylation27kmanifest))
save(IlluminaHumanMethylation27kmanifest, compress = "xz",
     file = "../../data/IlluminaHumanMethylation27kmanifest.rda")
