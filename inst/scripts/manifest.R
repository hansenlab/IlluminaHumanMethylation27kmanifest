library(minfi)

manifest <- "~/Desktop/humanmethylation27_270596_v1-2.csv"
stopifnot(file.exists(manifest))
maniTmp <- minfi:::read.manifest.27k(manifest)

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





