cat("Retrieving file list\n")
data.file.dir = "../rawdata/BRCA/Methylation/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/"
normal.file.names <- list.files(data.file.dir, pattern=".*11A-[A-Z0-9]{3}-[A-Z0-9]{4}-[0-9]{2}.*")
cancer.file.names <- list.files(data.file.dir, pattern=".*01A-[A-Z0-9]{3}-[A-Z0-9]{4}-[0-9]{2}.*")

normal.barcodes <- sub(".*TCGA", "TCGA", sub("-11A-[A-Z0-9]{3}-[A-Z0-9]{4}-[0-9]{2}.*", "", normal.file.names))
cancer.barcodes <- sub(".*TCGA", "TCGA", sub("-01A-[A-Z0-9]{3}-[A-Z0-9]{4}-[0-9]{2}.*", "", cancer.file.names))

cancer.barcodes <- cancer.barcodes[!duplicated(cancer.barcodes)]

matched.barcodes   <- cancer.barcodes[  cancer.barcodes %in% normal.barcodes ]
unmatched.barcodes <- cancer.barcodes[!(cancer.barcodes %in% normal.barcodes)]

matched.normal.file.names <- normal.file.names[sub(".*TCGA", "TCGA",
                                              sub("-11A-[A-Z0-9]{3}-[A-Z0-9]{4}-[0-9]{2}.*", "",
                                              normal.file.names)) %in% matched.barcodes]
matched.cancer.file.names <- cancer.file.names[sub(".*TCGA", "TCGA",
                                                      sub("-01A-[A-Z0-9]{3}-[A-Z0-9]{4}-[0-9]{2}.*", "",
                                                      cancer.file.names)) %in% matched.barcodes]
unmatched.cancer.file.names <- cancer.file.names[sub(".*TCGA", "TCGA",
                                                        sub("-01A-[A-Z0-9]{3}-[A-Z0-9]{4}-[0-9]{2}.*", "",
                                                        cancer.file.names)) %in% unmatched.barcodes]

save(matched.normal.file.names,   file="../Rdata/BRCA-new/methylation/file-lists/matched-normal-file-names.Rdata")
save(matched.cancer.file.names,   file="../Rdata/BRCA-new/methylation/file-lists/matched-cancer-file-names.Rdata")
save(unmatched.cancer.file.names, file="../Rdata/BRCA-new/methylation/file-lists/unmatched-cancer-file-names.Rdata")
