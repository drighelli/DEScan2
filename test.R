library(devtools)
install_github(repo = "jnkoberstein/DEScan", 
               auth_token = "293d0010ce53d248175cee2d67ff248395f1a0a2")
library(DEScan)
library(RUVSeq)
library(EDASeq)

bed.files <- system.file("extdata", 
                         "P43615_Sample_FC1_Input_fwd_chr19_Smartfilter.bed.zip", 
                         package = "DEScan")
print(bed.files)

peaks <- findPeaks(bed.files, chr = 19, filetype = "bed", fraglen = 200,
                   rlen = 100, min_bin = 50, max_win = 20, blocksize = 10000,
                   zthresh = 5, min_count = 0.1, verbose = FALSE, save = FALSE) 


peak.files <- system.file("extdata",
                          "")
regions <- finalRegions(peak_path = , zthresh = 20,  
                        min_carriers = 3, save_file = "bed", verbose = FALSE)

