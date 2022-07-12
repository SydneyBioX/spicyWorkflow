library(EBImage)

# Following instructions from the RISOM 2022 paper.
# https://doi.org/10.1016/j.cell.2021.12.023
# Images are downloaded from https://data.mendeley.com/datasets/d87vg86zd8

dir.create("Data")
pathToImages <- "/dski/nobackup/biostat/datasets/spatial/Risom2022_BreastCancer_IMC/data/Imaged_Data/raw_tifs_masks"
patients <- dir(pathToImages)
patients <- patients[-grep("Normal",patients)]
patients <- patients[-grep("Tonsil",patients)]

tmp <- sapply(patients, function(patient){
files <- list.files(paste0(pathToImages,"/",patient,"/TIFs"))
img <- EBImage::readImage(paste0(pathToImages,"/",patient,"/TIFs/", files))
img <- img[,,!dimnames(img)[[3]]%in%c("Au", "Aumask", "Background", "C", "Ca40")]
writeImage(img, paste0("Data/images/", patient, ".tiff"), compression = "LZW")
})

file.rename("Data/images/Point6202_pt1027_20594.tiff", "Data/images/Point6202_pt1026_20594.tiff")






