# Package 'spicyWorkflow' Review
Thank you for submition your package to Bioconductor. The package passed check and build. But when build the vignettes, my system got crashed and restarted. My system is  macOS Big Sur M1 chip. I don't think this is acceptable. Please try to fix the issue and answer the following comments line by line when you are ready for a second review.
Code: Note: please condsider; Important: must be addressed.

## The DESCRIPTION file
- [ ] At least 3 sentences in `Description` field is required.
- [ ] update the R required version to 4.3.0
- [ ] NOTE: Consider adding the maintainer's ORCID iD in 'Authors@R'
      with 'comment=c(ORCID="...")'
 
## R code:
- [ ] Important: fix the hard coding:
        * rmd file vignettes/spicyWorkflow.Rmd
            - line 111: pathToImages <- "../inst/extdata/images" # suggest `?system.file`
            - line 128: dir.create("../inst/extdata/h5Files") # change to a tmp file
            - line 132: h5FilesPath = "../inst/extdata/h5Files",
            - line 147: clinical <- read.csv("../inst/extdata/1-s2.0-S0092867421014860-mmc1.csv")
- [ ] Important: Use `vapply` to replace `sapply`
        * rmd file vignettes/spicyWorkflow.Rmd
            - line 118: files <- sapply(imageDirs, list.files, pattern = "tif", full.names = TRUE, simplify = FALSE)
 - [ ] Auto detect available cores and minimize the used cores for the vignette at
         * rmd file vignettes/spicyWorkflow.Rmd
             - line 91: nCores <- 5
 
## Documentation
- [ ] Important: Vignette should have an *Installation* section.
	* rmd file vignettes/spicyWorkflow.Rmd
- [ ] Note: Vignette includes `motivation for submitting to Bioconductor` as part of the abstract/intro of the main vignette.
	* rmd file vignettes/spicyWorkflow.Rmd
- [ ] Important: Cause system crash by build the  vignettes/spicyWorkflow.Rmd.
- [ ] Important: move the `docs` folder and the `_pkgdown.yml` to a separate branch like `gh-pages`.

