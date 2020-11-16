library(MSnbase)
library(ggplot2)
library(stringr)
library(miscTools)
library(dplyr)

isotopeEnrichment <- function (
  mqdata, rawfilepath, special_residues = c(),
  ppmTol = 5, rtWindow = 0.5, base_isotopes = 4,
  fit_max_iter = 200, plotOutputs = FALSE ) {

  # read mq data
  print('Reading MQ data...')
  mqdata <- read.csv(mqdata, sep='\t')

  # remove unidentified entries
  print('Filtering MS/MS scan data...')
  mqdata <- subset(mqdata, Potential.contaminant != '+')
  mqdata <- subset(mqdata, Reverse != '+')
  
  # select smaller subset of cols
  cols <- c('Sequence', 'Modifications', 'Retention.time', 'Mass', 'Charges', 'Intensity')
  mqdata <- mqdata[ , cols]

  # Charges column can have multiple values - '2;3;4'
  # split these and duplicate the row for each charge
  mqdataUnwrapped <- NULL
  for (rowi in 1:nrow(mqdata) ) {
    for (charge in strsplit(mqdata[rowi, 'Charges'],';',fixed=TRUE)[[1]]) {
      newRow <- mqdata[rowi,]
      newRow$Charge <- as.numeric(charge)        
      mqdataUnwrapped <- rbind(mqdataUnwrapped, newRow)
    }
  }

  mqdata <- mqdataUnwrapped

  # calculated m/z ratios
  mqdata$m.z <- (mqdata$Mass + mqdata$Charge * 1.007276 ) / mqdata$Charge
  
  # add ll/hl tol limits
  mqdata$ll <- mqdata$m.z - ppmTol / 1000000 * mqdata$m.z
  mqdata$hl <- mqdata$m.z + ppmTol / 1000000 * mqdata$m.z

  mqdata$rt_ll <- mqdata$Retention.time - rtWindow
  mqdata$rt_hl <- mqdata$Retention.time + rtWindow
  
  # read mzml file
  msfiles <- dir(rawfilepath, full.names = TRUE, pattern = 'mzML$')
  rawdata <- readMSData(msfiles, msLevel = 1, mode='onDisk')

  # subset out rt ranges - seems that conversion to seconds is needed
  rtr <- mqdata[, c('rt_ll', 'rt_hl')] *60

  # subset out mz ranges
  mzr <- mqdata[, c('ll', 'hl')]

  # add column for special residue count
  mqdata$special_residues <- 0
  special_residues <- c('S', 'G')
  for (sr in special_residues) {
    mqdata$special_residues <- mqdata$special_residues + str_count(mqdata$Sequence, sr)
  }

  # plot chromatograms
  # https://rdrr.io/bioc/MSnbase/man/Chromatogram-class.html
  print('Getting chromatogram data...')
  Chrom <- chromatogram(rawdata, rt = rtr, mz = mzr, missing = 0)
  outputTable <- find_isotope_intensities( rawdata, Chrom, mqdata, base_isotopes, ppmTol, fit_max_iter, plotOutputs )
  return(outputTable)
}

find_isotope_intensities <- function ( rawdata, Chrom, mqdata, base_isotopes, ppmTol, fit_max_iter, plotOutputs) {
  fp <- getwd()
  
  # get vector of sample names 
  msFiles <- phenoData(rawdata)$sampleNames
  
  # the default ggplot theme isn't nice
  mytheme <- theme_bw(base_size = 12)
  mytheme <- mytheme + theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.text = element_text(size=8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 10, color = 'black'),
    axis.title = element_text(size = 10, color = 'black'),
    panel.border = element_rect(linetype = "solid", colour = "black", size=1),
  )
  
  outputTable <- NULL
  
  for (chromi in 1:nrow(mqdata)) {
    
    plotData <- NULL
    fitData <- NULL
    
    for (filei in seq_along(msFiles)) {
      
      # get chromatogram data
      chromData <- as.data.frame(Chrom[chromi, filei]) 
      
      # retrieve matching peptide from MQ table
      mq_peptide <- mqdata[chromi, ]
      
      # add filename to peptide object - needed for output table
      mq_peptide$file <- msFiles[filei]
      
      # get parameters for fit guess
      max_int_row <- chromData[which.max(chromData$intensity),]
      mu <- max_int_row$rtime
      k <- max_int_row$intensity
      
      try ( {
        # fit gaussian to chromData
        # https://stats.stackexchange.com/questions/83022/how-to-fit-data-that-looks-like-a-gaussian
        
        res <- nls(
          intensity ~ k*exp(-1/2*(rtime-mu)^2/sigma^2), start=c(mu=mu,sigma=3,k=k),
          data = chromData,
          control=nls.control(maxiter = fit_max_iter, warnOnly=TRUE)
        )
      }, silent = TRUE )

      # place these as defaults 
      chromData$rtfitints <- chromData$rtime       
      v <- c(mu, 3, k)
      
      try ({
        # vector of fitting parameters
        v <- summary(res)$parameters[,"Estimate"]
        
        # compute new gaussian values and add to dataframe
        # sometimes this still fails with the error "element (2, 2) is zero, so the inverse cannot be computed "
        chromData$rtfitints <- unlist(lapply(chromData$rtime, function(x) v[3]*exp(-1/2*(x-v[1])^2/v[2]^2 )))             
      })
      
      # get isotope intensities from spectra
      isotopic_abundances <- get_ms_intensities(chromData, mq_peptide, rawdata, base_isotopes, ppmTol)
      
      # this returns a numeric vector - transform to dataframe 
      isotopic_abundances <- data.frame(t(isotopic_abundances))
      
      # merge peptide and intensity data
      entityData <- merge(mq_peptide, isotopic_abundances)
      
      # add tou output data frame
      if (is.null(outputTable)) {
        outputTable <- entityData
      } else {
        outputTable <- dplyr::bind_rows(outputTable, entityData)
      }
      
      if (plotOutputs == TRUE) {
        # add file name to dataframe
        chromData$fname <- msFiles[filei]
        plotData <- rbind(plotData, chromData)
        
        # calculate the gaussian fit value at each x point and plot a line
        fit <- data.frame('rtime' = seq(from = min(chromData$rtime), to = max(chromData$rtime), length.out = nrow(chromData) * 10))
        fit$intensity <- v[3]*exp(-1/2*(fit$rtime-v[1])^2/v[2]^2)
        fit$fname <- msFiles[filei]
        fitData <- rbind(fitData, fit)
   
      }
    }
    
    if (plotOutputs == TRUE) {
    
      p <- ggplot() + mytheme + scale_y_continuous(labels = scales::scientific) +
      
      # plot the experimental points 
      geom_point(
        aes(rtime, intensity, color = fname),
        data = plotData, shape = 1, size = 2
      ) + 
      
      # add gaussian fit
      geom_line(aes(rtime, intensity, color = fname), data = fitData) +

      xlab('Retention time (s)') +
      ylab('Intensity') +
      ggtitle(paste(c('m/z', format(round(mqdata[chromi, 'm.z'], 4), nsmall = 4)), collapse = ' '))
      
      
      fname <- paste(c(chromi, '.png'), collapse = '')
      f <- file.path(fp, fname)
      ggsave(plot = p, filename = f, dpi = 300)
    }
  }
  return (outputTable)
}

get_ms_intensities <- function (chromData, mq_peptide, rawdata, base_isotopes, ppmTol) {
  
  # find point of max intensity in gaussian fit
  # max_fit_row_index <- which.max(chromData$rtfitints) --- testing
  max_fit_row_index <- which.max(chromData$intensity)
  
  # take 1 spectra either side of this point
  hl <- max_fit_row_index + 1
  ll <- max_fit_row_index - 1

  # retrieve relevant spectra
  spectraToAverage <- rownames(na.omit(chromData[ll:hl,]))
  
  # create new data frame to hold isotopic abundances
  num_isotopes <- base_isotopes + mq_peptide$special_residues

  isotopicAbundanceDF <- data.frame(matrix(ncol = num_isotopes, nrow = 0))
  
  for (speci in seq_along(spectraToAverage)) {

    # get spectrum
    s <- spectraToAverage[speci]
    sp <- as.data.frame(rawdata[[s]])
    
    # for a given target, get isotopic abundances
    for (isoi in 1:num_isotopes) {
      
      # calculate extraction boundaries
      # isoi - 1 so that the first isotope is 0 - i.e. monoisotopic
      iso_target <- mq_peptide$m.z + (isoi-1) * 1.008664 / mq_peptide$Charge
      iso_ll <- iso_target - ppmTol / 1000000 * iso_target
      iso_hl <- iso_target + ppmTol / 1000000 * iso_target

      # subset and sum spectral intensities
      window <- subset(sp, mz > iso_ll & mz < iso_hl)
      iso_intensity <- sum(window$i)
      
      # add to dataframe
      isotopicAbundanceDF[speci, isoi] = iso_intensity
      
    }
  }
  return (colSums(isotopicAbundanceDF))
}

msmsScansFile <- 'D:/path/to/modificationSpecificPeptides.txt'
rawfilepath <- 'D:/path/to/mzmlfiles/'
special_residues <- c('S', 'G')

outputTable <- isotopeEnrichment( msmsScansFile, rawfilepath, special_residues, plotOutputs = FALSE)

