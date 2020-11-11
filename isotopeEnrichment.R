

isotopeEnrichment <- function(PyResultsDir,
                              returnCSV = T,
                              verbose = F) {
  
  #### get data
  
  data <- read.table(file = PyResultsDir, header = T, sep = "\t")
  
  ### removing peptides without G or S
  
  data = data[-c(which(data$Special.Residue.Count == 0)),]
  rownames(data) <- paste0(rep("X", length(rownames(data))), rownames(data))
  
  ### Treatments
  
  Treatments <- as.character(unique(data$File))
  
  test_Treatments <- c()
  
  for (i in 1:length(Treatments)) {
    
    test_Treatments[i] <- length(strsplit(Treatments, split = "")[[i]]) > 0
  }
  
  Treatments = Treatments[grep(T, test_Treatments)]
  
  TreatmentNo = length(Treatments)
  
  #### functions needed
  
  list2df <- function(x) 
  { 
    MAX.LEN <- max(sapply(x, length), na.rm = TRUE) 
    DF <- data.frame(lapply(x, function(x) c(x, rep(NA, MAX.LEN - length(x))))) 
    colnames(DF) <- paste("V", seq(ncol(DF)), sep = "")   
    DF 
  }
  
  get_isotopologues <- function(singlesingle_pep,
                                Treatments) {
    
    TreatmentNo = length(Treatments)
    
    dim_test <- na.omit(t(singlesingle_pep[1,10:dim(singlesingle_pep)[2]]))
    
    out_mat <- matrix(NA, nrow = dim(dim_test)[1], ncol = TreatmentNo)
    
    colnames(out_mat) <- Treatments
    
    for (i in 1:dim(singlesingle_pep)[1]) {
      
      df_single_row <- na.omit(t(singlesingle_pep[i,10:dim(singlesingle_pep)[2]]))
      
      alocation <- which(singlesingle_pep$File[i] == colnames(out_mat))
      
      out_mat[,alocation] <- as.numeric(df_single_row[1:dim(dim_test)[1],1])
      
    }
    
    IsoIDs <- paste0(as.character(rep(unique(singlesingle_pep$Sequence), dim(dim_test)[1])),
                     "_",
                     as.character(rep("Z", dim(dim_test)[1])),
                     as.character(rep(unique(singlesingle_pep$Charge), dim(dim_test)[1])),
                     "_",
                     apply(list2df(strsplit(rownames(df_single_row), "\\.")), 2, as.character)[2,])
    
    return_mat <- cbind("Measurements/Samples" = IsoIDs, out_mat)
    
  }
  
  #### Main
  
  ##### Measurement file
  
  ###### First step must be subset the matrix to group of identical peptides in bins
  
  peptides <- unique(data$Sequence)
  
  cat(paste0(length(peptides), "Peptides found"))
  
  MeasurementFile <- matrix(NA, nrow = 0, ncol = (TreatmentNo + 1))
  
  Molecule <- c()
  
  for (i in 1:length(peptides)) {
    
    test_sequence <- length(strsplit(as.character(peptides[i]), "")[[1]])
    
    if (test_sequence > 0){
      
      single_pep <- data[grep(paste0("\\b",peptides[i], "\\b"), data$Sequence),]
      
      ## filtering out repeated peptide peaks
      
      single_pep = single_pep[!duplicated(single_pep$Intensity.0),]
      
      ## grabbing identical peptides but differently charged separately
      
      test_charges <- unique(single_pep$Charge)
      
      if (verbose == T) {
        
        cat("...", "\n")
        cat(paste0("Peptide # ", i, " : ", peptides[i]))
        cat("\n")
        cat(paste0("Charge State(s): ", length(test_charges)))
        cat("...", "\n")
        
      }
      
      if (length(test_charges) > 1){
        
        MeasurMat <- matrix(NA, nrow = 0, ncol = (TreatmentNo + 1))
        
        for (j in 1:length(test_charges)) {
          
          singlesingle_pep <- single_pep[grep(test_charges[j], single_pep$Charge),]
          
          ## grabbing molecular formula and potentially labelled residue number
          
          runner <- paste0(gsub(" ", x = unique(singlesingle_pep$Formula),
                                replacement = "", fixed = T),
                           "LabN",
                           unique(singlesingle_pep$Special.Residue.Count),
                           collapse = "")
          
          Molecule <- c(Molecule, runner)
          
          ## grabbing the isotopologues
          
          MeasurMat <- rbind(MeasurMat,
                             get_isotopologues(singlesingle_pep,
                                               Treatments = as.character(unique(data$File)[1:TreatmentNo])))  
        }
        
      } else {
        
        ## peptides with only one charge
        
        singlesingle_pep <- single_pep
        
        ## grabbing molecular formula and potentially labelled residue number
        
        runner <- paste0(gsub(" ", x = unique(singlesingle_pep$Formula),
                              replacement = "", fixed = T),
                         "LabN",
                         unique(singlesingle_pep$Special.Residue.Count),
                         collapse = "")
        
        Molecule <- c(Molecule, runner)
        
        ## grabbing the isotopologues
        
        MeasurMat <- get_isotopologues(singlesingle_pep,
                                       Treatments = as.character(unique(data$File)[1:TreatmentNo]))
        
      } 
    }
    MeasurementFile <- rbind(MeasurementFile, MeasurMat)
  }
  
  ## identifying and dealing with duplicated peptides (Enrichment does not allow duplicated names)
  
  MeasurementFile = MeasurementFile[-c(which(duplicated(x = MeasurementFile, MARGIN = 1) == T)),]
  
  ##### Molecule file
  
  Molecule_names <- list2df(strsplit(MeasurementFile[,1], "_"))
  
  Sequence_part <- apply(Molecule_names[1,], 2, as.character)
  
  charges_part <- apply(Molecule_names[2,], 2, as.character)
  
  MoleculeFile <- as.data.frame(cbind(Molecule = unique(paste0(Sequence_part, "_", charges_part)),
                                      "MS ion or MS/MS product ion" = Molecule,
                                      "MS/MS neutral loss" = NA))
  
  if (returnCSV == T) {
    
    write.csv(x = MoleculeFile, file = "MoleculeFile.csv", quote = F, row.names = F)
    write.csv(x = MeasurementFile, file = "MeasurementFile.csv", quote = F, row.names = F)
    
  }
  return(list(MeasurementFile,
              MoleculeFile))
}