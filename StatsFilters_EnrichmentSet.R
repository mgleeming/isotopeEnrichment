



StatsFilters_EnrichmentSet <- function(EnrichmentFileDir,
                                       Treatment,
                                       LabelFactor) {
  
  
  ## Statistical filters to reduce the full list to significantly labelled peptides 
  
  ### Statistical filters 1 - data reduction based on quality
  
  enrichment_table <- read.csv(file = EnrichmentFileDir, header = T, row.names = 1)
  
  #### removing features with null standard deviation
  ##### (Levene´s test produces errors when there is no variance between treatments)
  
  enrichment_table = as.matrix(enrichment_table[-c(which(rowSds(as.matrix(enrichment_table)) == 0)),])
  
  #### tranforming zeros to VERY small values
  
  enrichment_table[which(enrichment_table == 0)] <- runif(1, min=0, max=0.00001)
  
  #### removing features with null variance in individual treatments
  ##### (Levene´s test produces errors when there is no variance between treatments)
  
  Treatment_vector <- Treatment
  
  mean_treatment_sd <- c()
  
  for (i in 1:dim(enrichment_table)[1]) {
    
    mean_treatment_sd[i] <- mean(aggregate(as.numeric(enrichment_table[i,]) ~ Treatment_vector,
                                           FUN = sd)[,2])
  }
  
  null_sd_features <- rownames(enrichment_table[which(mean_treatment_sd == 0),])
  
  if(length(null_sd_features) > 0) {
    
    enrichment_table = enrichment_table[-c(which(mean_treatment_sd == 0)),]
    
  }
  
  ### Statistical filters 2 - RandoDiStats data reduction based on class comparison
  
  
  test_stats_Lab <- OmicsUnivariateStats(class_comparison_mat = t(enrichment_table),
                                         Factor1 = LabelFactor,
                                         Contrast = F)
  
  test_stats_sf <- OmicsUnivariateStats(class_comparison_mat = t(enrichment_table),
                                        Factor1 = Treatment_vector,
                                        Contrast = F)
  
  #### grabbing the significant peptides from the tests
  
  print(colnames(test_stats_sf))
  
  return_colNr <- readline(prompt="which column indexes contain the labelled treatment pvalues? x1,y1,z1....: ")
  
  returnNum <- as.numeric(strsplit(return_colNr, ",")[[1]])
  
  sig_labTreatment <- list()
  
  for (i in 1:length(returnNum)) {
    
    runner <- rownames(test_stats[which(test_stats[,returnNum[i]] < 0.05),])
    
    sig_labTreatment[[i]] <- runner
    
  }

  sig_Treatments <- Reduce(union, sig_labTreatment)
  
  sig_Lab <- rownames(test_stats_Lab[which(test_stats_Lab[,"Factor1Labelled_P values"] < 0.05),])
  
  
  #### all significant peptides labeled during any treatment
  
  sig_ALL <- Reduce(union, list(sig_Treatments, sig_Lab))
  
  #### subset enrichment table
  
  enrichment_table_sig <- enrichment_table[sig_ALL,]
  
  
  ### Statistical filters 3 - data visualization based on class discovery 
  ### (quality check of the processing)
  
  data_mat <- as.matrix(enrichment_table_sig)
  
  df <- t(aggregate(t(data_mat) ~ Treatment_vector, FUN = mean))
  
  data_mat_means <- apply(df[2:dim(df)[1],], 2, as.numeric)
  colnames(data_mat_means) <- df[1,]
  
  type = gsub("s\\d+_", "", lapply(unique(Treatment_vector), FUN = as.character))
  
  ha = HeatmapAnnotation(df = data.frame(type = type))
  
  ### scaling the mat for resolution
  
  scaled_mat <- (data_mat_means - rowMeans(data_mat_means)) / rowSds(data_mat_means)
  
  Scal_Heatmap <- Heatmap(matrix = scaled_mat,
                          clustering_distance_columns = "pearson",
                          clustering_method_columns = "average"  ,
                          name = "Intensities",
                          km = 5, 
                          row_km_repeats = 1000, ## repeats to get consensus
                          col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1),
                                           c("white",
                                             "yellow",
                                             "darkgoldenrod1",
                                             "violet",
                                             "purple")),
                          top_annotation = ha , show_column_names = F, show_row_names = F,
                          cluster_columns = F,
                          cluster_rows = T,
                          clustering_method_rows = "average",
                          clustering_distance_rows = "pearson")
  
  print(Scal_Heatmap)
  
  return_decision <- readline(prompt="Do you want me to remove false enrichments? (Y/N): ")
  
  if (return_decision == "Y") {
    
    return_decision2 <- readline(prompt="Above which percentage do you consider false enrichments? (0-100%): ")
    
    rm_percent = as.numeric(return_decision2)/100
    
    ### Statistical filters 4 - filtering: in the significant set,
    ### i.e., enrichment_table_sig, 
    ### filter out things that are falsely "Enriched" in the controls above 5%
    
    rm_falseP_sig<-Reduce(union,
                          list(as.numeric(which(enrichment_table_sig[,"X201112_Federico_01"] > rm_percent)),
                               as.numeric(which(enrichment_table_sig[,"X201112_Federico_02"] > rm_percent)),
                               as.numeric(which(enrichment_table_sig[,"X201112_Federico_03"] > rm_percent)),
                               as.numeric(which(enrichment_table_sig[,"X201112_Federico_07"] > rm_percent)),
                               as.numeric(which(enrichment_table_sig[,"X201112_Federico_08"] > rm_percent)),
                               as.numeric(which(enrichment_table_sig[,"X201112_Federico_09"] > rm_percent))))
    
    enrichment_table_sig = enrichment_table_sig[-c(rm_falseP_sig),]
    
  }
  
  return(enrichment_table_sig)
  
}
