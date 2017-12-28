checkFileExistence <- function(input_file_name){
  if(! file.exists(input_file_name)){
    stop(paste("The following file does not exist: ", input_file_name, sep=""))
  }
}
  
  
#' Perform analysis on a WT - KO embryo pipeline: intensity, relative Jacobians, absolute Jacobians and struture analysis
#' 
#' @param csv_determinants -- the determinants.csv file that is produced by the Pydpiper pipeline
#' @param pipelinedir      -- top level directory of the pipeline. I.e., the one that contains the _processed directory
#' @param final_nlin_mask  -- a mask for the final non-linear average of the pipeline. Used to limit the analysis to areas within the embryo.
#' @param final_nlin_segmentation -- segmented labels for the final non linear embryo atlas
#' @param mapping_KO_WT    -- a csv file that maps Image.Name to Genotype (WT or KO)
#' @param fwhm             -- blur level to use for the Jacobian determinants
#' @param defs             -- label definitions 
#' 
#' @examples
#' EmbryoAnalysis("/hpf/largeprojects/MICe/stsatski/Satb2/embryo-pipe/determinants.csv", "/hpf/largeprojects/MICe/stsatski/Satb2/embryo-pipe/", "mask.mnc", "segmentation.mnc", "genotype_mapping.csv")
EmbryoAnalysis <- function(csv_determinants, pipelinedir, final_nlin_mask, final_nlin_segmentation, mapping_KO_WT, fwhm=0.2, defs="/axiom2/projects/software/embryo_atlas/version_2/E15.5_labels_version_2_mapping.csv"){
  library(RMINC)
  
  # first do some error checking
  # 1) first argument should be a csv file that is produced by the pipeline:
  checkFileExistence(csv_determinants)
  pipeline_csv <- read.csv(csv_determinants)
  # get all basenames
  pipeline_csv$basename <- as.character(pipeline_csv$xfm)
  for(i in 1:dim(pipeline_csv)[1]){
    pipeline_csv$basename[i] <- as.character(sub("_lsq6_lsq12_and_nlin_inverted.xfm", "", basename(as.character(pipeline_csv$xfm[i]))))
  }
  pipeline_csv_no_blur <- subset(pipeline_csv, pipeline_csv$fwhm == 0.0)
  #pipeline_csv_fwhm0.2 <- subset(pipeline_csv, pipeline_csv$fwhm == 0.2)
  
  if(is.null(pipeline_csv$fwhm)){
    stop(paste("The csv from the pipeline does not contain the column fwhm. Is this the correct csv file: ", csv_determinants, sep=""))
  }
  
  # 2) correct pipeline top level dir? We can simlpy check that one of the files from the csv exists:
  if(! file.exists(paste(pipelinedir, "/", pipeline_csv$log_full_det[1], sep=""))){
    # hmmm... this file does not exist, is the argument given actually a directory?
    if(! dir.exists(pipelinedir)){
      stop(paste("The provided argument for the top level directory of the pipeline does not exist: ", pipelinedir, sep-""))
    }
    # so the provided directory exists, but we can not find one of the input log determinant files....
    stop(paste("Can not find the following file:", pipelinedir, "/", pipeline_csv$log_full_det[1], sep=""))
  }
  
  # 3) is the mask a mask?
  checkFileExistence(final_nlin_mask)
  mask_file <- mincGetVolume(final_nlin_mask)
  if(max(mask_file) != 1){
    stop(paste("The maximum value in the mask file is not 1. Please provide a mask file with 0s and 1s. Now provided: ",  final_nlin_mask, sep=""))
  }
  
  # 4) is the segmentation file a readable file?
  checkFileExistence(final_nlin_segmentation)
  segmentation <- mincGetVolume(final_nlin_segmentation)
  
  # 5) can we find all the files from the KO/WT mapping in the determinant file?
  checkFileExistence(mapping_KO_WT)
  file_to_genotype_mapping <- read.csv(mapping_KO_WT)
  # we need to have both the Image.Name and the Genotype column:
  if(is.null(file_to_genotype_mapping$Image.Name)){
    stop(paste("The csv file mapping Genotype to Image.Name does not have the column Image.Name: ", mapping_KO_WT, sep=""))
  }
  if(is.null(file_to_genotype_mapping$Genotype)){
    stop(paste("The csv file mapping Genotype to Image.Name does not have the column Genotype: ", mapping_KO_WT, sep=""))
  }
  # we should have exactly the same number of samples in both input csv-s:
  if(dim(file_to_genotype_mapping)[1] != dim(pipeline_csv_no_blur)[1]){
    stop("The CSV file from the Pydpiper pipeline and the CSV file mapping filenames to genotype do not have the same number of samples...")
  }
  
  # add the mapping to the dataframe with the stats files information
  pipeline_csv_no_blur$Genotype <- factor(file_to_genotype_mapping$Genotype, levels=c("KO", "WT"))
  for(i in 1:dim(file_to_genotype_mapping)[1]){
    base_to_look_for <- sub(".mnc", "", as.character(file_to_genotype_mapping$Image.Name[i]))
    if(! base_to_look_for %in% pipeline_csv_no_blur$basename){
      stop(paste("Can not file the file: ", as.character(file_to_genotype_mapping$Image.Name[i]), " with base: ", base_to_look_for, " in the csv file that belongs to the Pydpiper pipeline.", sep=""))
    }
    pipeline_csv_no_blur$Genotype[which(pipeline_csv_no_blur$basename == base_to_look_for)] <- file_to_genotype_mapping$Genotype[i]
  }
  pipeline_csv_no_blur$Genotype <-relevel(pipeline_csv_no_blur$Genotype,ref="WT")
  
  message("\nWe are using the following mapping for the analysis, please verify:")
  print(cbind(pipeline_csv_no_blur$basename, as.character(pipeline_csv_no_blur$Genotype)))
  message("\nNumber of samples in each group:")
  print(table(pipeline_csv_no_blur$Genotype))
  message("\nNote that the reference level for the analysis is the WT group (same as what Michael's script did).\n")
  
  name_processed_dir <- dirname(dirname(dirname(as.character(pipeline_csv_no_blur$xfm[1]))))
  
  ###
  ### Intensity analysis
  ###
  message("Starting analysis (linear model) of the normalized and blurred intensities")
  pipeline_csv_no_blur$normalized_blur <- Sys.glob(paste(pipelinedir, "/", name_processed_dir, "/", pipeline_csv_no_blur$basename, "/resampled/", pipeline_csv_no_blur$basename, "*normalized_blur.mnc", sep=""))
  vs_intensities <- mincLm(pipeline_csv_no_blur$normalized_blur ~ pipeline_csv_no_blur$Genotype, mask=final_nlin_mask)
  qvals_intensities <- mincFDR(vs_intensities, mask=final_nlin_mask)
  message("FDR thresholds comparing image intensities (normalized_blur) between the genotypes:")
  print(qvals_intensities)
  mincWriteVolume(vs_intensities, "t_stats_intensities.mnc", "tvalue-pipeline_csv_no_blur$GenotypeKO")
  rm(vs_intensities, qvals_intensities)
  invisible(gc())
  
  ###
  ### Relative Jacobians -- (default: fwhm 0.2)
  ###
  message(paste("Starting analysis (linear model) of the relative Jacobians at a FWHM of: ", as.character(fwhm), sep=""))
  fwhm_for_test <- as.numeric(fwhm)
  pipeline_csv_fwhm_requested <- subset(pipeline_csv, pipeline_csv$fwhm == fwhm_for_test)
  if(dim(pipeline_csv_fwhm_requested)[1] != dim(pipeline_csv_no_blur)[1]){
    stop(paste("Can not find the requested level of blurring (fwhm) in the Pydpiper CSV. Value provided: ", as.character(fwhm), sep=""))
  }
  pipeline_csv_fwhm_requested$rel_jac_full_path <- Sys.glob(paste(pipelinedir, "/", pipeline_csv_fwhm_requested$log_nlin_det, sep=""))
  vs_relative <- mincLm(pipeline_csv_fwhm_requested$rel_jac_full_path ~ pipeline_csv_no_blur$Genotype, mask=final_nlin_mask)
  qvals_relative <- mincFDR(vs_relative, mask=final_nlin_mask)
  message(paste("FDR thresholds comparing relative Jacobians between the genotypes at FWHM: ", as.character(fwhm), sep=""))
  print(qvals_relative)
  mincWriteVolume(vs_relative, paste("t_stats_relative_Jacobians_fwhm", as.character(fwhm), ".mnc", sep=""), "tvalue-pipeline_csv_no_blur$GenotypeKO")
  rm(vs_relative, qvals_relative)
  invisible(gc())
  
  ###
  ### Absolute Jacobians -- (default: fwhm 0.2)
  ###
  message(paste("Starting analysis (linear model) of the absolute Jacobians at a FWHM of: ", as.character(fwhm), sep=""))
  pipeline_csv_fwhm_requested$abs_jac_full_path <- Sys.glob(paste(pipelinedir, "/", pipeline_csv_fwhm_requested$log_full_det, sep=""))
  vs_absolute <- mincLm(pipeline_csv_fwhm_requested$abs_jac_full_path ~ pipeline_csv_no_blur$Genotype, mask=final_nlin_mask)
  qvals_absolute <- mincFDR(vs_absolute, mask=final_nlin_mask)
  message(paste("FDR thresholds comparing absolute Jacobians between the genotypes at FWHM: ", as.character(fwhm), sep=""))
  print(qvals_absolute)
  mincWriteVolume(vs_absolute, paste("t_stats_absolute_Jacobians_fwhm", as.character(fwhm), ".mnc", sep=""), "tvalue-pipeline_csv_no_blur$GenotypeKO")
  rm(vs_absolute, qvals_absolute)
  invisible(gc())
  
  ###
  ### Atlas analysis
  ###
  message("Starting analysis (linear model) of the labels in the provided atlas.")
  pipeline_csv_no_blur$full_path_log_full_det <- Sys.glob(paste(pipelinedir, "/", pipeline_csv_no_blur$log_full_det, sep=""))
  vols <- anatGetAll(pipeline_csv_no_blur$full_path_log_full_det, atlas=final_nlin_segmentation, defs=defs)
  embryo_vols <- anatGetAll(pipeline_csv_no_blur$full_path_log_full_det, atlas=final_nlin_mask, defs="/projects/moush/mwong/E15.5_mouse_model/wholebody.csv")
  rownames(vols) <- as.character(pipeline_csv_no_blur$basename)
  vols_with_whole_body <- cbind(vols, embryo_vols)
  anatLm_out <- anatLm(~ Genotype, pipeline_csv_no_blur, vols_with_whole_body)
  write.csv(anatLm_out, "t-values-structure_analysis_absolute_volume.csv")
  anatFDR_out <- anatFDR(anatLm_out)
  message("FDR thresholds absolute volumes (atlas):")
  print(anatFDR_out)
  write.csv(anatFDR_out, "FDR-values-structure_analysis_absolute_volume.csv")
  rm(anatLm_out, anatFDR_out)
  invisible(gc())

  vols_normalized <- vols / as.numeric(embryo_vols)
  anatLm_out_rel <- anatLm(~ Genotype, pipeline_csv_no_blur, vols_normalized)
  write.csv(anatLm_out_rel, "t-values-structure_analysis_relative_volume.csv")
  anatFDR_out_rel <- anatFDR(anatLm_out_rel)
  message("FDR thresholds relative volumes (atlas):")
  print(anatFDR_out_rel)
  write.csv(anatFDR_out_rel, "FDR-values-structure_analysis_relative_volume.csv")
  rm(anatLm_out_rel, anatFDR_out_rel)
  invisible(gc())
}


