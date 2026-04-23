library('targets')
library('tarchetypes')
library(crew)
tar_option_set(controller = crew_controller_local(workers = 10)) 
Sys.setenv(VROOM_CONNECTION_SIZE = as.character(10 * 1024 * 1024)) #For any large datasets

####This runs a workflow for cleaning the MESA TOPMed multi-omics project proteomics data (O-link), 
#### formatting it into the following format: wide for metabolites & long for exams, 
#### and providing basic quality control (QC) metrics



tar_option_set(packages = c("dplyr", "tidyr", "tibble", "readr", "data.table", "bit64", 
                            "foreign", "quarto", "rlang", "purrr", "rcompanion", "knitr", "gtsummary", "labelled",
                            "kableExtra", "gt", "cli", "quarto", "sandwich", "lme4", "lmerTest", "ggplot2", "glmnet", "fastDummies", "doParallel"))

tar_source("/media/Analyses/MESA-MIND-Single-Time-Point-Proteins-SVD/R")

list(
  
  #---------------------------------------------------------------------------------------#
  #--------------------------------1. Build proteomics table------------------------------#
  #---------------------------------------------------------------------------------------#
  
  #--------------------------------1A. Initial table
  
  #Protein intensities
  tar_target(path_proteins, "/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Proteomics/Olink/SMP_IntensityNormalized_20251005.csv", format = "file"),
  #Mapping info
  tar_target(path_protein_info, "/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Proteomics/Olink/Mapping_SMP_Plate_20251005.csv", format = "file"),
  #Protein keys
  tar_target(path_protein_keys, "/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Proteomics/Olink/MESAOlink3k_proteinKeys_03292023.csv", format = "file"),
  
  #Bridging file
  tar_target(path_bridge,"/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESA-SHARE_IDList_Labeled.csv", format = "file"),
  
  #Run function
  tar_target(build_proteins_out,   
             build_protein_table_function(
               path_proteins = path_proteins,
               path_protein_info = path_protein_info,
               path_protein_keys = path_protein_keys,
               path_bridge = path_bridge)
  ),
  
  #Save outputs
  #QC info
  tar_target(build_proteins_QC, build_proteins_out$QC_info_out),
  #Proteins
  tar_target(Proteins_long, build_proteins_out$Formatted_proteins_out),
  #N in formatted file by exam
  tar_target(Proteins_long_N, build_proteins_out$N_by_exam),
  
  #--------------------------------1B. Mapping file
  
  #QC info
  tar_target(protein_info, QC_proteins_function(path_proteins = path_proteins,
                                                path_protein_info = path_protein_info,
                                                path_protein_keys = path_protein_keys,
                                                Proteins_long = Proteins_long)),
  tar_target(Protein_mapping_file, protein_info$final_proteins_mapping),
  
  #--------------------------------1C. Final clean protein table
  
  #Build proteins
  tar_target(proteins_clean, final_proteins_function(Proteins_long = Proteins_long, 
                                                     QC_file = Protein_mapping_file)
  ),
  
  
  #Save proteins
  tar_target(Proteins_long_clean, proteins_clean$Final_proteins),
  #N in cleaned & formatted file by exam
  tar_target(Proteins_clean_N, proteins_clean$N_by_exam),
  
  
  #---------------------------------------------------------------------------------------#
  #--------------------------------2. Build traits table----------------------------------#
  #---------------------------------------------------------------------------------------#
  
  #E1_covs
  tar_target(path_E1_covs, "/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAe1FinalLabel02092016.dta", format = "file"),
  #E5_covs
  tar_target(path_E5_covs, "/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAe5_FinalLabel_20140613.dta", format = "file"),
  #E6_covs
  tar_target(path_E6_covs, "/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAe6_FinalLabel_20220513.dta", format = "file"),
  #afib
  tar_target(path_afib, "/media/RawData/MESA/MESA-Phenotypes/MESA-SHARe-Phenos/SHARe_MesaEventsThruYear2020_AF_DS.txt", format = "file"),
  tar_target(path_apoe, "/media/RawData/MESA/MESA-Phenotypes/MESA-MIND/MESA_ApoE_03102014.sas7bdat", format = "file"),
  tar_target(path_cvd, "/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAEvThru2020AllCohort_20241120.dta", format = "file"),
  tar_target(path_mb, "/media/RawData/MESA/MESA-Phenotypes/MESA-MIND/MESAe6as253as301_BMRICMB_08052025.csv", format="file"),
  tar_target(path_evps, "/media/RawData/MESA/MESA-Phenotypes/MESA-MIND/MESAe6as253as301_BMRIPVS_20250310.csv", format="file"),
  tar_target(path_wmh, "/media/RawData/MESA/MESA-Phenotypes/MESA-MIND/MESAe6anyFIRST_BMRIWMHVol_20240422.csv", format="file"),
  tar_target(path_icv, "/media/RawData/MESA/MESA-Phenotypes/MESA-MIND/Archive/SHARe_AncilMesaAF_BMRIROIVol_DS.txt", format="file"),
  tar_target(path_wmfa, "/media/RawData/MESA/MESA-Phenotypes/MESA-MIND/mesae6anyfirst_bmriTotalFAMUSE_20250828.csv", format="file"),
  
  tar_target(build_traits, build_traits_function(path_E1_covs = path_E1_covs,
                                                 path_E5_covs = path_E5_covs,
                                                 path_E6_covs = path_E6_covs,
                                                 path_afib = path_afib,
                                                 path_bridge = path_bridge,
                                                 path_cvd = path_cvd,
                                                 path_apoe = path_apoe,
                                                 path_mb = path_mb,
                                                 path_evps = path_evps,
                                                 path_wmh = path_wmh,
                                                 path_icv = path_icv,
                                                 path_wmfa = path_wmfa,
                                                 cleaned_proteins = Proteins_long)),
  
  tar_target(traits_QC_info, build_traits$QC_info_out),
  tar_target(traits_db, build_traits$Traits_table),
  tar_target(mind_ids, build_traits$mind_ids),
  tar_target(protein_and_MIND_ids, build_traits$protein_and_MIND_ids),
  #tar_target(protein_and_MIND_and_cov_ids, build_traits$protein_and_MIND_and_cov_ids),
  tar_target(missing_cov_info, build_traits$missing_cov_info),


#--------------------------------------------------------------------------------------------------------#
#--------------------------------3. Run Cross-sectional PWAS -------------------------------------------#
#-------------------------------------------------------------------------------------------------------#

##-------------White matter hyperintensities-----------------#
tar_target(wmh_numeric_covs, c("icv", "age", "egfr", "BMI", "sbp", "ldl")),

tar_target(wmh_factor_covs, c("gender", "race", "edu", "htnmeds", "smoking", "E4",
                              "AFprevalent", "diabetes", "MIprevalent", "CHFprevalent")),


tar_target(WMH_E6_PWAS, cross_sectional_PWAS_function(cleaned_proteins = Proteins_long_clean, 
                                                      protein_mapping = Protein_mapping_file,
                                                      traits_db = traits_db, 
                                                      outcome = "wmh",
                                                      numeric_covariates = wmh_numeric_covs,
                                                      factor_covariates = wmh_factor_covs,
                                                      chosen_exam = 6)),

#Save as csv file
tar_target(save_WMH_E6_PWAS,
           {
             out_dir  <- "outputs"
             dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
             out_path <- file.path(out_dir, "WMH_E6_PWAS.csv")
             readr::write_csv(WMH_E6_PWAS, out_path)
             out_path
           },
           format = "file"
),



tar_target(WMH_E1_PWAS, cross_sectional_PWAS_function(cleaned_proteins = Proteins_long_clean, 
                                                      protein_mapping = Protein_mapping_file,
                                                      traits_db = traits_db, 
                                                      outcome = "wmh",
                                                      numeric_covariates = wmh_numeric_covs,
                                                      factor_covariates = wmh_factor_covs,
                                                      chosen_exam = 1)),

#Save as csv file
tar_target(save_WMH_E1_PWAS,
           {
             out_dir  <- "outputs"
             dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
             out_path <- file.path(out_dir, "WMH_E1_PWAS.csv")
             readr::write_csv(WMH_E1_PWAS, out_path)
             out_path
           },
           format = "file"
),

##-------------Perivascular spaces. ----------------#
tar_target(epvs_numeric_covs, c("icv", "age", "egfr", "BMI", "sbp", "ldl")),
tar_target(epvs_factor_covs, c("gender", "race", "edu", "htnmeds", "smoking", "E4",
                               "AFprevalent", "diabetes", "MIprevalent", "CHFprevalent")),


tar_target(EPVS_thalamus_E6_PWAS, cross_sectional_PWAS_function(cleaned_proteins = Proteins_long_clean, 
                                                                protein_mapping = Protein_mapping_file,
                                                                traits_db = traits_db, 
                                                                outcome = "epvs_thalamus",
                                                                numeric_covariates = epvs_numeric_covs,
                                                                factor_covariates = epvs_factor_covs,
                                                                chosen_exam = 6)),

#Save as csv file
tar_target(save_EPVS_thalamus_E6_PWAS,
           {
             out_dir  <- "outputs"
             dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
             out_path <- file.path(out_dir, "EPVS_thalamus_E6_PWAS.csv")
             readr::write_csv(EPVS_thalamus_E6_PWAS, out_path)
             out_path
           },
           format = "file"
),

tar_target(EPVS_thalamus_E1_PWAS, cross_sectional_PWAS_function(cleaned_proteins = Proteins_long_clean, 
                                                                protein_mapping = Protein_mapping_file,
                                                                traits_db = traits_db, 
                                                                outcome = "epvs_thalamus",
                                                                numeric_covariates = epvs_numeric_covs,
                                                                factor_covariates = epvs_factor_covs,
                                                                chosen_exam = 1)),


#Save as csv file
tar_target(save_EPVS_thalamus_E1_PWAS,
           {
             out_dir  <- "outputs"
             dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
             out_path <- file.path(out_dir, "EPVS_thalamus_E1_PWAS.csv")
             readr::write_csv(EPVS_thalamus_E1_PWAS, out_path)
             out_path
           },
           format = "file"
),

tar_target(EPVS_basalganglia_E6_PWAS, cross_sectional_PWAS_function(cleaned_proteins = Proteins_long_clean, 
                                                                    protein_mapping = Protein_mapping_file,
                                                                    traits_db = traits_db, 
                                                                    outcome = "epvs_basalganglia",
                                                                    numeric_covariates = epvs_numeric_covs,
                                                                    factor_covariates = epvs_factor_covs,
                                                                    chosen_exam = 6)),

#Save as csv file
tar_target(save_EPVS_basalganglia_E6_PWAS,
           {
             out_dir  <- "outputs"
             dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
             out_path <- file.path(out_dir, "EPVS_basalganglia_E6_PWAS.csv")
             readr::write_csv(EPVS_basalganglia_E6_PWAS, out_path)
             out_path
           },
           format = "file"
),

tar_target(EPVS_basalganglia_E1_PWAS, cross_sectional_PWAS_function(cleaned_proteins = Proteins_long_clean, 
                                                                    protein_mapping = Protein_mapping_file,
                                                                    traits_db = traits_db, 
                                                                    outcome = "epvs_basalganglia",
                                                                    numeric_covariates = epvs_numeric_covs,
                                                                    factor_covariates = epvs_factor_covs,
                                                                    chosen_exam = 1)),


#Save as csv file
tar_target(save_EPVS_basalganglia_E1_PWAS,
           {
             out_dir  <- "outputs"
             dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
             out_path <- file.path(out_dir, "EPVS_basalganglia_E1_PWAS.csv")
             readr::write_csv(EPVS_basalganglia_E1_PWAS, out_path)
             out_path
           },
           format = "file"
),

tar_target(EPVS_insula_E6_PWAS, cross_sectional_PWAS_function(cleaned_proteins = Proteins_long_clean, 
                                                              protein_mapping = Protein_mapping_file,
                                                              traits_db = traits_db, 
                                                              outcome = "epvs_insula",
                                                              numeric_covariates = epvs_numeric_covs,
                                                              factor_covariates = epvs_factor_covs,
                                                              chosen_exam = 6)),

#Save as csv file
tar_target(save_EPVS_insula_E6_PWAS,
           {
             out_dir  <- "outputs"
             dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
             out_path <- file.path(out_dir, "EPVS_insula_E6_PWAS.csv")
             readr::write_csv(EPVS_insula_E6_PWAS, out_path)
             out_path
           },
           format = "file"
),

tar_target(EPVS_insula_E1_PWAS, cross_sectional_PWAS_function(cleaned_proteins = Proteins_long_clean, 
                                                              protein_mapping = Protein_mapping_file,
                                                              traits_db = traits_db, 
                                                              outcome = "epvs_insula",
                                                              numeric_covariates = epvs_numeric_covs,
                                                              factor_covariates = epvs_factor_covs,
                                                              chosen_exam = 1)),


#Save as csv file
tar_target(save_EPVS_insula_E1_PWAS,
           {
             out_dir  <- "outputs"
             dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
             out_path <- file.path(out_dir, "EPVS_insula_E1_PWAS.csv")
             readr::write_csv(EPVS_insula_E1_PWAS, out_path)
             out_path
           },
           format = "file"
),
##-------------Fractional ansiotropy ----------------#
tar_target(fa_numeric_covs, c("age", "egfr", "BMI", "sbp", "ldl")),
tar_target(fa_factor_covs, c("gender", "race", "edu", "htnmeds", "smoking", "E4",
                             "AFprevalent", "diabetes", "MIprevalent", "CHFprevalent")),


tar_target(FA_E6_PWAS, cross_sectional_PWAS_function(cleaned_proteins = Proteins_long_clean, 
                                                     protein_mapping = Protein_mapping_file,
                                                     traits_db = traits_db, 
                                                     outcome = "fa",
                                                     numeric_covariates = fa_numeric_covs,
                                                     factor_covariates = fa_factor_covs,
                                                     chosen_exam = 6)),

#Save as csv file
tar_target(save_FA_E6_PWAS,
           {
             out_dir  <- "outputs"
             dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
             out_path <- file.path(out_dir, "FA_E6_PWAS.csv")
             readr::write_csv(FA_E6_PWAS, out_path)
             out_path
           },
           format = "file"
),

tar_target(FA_E1_PWAS, cross_sectional_PWAS_function(cleaned_proteins = Proteins_long_clean, 
                                                     protein_mapping = Protein_mapping_file,
                                                     traits_db = traits_db, 
                                                     outcome = "fa",
                                                     numeric_covariates = fa_numeric_covs,
                                                     factor_covariates = fa_factor_covs,
                                                     chosen_exam = 1)),

#Save as csv file
tar_target(save_FA_E1_PWAS,
           {
             out_dir  <- "outputs"
             dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
             out_path <- file.path(out_dir, "FA_E1_PWAS.csv")
             readr::write_csv(FA_E1_PWAS, out_path)
             out_path
           },
           format = "file"
),

##-------------LASSO ----------------#

tar_target(FA_E1_LASSO, cross_sectional_LASSO_function(cleaned_proteins  = Proteins_long_clean,
                                                       traits_db = traits_db,
                                                       numeric_covariates = fa_numeric_covs, 
                                                       factor_covariates = fa_factor_covs, 
                                                       outcome = "fa", 
                                                       chosen_exam = 1, 
                                                       MWASres = FA_E1_PWAS, 
                                                       threshold = 0.05)),

tar_target(FA_E1_LASSO_res, FA_E1_LASSO$LASSO_output), #LASSO for Quarto
tar_target(FA_E1_LASSO_forcsv, FA_E1_LASSO$LASSO_output_all_vars), #Full LASSO for csv 

#Save full LASSO as csv file
tar_target(save_FA_E1_LASSO_forcsv,
           {
             out_dir  <- "outputs"
             dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
             out_path <- file.path(out_dir, "FA_E1_PWAS_fullLASSO.csv")
             readr::write_csv(FA_E1_LASSO_forcsv, out_path)
             out_path
           },
           format = "file"
),




##-------------Microbleeds ----------------#
tar_target(mb_numeric_covs, c("age", "egfr", "BMI", "sbp", "ldl")),
tar_target(mb_factor_covs, c("gender", "race", "edu", "htnmeds", "smoking", "E4",
                             "AFprevalent", "diabetes", "MIprevalent", "CHFprevalent")),


tar_target(MB_E6_PWAS, cross_sectional_logistic_PWAS_function(cleaned_proteins = Proteins_long_clean, 
                                                              protein_mapping = Protein_mapping_file,
                                                              traits_db = traits_db, 
                                                              outcome = "mb_present",
                                                              numeric_covariates = mb_numeric_covs,
                                                              factor_covariates = mb_factor_covs,
                                                              chosen_exam = 6)),

#Save as csv file
tar_target(save_MB_E6_PWAS,
           {
             out_dir  <- "outputs"
             dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
             out_path <- file.path(out_dir, "MB_E6_PWAS.csv")
             readr::write_csv(MB_E6_PWAS, out_path)
             out_path
           },
           format = "file"
),

tar_target(MB_E1_PWAS, cross_sectional_logistic_PWAS_function(cleaned_proteins = Proteins_long_clean, 
                                                              protein_mapping = Protein_mapping_file,
                                                              traits_db = traits_db, 
                                                              outcome = "mb_present",
                                                              numeric_covariates = mb_numeric_covs,
                                                              factor_covariates = mb_factor_covs,
                                                              chosen_exam = 1)),

#Save as csv file
tar_target(save_MB_E1_PWAS,
           {
             out_dir  <- "outputs"
             dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
             out_path <- file.path(out_dir, "MB_E1_PWAS.csv")
             readr::write_csv(MB_E1_PWAS, out_path)
             out_path
           },
           format = "file"
),


#---------------------------------------------------------------------------------------#
#--------------------------------Quarto output------------------------------------------#
#---------------------------------------------------------------------------------------#

tarchetypes::tar_quarto(
  build_quarto,
  path = "/media/Analyses/MESA-MIND-Single-Time-Point-Proteins-SVD/MESA-MIND-Single-Time-Point-Proteins-SVD.qmd",
  quiet = FALSE
)



)


