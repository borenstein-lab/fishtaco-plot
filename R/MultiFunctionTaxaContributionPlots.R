#' A function to plot different versions of the FishTaco plot for multiple functions
#'
#' @param input_dir The directory where all the FiShTaCo input files can be found (default: NULL).
#' @param input_prefix The prefix for all input files (e.g., pathway_with_t2f, default: NULL).
#' @param input_suffix The suffix for all input files (default: ".tab").
#' @param input_score The calculated shift scoring metric (default: wilcoxon).
#' @param input_contribution The taxon-level contribution profile (default: "taxa_contributions").
#' @param input_original The original (metagenome-based) calculated shift (default: "original_value").
#' @param input_permutation The permutation type used (default: "multi_taxa").
#' @param input_predicted_da The predicted (taxa-based) shift score (default: "predicted_DA_value").
#' @param input_predicted_function The predicted (taxa-based) functional abundance (default: "predicted_function_abundance").
#' @param input_predicted_residual The residual functional abundance (default: "residual_function_abundance").
#' @param input_predicted_abundance_agreement The agreement between the taxa- and metagenome based functional abundance (default: "predicted_function_agreement").
#' @param input_inference_copy_number The inferred copy number for each function (default: "taxa_learned_copy_num").
#' @param input_taxa_da The calculated taxonomic composition shift score (default: "DA_taxa").
#' @param input_function_abundance The metagenome-based functional abundance (default: NULL).
#' @param input_function_meta Metadata on the given functions (default: NULL).
#' @param input_taxa_taxonomy Phylogenetic assignment for each taxon (default: NULL).
#' @param input_taxa_vs_function The genomic content of each taxon (default: NULL).
#' @param input_function_counts The median count of functions across samples (default: NULL).
#' @param input_function_stats Statistics about each function across the taxa (default: NULL).
#' @param input_function_filter_list A list of specific functions to plot (default: NULL).
#' @param input_function_filter_file A file containing a list of specific functions to plot (default: NULL).
#' @param plot_type The type of plot to be plotted (default: "bars").
#' @param sort_by How to sort the functions in the plot (default: "predicted_da").
#' @param min_contribution Minimum taxon contribution to plot (default: 0).
#' @param min_cont_as_separate Minimum contribution to plot as separate species (default: 0.2).
#' @param max_function_to_show Maximum number of functions to show (default: NULL).
#' @param show_only_diff_abun_taxa Show only taxa that shift in abundance (default: FALSE).
#' @param show_only_pos Show only taxa that contribute positively (drive) to the shift (default: FALSE).
#' @param show_only_neg Show only taxa that contribute negatively (attenuate) to the shift (default: FALSE).
#' @param add_original_da_markers Add diamond markers for the metagenome-based shift score (default: FALSE).
#' @param show_only_taxa_with_function Show only taxa that have the function in their genome (default: FALSE).
#' @param show_only_enriched_taxa Show only taxa that are case-associated (default: FALSE).
#' @param show_only_enriched_functions Show only functions that are case-enriched (default: FALSE).
#' @param show_only_depleted_taxa Show only taxa that are control-associated (default: FALSE).
#' @param show_only_depleted_functions Show only functions that are control-enriched (default: FALSE).
#' @param scale_pos_by_original_pos_sum Scale driving contribution by their sum (default: FALSE)
#' @param scale_neg_by_original_neg_sum Scale attenuating contribution by their sum (default: FALSE)
#' @param add_facet_labels Add labels to each function (default: FALSE).
#' @param add_names_in_bars Add taxa names in the bars (default: FALSE).
#' @param add_taxa_da Add a top bar of taxonomic shifts (default: FALSE).
#' @param add_case_control_line Add a line separating case- and control-enriched functions (default: FALSE).
#' @param add_predicted_da_markers Add diamond markers for the taxa-based shift score (default: TRUE).
#' @param input_bar_width Defines the width of the bars (default: NULL).
#' @param verbose Run with verbose (default: FALSE).
#' @param flip_coord Flip axis to show horizontal bars (default: TRUE).
#' @param split_function_names_2_lines Splits function names into 2 lies (default: TRUE).
#' @param separate_enriched_depleted_taxa Separate case- and control-associated taxa (default: TRUE).
#' @param use_facets_for_separation Use facets for separation of functions (default: FALSE).
#' @param flip_case_control Switch between cases and controls (default: FALSE).
#'
#'
#'
#' @return a handle to the resulting plot.
#' @export
MultiFunctionTaxaContributionPlots <- function(input_dir=NULL,
                                               input_prefix=NULL,
                                               input_suffix=".tab",
                                               input_score="wilcoxon",
                                               input_contribution="taxa_contributions",
                                               input_original="original_value",
                                               input_permutation="multi_taxa",
                                               input_predicted_da="predicted_DA_value",
                                               input_predicted_function="predicted_function_abundance",
                                               input_predicted_residual="residual_function_abundance",
                                               input_predicted_abundance_agreement="predicted_function_agreement",
                                               input_inference_copy_number="taxa_learned_copy_num",
                                               input_taxa_da="DA_taxa",
                                               input_function_abundance=NULL,
                                               input_function_meta=NULL,
                                               input_taxa_taxonomy=NULL,
                                               input_taxa_vs_function=NULL,
                                               input_function_counts=NULL,
                                               input_function_stats=NULL,
                                               input_function_filter_list=NULL,
                                               input_function_filter_file=NULL,
                                               plot_type="bars",
                                               sort_by="predicted_da",
                                               min_contribution=0,
                                               min_cont_as_separate=0.025,
                                               min_function_counts=0,
                                               min_shift_explained=0,
                                               max_function_to_show=NULL,
                                               show_only_diff_abun_taxa=FALSE,
                                               show_only_pos=FALSE,
                                               show_only_neg=FALSE,
                                               show_only_taxa_with_function=FALSE,
                                               show_only_enriched_taxa=FALSE,
                                               show_only_enriched_functions=FALSE,
                                               show_only_depleted_taxa=FALSE,
                                               show_only_depleted_functions=FALSE,
                                               scale_pos_by_original_pos_sum=FALSE,
                                               scale_neg_by_original_neg_sum=FALSE,
                                               add_facet_labels=FALSE,
                                               add_names_in_bars=FALSE,
                                               add_taxa_da=FALSE,
                                               add_case_control_line=FALSE,
                                               add_predicted_da_markers=TRUE,
                                               add_original_da_markers=FALSE,
                                               add_phyla_names_to_species_label=FALSE,
                                               color_small_cont_by_phyla=FALSE,
                                               input_bar_width=NULL,
                                               verbose=FALSE,
                                               flip_coord=TRUE,
                                               split_function_names_2_lines=TRUE,
                                               separate_enriched_depleted_taxa=TRUE,
                                               use_facets_for_separation=FALSE,
                                               flip_case_control=FALSE) {




  #####################
  # example input
  #####################
  use_example_input = FALSE
  
  if (use_example_input) { 
    input_dir = "/Volumes/ohadm/OhadM/METAFIT/HMP_DATA/TONGUE_DORSUM_vs_BUCCAL_MUCOSA/Output_METAPHLAN_REMOVE_RESIDUAL_SCALE_PERMUTED_MAPPING"
    # "/Volumes/ohadm/OhadM/METAFIT/HMP_DATA/TONGUE_DORSUM_vs_BUCCAL_MUCOSA/Output_METAPHLAN_REMOVE_RESIDUAL_SCALE_PERMUTED_MAPPING"
    # "/Volumes/ohadm/OhadM/METAFIT/IGC/QIN_T2D/Output_Genus_REMOVE_RESIDUAL_SCALE_PERMUTED"

    input_prefix = "pathway_with_t2f"# "module_with_t2f"# "class_with_t2f" # "ko_with_t2f" # "pathway_with_t2f"

    input_score = "wilcoxon" # "wilcoxon"
    input_contribution = "taxa_contributions"
    input_taxa_da = "DA_taxa"
    input_predicted_da = "predicted_DA_value"
    input_inference_copy_number="taxa_learned_copy_num"
    input_predicted_function = "predicted_function_abundance"
    input_predicted_residual="residual_function_abundance"
    input_predicted_abundance_agreement="predicted_function_agreement"
    input_original = "original_value"
    input_permutation = "multi_taxa" #"single_taxa"
    input_suffix = ".tab"

    plot_type= "function_DA_bars"
    # return_final_taxa_list
    # percent_unknown_no_negative function_DA_bars predicted_DA_bars percentage_bars percent_predicted
    # percent_abundance_explained_pearson percent_abundance_explained
    # bars sum_vs_original sum_vs_predicted percent_predicted_diamonds

    input_function_meta= "/Volumes/ohadm/OhadM/MUSiCC/Matrices/PATHWAYvsNAME_BACTERIAL_KEGG_2013_07_15.lst"
    #"/Volumes/ohadm/OhadM/MUSiCC/Matrices/KOvsNAME_KEGG_2013_07_15.lst"
    #"/Volumes/ohadm/OhadM/MUSiCC/Matrices/MODULEvsNAME_BACTERIAL_KEGG_2013_07_15.lst"
    #"/Volumes/ohadm/OhadM/MUSiCC/Matrices/PATHWAYvsNAME_BACTERIAL_KEGG_2013_07_15.lst"
    # NULL

    input_function_counts = "/Volumes/ohadm/OhadM/METAFIT/HMP_DATA/TONGUE_DORSUM_vs_BUCCAL_MUCOSA/median_PATHWAY_counts.tab"
    # "/Volumes/ohadm/OhadM/METAFIT/IGC/QIN_T2D/median_KO_counts.tab"
    # "/Volumes/ohadm/OhadM/METAFIT/IGC/QIN_T2D/median_MODULE_counts.tab"
    # "/Volumes/ohadm/OhadM/METAFIT/HMP_DATA/TONGUE_DORSUM_vs_BUCCAL_MUCOSA/median_PATHWAY_counts.tab"

    input_function_stats = NULL
    # "/Volumes/ohadm/OhadM/METAFIT/HMP_DATA/TONGUE_DORSUM_vs_BUCCAL_MUCOSA/PATHWAY_stats.tab"
    # NULL
    input_function_abundance=NULL
    # "/Volumes/ohadm/OhadM/METAFIT/MetaHIT/Data/WGS_PATHWAY_vs_SAMPLE_MUSiCC.tab"
    # NULL
    input_taxa_vs_function=NULL#"/Volumes/ohadm/OhadM/METAFIT/HMP_DATA/METAPHLAN_taxa_vs_PATHWAY.tab"# NULL

    input_taxa_taxonomy = "/Volumes/ohadm/OhadM/METAFIT/HMP_DATA/METAPHLAN_taxa_vs_TAXONOMY.tab"
    # "/Volumes/ohadm/OhadM/METAFIT/IGC/Data/genus_vs_TAXONOMY.tab"
    # "/Volumes/ohadm/OhadM/METAFIT/HMP_DATA/METAPHLAN_taxa_vs_TAXONOMY.tab"

    input_function_filter_file = "/Volumes/ohadm/OhadM/MUSiCC/Matrices/PATHWAY_BACTERIAL_KEGG_2013_07_15.lst"
    #"/Volumes/ohadm/OhadM/MUSiCC/Matrices/MODULE_BACTERIAL_KEGG_2013_07_15.lst"
    # "/Volumes/ohadm/OhadM/MUSiCC/Matrices/KO_BACTERIAL_KEGG_2013_07_15.lst"
    # "/Volumes/ohadm/OhadM/MUSiCC/Matrices/PATHWAY_BACTERIAL_KEGG_2013_07_15.lst"

    input_function_filter_list = c("ko02040", "ko00540","ko00020")
    # c("ko02040", "ko00540","ko00020")
    # c("ko02040", "ko00540", "ko00020") # NULL# c("M00217", "M00216", "M00356") c("K01426","K00128","K01464")

    min_contribution = 0
    min_cont_as_separate = 0.2# 0.5 0.2 0.025
    min_function_counts = 20#0#3#5#20
    min_shift_explained=0
    sort_by = "list" # "da"# "strongest" # da "list" "unknown" "predicted_da"
    max_function_to_show = 200
    add_names_in_bars = TRUE
    add_taxa_da = FALSE
    input_bar_width = 0.6
    verbose = TRUE
    flip_coord = FALSE

    flip_case_control=TRUE#FALSE#TRUE

    add_original_da_markers=TRUE
    show_only_taxa_with_function=FALSE
    show_only_diff_abun_taxa=TRUE

    show_only_pos=FALSE
    scale_pos_by_original_pos_sum=FALSE
    show_only_enriched_taxa=FALSE
    show_only_enriched_functions=FALSE

    show_only_neg=FALSE
    scale_neg_by_original_neg_sum=FALSE
    show_only_depleted_taxa=FALSE
    show_only_depleted_functions=FALSE
    add_case_control_line=FALSE
    add_predicted_da_markers=FALSE

    separate_enriched_depleted_taxa=TRUE
    use_facets_for_separation=FALSE
    add_facet_labels=FALSE
    color_small_cont_by_phyla=FALSE
    add_phyla_names_to_species_label=FALSE
    split_function_names_2_lines=TRUE
  }
  ####################################################################################
  # END of example input
  ####################################################################################


  ####################################################################################
  # CHECK INPUT FOR INCONSISTENCY
  ####################################################################################

  if (show_only_enriched_taxa && show_only_depleted_taxa) {
    print(paste("Error: show_only_enriched_taxa && show_only_depleted_taxa"))
    return(NULL)
  }

  if (show_only_enriched_functions && show_only_depleted_functions) {
    print(paste("Error: show_only_enriched_functions && show_only_depleted_functions"))
    return(NULL)
  }

  if (show_only_pos && show_only_neg) {
    print(paste("Error: show_only_pos && show_only_neg"))
    return(NULL)
  }

  ####################################################################################
  # READ INPUT FILES
  ####################################################################################

  #-----------------------
  # read taxa contribution
  #-----------------------
  input_file_name = paste(input_dir,"/",input_prefix,"_STAT_",input_contribution,"_SCORE_",
                          input_score,"_ASSESSMENT_",input_permutation,input_suffix, sep="")

  if (!file.exists(input_file_name)) {
    print(paste("Error: input file",  input_file_name,"does not exist"))
    return(NULL)
  }

  if (verbose) {
    print(paste("reading file:",input_file_name))
  }

  df = read.table(input_file_name, sep = "\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, quote="")

  if (flip_case_control) {
    df[, -1] = -1 * df[, -1, drop=FALSE]
  }

  functions = names(df)[2:length(names(df))]

  if (verbose) {
    print("functions:")
    print(paste(functions))
  }

  #-------------------------------
  # save the original positive and
  # negative contribution sums for
  # each function
  #-------------------------------
  tmp_df = df[, -1]
  tmp_df[tmp_df < 0] = 0
  original_sum_of_pos_contributions = data.frame(apply(tmp_df, 2, sum))
  names(original_sum_of_pos_contributions) = "PosSum"

  tmp_df = df[, -1]
  tmp_df[tmp_df > 0] = 0
  original_sum_of_neg_contributions = data.frame(apply(tmp_df, 2, sum))
  names(original_sum_of_neg_contributions) = "NegSum"

  #-------------------------------
  # read function filter file/list
  #-------------------------------
  if (!is.null(input_function_filter_file)) {
    function_filter_list_from_file = read.table(input_function_filter_file, sep = "\t", header=FALSE, stringsAsFactors=FALSE, quote="")
    functions = functions[functions %in% function_filter_list_from_file$V1]
    if (length(functions) == 0){
      print(paste("Error: the filter file:",input_function_filter_file,
                  "did not have any intersection with the functions in the input file:",
                  input_file_name))
      return(NULL)
    }
    df = df[c("Taxa",functions)]
  }

  if (!is.null(input_function_filter_list)) {
    function_filter_list_present = input_function_filter_list[input_function_filter_list %in% functions]
    functions = function_filter_list_present
    df = df[c("Taxa",functions)]
  }

  if (verbose) {
    print("functions after filter file/list:")
    print(paste(functions))
  }

  #-------------------------------
  # read taxonomy file
  #-------------------------------
  if (!is.null(input_taxa_taxonomy)) {
    taxonomy = read.table(input_taxa_taxonomy, sep = "\t", header=FALSE, stringsAsFactors=FALSE, quote="")
    rownames(taxonomy) = taxonomy[,1]
    taxonomy[,1] = NULL
    names(taxonomy) = c('kingdom','phylum','class','order','family','genus','species')
  }

  #-------------------------------
  # read function counts file
  # filter out low-counts functions
  #-------------------------------
  if (!is.null(input_function_counts)) {
    function_counts = read.table(input_function_counts, sep = "\t", header=TRUE, stringsAsFactors=FALSE, quote="")
    function_counts = function_counts[function_counts$Stats %in% functions,]
    function_counts = function_counts[function_counts$Median >= min_function_counts,]
    functions = functions[functions %in% function_counts$Stats]
    #print(functions)
    df = df[c("Taxa",functions)]
  }

  #-------------------------------
  # read function stats file
  #-------------------------------
  if (!is.null(input_function_stats)) {
    function_stats = read.table(input_function_stats, sep = "\t", header=TRUE, stringsAsFactors=FALSE, quote="")
    function_stats = function_stats[function_stats$Stats %in% functions,]
    row.names(function_stats) = function_stats[,1]
    function_stats[,1] = NULL
  }

  #-------------------------------
  # read taxa-to-function file
  # if we want to show only the
  # taxa containing the function,
  # zero out taxa not containing the
  # function (for each function)
  #-------------------------------
  if (!is.null(input_taxa_vs_function)) {
    taxa_to_function = read.table(input_taxa_vs_function, sep = "\t", header=TRUE, stringsAsFactors=FALSE, quote="")
    row.names(taxa_to_function) = taxa_to_function[,1]
    taxa_to_function[,1] = NULL
    taxa_to_function = taxa_to_function[, functions]

    if (show_only_taxa_with_function) {
      for (i in 1:length(functions)) {
        non_containing_taxa = row.names(taxa_to_function[taxa_to_function[, functions[i]] == 0, ])
        df[df$Taxa %in% non_containing_taxa, functions[i]] = 0
      }
    }
  }

  #-------------------------------
  # read predicted abundance
  # agreement file
  #-------------------------------
  if (!is.null(input_predicted_abundance_agreement)) {
    predicted_abundance_agreement = read.table(paste(input_dir,"/",input_prefix,"_STAT_",
                                                     input_predicted_abundance_agreement,"_SCORE_",
                                                     input_score,"_ASSESSMENT_",input_permutation,input_suffix, sep=""),
                                               sep = "\t", header=TRUE, stringsAsFactors=FALSE, quote="")

    predicted_abundance_agreement = predicted_abundance_agreement[predicted_abundance_agreement$KO %in% functions,]
  }

  #-------------------------------
  # read predicted DA/shift file
  #-------------------------------
  if (!is.null(input_predicted_da)) { #
    predicted_da = read.table(paste(input_dir,"/",input_prefix,"_STAT_",input_predicted_da,"_SCORE_",
                                    input_score,"_ASSESSMENT_",input_permutation,input_suffix, sep=""),
                              sep = "\t", header=TRUE, stringsAsFactors=FALSE, quote="")

    predicted_da = predicted_da[predicted_da$KO %in% functions,]
    if (flip_case_control) {
      predicted_da[, 2] = -1 * predicted_da[, 2]
    }
  }

  #-------------------------------
  # read original DA/shift data
  # for functions and sort accordingly
  #-------------------------------
  if (!is.null(input_original)) {
    function_da = read.table(paste(input_dir,"/",input_prefix,"_STAT_",input_original,"_SCORE_",
                                   input_score,"_ASSESSMENT_",input_permutation,input_suffix, sep=""),
                             sep = "\t", header=TRUE, stringsAsFactors=FALSE, quote="")

    function_da = function_da[function_da$KO %in% functions,]

    if (flip_case_control) {
      function_da[, 2] = -1 * function_da[, 2]
    }

    if (sort_by == "da") {
      sort_index = sort(function_da[,input_score], index.return=TRUE, decreasing=TRUE)$ix
    } else if (sort_by == "predicted_da"){
      sort_index = sort(predicted_da[,input_score], index.return=TRUE, decreasing=TRUE)$ix
    } else if (sort_by == "list"){
      tmp_da = function_da
      row.names(tmp_da) = tmp_da[, 1]
      tmp_da[, 1] = NULL
      tmp_da$index = 1:length(functions)
      tmp_da[, 1] = NULL
      sort_index = tmp_da[functions,] # sort by list, with repsect to location in DA dataframe
    } else if (sort_by == "strongest") {
      sort_index = 1:length(functions)
    }

    functions = function_da[sort_index, 1]

    # if we want to show only enriched functions, remove non
    if (show_only_enriched_functions) {
      functions = functions[functions %in% function_da$KO[function_da[,input_score] > 0]]
    }
    # if we want to show only depleted functions, remove non
    if (show_only_depleted_functions) {
      functions = functions[functions %in% function_da$KO[function_da[,input_score] < 0]]
    }

    if (!is.null(input_predicted_da)) {
      predicted_da = predicted_da[predicted_da$KO %in% functions,]
    }
    if (!is.null(input_predicted_abundance_agreement)) {
      predicted_abundance_agreement = predicted_abundance_agreement[predicted_abundance_agreement$KO %in% functions,]
    }

    df = df[c("Taxa",functions)]
    original_sum_of_pos_contributions = original_sum_of_pos_contributions[functions, "PosSum", drop=FALSE]
    original_sum_of_neg_contributions = original_sum_of_neg_contributions[functions, "NegSum", drop=FALSE]
  }

  #-------------------------------
  # compute the percent of DA/shift
  # explained for each function and
  # filter out functions with
  # badly predicted shifts
  #-------------------------------
  if (!is.null(min_shift_explained) && !is.null(input_original) && !is.null(input_predicted_da)) {
    merged_predicted_original_DA = merge(predicted_da, function_da, by="KO")
    merged_predicted_original_DA$explained = merged_predicted_original_DA[, 2] / merged_predicted_original_DA[, 3]
    filtered_functions = merged_predicted_original_DA$KO[merged_predicted_original_DA$explained > min_shift_explained]
    # update functions to filtered ones
    functions = functions[functions %in% filtered_functions]
    df = df[c("Taxa",functions)]
    function_da = function_da[function_da$KO %in% functions, ]
    predicted_da = predicted_da[predicted_da$KO %in% functions,]
    if (!is.null(input_predicted_abundance_agreement)) {
      predicted_abundance_agreement = predicted_abundance_agreement[predicted_abundance_agreement$KO %in% functions,]
    }
    if (!is.null(input_function_stats)) {
      function_stats = function_stats[row.names(function_stats) %in% functions,]
    }
    original_sum_of_pos_contributions = original_sum_of_pos_contributions[functions, "PosSum", drop=FALSE]
    original_sum_of_neg_contributions = original_sum_of_neg_contributions[functions, "NegSum", drop=FALSE]
  }


  #-------------------------------
  # if we have more functions then the maximum
  # remove the least DA functions
  #-------------------------------
  if (!is.null(max_function_to_show) && (length(functions) > max_function_to_show)) {
    if (show_only_enriched_functions) {
      functions = functions[(length(functions)-max_function_to_show+1):length(functions)]
    } else if (show_only_depleted_functions) {
      functions = functions[1:max_function_to_show]
    } else {
      cases_to_show = ceiling(max_function_to_show/2)
      controls_to_show = floor(max_function_to_show/2)
      functions = functions[c(1:controls_to_show, (length(functions)-cases_to_show+1):length(functions))]
    }


    df = df[c("Taxa",functions)]
    function_da = function_da[function_da$KO %in% functions, ]
    predicted_da = predicted_da[predicted_da$KO %in% functions,]
    if (!is.null(input_predicted_abundance_agreement)) {
      predicted_abundance_agreement = predicted_abundance_agreement[predicted_abundance_agreement$KO %in% functions,]
    }
    if (!is.null(input_function_stats)) {
      function_stats = function_stats[row.names(function_stats) %in% functions,]
    }
    original_sum_of_pos_contributions = original_sum_of_pos_contributions[functions, "PosSum", drop=FALSE]
    original_sum_of_neg_contributions = original_sum_of_neg_contributions[functions, "NegSum", drop=FALSE]
  }

  #-------------------------------
  # read taxa DA/shift data
  #-------------------------------
  if (!is.null(input_taxa_da)) { # we have DA information for taxa
    taxa_da = read.table(paste(input_dir,"/",input_prefix,"_STAT_",input_taxa_da,"_SCORE_",
                               input_score,"_ASSESSMENT_",input_permutation,input_suffix, sep=""),
                         sep = "\t", header=TRUE, stringsAsFactors=FALSE, quote="")

    taxa_da_original_names = taxa_da

    # add unkown to taxa_da
    if (sum(df$Taxa == "Unknown") > 0) {
      taxa_da = rbind(taxa_da, c("Unknown",0,0,0,1,0,0,0,0,0,"Unknown"))
      taxa_da$StatValue = as.numeric(taxa_da$StatValue)
    }

    if (flip_case_control) {
      taxa_da$StatValue = -1 * taxa_da$StatValue
    }

    # if we want to show only enriched/depleted taxa, then zero out all
    # contributions from other taxa
    if (show_only_enriched_taxa){
      df[df$Taxa %in% taxa_da$Taxa[taxa_da$StatValue <= 0], -1] = 0
    }
    if (show_only_depleted_taxa){
      df[df$Taxa %in% taxa_da$Taxa[taxa_da$StatValue >= 0], -1] = 0
    }
    # if we want to remove taxa with 0 shift in any direction
    if (show_only_diff_abun_taxa) {
      df[df$Taxa %in% taxa_da$Taxa[taxa_da$StatValue == 0], -1] = 0
    }

  }

  #-------------------------------
  # remove taxa with very low
  # contributions across all functions
  #-------------------------------
  if (!is.null(min_contribution)) {
    max_val = apply(abs(df[, -1, FALSE]), 1, max)
    taxa_with_contribution = which(max_val > min_contribution)
    df = df[taxa_with_contribution,]

    if (!is.null(input_taxa_da)) { # we have DA information for taxa
      taxa_da = taxa_da[taxa_with_contribution,]
      if(sum(df$Taxa != taxa_da$Taxa)) {
        stop("df$Taxa != taxa_da$Taxa")
      }
    }
  }

  #-------------------------------
  # read inferred copy number
  #-------------------------------
  if (plot_type == "inferred_copy_number_heatmap") {
    if (!is.null(input_inference_copy_number)) { #
      inferred_copy_number = read.table(paste(input_dir,"/",input_prefix,"_STAT_",input_inference_copy_number,"_SCORE_",
                                              input_score,"_ASSESSMENT_",input_permutation,input_suffix, sep=""),
                                        sep = "\t", header=TRUE, stringsAsFactors=FALSE, quote="")

      inferred_copy_number = inferred_copy_number[inferred_copy_number$Taxa %in% df$Taxa, c("Taxa", functions)]
    }
  }

  #-------------------------------
  # read function metadata and
  # rename functions by meta-names
  #-------------------------------
  if (!is.null(input_function_meta)) {

    if (!file.exists(input_function_meta)) {
      print(paste("Error: input file",input_function_meta,"does not exist"))
      return(NULL)
    }

    meta = read.table(input_function_meta, sep = "\t", header=FALSE, stringsAsFactors=FALSE,  quote="")
    rownames(meta) = meta[,1]
    meta[,1] = NULL
    # if we have non-unique meta-names for the functions, add the
    # function identifier to the meta-names
    if (length(unique(meta$V2)) < length(meta$V2)) {
      meta$V2 = paste0(row.names(meta),": ",meta$V2)
    }

    if (split_function_names_2_lines) {
      location_to_split = ceiling((sapply(gregexpr("\\W+", meta$V2), length)) / 2)
      for (i in 1:length(location_to_split)) {
        index_to_replace = gregexpr("\\W+", meta$V2[i])[[1]][location_to_split[i]]
        if (index_to_replace > 0) {
          substr(meta$V2[i], index_to_replace, index_to_replace) = "\n"
        }
      }
    }


    #print(head(meta))
    old_functions = functions
    functions = meta[functions,]
    missing_functions_from_meta = data.frame(old_functions[is.na(functions)])
    row.names(missing_functions_from_meta) = old_functions[is.na(functions)]
    names(missing_functions_from_meta) = "V2"
    meta = rbind(meta, missing_functions_from_meta)


    functions[is.na(functions)] = old_functions[is.na(functions)]
    #print(functions)
    names(df) = c('Taxa',functions)

    if (!is.null(input_function_stats)) {
      function_stats = merge(meta, function_stats, by="row.names")
      row.names(function_stats) = function_stats[,2]
      function_stats[,1] = NULL
      function_stats[,1] = NULL
    }
  }

  #-------------------------------
  # read real function abundance data
  #-------------------------------
  if (!is.null(input_function_abundance)) { # we have function abundance data
    if (!file.exists(input_function_abundance)) {
      print(paste("Error: input file",  input_function_abundance, "does not exist"))
      return(NULL)
    }
    real_function_abun = read.table(input_function_abundance, sep = "\t", header=TRUE,
                                    stringsAsFactors=FALSE, quote="")

    names(real_function_abun) = c("KO", names(real_function_abun)[2:length(names(real_function_abun))])
  }

  #-------------------------------
  # read predicted function abundance
  #-------------------------------
  if (!is.null(input_predicted_function)) { # we have predicted function abundance data
    if (!file.exists(paste0(input_dir,"/",input_prefix,"_STAT_",
                            input_predicted_function,"_SCORE_",
                            input_score,"_ASSESSMENT_",
                            input_permutation,input_suffix))) {
      print(paste("Error: input file",  paste0(input_dir,"/",input_prefix,"_STAT_",
                                               input_predicted_function,"_SCORE_",
                                               input_score,"_ASSESSMENT_",
                                               input_permutation,input_suffix), "does not exist"))
      return(NULL)
    }
    predicted_function_abun = read.table(paste0(input_dir,"/",input_prefix,"_STAT_",
                                                input_predicted_function,"_SCORE_",
                                                input_score,"_ASSESSMENT_",
                                                input_permutation,input_suffix), sep = "\t", header=TRUE,
                                         stringsAsFactors=FALSE, quote="")

    names(predicted_function_abun) = c("KO", names(predicted_function_abun)[2:length(names(predicted_function_abun))])
  }

  #-------------------------------
  # read residual of predicted
  # function abundance
  #-------------------------------
  if (!is.null(input_predicted_residual)) { # we have residual function abundance data
    if (!file.exists(paste0(input_dir,"/",input_prefix,"_STAT_",
                            input_predicted_residual,"_SCORE_",
                            input_score,"_ASSESSMENT_",
                            input_permutation,input_suffix))) {
      print(paste("Error: input file",  paste0(input_dir,"/",input_prefix,"_STAT_",
                                               input_predicted_residual,"_SCORE_",
                                               input_score,"_ASSESSMENT_",
                                               input_permutation,input_suffix), "does not exist"))
      return(NULL)
    }

    residual_function_abun = read.table(paste0(input_dir,"/",input_prefix,"_STAT_",
                                               input_predicted_residual,"_SCORE_",
                                               input_score,"_ASSESSMENT_",
                                               input_permutation,input_suffix), sep = "\t", header=TRUE,
                                        stringsAsFactors=FALSE, quote="")

    names(residual_function_abun) = c("KO", names(residual_function_abun)[2:length(names(residual_function_abun))])
  }

  num_of_functions = length(functions)
  num_of_taxa = length(df$Taxa)

  #-------------------------------
  # set bar width
  #-------------------------------
  if (!is.null(input_bar_width)) {
    bar_width=input_bar_width
  } else {
    bar_width = 0.3
  }

  #-------------------------------
  # print final functions/taxa
  #-------------------------------
  if (verbose) {
    print(functions)
    print(paste("#Taxa:",num_of_taxa,"#Functions:",num_of_functions))
  }

  #-------------------------------
  # make sure there are functions
  #-------------------------------
  if (num_of_functions < 1) {
    print("No functions to plot, exiting...")
    return(NULL)
  }
  
  ########################################################################################
  # FINISHED READING INPUT FILES
  ########################################################################################










  ############################################
  # PREPARE DATA FOR PLOT
  ############################################
  if (verbose) {
    print("Preparing data for plots")
  }

  # separate each function to pos and neg
  df_pos = data.frame(df['Taxa'])
  df_neg = data.frame(df['Taxa'])
  for (i in 1:num_of_functions){
    #print(i)
    #print(functions[i])
    #print(df[[functions[i]]])
    pos = df[[functions[i]]]
    pos[pos < 0] = 0
    pos_name = paste(functions[i],"-pos",sep="")
    df_pos[functions[i]] = pos

    neg = df[[functions[i]]]
    neg[neg > 0] = 0
    neg_name = paste(functions[i],"-neg",sep="")
    df_neg[functions[i]] = neg
  }

  if (!is.null(input_taxa_da)) { # we have DA information for taxa
    if(sum(df_pos$Taxa != taxa_da$Taxa)) {
      stop("df_pos$Taxa != taxa_da$Taxa")
    }
  }

  if (!is.null(input_taxa_taxonomy)) {
    # find the strong taxa that contribute the most and leave them intact,
    # but try to find their taxonomy name
    cont_fraction_pos = sweep(df_pos[,-1,FALSE], 2, colSums(df_pos[,-1,FALSE]), FUN='/')
    cont_fraction_pos[apply(cont_fraction_pos, 2, is.nan)] = 0

    cont_fraction_neg = sweep(df_neg[,-1,FALSE], 2, colSums(df_neg[,-1,FALSE]), FUN='/')
    cont_fraction_neg[apply(cont_fraction_neg, 2, is.nan)] = 0

    # list the strong contributors. If we have unknown, we want to show it even if 
    # it has a very slight contribution
    if (show_only_pos) {
      strong_contributors = sort(unique(c(which(apply(cont_fraction_pos, 1, max) > min_cont_as_separate),
                                          which(df$Taxa == "Unknown"))))
    } else if (show_only_neg) {
      strong_contributors = sort(unique(c(which(apply(cont_fraction_neg, 1, max) > min_cont_as_separate),
                                          which(df$Taxa == "Unknown"))))
    } else {
      strong_contributors = sort(unique(c(which(apply(cont_fraction_pos, 1, max) > min_cont_as_separate),
                                          which(apply(cont_fraction_neg, 1, max) > min_cont_as_separate),
                                          which(df$Taxa == "Unknown"))))
    }

    # get the sorted list of sum of contributions, if we want to sort by the strongest contributor
    sorted_list_of_contributors = sort(rowSums(df_pos[,-1,FALSE]), decreasing = TRUE,index.return = TRUE)$ix
    sorted_list_of_contributors = sorted_list_of_contributors[sorted_list_of_contributors %in% strong_contributors]

    #print(strong_contributors)
    #print(sorted_list_of_contributors)
    #print(df_pos$Taxa)
    #print(taxa_da$Taxa)

    if (length(strong_contributors) > 0) {
      for (i in 1:length(strong_contributors)) {
        curr_taxa = df_pos$Taxa[strong_contributors[i]]
        #print(strong_contributors[i])
        #print(curr_taxa)
        #print(taxonomy[curr_taxa,"species"])
        #print(taxonomy[curr_taxa,"genus"])
        if (!is.na(taxonomy[curr_taxa,"species"])){
          #print(paste(taxonomy[curr_taxa,"genus"],taxonomy[curr_taxa,"species"]))
          if (taxonomy[curr_taxa,"species"] != "s__") {
            curr_taxonomy_names = paste(taxonomy[curr_taxa,"phylum"], taxonomy[curr_taxa,"genus"], taxonomy[curr_taxa,"species"])
          } else if (taxonomy[curr_taxa,"genus"] != "g__") {
            if (!is.na(as.numeric(curr_taxa))) {
              curr_taxonomy_names = paste(taxonomy[curr_taxa,"phylum"], taxonomy[curr_taxa,"genus"], paste("s__",curr_taxa,sep=""))
            } else {
              curr_taxonomy_names = paste(taxonomy[curr_taxa,"phylum"], taxonomy[curr_taxa,"genus"])
            }
          } else if (taxonomy[curr_taxa,"family"] != "f__") {
            curr_taxonomy_names = paste(taxonomy[curr_taxa,"phylum"], taxonomy[curr_taxa,"family"], paste("s__",curr_taxa,sep=""))
          } else if (taxonomy[curr_taxa,"order"] != "o__") {
            curr_taxonomy_names = paste(taxonomy[curr_taxa,"phylum"], taxonomy[curr_taxa,"order"], paste("s__",curr_taxa,sep=""))
          } else if (taxonomy[curr_taxa,"order"] != "c__") {
            curr_taxonomy_names = paste(taxonomy[curr_taxa,"phylum"], taxonomy[curr_taxa,"class"], paste("s__",curr_taxa,sep=""))
          } else {
            curr_taxonomy_names = paste(taxonomy[curr_taxa,"kingdom"], taxonomy[curr_taxa,"phylum"], curr_taxa)
          }
        } else if (curr_taxa == "Unknown"){
          curr_taxonomy_names = curr_taxa
        } else {
          print(paste("Error: taxa",curr_taxa,"is missing a taxonomy assignment in the taxonomy file"))
          return(NULL)
        }

        if (strong_contributors[i] == sorted_list_of_contributors[1]){
          df_pos_name_of_strongest = curr_taxonomy_names
          #print(df_pos_name_of_strongest)
        }

        # if we want to separate enricehd/depleted taxa, add the information
        # to the taxonomy name
        if (separate_enriched_depleted_taxa) {
          if (taxa_da[taxa_da$Taxa == curr_taxa, "StatValue"] >= 0) {
            curr_taxonomy_names = paste0(curr_taxonomy_names,":ENRICHED")
          } else {
            curr_taxonomy_names = paste0(curr_taxonomy_names,":DEPLETED")
          }
        }

        df_pos$Taxa[strong_contributors[i]] = curr_taxonomy_names
        df_neg$Taxa[strong_contributors[i]] = curr_taxonomy_names

        if (!is.null(input_taxa_da)) { # we have DA information for taxa
          taxa_da$Taxa[strong_contributors[i]] = curr_taxonomy_names
        }
      }
    }

    #print(df_pos$Taxa)
    #print(taxa_da$Taxa)

    # for all the other taxa (that contribute less), merge them by pyhlum level
    for (i in 1:length(df_pos$Taxa)) {
      if (i %in% strong_contributors){ # skip strong contributors as they already have taxonomy assigned
        next
      }
      curr_taxa = df_pos$Taxa[i]
      #print(curr_taxa)
      if ( !is.na(taxonomy[curr_taxa,"phylum"]) && color_small_cont_by_phyla){

        if (separate_enriched_depleted_taxa) {
          if (taxa_da[taxa_da$Taxa == curr_taxa, "StatValue"] > 0) {
            df_pos$Taxa[i] = paste0(taxonomy[curr_taxa,"phylum"],":ENRICHED")
            df_neg$Taxa[i] = paste0(taxonomy[curr_taxa,"phylum"],":ENRICHED")
          } else {
            df_pos$Taxa[i] = paste0(taxonomy[curr_taxa,"phylum"],":DEPLETED")
            df_neg$Taxa[i] = paste0(taxonomy[curr_taxa,"phylum"],":DEPLETED")
          }
        } else {
          df_pos$Taxa[i] = paste(taxonomy[curr_taxa,"phylum"])
          df_neg$Taxa[i] = paste(taxonomy[curr_taxa,"phylum"])
        }


        if (!is.null(input_taxa_da)) { # we have DA information for taxa
          taxa_da$Taxa[i] = paste(taxonomy[curr_taxa,"phylum"])
        }
        
      } else {
        
        if (separate_enriched_depleted_taxa) {
          if (taxa_da[taxa_da$Taxa == curr_taxa, "StatValue"] > 0) {
            df_pos$Taxa[i] = paste0("ZOther taxa",":ENRICHED")
            df_neg$Taxa[i] = paste0("ZOther taxa",":ENRICHED")
          } else {
            df_pos$Taxa[i] = paste0("ZOther taxa",":DEPLETED")
            df_neg$Taxa[i] = paste0("ZOther taxa",":DEPLETED")
          }
        } else {
          df_pos$Taxa[i] = "ZOther taxa"
          df_neg$Taxa[i] = "ZOther taxa"
        }
        
        if (!is.null(input_taxa_da)) { # we have DA information for taxa
          taxa_da$Taxa[i] = "ZOther taxa"
        }
      }
    }

    #print(df_pos$Taxa)
    #print(taxa_da$Taxa)

    #-----------------------
    # aggregate similar taxa
    # to their sum and add
    # taxa counts for phyla
    #-----------------------

    if (separate_enriched_depleted_taxa) {
      tmp_df = df_pos[, 1, drop=FALSE]
      tmp_df$count = 1
      tmp_df$Taxa = sub(":DEPLETED","",(sub(":ENRICHED","",tmp_df$Taxa)))
      taxa_counts = aggregate(tmp_df$count, by=list(tmp_df$Taxa), sum)
      taxa_counts$left = "("
      taxa_counts$right = ")"
      taxa_counts$pasted = do.call(paste0, taxa_counts[c(1, 3, 2, 4)])
      taxa_counts$x = NULL
      taxa_counts$left = NULL
      taxa_counts$right = NULL
      names(taxa_counts) = c("Taxa","Full")
      # remove the (1) where only 1 species is present from species or genus (but not from pyhla)
      remove_indices = grepl("(1)", taxa_counts$Full, fixed=TRUE) & grepl("s__", taxa_counts$Full, fixed=TRUE)
      taxa_counts$Full[remove_indices] = taxa_counts$Taxa[remove_indices]
      remove_indices = grepl("(1)", taxa_counts$Full, fixed=TRUE) & grepl("g__", taxa_counts$Full, fixed=TRUE)
      taxa_counts$Full[remove_indices] = taxa_counts$Taxa[remove_indices]
      # remove the (1) from unknown if present
      remove_indices = grepl("(1)", taxa_counts$Full, fixed=TRUE) & grepl("Unknown", taxa_counts$Full, fixed=TRUE)
      taxa_counts$Full[remove_indices] = taxa_counts$Taxa[remove_indices]
    } else {
      df_pos$count = 1
      taxa_counts = aggregate(df_pos$count, by=list(df_pos$Taxa), sum)
      taxa_counts$left = "("
      taxa_counts$right = ")"
      taxa_counts$pasted = do.call(paste0, taxa_counts[c(1, 3, 2, 4)])
      taxa_counts$x = NULL
      taxa_counts$left = NULL
      taxa_counts$right = NULL
      names(taxa_counts) = c("Taxa","Full")
      df_pos$count = NULL
      # remove the (1) where only 1 species is present
      remove_indices = grepl("(1)", taxa_counts$Full, fixed=TRUE) & grepl("s__", taxa_counts$Full, fixed=TRUE)
      taxa_counts$Full[remove_indices] = taxa_counts$Taxa[remove_indices]
      remove_indices = grepl("(1)", taxa_counts$Full, fixed=TRUE) & grepl("g__", taxa_counts$Full, fixed=TRUE)
      taxa_counts$Full[remove_indices] = taxa_counts$Taxa[remove_indices]
      # remove the (1) from unknown if present
      remove_indices = grepl("(1)", taxa_counts$Full, fixed=TRUE) & grepl("Unknown", taxa_counts$Full, fixed=TRUE)
      taxa_counts$Full[remove_indices] = taxa_counts$Taxa[remove_indices]
    }

    # find the full name of the strongest contributor
    #print(taxa_counts$Taxa)
    #print(df_pos$Taxa)
    #print(strong_contributors)
    #print(df_pos_name_of_strongest)
    #print(which(taxa_counts$Taxa %in% df_pos_name_of_strongest))
    if (length(strong_contributors) > 0) {
      final_name_of_strong = taxa_counts$Full[which(taxa_counts$Taxa %in% df_pos_name_of_strongest)]
      #print(final_name_of_strong)
    }


    original_names = names(df_pos)
    df_pos = aggregate(df_pos[,-1,drop=FALSE], by=list(df_pos$Taxa), sum)
    names(df_pos) = original_names
    df_neg = aggregate(df_neg[,-1,drop=FALSE], by=list(df_neg$Taxa), sum)
    names(df_neg) = original_names

    if (separate_enriched_depleted_taxa) {
      df_pos$Shift = sub(".*:", "", df_pos$Taxa)
      df_pos$Taxa = sub(":ENRICHED", "", sub(":DEPLETED", "", df_pos$Taxa))
      df_pos = df_pos[, c("Taxa","Shift",functions)]
      original_names = names(df_pos)
      df_pos = merge(taxa_counts, df_pos, by="Taxa")
      df_pos$Taxa = NULL
      names(df_pos) = original_names

      df_neg$Shift = sub(".*:", "", df_neg$Taxa)
      df_neg$Taxa = sub(":ENRICHED", "", sub(":DEPLETED", "", df_neg$Taxa))
      df_neg = df_neg[, c("Taxa","Shift",functions)]
      df_neg = merge(taxa_counts, df_neg, by="Taxa")
      df_neg$Taxa = NULL
      names(df_neg) = original_names

    } else {
      df_pos = merge(taxa_counts, df_pos, by="Taxa")
      df_pos$Taxa = NULL
      names(df_pos) = original_names

      df_neg = merge(taxa_counts, df_neg, by="Taxa")
      df_neg$Taxa = NULL
      names(df_neg) = original_names
    }


  }

  # if we want to sort by the strongest, then sort cases and controls separately, according to it
  if ((length(strong_contributors) > 0) && (sort_by == "strongest")) {
    #print(names(df_pos))
    num_cases = sum(function_da[,2] > 0)
    sorted_names_by_strongest_for_cases = names(sort(df_pos[which(df_pos$Taxa %in% final_name_of_strong),
                                                            2:(num_cases+1)], decreasing=TRUE))

    if (num_cases+2 < length(names(df_pos))) {
      sorted_names_by_strongest_for_controls = names(sort(df_pos[which(df_pos$Taxa %in% final_name_of_strong),
                                                                 (num_cases+2):length(names(df_pos))], decreasing=TRUE))
    } else {
      sorted_names_by_strongest_for_controls = names(df_pos)[num_cases+2]
    }
    #print(sorted_names_by_strongest_for_cases)
    #print(sorted_names_by_strongest_for_controls)

    sorted_names_by_strongest_joined = c(sorted_names_by_strongest_for_cases, sorted_names_by_strongest_for_controls)

    #print(functions)
    #print(sorted_names_by_strongest_joined)
    functions = sorted_names_by_strongest_joined
    #print("sorted by strongest:")
    #print(functions)

    # reorder the columns
    #print(names(df_pos))
    #print(functions)
    df_pos = df_pos[, c("Taxa", functions)]
    df_neg = df_neg[, c("Taxa", functions)]
  }

  if (plot_type == "percentage_bars") {
    if (separate_enriched_depleted_taxa) {
      if (scale_pos_by_original_pos_sum) {
        df_pos[,-2] = sweep(df_pos[,-2,drop=FALSE], 2, original_sum_of_pos_contributions$PosSum, FUN='/')
      } else {
        df_pos[,-2] = sweep(df_pos[,-2,drop=FALSE], 2, colSums(df_pos[,-2,FALSE]), FUN='/')
      }
      if (scale_neg_by_original_neg_sum) {
        df_neg[,-2] = sweep(df_neg[,-2,drop=FALSE], 2, original_sum_of_neg_contributions$NegSum, FUN='/')
      } else {
        df_neg[,-2] = -1 * sweep(df_neg[,-2,drop=FALSE], 2, colSums(df_neg[,-2,FALSE]), FUN='/')
      }
    } else {
      if (scale_pos_by_original_pos_sum) {
        df_pos[,-1] = sweep(df_pos[,-1,drop=FALSE], 2, original_sum_of_pos_contributions$PosSum, FUN='/')
      } else {
        df_pos[,-1] = sweep(df_pos[,-1,drop=FALSE], 2, colSums(df_pos[,-1,FALSE]), FUN='/')
      }
      if (scale_neg_by_original_neg_sum) {
        df_neg[,-1] = sweep(df_neg[,-1,drop=FALSE], 2, original_sum_of_neg_contributions$NegSum, FUN='/')
      } else {
        df_neg[,-1] = -1 * sweep(df_neg[,-1,drop=FALSE], 2, colSums(df_neg[,-1,FALSE]), FUN='/')
      }
    }

  }



  ########################
  # COLORS!!!
  # set the color pallete
  ########################
  if (color_small_cont_by_phyla) {
    num_of_final_phyla = sum(!grepl('\\s',unique(df_pos$Taxa), perl=TRUE))
  } else {
    num_of_final_phyla = length(unique(gsub(" g__.*","",grep("p__.* g__",df_pos$Taxa, value=TRUE))))
  }
  num_of_final_taxa = length(unique(df_pos$Taxa))
  #print(num_of_final_phyla)
  #print(num_of_final_taxa)
  cPalette = data.frame(unique(df_pos$Taxa))
  names(cPalette) = "Taxa"
  cPalette$h = 0
  cPalette$s = 0
  cPalette$v = 0

  # the 5 major phyla in our data
  num_Actinobacteria = length(grep("Actinobacteria", cPalette$Taxa))
  num_Bacteroidetes = length(grep("Bacteroidetes", cPalette$Taxa))
  num_Firmicutes = length(grep("Firmicutes", cPalette$Taxa))
  num_Proteobacteria = length(grep("Proteobacteria", cPalette$Taxa))
  num_Fusobacteria = length(grep("Fusobacteria", cPalette$Taxa))

  number_of_other_phyla = num_of_final_phyla - ((num_Actinobacteria > 0) + (num_Bacteroidetes > 0) + (num_Firmicutes > 0) +
    (num_Proteobacteria > 0) + (num_Fusobacteria > 0))

  # Actinobacteria = green_pallete
  cPalette[grep("Actinobacteria", cPalette$Taxa), -1] = expand.grid(h=0.4, s=seq(0.3,1,length.out=num_Actinobacteria), v=0.9)

  # Fusobacteria = orange_pallete
  cPalette[grep("Fusobacteria", cPalette$Taxa), -1] = expand.grid(h=0.2, s=seq(0.3,1,length.out=num_Fusobacteria), v=0.9)

  # Bacteroidetes = purple_pallete
  cPalette[grep("Bacteroidetes", cPalette$Taxa), -1] = expand.grid(h=0.8, s=seq(0.3,1,length.out=num_Bacteroidetes), v=0.9)

  # Firmicutes = blue_pallete
  cPalette[grep("Firmicutes", cPalette$Taxa), -1] = expand.grid(h=0.6, s=seq(0.3,1,length.out=num_Firmicutes), v=0.9)

  # Proteobacteria = red_pallete
  cPalette[grep("Proteobacteria", cPalette$Taxa), -1] = expand.grid(h=0, s=seq(0.3,1,length.out=num_Proteobacteria), v=0.9)

  #print(cPalette)
  #print(number_of_other_phyla)
  if (number_of_other_phyla > 0){
    # get list of other pyhla
    taxa_from_other_phyla = cPalette$Taxa[!grepl("Actinobacteria", cPalette$Taxa) & !grepl("Fusobacteria", cPalette$Taxa) & !grepl("Bacteroidetes", cPalette$Taxa) & !grepl("Firmicutes", cPalette$Taxa) & !grepl("Proteobacteria", cPalette$Taxa)]
    #print(taxa_from_other_phyla)
    # extract only phyla names from the taxa names
    other_phyla_names_uniq = unique(gsub(" .*", "",  gsub('\\d+',"",  gsub('\\)', "", gsub('\\(', "", taxa_from_other_phyla)))))
    # set their colors....
    if (number_of_other_phyla <= 5) { # if < 5, then we can fit them in nicely in the color wheel
      other_hues = seq(from = 0.1, to=0.9, by=0.2)
      for (i in 1:length(other_phyla_names_uniq)) {
        curr_num = length(grep(other_phyla_names_uniq[i], cPalette$Taxa))
        cPalette[grep(other_phyla_names_uniq[i], cPalette$Taxa), -1] =
          expand.grid(h=other_hues[i], s=seq(0.3,1,length.out=curr_num), v=0.9)
      }
    } else if (number_of_other_phyla <= 15) { # if < 5, then we can fit them in nicely in the color wheel
      other_hues = seq(from = 0, to=1, by=0.05)[c(2,3,4,6,7,8,10,11,12,14,15,16,18,19,20)]
      for (i in 1:length(other_phyla_names_uniq)) {
        curr_num = length(grep(other_phyla_names_uniq[i], cPalette$Taxa))
        cPalette[grep(other_phyla_names_uniq[i], cPalette$Taxa), -1] =
          expand.grid(h=other_hues[i], s=seq(0.3,1,length.out=curr_num), v=0.9)
      }
    } else {
      other_hues = seq(from = 0, to=1, length.out=number_of_other_phyla)
      for (i in 1:length(other_phyla_names_uniq)) {
        curr_num = length(grep(other_phyla_names_uniq[i], cPalette$Taxa))
        cPalette[grep(other_phyla_names_uniq[i], cPalette$Taxa), -1] =
          expand.grid(h=other_hues[i], s=seq(0.3,1,length.out=curr_num), v=0.9)
      }
    }
  }

  is_other = length(grep("ZOther taxa", cPalette$Taxa))
  is_unknown = length(grep("Unknown", cPalette$Taxa))
  if (is_other > 0) {
    cPalette[grep("ZOther taxa", cPalette$Taxa), -1] = c(0,0,0.5)
  }
  if (is_unknown > 0) {
    cPalette[grep("Unknown", cPalette$Taxa), -1] = c(0,0,0)
  }

  #print(cPalette)
  #print(df_pos)

  # if we have DA information for taxa, than add this as another bar at the bottom
  if (add_taxa_da && !is.null(input_taxa_da) && (plot_type != "percent_unknown")
      && (plot_type != "percent_abundance_explained")
      && (plot_type != "percent_abundance_explained_pearson")
      && (plot_type != "percent_unknown_no_negative")
      && (plot_type != "percent_predicted")
      && (plot_type != "percent_predicted_diamonds")
      && (plot_type != "percent_unknown_distribution")
      && (plot_type != "percent_unknown_cumulative")
      && (plot_type != "function_DA_bars")
      && (plot_type != "predicted_DA_bars")
      && (plot_type != "predicted_and_function_DA_bars")
      && (plot_type != "percentage_bars")) {
    #print(taxa_da[,1:2])
    taxa_da_bar = aggregate(taxa_da[,4], by=list(taxa_da$Taxa), sum)
    #print(taxa_da_bar[,1:2])
    names(taxa_da_bar) = c("Taxa","StatValue")
    taxa_da_bar$StatValue[which(taxa_da_bar$Taxa == "ZOther taxa")] = 0
    taxa_da_bar = merge(taxa_counts, taxa_da_bar, by="Taxa")
    #print(taxa_da_bar[,1:2])
    taxa_da_bar$Taxa = NULL
    names(taxa_da_bar) = c("Taxa","StatValue")
    pos_ind = which(taxa_da_bar$StatValue > 0)
    neg_ind = which(taxa_da_bar$StatValue < 0)
    taxa_da_bar_pos = taxa_da_bar
    taxa_da_bar_neg = taxa_da_bar
    taxa_da_bar_pos$StatValue[neg_ind] = 0
    taxa_da_bar_neg$StatValue[pos_ind] = 0

    #print(taxa_da_bar_pos[,1:2])
    #print(df_pos[,1:2])
    df_pos$Taxa_Differential_abundance = c(taxa_da_bar_pos$StatValue / 10)
    df_neg$Taxa_Differential_abundance = c(taxa_da_bar_neg$StatValue / 10)
  }

  if (separate_enriched_depleted_taxa) {
    if (!use_facets_for_separation) { 
      df_pos.m = reshape2::melt(df_pos, id.vars=c(1,2))
      df_pos.m$Taxa = factor(df_pos.m$Taxa, levels=unique(df_pos.m$Taxa))
      # reverse the order of the functions because faceting reverses the order
      df_pos.m$variable = factor(df_pos.m$variable, levels=(levels(df_pos.m$variable)))
  
      df_neg.m = reshape2::melt(df_neg, id.vars=c(1,2))
      df_neg.m$Taxa = factor(df_neg.m$Taxa, levels=unique(df_neg.m$Taxa))
      # reverse the order of the functions because faceting reverses the order
      df_neg.m$variable = factor(df_neg.m$variable, levels=(levels(df_neg.m$variable)))
    } else {
      df_pos.m = reshape2::melt(df_pos, id.vars=c(1,2))
      df_pos.m$Taxa = factor(df_pos.m$Taxa, levels=unique(df_pos.m$Taxa))
      # reverse the order of the functions because faceting reverses the order
      df_pos.m$variable = factor(df_pos.m$variable, levels=rev(levels(df_pos.m$variable)))
      
      df_neg.m = reshape2::melt(df_neg, id.vars=c(1,2))
      df_neg.m$Taxa = factor(df_neg.m$Taxa, levels=unique(df_neg.m$Taxa))
      # reverse the order of the functions because faceting reverses the order
      df_neg.m$variable = factor(df_neg.m$variable, levels=rev(levels(df_neg.m$variable)))
    }

  } else {
    df_pos.m = reshape2::melt(df_pos)
    df_pos.m$Taxa = factor(df_pos.m$Taxa, levels=unique(df_pos.m$Taxa))

    df_neg.m = reshape2::melt(df_neg)
    df_neg.m$Taxa = factor(df_neg.m$Taxa, levels=unique(df_neg.m$Taxa))
  }

  if (verbose) {
    print("Finished Preparing data for plots")
    print(head(df_pos.m))
    print(head(df_neg.m))
    print(cPalette)
  }
  
  # if needed, clean the taxa names a bit
  if (!add_phyla_names_to_species_label) {
    df_pos.m$Taxa = as.character(df_pos.m$Taxa)
    df_pos.m$Taxa[grepl("g__", df_pos.m$Taxa)] = gsub("s__","",
                                                      gsub("g__","",
                                                           gsub("p__.*g__","g__",
                                                                df_pos.m$Taxa[grepl("g__", df_pos.m$Taxa)])))
    
    df_pos.m$Taxa = factor(df_pos.m$Taxa, levels=unique(df_pos.m$Taxa))
    
    df_neg.m$Taxa = as.character(df_neg.m$Taxa)
    df_neg.m$Taxa[grepl("g__", df_neg.m$Taxa)] = gsub("s__","",
                                                      gsub("g__","",
                                                           gsub("p__.*g__","g__",
                                                                df_neg.m$Taxa[grepl("g__", df_neg.m$Taxa)])))
    
    df_neg.m$Taxa = factor(df_neg.m$Taxa, levels=unique(df_neg.m$Taxa))
    
  }
  
  # remove the Z from the "ZOther Taxa" (used for ordering as last)
  if (sum(grepl("ZOther", df_pos.m$Taxa)) > 0) {
    df_pos.m$Taxa = as.character(df_pos.m$Taxa)
    df_pos.m$Taxa[grepl("ZOther", df_pos.m$Taxa)] = gsub("Z","",df_pos.m$Taxa[grepl("ZOther", df_pos.m$Taxa)])  
    df_pos.m$Taxa = factor(df_pos.m$Taxa, levels=unique(df_pos.m$Taxa))
    
    df_neg.m$Taxa = as.character(df_neg.m$Taxa)
    df_neg.m$Taxa[grepl("ZOther", df_neg.m$Taxa)] = gsub("Z","",df_neg.m$Taxa[grepl("ZOther", df_neg.m$Taxa)])    
    df_neg.m$Taxa = factor(df_neg.m$Taxa, levels=unique(df_neg.m$Taxa))
  }
  
  final_list_of_taxa_to_plot = unique(c(as.character(df_pos.m$Taxa), as.character(df_neg.m$Taxa)))
  
  
  ############################################
  # FINISHED PREPARING DATA FOR PLOT
  ############################################


  
  
  
  
  
  

  ############################################
  # BAR PLOT
  ############################################
  if (verbose) {
    print("Plotting...")
  }

  BarPlot = ggplot()

  if (separate_enriched_depleted_taxa) {
    
    if (!use_facets_for_separation) { 
      no_facet_df_pos.m = df_pos.m
      no_facet_df_pos.m$x_val = 0
      for (i in 1:length(levels(df_pos.m$variable))) {
        no_facet_df_pos.m$x_val[no_facet_df_pos.m$variable == levels(df_pos.m$variable)[i]] = i
      }
      # add/subtract 0.5 for enriched/depleted
      no_facet_df_pos.m$x_val[no_facet_df_pos.m$Shift == "ENRICHED"] = 
        no_facet_df_pos.m$x_val[no_facet_df_pos.m$Shift == "ENRICHED"] + 0.175
      
      no_facet_df_pos.m$x_val[no_facet_df_pos.m$Shift == "DEPLETED"] = 
        no_facet_df_pos.m$x_val[no_facet_df_pos.m$Shift == "DEPLETED"] - 0.175
  
      no_facet_df_neg.m = df_neg.m
      no_facet_df_neg.m$x_val = 0
      for (i in 1:length(levels(df_neg.m$variable))) {
        no_facet_df_neg.m$x_val[no_facet_df_neg.m$variable == levels(df_neg.m$variable)[i]] = i
      }
      # add/subtract 0.5 for enriched/depleted
      no_facet_df_neg.m$x_val[no_facet_df_neg.m$Shift == "ENRICHED"] = 
        no_facet_df_neg.m$x_val[no_facet_df_neg.m$Shift == "ENRICHED"] + 0.175
      
      no_facet_df_neg.m$x_val[no_facet_df_neg.m$Shift == "DEPLETED"] = 
        no_facet_df_neg.m$x_val[no_facet_df_neg.m$Shift == "DEPLETED"] - 0.175
  
      if (!show_only_neg) {
        BarPlot = BarPlot +
          # A hack to hide the slashes: first graph the bars with no outline and add the legend,
          # then graph the bars again with outline, but with a blank legend.
          geom_bar(data=no_facet_df_pos.m, aes(x=x_val, fill=Taxa, y=value), stat = "identity", width=bar_width) +
          geom_bar(data=no_facet_df_pos.m, aes(x=x_val, fill=Taxa, y=value), stat = "identity",
                   colour="black", width=bar_width, show_guide=FALSE)
      }
      
      if (!show_only_pos) {
        BarPlot = BarPlot +
          geom_bar(data=no_facet_df_neg.m, aes(x=x_val, fill=Taxa, y=value), stat = "identity", width=bar_width) +
          geom_bar(data=no_facet_df_neg.m, aes(x=x_val, fill=Taxa, y=value), stat = "identity",
                   colour="black", width=bar_width, show_guide=FALSE)
      }
      
      # set the x-axis labels
      BarPlot = BarPlot + scale_x_continuous(breaks=1:length(levels(df_pos.m$variable)), labels=levels(df_pos.m$variable))
    
    } else { # use facets
      
      if (!show_only_neg) {
        BarPlot = BarPlot +
          # A hack to hide the slashes: first graph the bars with no outline and add the legend,
          # then graph the bars again with outline, but with a blank legend.
          geom_bar(data=df_pos.m, aes(x=Shift, fill=Taxa, y=value), stat = "identity", width=bar_width) +
          geom_bar(data=df_pos.m, aes(x=Shift, fill=Taxa, y=value), stat = "identity",
                   colour="black", width=bar_width, show_guide=FALSE)
      }
      
      if (!show_only_pos) {
        BarPlot = BarPlot +
          geom_bar(data=df_neg.m, aes(x=Shift, fill=Taxa, y=value), stat = "identity", width=bar_width) +
          geom_bar(data=df_neg.m, aes(x=Shift, fill=Taxa, y=value), stat = "identity",
                   colour="black", width=bar_width, show_guide=FALSE)
      }
      
      BarPlot = BarPlot +
        geom_abline(data=df_pos.m, intercept = 0, slope = 0, size = 1) +
        scale_y_continuous(breaks = round(seq(min(original_sum_of_neg_contributions$NegSum),
                                              max(original_sum_of_pos_contributions$PosSum), by = 1),1)) +
        facet_wrap( ~ variable, ncol=1) +
        theme(panel.margin=unit(0, "cm"),
              plot.title = element_text(face="bold", size=20),
              axis.title = element_text(face="bold", size=20),
              axis.text.x = element_text(face="bold", size=15),
              panel.grid=element_blank(),
              panel.border=element_blank(),
              axis.line = element_line(colour = 'black', size = 1))
      if (!add_facet_labels) {
        BarPlot = BarPlot + theme(strip.text = element_blank(),
                                  strip.background = element_rect(fill = "transparent",colour = NA))
      } else { 
        BarPlot = BarPlot + theme(strip.text = element_text(size=20))
      }
    }
    
  } else { # don't separate enriched/depleted taxa

    if (!show_only_neg) {
      BarPlot = BarPlot +
        # A hack to hide the slashes: first graph the bars with no outline and add the legend,
        # then graph the bars again with outline, but with a blank legend.
        geom_bar(data=df_pos.m, aes(x=variable, fill=Taxa, y=value), stat = "identity", width=bar_width) +
        geom_bar(data=df_pos.m, aes(x=variable, fill=Taxa, y=value), stat = "identity",
                 colour="black", width=bar_width, show_guide=FALSE)
    }

    if (!show_only_pos) {
      BarPlot = BarPlot +
        geom_bar(data=df_neg.m, aes(x=variable, fill=Taxa, y=value), stat = "identity", width=bar_width) +
        geom_bar(data=df_neg.m, aes(x=variable, fill=Taxa, y=value), stat = "identity",
                 colour="black", width=bar_width, show_guide=FALSE)
    }

    BarPlot = BarPlot + geom_abline(intercept = 0, slope = 0, size = 1) +
      theme(plot.title = element_text(face="bold", size=20)) +
      theme(axis.title = element_text(face="bold", size=20)) +
      theme(axis.text.x = element_text(face="bold", size=15))
  }

  if (verbose) {
    print("Made basic plot, now flipping coordinates if needed...")
  }

  # flip coordinates
  if (flip_coord) {
    if (plot_type == "percentage_bars") {
      if (show_only_pos) {
        BarPlot = BarPlot + coord_flip(ylim = c(0, 1), xlim=c(0.0, num_of_functions+0.6))
      } else if (show_only_neg) {
        BarPlot = BarPlot + coord_flip(ylim = c(-1, 0))
      } else {
        BarPlot = BarPlot + coord_flip(ylim = c(-1, 1))
      }
    } else { # normal bars
      if (separate_enriched_depleted_taxa && !use_facets_for_separation) {
        BarPlot = BarPlot + coord_flip(xlim=c(0.5,num_of_functions+0.5))
      } else {
        BarPlot = BarPlot + coord_flip()
      }    
    }
   
  } 

  # add line emphesizing the zero
  BarPlot = BarPlot + geom_abline(intercept = 0, slope = 0, size = 1, colour="black")
  
  if (verbose) {
    print("Changing colors to my pallete...")
  }

  BarPlot = BarPlot +
    ylab(paste(input_score, "test statistic")) +
    xlab("Function") +
    scale_fill_manual(values=apply(cPalette[,-1], 1, function(x) hsv(x[1],x[2],x[3]))) +
    ggtitle(paste("#Taxa:",num_of_taxa,"#Func:",num_of_functions,"Taxa contribution to",input_score,"when assessing with",input_permutation)) 


  # add line separating cases and controls
  num_cases = sum(function_da[,2] > 0)
  if (num_cases > 0) {
    if (add_case_control_line) {
      if (verbose) {
        print("Adding case/control line...")
      }
      BarPlot = BarPlot + geom_vline(xintercept = num_cases + 0.5, size = 1, linetype = "dashed")
    }
  }

  if (plot_type == "percentage_bars") {
    if (verbose) {
      print("Changing to percentage bars...")
    }
    BarPlot = BarPlot +
      ylab(paste(input_score, "test statistic fraction")) +
      scale_y_continuous(breaks=seq(-1,1,0.1))
  }

  # if we have DA information for functions, add markers (red diamonds)
  if (!is.null(input_original) && (plot_type != "percentage_bars")) {
    if (verbose) {
      print("Adding original DA markers...")
    }
    # add markers
    rownames(function_da) = function_da$KO
    function_da$KO = NULL
    if (!is.null(input_predicted_da)) {
      rownames(predicted_da) = predicted_da$KO
      predicted_da$KO = NULL
    }
    if (!is.null(input_predicted_abundance_agreement)) {
      rownames(predicted_abundance_agreement) = predicted_abundance_agreement$KO
      predicted_abundance_agreement$KO = NULL
    }

    if (!is.null(input_function_meta)) {
      function_da = merge(meta, function_da, by="row.names")[, 2:3]
      rownames(function_da) = function_da$V2
      function_da$V2 = NULL
      if (!is.null(input_predicted_da)) {
        predicted_da = merge(meta, predicted_da, by="row.names")[, 2:3]
        rownames(predicted_da) = predicted_da$V2
        predicted_da$V2 = NULL
      }
      if (!is.null(input_predicted_abundance_agreement)) {
        predicted_abundance_agreement = merge(meta, predicted_abundance_agreement, by="row.names")[, -1]
        rownames(predicted_abundance_agreement) = predicted_abundance_agreement$V2
        predicted_abundance_agreement$V2 = NULL
      }
    }
    #print(functions)
    #print(function_da)
    function_da$y_vals = function_da[functions, input_score]
    #print(function_da)
    function_da$x_val = 1:length(function_da$y_vals)
    function_da[input_score] = NULL
    rownames(function_da) = functions
    if (!is.null(input_predicted_da)) {
      predicted_da$x_vals = predicted_da[functions, input_score]
      predicted_da[input_score] = NULL
      rownames(predicted_da) = functions
    }

    if (add_taxa_da && !is.null(input_taxa_da)
        && (plot_type != "percent_unknown")
        && (plot_type != "percent_abundance_explained")
        && (plot_type != "percent_abundance_explained_pearson")
        && (plot_type != "percent_predicted")
        && (plot_type != "percent_predicted_diamonds")
        && (plot_type != "percent_unknown_no_negative")
        && (plot_type != "percent_unknown_distribution")
        && (plot_type != "percent_unknown_cumulative")
        && (plot_type != "taxa_corr_heatmap")
        && (plot_type != "function_DA_bars")
        && (plot_type != "predicted_and_function_DA_bars")
        && (plot_type != "predicted_DA_bars")) {
      function_da = rbind(function_da, c(0, length(function_da$y_vals)+1))
      BarPlot = BarPlot + geom_vline(xintercept = length(function_da$y_vals) - 0.5, size = 1, linetype = "longdash")
    }

    if (add_original_da_markers) {
      BarPlot = BarPlot + geom_point(data=function_da, aes(x_val, y_vals), colour="black", fill="red", size=4, shape=23)
    }
  }

  # add markers of the predicted DA of the pseudo-metagenome
  # which for shapley will be equal to the sumof positive and negative
  # white diamonds
  if (plot_type != "percentage_bars" && add_predicted_da_markers) {
    if (verbose) {
      print("Adding predicted DA markers...")
    }
    
    if (separate_enriched_depleted_taxa) {
      if (!use_facets_for_separation) {
        markers_of_predicted_da = predicted_da
        names(markers_of_predicted_da) = "y_vals"
        markers_of_predicted_da$x_val = 0
        for (i in 1:length(levels(df_pos.m$variable))) {
          markers_of_predicted_da$x_val[rownames(markers_of_predicted_da) == levels(df_neg.m$variable)[i]] = i
        }      
      } else { # use facets
        markers_of_predicted_da = predicted_da[levels(df_pos.m$variable), ,drop=FALSE]
        markers_of_predicted_da$variable = levels(df_pos.m$variable)
        names(markers_of_predicted_da) = c("y_vals", "variable")
        markers_of_predicted_da$x_val = 1.5
      }
    
    } else { # don't separate enriched/depleted taxa 
      markers_of_predicted_da = predicted_da
      names(markers_of_predicted_da) = "y_vals"
      markers_of_predicted_da$x_val = row.names(markers_of_predicted_da)
    }

    BarPlot = BarPlot +
      geom_point(data=markers_of_predicted_da, aes(x_val, y_vals), colour="black", fill="white", size=4, shape=23)
  }

  # if needed, add the taxa names in bars (only the last name in each taxa)
  if (add_names_in_bars) {
    if (separate_enriched_depleted_taxa) {
      
      if (!use_facets_for_separation) {
        
        if (!show_only_neg) {
          pos_labels = ggplot_build(BarPlot)$data[[1]]
          pos_labels$label = levels(df_pos.m$Taxa)[pos_labels$group]
          pos_labels$label_x = pos_labels$x
          pos_labels$label_y = (pos_labels$ymax + pos_labels$ymin) / 2
          pos_labels$plot_label = pos_labels$ymax-pos_labels$ymin
          pos_labels = pos_labels[pos_labels$plot_label > 0.01, c("label","label_x","label_y")]
          pos_labels$label = gsub(".* ","", pos_labels$label)
          pos_labels = pos_labels[!grepl("taxa",pos_labels$label),]
          BarPlot = BarPlot + geom_text(data=pos_labels, aes(x=label_x, y=label_y, label=label),
                                        stat = "identity", angle = 90, size = 3)
        }
        
        if (!show_only_pos) {
          neg_labels = ggplot_build(BarPlot)$data[[3]]
          neg_labels$label = levels(df_neg.m$Taxa)[neg_labels$group]
          neg_labels$label_x = neg_labels$x
          neg_labels$label_y = (neg_labels$ymax + neg_labels$ymin) / 2
          neg_labels$plot_label = neg_labels$ymax-neg_labels$ymin
          neg_labels = neg_labels[neg_labels$plot_label < -0.01, c("label","label_x","label_y")]
          neg_labels$label = gsub(".* ","", neg_labels$label)
          neg_labels = neg_labels[!grepl("taxa",neg_labels$label),]
          BarPlot = BarPlot + geom_text(data=neg_labels, aes(x=label_x, y=label_y, label=label),
                                        stat = "identity", angle = 90, size = 3)
        }       
    
      } else { # use facets
        
        if (!show_only_neg) {
          pos_labels = ggplot_build(BarPlot)$data[[1]]
          pos_labels$variable = rev(df_pos.m$variable)
          taxa_vs_color = cbind(cPalette, apply(cPalette[,-1], 1, function(x) hsv(x[1],x[2],x[3])))[,c(1,5)]
          names(taxa_vs_color) = c("Taxa","fill")
          pos_labels$label_x = pos_labels$x
          pos_labels$label_y = (pos_labels$ymax + pos_labels$ymin) / 2
          pos_labels$plot_label = pos_labels$ymax-pos_labels$ymin
          pos_labels = merge(taxa_vs_color, pos_labels, by="fill")
          pos_labels = pos_labels[pos_labels$plot_label > 0.01, c("Taxa","label_x","label_y","variable")]
          pos_labels$Taxa = gsub(".* ","", pos_labels$Taxa)
          BarPlot = BarPlot + geom_text(data=pos_labels, aes(x=label_x, y=label_y, label=Taxa),
                                        stat = "identity", angle = 90, size = 3)
        } 

        if (!show_only_pos) {
          neg_labels = ggplot_build(BarPlot)$data[[3]]
          neg_labels$variable = rev(df_neg.m$variable)
          taxa_vs_color = cbind(cPalette, apply(cPalette[,-1], 1, function(x) hsv(x[1],x[2],x[3])))[,c(1,5)]
          names(taxa_vs_color) = c("Taxa","fill")
          neg_labels$label_x = neg_labels$x
          neg_labels$label_y = (neg_labels$ymax + neg_labels$ymin) / 2
          neg_labels$plot_label = neg_labels$ymax-neg_labels$ymin
          neg_labels = merge(taxa_vs_color, neg_labels, by="fill")
          neg_labels = neg_labels[neg_labels$plot_label < -0.01, c("Taxa","label_x","label_y","variable")]
          neg_labels$Taxa = gsub(".* ","", neg_labels$Taxa)
          BarPlot = BarPlot + geom_text(data=neg_labels, aes(x=label_x, y=label_y, label=Taxa),
                                        stat = "identity", angle = 90, size = 3)
  
        }
    
      }
      
    } else {  # don't separate enriched/depleted taxa
      if (!show_only_neg) {
        pos_labels = ggplot_build(BarPlot)$data[[1]]
        pos_labels$label = df_pos$Taxa
        pos_labels$label_x = pos_labels$x
        pos_labels$label_y = (pos_labels$ymax + pos_labels$ymin) / 2
        pos_labels$plot_label = pos_labels$ymax-pos_labels$ymin
        pos_labels = pos_labels[pos_labels$plot_label > 0.01, c("label","label_x","label_y")]
        pos_labels$label = gsub(".* ","", pos_labels$label)
        BarPlot = BarPlot + geom_text(data=pos_labels, aes(x=label_x, y=label_y, label=label),
                                      stat = "identity", angle = 90, size = 3)
      }
      if (!show_only_pos) {
        neg_labels = ggplot_build(BarPlot)$data[[2]]
        neg_labels$label = df_pos$Taxa
        neg_labels$label_x = neg_labels$x
        neg_labels$label_y = (neg_labels$ymax + neg_labels$ymin) / 2
        neg_labels$plot_label = neg_labels$ymax-neg_labels$ymin
        neg_labels = neg_labels[neg_labels$plot_label < -0.01, c("label","label_x","label_y")]
        neg_labels$label = gsub(".* ","", neg_labels$label)
        BarPlot = BarPlot + geom_text(data=neg_labels, aes(x=label_x, y=label_y, label=label),
                                      stat = "identity", angle = 90, size = 3)
      }
    }

  }

  ############################################
  # END OF BAR PLOT
  ############################################

  
  
  
  
  
  
  

  #################################################################################################
  #################################################################################################
  # OTHER PLOTS
  #################################################################################################
  #################################################################################################


  #################################################################################################
  # plot the ratios of pos sum to predicted DA, and also the ratio between the enriched pos to
  # depelted pos
  #################################################################################################
  if (plot_type == "pos_neg_ratios") {
    if (verbose) {
      print("Plotting pos_neg_ratios plot...")
    }
    if (separate_enriched_depleted_taxa) {

      # POS for ENRICHED FUNCTIONS
      ratio_pos_predicted = merge(data.frame(apply(df_pos[,c(-1, -2)], 2, sum)),
                                  predicted_da[predicted_da[,1] > 0, ,drop=FALSE], by="row.names")
      row.names(ratio_pos_predicted) = ratio_pos_predicted[, 1]
      ratio_pos_predicted[, 1] = NULL
      ratio_pos_predicted$ratio = ratio_pos_predicted[, 1] / ratio_pos_predicted[, 2]
      names(ratio_pos_predicted) = c("PosSum","Predicted","RatioPosSumPredicted")

      ratio_pos_enriched_depleted = merge(data.frame(apply(df_pos[df_pos$Shift == "ENRICHED", c(-1, -2)], 2, sum)),
            data.frame(apply(df_pos[df_pos$Shift == "DEPLETED", c(-1, -2)], 2, sum)), by="row.names")
      row.names(ratio_pos_enriched_depleted) = ratio_pos_enriched_depleted[, 1]
      ratio_pos_enriched_depleted[, 1] = NULL
      ratio_pos_enriched_depleted$ratio = ratio_pos_enriched_depleted[, 1] / ratio_pos_enriched_depleted[, 2]
      names(ratio_pos_enriched_depleted) = c("PosSumEnriched","PosSumDepleted","RatioPosEnrichedDepleted")
      ratio_pos_enriched_depleted = merge(ratio_pos_enriched_depleted, predicted_da[predicted_da[,1] > 0, ,drop=FALSE], by="row.names")
      row.names(ratio_pos_enriched_depleted) = ratio_pos_enriched_depleted[, 1]
      ratio_pos_enriched_depleted[, 1] = NULL

      if (!is.null(input_function_stats)) {
        ratio_pos_enriched_depleted = merge(ratio_pos_enriched_depleted, function_stats, by="row.names")
        row.names(ratio_pos_enriched_depleted) = ratio_pos_enriched_depleted[, 1]
        ratio_pos_enriched_depleted[, 1] = NULL
      }

      ratio_pos_data = merge(ratio_pos_predicted, ratio_pos_enriched_depleted, by="row.names")
      row.names(ratio_pos_data) = ratio_pos_data[, 1]
      ratio_pos_data[, 1] = NULL

      ScatterPlot1 = ggplot(ratio_pos_predicted, aes(x=RatioPosSumPredicted, y=Predicted)) +
        geom_point(size=5) +
        ggtitle(paste("R =", cor(ratio_pos_predicted$Predicted, ratio_pos_predicted$RatioPosSumPredicted)))

      ScatterPlot2 = ggplot(ratio_pos_data, aes(x=RatioPosEnrichedDepleted, y=Predicted)) +
        geom_point(size=5) +
        geom_vline(xintercept=1, size=2, linetype="dashed") +
        ggtitle(paste("R =", cor(ratio_pos_data$Predicted, ratio_pos_data$RatioPosEnrichedDepleted)))

      # NEG for DEPLETED FUNCTIONS
      ratio_neg_predicted = merge(data.frame(apply(df_neg[,c(-1, -2)], 2, sum)),
                                  predicted_da[predicted_da[,1] < 0, ,drop=FALSE], by="row.names")
      row.names(ratio_neg_predicted) = ratio_neg_predicted[, 1]
      ratio_neg_predicted[, 1] = NULL
      ratio_neg_predicted$ratio = ratio_neg_predicted[, 1] / ratio_neg_predicted[, 2]
      names(ratio_neg_predicted) = c("NegSum","Predicted","RatioNegSumPredicted")

      ratio_neg_enriched_depleted = merge(data.frame(apply(df_neg[df_neg$Shift == "ENRICHED", c(-1, -2)], 2, sum)),
                                          data.frame(apply(df_neg[df_neg$Shift == "DEPLETED", c(-1, -2)], 2, sum)), by="row.names")
      row.names(ratio_neg_enriched_depleted) = ratio_neg_enriched_depleted[, 1]
      ratio_neg_enriched_depleted[, 1] = NULL
      ratio_neg_enriched_depleted$ratio = ratio_neg_enriched_depleted[, 2] / ratio_neg_enriched_depleted[, 1]
      names(ratio_neg_enriched_depleted) = c("NegSumEnriched","NegSumDepleted","RatioNegDepletedEnriched")
      ratio_neg_enriched_depleted = merge(ratio_neg_enriched_depleted, predicted_da[predicted_da[,1] < 0, ,drop=FALSE], by="row.names")
      row.names(ratio_neg_enriched_depleted) = ratio_neg_enriched_depleted[, 1]
      ratio_neg_enriched_depleted[, 1] = NULL

      if (!is.null(input_function_stats)) {
        ratio_neg_enriched_depleted = merge(ratio_neg_enriched_depleted, function_stats, by="row.names")
        row.names(ratio_neg_enriched_depleted) = ratio_neg_enriched_depleted[, 1]
        ratio_neg_enriched_depleted[, 1] = NULL
      }

      ratio_neg_data = merge(ratio_neg_predicted, ratio_neg_enriched_depleted, by="row.names")
      row.names(ratio_neg_data) = ratio_neg_data[, 1]
      ratio_neg_data[, 1] = NULL

      ScatterPlot3 = ggplot(ratio_neg_predicted, aes(x=RatioNegSumPredicted, y=Predicted)) +
        geom_point(size=5) +
        ggtitle(paste("R =", cor(ratio_neg_predicted$Predicted, ratio_neg_predicted$RatioNegSumPredicted)))

      ScatterPlot4 = ggplot(ratio_neg_data, aes(x=RatioNegDepletedEnriched, y=Predicted)) +
        geom_point(size=5) +
        geom_vline(xintercept=1, size=2, linetype="dashed") +
        ggtitle(paste("R =", cor(ratio_neg_data$Predicted, ratio_neg_data$RatioNegDepletedEnriched)))

      ScatterPlot = arrangeGrob(ScatterPlot1, ScatterPlot2, ScatterPlot3, ScatterPlot4, ncol = 2)

    } else {
      print("Error: pos_neg_ratios can only be used with separate_enriched_depleted_taxa=TRUE")
      return(NULL)
    }
  }

  #################################################################################################
  # plot only the function DA scores as bars
  #################################################################################################
  if (plot_type == "function_DA_bars" || plot_type == "predicted_DA_bars") {
    function_da$func_name = row.names(function_da)
    predicted_da$x_val = function_da$x_val
    predicted_da$y_vals = predicted_da$x_vals
    if (plot_type == "function_DA_bars") {
      data_to_plot =  function_da 
    } else {
      data_to_plot =  predicted_da 
    }
    
    DABarPlot = ggplot() +
      geom_bar(data=data_to_plot, aes(x=x_val, y=y_vals), stat = "identity", colour="black", width=bar_width)
    
    if (add_predicted_da_markers) {
      DABarPlot = DABarPlot +   
        geom_point(data=predicted_da, aes(x=x_val, y=y_vals), colour="black", fill="white", size=4, shape=23)
    }

    if (add_original_da_markers) {
      DABarPlot = DABarPlot +   
        geom_point(data=function_da, aes(x=x_val, y=y_vals), colour="black", fill="red", size=4, shape=23)
    }
        
    DABarPlot = DABarPlot +  
      geom_abline(intercept = 0, slope=0, size = 1)  +
      ggtitle(paste("#Taxa:",num_of_taxa,"#Func:",num_of_functions,"Taxa contribution to",input_score,"when assessing with",input_permutation)) +
      ylab(paste(input_score, "test statistic")) +
      xlab("Function") +
      scale_x_discrete(labels=function_da$func_name, breaks=function_da$x_val, limits=function_da$x_val)
    
    
    if (flip_coord) {
      DABarPlot = DABarPlot + 
        scale_y_continuous(labels=ggplot_build(BarPlot)$panel$ranges[[1]]$x.labels,
                           breaks=ggplot_build(BarPlot)$panel$ranges[[1]]$x.major_source,
                           limits=ggplot_build(BarPlot)$panel$ranges[[1]]$x.range) + 
        coord_flip()
    } else {
      DABarPlot = DABarPlot +
        scale_y_continuous(labels=ggplot_build(BarPlot)$panel$ranges[[1]]$y.labels,
                           breaks=ggplot_build(BarPlot)$panel$ranges[[1]]$y.major_source,
                           limits=ggplot_build(BarPlot)$panel$ranges[[1]]$y.range)
    }
    

  }

  #################################################################################################
  # Spearman correlation heatmap of taxa contribution across functions
  #################################################################################################
  if (plot_type == "taxa_corr_heatmap") {
    spearman_corr_matrix = cor(df[,-1], method="spearman")
    strict_pairwise_corr_vector = spearman_corr_matrix[upper.tri(spearman_corr_matrix)]
    mean_pairwise_corr = mean(strict_pairwise_corr_vector)
    median_pairwise_corr = median(strict_pairwise_corr_vector)

    #print(function_da)
    pos_functions = row.names(function_da[1:num_cases,])
    neg_functions = row.names(function_da[(num_cases+1):dim(function_da)[1],])
    #print(pos_functions)
    #print(neg_functions)
    if (length(pos_functions) > 2) {
      spearman_corr_matrix_pos = cor(df[,pos_functions], method="spearman")
      strict_pairwise_corr_vector_pos = spearman_corr_matrix_pos[upper.tri(spearman_corr_matrix_pos)]
      mean_pairwise_corr_pos = mean(strict_pairwise_corr_vector_pos)
      median_pairwise_corr_pos = median(strict_pairwise_corr_vector_pos)
    } else {
      mean_pairwise_corr_pos = 0
      median_pairwise_corr_pos = 0
    }
    if (length(neg_functions) > 2) {
      spearman_corr_matrix_neg = cor(df[,neg_functions], method="spearman")
      strict_pairwise_corr_vector_neg = spearman_corr_matrix_neg[upper.tri(spearman_corr_matrix_neg)]
      mean_pairwise_corr_neg = mean(strict_pairwise_corr_vector_neg)
      median_pairwise_corr_neg = median(strict_pairwise_corr_vector_neg)
    } else {
      mean_pairwise_corr_neg = 0
      median_pairwise_corr_neg = 0
    }

    #print(paste("Mean/Median Spearman:",format(mean_pairwise_corr, digits=2, nsmall=2),"/",format(median_pairwise_corr, digits=2, nsmall=2)))
    #print(paste("Mean/Median Spearman POS:",format(mean_pairwise_corr_pos, digits=2, nsmall=2),"/",format(median_pairwise_corr_pos, digits=2, nsmall=2)))
    #print(paste("Mean/Median Spearman NEG:",format(mean_pairwise_corr_neg, digits=2, nsmall=2),"/",format(median_pairwise_corr_neg, digits=2, nsmall=2)))

    spearman_corr_matrix.m = reshape2::melt(spearman_corr_matrix)
    names(spearman_corr_matrix.m) = c("Function", "Function2", "Spearman")

    HeatmapPlot = ggplot(spearman_corr_matrix.m, aes(Function2, Function)) +
      geom_tile(aes(fill = Spearman), colour = "white") +
      scale_fill_gradient2("Spearman", low=muted("blue"), high=muted("red")) +
      theme(axis.title.x = element_blank()) +
      theme(axis.text.x = element_blank())

    if (add_case_control_line) {
      HeatmapPlot = HeatmapPlot +
        geom_vline(xintercept = num_cases + 0.5, size = 1, linetype = "dashed") +
        geom_abline(intercept = num_cases + 0.5, slope = 0, size = 2, linetype = "longdash")
    }

    HeatmapPlot = HeatmapPlot +
      ggtitle(paste("#Taxa:",num_of_taxa,"#Func:",num_of_functions,
                    "Mean Spearman All:", format(mean_pairwise_corr, digits=2, nsmall=2),
                    "Pos:", format(mean_pairwise_corr_pos, digits=2, nsmall=2),
                    "Neg:", format(mean_pairwise_corr_neg, digits=2, nsmall=2),
                    "score:",input_score,"assessing with:",input_permutation))

  }

  #################################################################################################
  # PERCENT OF FUNCTION ABUNDANCE EXPLAINED
  #################################################################################################
  if (plot_type == "percent_abundance_explained" || plot_type == "percent_abundance_explained_pearson") {
    df_for_percent_abun_explained = merge(predicted_abundance_agreement, function_da, by="row.names")
    if (plot_type == "percent_abundance_explained") {
      df_for_percent_abun_explained = df_for_percent_abun_explained[, c("Row.names","R.2","x_val")]
    } else { # percent_abundance_explained_pearson
      df_for_percent_abun_explained = df_for_percent_abun_explained[, c("Row.names","PearsonCorr","x_val")]
    }

    row.names(df_for_percent_abun_explained) = df_for_percent_abun_explained[,1]

    # sort functions
    if (sort_by == "unknown") {
      if (flip_coord) {
        sorted_functions_by_unknown = df_for_percent_abun_explained$id[sort(df_for_percent_abun_explained$value, decreasing=TRUE, index.return=TRUE)$ix]
      } else {
        sorted_functions_by_unknown = df_for_percent_abun_explained$id[sort(df_for_percent_abun_explained$value, decreasing=FALSE, index.return=TRUE)$ix]
      }
    } else if (sort_by == "da" || sort_by == "predicted_da" || sort_by == "list") {
      sorted_functions_by_unknown = functions
    }

    df_for_percent_abun_explained = df_for_percent_abun_explained[sorted_functions_by_unknown, c(1,2)]
    names(df_for_percent_abun_explained) = c("id", "value")
    df_for_percent_abun_explained.m = reshape2::melt(df_for_percent_abun_explained)
    df_for_percent_abun_explained.m$id = factor(df_for_percent_abun_explained.m$id, levels=df_for_percent_abun_explained.m$id)
    if (plot_type == "percent_abundance_explained") {
      df_for_percent_abun_explained.m$value = 100 * df_for_percent_abun_explained.m$value
    }

    UnknownPlot = ggplot() +
      geom_bar(data=df_for_percent_abun_explained.m, aes(x=id, y=value), stat="identity",
               width=0.8, fill="gray", colour="black")

    if (add_case_control_line) {
      UnknownPlot = UnknownPlot +
        geom_vline(xintercept = num_cases + 0.5, size = 1, linetype = "dashed")
    }

    if (plot_type == "percent_abundance_explained") {
      scale_y_continuous(breaks = c(0,50,100))
      coord_cartesian(ylim = c(0, 1))
    }

    if (flip_coord) {
      UnknownPlot = UnknownPlot + coord_flip()
    }

    UnknownPlot = UnknownPlot +
      geom_vline(xintercept=0) +
      geom_abline(slope=0, intercept=0) +
      xlab("Percent of functions") +
      ylab("Percent of explained functional shift") +
      theme(plot.title = element_text(face="bold", size=20)) +
      theme(axis.title = element_text(face="bold", size=20))

  }




  #################################################################################################
  # PERCENT UNKNOWN
  #################################################################################################
  if (plot_type == "percent_unknown" || plot_type == "percent_unknown_no_negative"
      || plot_type == "percent_unknown_distribution" || plot_type == "percent_unknown_cumulative"
      || plot_type == "percent_predicted" || plot_type == "percent_predicted_diamonds") {
    df_for_percent_unknown = merge(predicted_da, function_da, by="row.names")
    df_for_percent_unknown = df_for_percent_unknown[,1:3]

    if (plot_type == "percent_predicted_diamonds") {
      df_for_percent_unknown = df_for_percent_unknown[, 1:3]
      row.names(df_for_percent_unknown) = df_for_percent_unknown[, 1]
      df_for_percent_unknown[, 1] = NULL
      names(df_for_percent_unknown) = c("Predicted DA", "Function DA")
      df_for_percent_unknown = t(df_for_percent_unknown)
      df_for_neg_predicted = df_for_percent_unknown

    } else {
      df_for_percent_unknown$agree = sign(df_for_percent_unknown$x_vals * df_for_percent_unknown$y_vals)
      df_for_percent_unknown$Explained = abs(df_for_percent_unknown$x_vals) / abs(df_for_percent_unknown$y_vals)
      df_for_percent_unknown$Unexplained = 1 - df_for_percent_unknown$Explained
      df_for_percent_unknown$Unexplained[df_for_percent_unknown$agree < 0] = 1
      df_for_percent_unknown$Unexplained[df_for_percent_unknown$Unexplained < 0] = 0
      df_for_percent_unknown$Unexplained[df_for_percent_unknown$Unexplained > 1] = 1
      df_for_percent_unknown$Explained[df_for_percent_unknown$Explained > 1] = 1

      df_for_neg_predicted = df_for_percent_unknown[,c("Row.names", "agree", "Explained")]
      df_for_neg_predicted$Explained[df_for_neg_predicted$agree < 0] = -1 * df_for_neg_predicted$Explained[df_for_neg_predicted$agree < 0]
      df_for_neg_predicted$Explained[df_for_neg_predicted$agree > 0] = 0
      df_for_neg_predicted = df_for_neg_predicted[, c(1,3)]
      row.names(df_for_neg_predicted) = df_for_neg_predicted[, 1]
      df_for_neg_predicted[, 1] = NULL
      names(df_for_neg_predicted) = "Predicted DA with opposite sign"

      df_for_percent_unknown$Explained[df_for_percent_unknown$agree < 0] = 0
      df_for_percent_unknown = df_for_percent_unknown[,c(1,5,6)]
      row.names(df_for_percent_unknown) = df_for_percent_unknown[, 1]
      df_for_percent_unknown[, 1] = NULL
      names(df_for_percent_unknown) = c("Predicted DA", "Unexplained DA")
      df_for_percent_unknown = t(df_for_percent_unknown)
      df_for_neg_predicted = t(df_for_neg_predicted)
    }

    average_percent_unknown = mean(df_for_percent_unknown[2, ])

    # sort functions
    if (sort_by == "unknown") {
      if (flip_coord) {
        sorted_functions_by_unknown = names(sort(df_for_percent_unknown[2, ] + abs(df_for_neg_predicted[1,]), decreasing=TRUE))
      } else {
        sorted_functions_by_unknown = names(sort(df_for_percent_unknown[2, ] + abs(df_for_neg_predicted[1,]), decreasing=FALSE))
      }
    } else if (sort_by == "da" || sort_by == "predicted_da" || sort_by == "list") {
      sorted_functions_by_unknown = functions
    }

    df_for_percent_unknown = df_for_percent_unknown[, sorted_functions_by_unknown, drop=FALSE]
    df_for_neg_predicted = df_for_neg_predicted[, sorted_functions_by_unknown, drop=FALSE]

    df_for_percent_unknown.m = reshape2::melt(df_for_percent_unknown)
    df_for_neg_predicted.m = reshape2::melt(df_for_neg_predicted)
    df_for_neg_predicted.m$Var1 = "Predicted DA with opposite sign"
    df_for_neg_predicted.m$Var2 = 1:length(row.names(df_for_neg_predicted.m))

    if (plot_type == "percent_unknown_distribution") {
      df_for_percent_unknown_distribution.m = reshape2::melt(df_for_percent_unknown[1,])
      UnknownPlot =
        ggplot(df_for_percent_unknown_distribution.m, aes(x=value, y=..count../sum(..count..))) +
        geom_freqpoly(size=2, binwidth=0.1, origin=-0.15) + ylab("Fraction of functions") +
        xlab("Percentage of explained functional shift") + coord_cartesian(xlim=c(0,1))

    } else if (plot_type == "percent_unknown_cumulative") {
      df_for_percent_unknown_distribution.m = reshape2::melt(df_for_percent_unknown[2,])
      UnknownPlot =
        ggplot(df_for_percent_unknown_distribution.m, aes(x=value)) +
        stat_ecdf(size=2) +
        ylab("Fraction of functions") + xlab("Percentage of explained functional shift") +
        coord_cartesian(xlim=c(-0.001,1.001), ylim=c(0,1))

    } else if (plot_type == "percent_predicted") {
      df_for_percent_unknown.m$value = 100 * df_for_percent_unknown.m$value
      UnknownPlot = ggplot() +
        geom_bar(data=df_for_percent_unknown.m[df_for_percent_unknown.m$Var1 == "Predicted DA", ],
                 aes(x=Var2, y=value, fill=Var1), stat = "identity", width = 1) +
        scale_y_continuous(breaks = c(0,50,100)) +
        coord_cartesian(ylim = c(0, 1))

    } else if (plot_type == "percent_predicted_diamonds") {
      UnknownPlot = ggplot() +
        geom_point(data=df_for_percent_unknown.m[df_for_percent_unknown.m$Var1 == "Predicted DA", ],
                 aes(x=Var2, y=value, fill=Var1), stat = "identity", colour="black", size=4, shape=23) +
        geom_point(data=df_for_percent_unknown.m[df_for_percent_unknown.m$Var1 == "Function DA", ],
                 aes(x=Var2, y=value, fill=Var1), stat = "identity", colour="black", size=4, shape=23) +
        scale_fill_manual(values = c("Function DA" = "red", "Predicted DA" = "white"),
                          labels = c("gene-centric\nfunctional shift", "taxa-based\nfunctional shift"))
        if (add_case_control_line) {
          UnknownPlot = UnknownPlot + geom_vline(xintercept = num_cases + 0.5, size = 1, linetype = "dashed")
        }

    } else {
      UnknownPlot = ggplot() +
        geom_bar(data=df_for_percent_unknown.m, aes(x=Var2, y=value, fill=Var1), stat = "identity", width = 1) +
        scale_y_continuous(labels = percent_format())
    }

    if (plot_type == "percent_unknown") {
      UnknownPlot = UnknownPlot +
        geom_bar(data=df_for_neg_predicted.m, aes(x=Var2, y=value, fill=Var1), stat = "identity", width = 1)
    }

    if (flip_coord) {
      UnknownPlot = UnknownPlot + coord_flip()
    }

    if (plot_type != "percent_unknown_distribution" && plot_type != "percent_unknown_cumulative") {
      UnknownPlot = UnknownPlot +
        geom_vline(xintercept=0) +
        geom_abline(slope=0, intercept=0) +
        xlab("Percent of functions") +
        ylab("Percent of explained functional shift") +
        theme(legend.title = element_blank()) +
        theme(plot.title = element_text(face="bold", size=20)) +
        theme(axis.title = element_text(face="bold", size=20)) +
        ggtitle(paste("#Taxa:",num_of_taxa,"#Func:",num_of_functions,"Average % unexplained:",
                      paste(format(100 * average_percent_unknown, digits=2, nsmall=2),"%", sep=""),
                    "score:",input_score,"assessing with:",input_permutation))

      if (flip_coord) {
        UnknownPlot = UnknownPlot + theme(axis.text.x = element_text(face="bold", size=15))
        if (num_of_functions < 30) {
          UnknownPlot = UnknownPlot +  theme(axis.text.y = element_text(face="bold", size=12))
        } else {
          UnknownPlot = UnknownPlot +  theme(axis.text.y = element_blank())
        }
      } else {
        UnknownPlot = UnknownPlot + theme(axis.text.y = element_text(face="bold", size=15))
        if (num_of_functions < 30) {
          UnknownPlot = UnknownPlot +  theme(axis.text.x = element_text(face="bold", size=12))
        } else {
          UnknownPlot = UnknownPlot +  theme(axis.text.x = element_blank())
        }
      }
    }

  }


  #################################################################################################
  # SCATTER OF TAXA MEAN ABUNDANCE IN CASES AND CONTROLS VS. ITS AVERAGE CONTRIBUTION
  #################################################################################################
  if (plot_type == "abun_cases_controls_vs_contribution") {

    data_abun_vs_cont = df[df$Taxa != "Unknown", 1:2]
    data_abun_vs_cont$AvgContribution = rowMeans((df[df$Taxa != "Unknown", -1]))
    data_abun_vs_cont = data_abun_vs_cont[,c(1,3)]
    data_abun_vs_cont = merge(data_abun_vs_cont, taxa_da_original_names,  by="Taxa")[1:4]

    top_avg_contributors = sort(abs(data_abun_vs_cont$AvgContribution),
                                index.return = TRUE, decreasing = TRUE)$ix[1:10]

    AbunVsContPlot =
      ggplot(data_abun_vs_cont, aes(x=meanCases, y=meanControls)) +
      geom_point(aes(fill=AvgContribution), colour="black",pch=21, size=5) +
      scale_fill_gradient2(low="red", high="blue") +
      ggtitle(paste("#Taxa:",num_of_taxa,"#Func:",num_of_functions, "score:",input_score,"assessing with:",input_permutation)) +
      geom_abline(slope=0, intercept=0) +
      geom_vline(xintercept=0) +
      theme(plot.title = element_text(face="bold", size=20)) +
      theme(axis.title = element_text(face="bold", size=20)) +
      theme(axis.text.x = element_text(face="bold", size=15)) +
      theme(axis.text.y = element_text(face="bold", size=15)) +
      geom_abline(slope=1, intercept=0) +
      geom_text(data=data_abun_vs_cont[top_avg_contributors,],
                aes(x=meanCases, y=meanControls, label=Taxa),hjust=0, vjust=0)

  }

  #################################################################################################
  # SCATTER OF TAXA DIFFERENCE IN MEAN ABUNDANCE IN CASES AND CONTROLS VS. ITS AVERAGE CONTRIBUTION
  #################################################################################################
  if (plot_type == "abun_difference_vs_contribution") {

    data_abun_vs_cont = df[df$Taxa != "Unknown", 1:2]
    data_abun_vs_cont$AvgContribution = rowMeans((df[df$Taxa != "Unknown", -1]))
    data_abun_vs_cont = data_abun_vs_cont[,c(1,3)]
    data_abun_vs_cont = merge(data_abun_vs_cont, taxa_da_original_names,  by="Taxa")[1:4]
    data_abun_vs_cont$diff = data_abun_vs_cont$meanCases - data_abun_vs_cont$meanControls
    data_abun_vs_cont = data_abun_vs_cont[,c(1,2,5)]

    top_avg_contributors = sort(abs(data_abun_vs_cont$AvgContribution),
                                index.return = TRUE, decreasing = TRUE)$ix[1:10]

    AbunVsContPlot =
      ggplot(data_abun_vs_cont, aes(x=diff, y=AvgContribution)) +
      geom_point(colour="black", size=5) +
      ggtitle(paste("#Taxa:",num_of_taxa,"#Func:",num_of_functions, "score:",input_score,"assessing with:",input_permutation)) +
      geom_abline(slope=0, intercept=0) +
      geom_vline(xintercept=0) +
      theme(plot.title = element_text(face="bold", size=20)) +
      theme(axis.title = element_text(face="bold", size=20)) +
      theme(axis.text.x = element_text(face="bold", size=15)) +
      theme(axis.text.y = element_text(face="bold", size=15)) +
      geom_text(data=data_abun_vs_cont[top_avg_contributors,],
                aes(x=diff, y=AvgContribution, label=Taxa),hjust=0, vjust=0)

  }

  #################################################################################################
  # SCATTER OF TAXA MEAN ABUNDANCE IN CASES AND CONTROLS VS. ITS AVERAGE CONTRIBUTION
  #################################################################################################
  if (plot_type == "taxa_avg_abun_vs_contribution") {

    data_abun_vs_cont = df[df$Taxa != "Unknown", 1:2]
    data_abun_vs_cont$AvgContribution = rowMeans((df[df$Taxa != "Unknown", -1]))
    data_abun_vs_cont = data_abun_vs_cont[,c(1,3)]
    data_abun_vs_cont = merge(data_abun_vs_cont, taxa_da_original_names,  by="Taxa")[1:4]
    data_abun_vs_cont$AvgAbun = (data_abun_vs_cont$meanCases + data_abun_vs_cont$meanControls) / 2
    data_abun_vs_cont = data_abun_vs_cont[,c(1,2,5)]

    top_avg_contributors = sort(abs(data_abun_vs_cont$AvgContribution),
                                index.return = TRUE, decreasing = TRUE)$ix[1:10]

    AbunVsContPlot =
      ggplot(data_abun_vs_cont, aes(x=AvgAbun, y=AvgContribution)) +
      geom_point(colour="black", size=5) +
      ggtitle(paste("#Taxa:",num_of_taxa,"#Func:",num_of_functions, "score:",input_score,"assessing with:",input_permutation)) +
      geom_abline(slope=0, intercept=0) +
      geom_vline(xintercept=0) +
      theme(plot.title = element_text(face="bold", size=20)) +
      theme(axis.title = element_text(face="bold", size=20)) +
      theme(axis.text.x = element_text(face="bold", size=15)) +
      theme(axis.text.y = element_text(face="bold", size=15)) +
      geom_text(data=data_abun_vs_cont[top_avg_contributors,],
                aes(x=AvgAbun, y=AvgContribution, label=Taxa),hjust=0, vjust=0)

  }

  #################################################################################################
  # SCATTER OF TAXA MEAN ABUNDANCE IN CASES AND CONTROLS VS. ITS AVERAGE ABSOLUTE CONTRIBUTION
  #################################################################################################
  if (plot_type == "taxa_avg_abun_vs_abs_contribution") {

    data_abun_vs_cont = df[df$Taxa != "Unknown", 1:2]
    data_abun_vs_cont$AvgAbsContribution = abs(rowMeans((df[df$Taxa != "Unknown", -1])))
    data_abun_vs_cont = data_abun_vs_cont[,c(1,3)]
    data_abun_vs_cont = merge(data_abun_vs_cont, taxa_da_original_names,  by="Taxa")[1:4]
    data_abun_vs_cont$AvgAbun = (data_abun_vs_cont$meanCases + data_abun_vs_cont$meanControls) / 2
    data_abun_vs_cont = data_abun_vs_cont[,c(1,2,5)]

    top_avg_contributors = sort(abs(data_abun_vs_cont$AvgAbsContribution),
                                index.return = TRUE, decreasing = TRUE)$ix[1:10]

    cor_value = cor(data_abun_vs_cont$AvgAbun, data_abun_vs_cont$AvgAbsContribution, method="pearson")

    AbunVsContPlot =
      ggplot(data_abun_vs_cont, aes(x=AvgAbun, y=AvgAbsContribution)) +
      geom_point(colour="black", size=5) +
      ggtitle(paste("corr =",format(cor_value, digits=2, nsmall=2),"#Taxa:",num_of_taxa,"#Func:",num_of_functions, "score:",input_score,"assessing with:",input_permutation)) +
      geom_abline(slope=0, intercept=0) +
      geom_vline(xintercept=0) +
      geom_smooth(method=lm, se=FALSE, color="black") +
      theme(plot.title = element_text(face="bold", size=20)) +
      theme(axis.title = element_text(face="bold", size=20)) +
      theme(axis.text.x = element_text(face="bold", size=15)) +
      theme(axis.text.y = element_text(face="bold", size=15)) +
      geom_text(data=data_abun_vs_cont[top_avg_contributors,],
                aes(x=AvgAbun, y=AvgAbsContribution, label=Taxa),hjust=1, vjust=0)

  }

  #################################################################################################
  # SCATTER OF TAXA DIFFERNTIAL ABUNDANCE IN CASES AND CONTROLS VS. ITS AVERAGE CONTRIBUTION
  #################################################################################################
  if (plot_type == "taxa_DA_vs_contribution") {

    data_abun_vs_cont = df[df$Taxa != "Unknown", 1:2]
    data_abun_vs_cont$AvgContribution = rowMeans((df[df$Taxa != "Unknown", -1]))
    data_abun_vs_cont = data_abun_vs_cont[,c(1,3)]
    data_abun_vs_cont = merge(data_abun_vs_cont, taxa_da_original_names,  by="Taxa")[c(1,2,5)]

    top_avg_contributors = sort(abs(data_abun_vs_cont$AvgContribution),
                                index.return = TRUE, decreasing = TRUE)$ix[1:10]

    AbunVsContPlot =
      ggplot(data_abun_vs_cont, aes(x=StatValue, y=AvgContribution)) +
      geom_point(colour="black", size=5) +
      ggtitle(paste("#Taxa:",num_of_taxa,"#Func:",num_of_functions, "score:",input_score,"assessing with:",input_permutation)) +
      geom_abline(slope=0, intercept=0) +
      geom_vline(xintercept=0) +
      theme(plot.title = element_text(face="bold", size=20)) +
      theme(axis.title = element_text(face="bold", size=20)) +
      theme(axis.text.x = element_text(face="bold", size=15)) +
      theme(axis.text.y = element_text(face="bold", size=15)) +
      geom_text(data=data_abun_vs_cont[top_avg_contributors,],
                aes(x=StatValue, y=AvgContribution, label=Taxa),hjust=0, vjust=0)

  }

  #################################################################################################
  # SCATTER OF SUM OF INDIVIDUAL CONTRIBUTIONS VS. ORIGINAL DIFFERENTIAL ABUNDANCE SCORE
  #################################################################################################
  if (plot_type == "sum_vs_original" || plot_type == "sum_vs_predicted") {

    if (plot_type == "sum_vs_original") {
      data_scatter = merge(data.frame(colSums(df[,-1,FALSE])), function_da,  by="row.names")[1:3]
    } else {
      data_scatter = merge(data.frame(colSums(df[,-1,FALSE])), predicted_da,  by="row.names")[1:3]
    }

    names(data_scatter) = c("function","pos_neg_sum","original")

    ScatterPlot =  ggplot() + geom_point(data=data_scatter, aes(x=pos_neg_sum, y=original), size=5)

    min_xy = min(ggplot_build(ScatterPlot)$panel$ranges[[1]]$y.range[1], ggplot_build(ScatterPlot)$panel$ranges[[1]]$x.range[1])
    max_xy = max(ggplot_build(ScatterPlot)$panel$ranges[[1]]$y.range[2], ggplot_build(ScatterPlot)$panel$ranges[[1]]$x.range[2])

    ScatterPlot <- ScatterPlot + xlim(min_xy, max_xy) + ylim(min_xy, max_xy) +
      geom_abline(slope=0, intercept=0) +
      ggtitle(paste("#Taxa:",num_of_taxa,"#Func:",num_of_functions, "score:",input_score,"assessing with:",input_permutation)) +
      geom_vline(xintercept=0) +
      theme(plot.title = element_text(face="bold", size=20)) +
      theme(axis.title = element_text(face="bold", size=20)) +
      theme(axis.text.x = element_text(face="bold", size=15)) +
      theme(axis.text.y = element_text(face="bold", size=15)) +
      geom_abline(slope=1, intercept=0)
  }

  #################################################################################################
  # CURVE OF SUM OF KNOWN INDIVIDUAL CONTRIBUTIONS VS. ORIGINAL DIFFERENTIAL ABUNDANCE SCORE
  #################################################################################################
  if (plot_type == "curve_sum_known_vs_original") {
    df_known = df[which(df$Taxa != "Unknown"),]
    data_scatter = merge(data.frame(colSums(df_known[,-1,FALSE])), function_da,  by="row.names")[1:3]
    names(data_scatter) = c("function","pos_neg_sum","original")
    # sort the data by the original diff. abun. value
    sort_index = sort(data_scatter$original, decreasing = TRUE, index.return = TRUE)$ix
    data_scatter = data_scatter[sort_index, ]
    data_scatter$x = 1:length(data_scatter$original)
    num_pos_sign = sum(data_scatter$original > 0)

    ScatterPlot =  ggplot() +
      geom_point(data=data_scatter, aes(x=x, y=original), colour="black", fill="red", size=4, shape=23) +
      geom_point(data=data_scatter, aes(x=x, y=pos_neg_sum), colour="black", fill="black", size=4)

    ScatterPlot <- ScatterPlot +
      ylab("Functional shift") +
      xlab("Function") +
      coord_cartesian(xlim = c(-1, length(data_scatter$original)+1)) +
      geom_abline(slope=0, intercept=0) +
      ggtitle(paste("#Taxa:",num_of_taxa,"#Func:",num_of_functions, "score:",input_score,"assessing with:",input_permutation)) +
      geom_vline(xintercept=(num_pos_sign+0.5)) +
      theme(plot.title = element_text(face="bold", size=20)) +
      theme(axis.title = element_text(face="bold", size=20)) +
      theme(axis.text.x = element_blank()) +
      theme(axis.text.y = element_text(face="bold", size=15))
  }

  #################################################################################################
  # CURVE OF PREDICTED DA GIVEN THE TAXA AND GENOMIC CONTENT VS. ORIGINAL DIFFERENTIAL ABUNDANCE SCORE
  #################################################################################################
  if (plot_type == "curve_predicted_da_vs_original") {
    data_scatter = merge(predicted_da, function_da,  by="row.names")[1:3]
    names(data_scatter) = c("function","expected_da","original")
    # sort the data by the original diff. abun. value
    sort_index = sort(data_scatter$original, decreasing = TRUE, index.return = TRUE)$ix
    data_scatter = data_scatter[sort_index, ]
    data_scatter$x = 1:length(data_scatter$original)
    num_pos_sign = sum(data_scatter$original > 0)

    ScatterPlot =  ggplot() +
      geom_point(data=data_scatter, aes(x=x, y=original), colour="black", fill="red", size=4, shape=23) +
      geom_point(data=data_scatter, aes(x=x, y=expected_da), colour="black", fill="black", size=4)

    ScatterPlot <- ScatterPlot +
      ylab("Functional shift") +
      xlab("Function") +
      coord_cartesian(xlim = c(-1, length(data_scatter$original)+1)) +
      geom_abline(slope=0, intercept=0) +
      ggtitle(paste("#Taxa:",num_of_taxa,"#Func:",num_of_functions, "score:",
                    input_score,"assessing with:",input_permutation)) +
      geom_vline(xintercept=(num_pos_sign+0.5)) +
      theme(plot.title = element_text(face="bold", size=20)) +
      theme(axis.title = element_text(face="bold", size=20)) +
      theme(axis.text.x = element_blank()) +
      theme(axis.text.y = element_text(face="bold", size=15))
  }

  #################################################################################################
  # BOXPLOTS OF THE PREDICTED FUNCTIONAL ABUNDANCE WITH THE REAL AND RESIDUALS
  #################################################################################################
  if (plot_type == "boxplot_predicted_abun_vs_real_vs_residual") {
    # sort the data by the original real abundacne
    sort_index = sort(rowMeans(real_function_abun[,-1]), decreasing = TRUE, index.return=TRUE)$ix
    real_function_abun.m = reshape2::melt(real_function_abun[sort_index, ])
    predicted_function_abun.m = reshape2::melt(predicted_function_abun)
    residual_function_abun.m = reshape2::melt(residual_function_abun)
    data_box = merge(merge(real_function_abun.m, predicted_function_abun.m, by = c("KO","variable")),
                     residual_function_abun.m, by = c("KO","variable"))
    data_box = data_box[,c(1,3:5)]
    names(data_box) = c("Function", "Real", "Predicted", "Residual")
    data_box.m = reshape2::melt(data_box)

    BoxPlot =  ggplot() +
      geom_boxplot(data=data_box.m, aes(x=reorder(Function, value, FUN=median), y=value, fill=variable)) +
      ylab("Function abundance") +
      xlab("Function") +
      theme(plot.title = element_text(face="bold", size=20)) +
      theme(axis.title = element_text(face="bold", size=20)) +
      theme(axis.text.x = element_blank()) +
      theme(axis.text.y = element_text(face="bold", size=15))
  }


  #################################################################################################
  # FUNCTION CONTRIBUTIONS STATS VS. FUNCTION INCIDENCE STATS
  #################################################################################################

  if (plot_type == "function_contributions_vs_stats") {
    # compute some statistics for each function regarding the distribution of contributions it has
    contMatrix = data.matrix(df[,-1])
    func_compare_stats = data.frame(colSds(contMatrix))
    func_compare_stats$ContMAD = colMads(contMatrix)
    func_compare_stats$ContMaxFrac = apply(sweep(abs(contMatrix), 2, colSums(abs(contMatrix)), FUN='/'),
                                           2, max)

    names(func_compare_stats) = c("ContSD","ContMAD","ContMaxFrac")

    # compare with statistics about the distribution of this function in genomes in general
    #print(func_compare_stats)
    #print(function_stats)
    func_compare_stats = merge(func_compare_stats, function_stats, by="row.names")
    row.names(func_compare_stats) = func_compare_stats[,1]
    func_compare_stats[,1] = NULL
    corr_mat = cor(data.matrix(func_compare_stats), data.matrix(func_compare_stats))

    ScatterPlot =  ggplot(func_compare_stats, aes(x=Entropy, y=ContMaxFrac)) +
      geom_point(size=5) +
      theme(plot.title = element_text(face="bold", size=20)) +
      theme(axis.title = element_text(face="bold", size=20)) +
      theme(axis.text.x = element_text(face="bold", size=15)) +
      theme(axis.text.y = element_text(face="bold", size=15)) +
      geom_smooth(method=lm) +
      ggtitle(paste("#Taxa:",num_of_taxa,"#Func:",num_of_functions, "Correlation:", format(corr_mat["ContMaxFrac","Entropy"], digits=2, nsmall=2),
                    "score:",input_score,"assessing with:",input_permutation))
  }


  #################################################################################################
  # FUNCTION CONTRIBUTIONS STATS VS. FUNCTION INCIDENCE STATS
  #################################################################################################

  if (plot_type == "inferred_copy_number_heatmap") {

    inferred_copy_number.m = reshape2::melt(inferred_copy_number)
    names(inferred_copy_number.m) = c("Taxa","Function","Copies")

    # if present normalize by the median count for this function in this dataset
    if (!is.null(input_function_counts)) {
      function_counts$Function  = function_counts$Stats
      inferred_copy_number.m = merge(inferred_copy_number.m, function_counts, by="Function")
      inferred_copy_number.m$Copies = inferred_copy_number.m$Copies / inferred_copy_number.m$Median
    }

    HeatmapPlot = ggplot(inferred_copy_number.m, aes(x=Function, y=Taxa, fill=Copies)) + geom_tile() +
      scale_fill_continuous(limits=c(0, 2), breaks=seq(0,2,by=0.25))

  }


  #################################################################################################
  # RETURN PLOTS
  #################################################################################################

  if (plot_type == "bars" || plot_type == "percentage_bars") {
    return(BarPlot)
  } else if (plot_type == "function_DA_bars" || plot_type == "predicted_DA_bars" 
             || plot_type == "predicted_and_function_DA_bars") {
    return(DABarPlot)
  } else if (plot_type == "sum_vs_original" || plot_type == "sum_vs_predicted"
             || plot_type == "function_contributions_vs_stats"
             || plot_type == "curve_sum_known_vs_original"
             || plot_type == "curve_predicted_da_vs_original"
             || plot_type == "pos_neg_ratios") {
    return(ScatterPlot)
  } else if (plot_type == "taxa_da") {
    return(TaxaBarPlot)
  } else if (plot_type == "percent_unknown" || plot_type == "percent_unknown_no_negative"
             || plot_type == "percent_unknown_distribution" || plot_type == "percent_unknown_cumulative"
             || plot_type == "percent_predicted" || plot_type == "percent_predicted_diamonds"
             || plot_type == "percent_abundance_explained" || plot_type == "percent_abundance_explained_pearson") {
    return(UnknownPlot)
  } else if (plot_type == "taxa_corr_heatmap"
             || plot_type == "inferred_copy_number_heatmap") {
    return(HeatmapPlot)
  } else if (plot_type == "abun_cases_controls_vs_contribution"
             || plot_type == "abun_difference_vs_contribution"
             || plot_type == "taxa_DA_vs_contribution"
             || plot_type == "taxa_avg_abun_vs_contribution"
             || plot_type == "taxa_avg_abun_vs_abs_contribution") {
    return(AbunVsContPlot)
  } else if (plot_type == "boxplot_predicted_abun_vs_real_vs_residual") {
    return(BoxPlot)
  } else if (plot_type == "return_final_taxa_list") {
    return(final_list_of_taxa_to_plot)
  }
  else {
    print("Error: unknown plot type")
    return(NULL)
  }

}

