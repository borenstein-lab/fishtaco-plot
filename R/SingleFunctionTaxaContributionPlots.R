#' Plots one of the following plots for a single function:
#' - Steps plot with mean and std contributions for all taxa
#' - Ring plot where the contributions are shown as a ring
#' - Bar plot where the contributions are shown as a bar
#'
#' @param input_dir The directory where all the FiShTaCo input files can be found (default: NULL).
#' @param input_prefix The prefix for all input files (e.g., pathway_with_t2f, default: NULL).
#' @param input_suffix The suffix for all input files (default: ".tab").
#' @param input_score The calculated shift score for all input files (default: wilcoxon).
#' @param input_matrix_mean The mean contribution stat name (default: "mean_stat").
#' @param input_matrix_std The contribution std stat name (default: "std_stat").
#' @param input_permutation The permutation type used for all input files (default: "multi_taxa").
#' @param input_original The original shift stat name (default: "original_value").
#' @param input_function The function to plot (default: NULL).
#' @param min_contribution The minimum taxon contribution to include in plot  (default: 0).
#' @param input_function_meta Metadata file for the functions to include (default: NULL).
#' @param input_taxa_taxonomy Taxonomy file for the given taxa (default: NULL).
#' @param plot_type defines the type of plot to be plotted (default: "steps").
#' @param flip_case_control flip the case and control contribution (default: FALSE).
#' @param class_label_1 label for the the controls (default: NULL).
#' @param class_label_2 label for the the cases (default: NULL).
#' @param add_phyla_names_to_species_label indicating whether to add the phylum to each species (default: FALSE).
#'
#' @return a handle to the resulting plot.
#' @export

SingleFunctionTaxaContributionPlots <- function(input_dir=NULL, 
                                                input_prefix=NULL, 
                                                input_suffix=".tab", 
                                                input_score="wilcoxon",
                                                input_matrix_mean="mean_stat", 
                                                input_matrix_std="std_stat",
                                                input_permutation="multi_taxa", 
                                                input_original="original_value",
                                                input_function=NULL,                                             
                                                min_contribution=0, 
                                                input_function_meta=NULL,
                                                input_taxa_taxonomy=NULL, 
                                                plot_type="steps",
                                                flip_case_control=FALSE, 
                                                class_label_1=NULL,
                                                class_label_2=NULL,
                                                add_phyla_names_to_species_label=FALSE) {

  
  use_example_input = FALSE
  
  if (use_example_input) { 
    print("Using example input for testing...")
    # example input
    input_dir="/Volumes/ohadm/OhadM/METAFIT/HMP_DATA/TONGUE_DORSUM_vs_BUCCAL_MUCOSA/Output_METAPHLAN_REMOVE_RESIDUAL_SCALE_PERMUTED"
    input_prefix="pathway_with_t2f"
    input_score="wilcoxon"
    input_matrix_mean="mean_stat"
    input_matrix_std="std_stat"
    input_original="original_value"
    input_permutation="multi_taxa"
    input_function="ko00020"
    input_suffix=".tab"
    plot_type="steps" # ring steps bar
    input_function_meta="/Volumes/ohadm/OhadM/MUSiCC/Matrices/PATHWAYvsNAME_BACTERIAL_KEGG_2013_07_15.lst"
    input_taxa_taxonomy = "/Volumes/ohadm/OhadM/METAFIT/HMP_DATA/METAPHLAN_taxa_vs_TAXONOMY.tab"
    min_contribution=0.025
    flip_case_control=TRUE
    class_label_1="TONGUE"
    class_label_2="BUCCAL"
    add_phyla_names_to_species_label=FALSE
  }


  #
  ##############
  # READ INPUT
  ##############
  mean_stat = read.table(paste(input_dir,'/',input_prefix,"_STAT_",input_matrix_mean,"_SCORE_",input_score,"_ASSESSMENT_",input_permutation,".tab", sep=""), sep = "\t", header=TRUE, stringsAsFactors=FALSE)
  std_stat = read.table(paste(input_dir,'/',input_prefix,"_STAT_",input_matrix_std,"_SCORE_",input_score,"_ASSESSMENT_",input_permutation,".tab", sep=""), sep = "\t", header=TRUE, stringsAsFactors=FALSE)
  original_stat = read.table(paste(input_dir,'/',input_prefix,"_STAT_",input_original,"_SCORE_",input_score,"_ASSESSMENT_",input_permutation,".tab", sep=""), sep = "\t", header=TRUE, stringsAsFactors=FALSE)
  original_stat[, input_score] = 0

  if (flip_case_control) {
    mean_stat[, -1] = -1 * mean_stat[, -1]
    std_stat[, -1] = -1 * std_stat[, -1]
    tmp_label = class_label_1
    class_label_1 = class_label_2
    class_label_2 = tmp_label
  }

  if (!is.null(input_function)){
    function_index = which(original_stat$KO == input_function)
  } else {
    function_index = 1
    input_function = original_stat$KO[function_index]
  }

  if (!is.null(input_function_meta)) {
    meta = read.table(input_function_meta, sep = "\t", header=FALSE, stringsAsFactors=FALSE,  quote="")
    rownames(meta) = meta[,1]
    meta[,1] = NULL
    function_name = meta[input_function,]
  } else {
    function_name = input_function
  }

  if (!is.null(input_taxa_taxonomy)) {
    print("Uploading taxonomy file")
    
    taxonomy = read.table(input_taxa_taxonomy, sep = "\t", header=FALSE, stringsAsFactors=FALSE)
    rownames(taxonomy) = taxonomy[,1]
    taxonomy[,1] = NULL
    names(taxonomy) = c('kingdom','phylum','class','order','family','genus','species')

    for (i in 1:length(mean_stat$Taxa)) {
      curr_taxa = mean_stat$Taxa[i]
      #print(curr_taxa)
      if (!is.na(taxonomy[curr_taxa,"species"])){
        #print(paste(taxonomy[curr_taxa,"genus"],taxonomy[curr_taxa,"species"]))
        if (taxonomy[curr_taxa,"species"] != " s__") {
          curr_taxonomy_names = paste(taxonomy[curr_taxa,"phylum"], taxonomy[curr_taxa,"genus"], taxonomy[curr_taxa,"species"])
        } else if (taxonomy[curr_taxa,"genus"] != " g__") {
          curr_taxonomy_names = paste(taxonomy[curr_taxa,"phylum"], taxonomy[curr_taxa,"genus"],paste("s__",curr_taxa,sep=""))
        } else if (taxonomy[curr_taxa,"family"] != " f__") {
          curr_taxonomy_names = paste(taxonomy[curr_taxa,"phylum"], taxonomy[curr_taxa,"family"],paste("s__",curr_taxa,sep=""))
        } else if (taxonomy[curr_taxa,"order"] != " o__") {
          curr_taxonomy_names = paste(taxonomy[curr_taxa,"phylum"], taxonomy[curr_taxa,"order"],paste("s__",curr_taxa,sep=""))
        } else if (taxonomy[curr_taxa,"order"] != " c__") {
          curr_taxonomy_names = paste(taxonomy[curr_taxa,"phylum"], taxonomy[curr_taxa,"class"],paste("s__",curr_taxa,sep=""))
        } else {
          curr_taxonomy_names = paste(taxonomy[curr_taxa,"kingdom"], taxonomy[curr_taxa,"phylum"], curr_taxa)
        }
      } else if (curr_taxa == "Unknown"){
        curr_taxonomy_names = curr_taxa
      } else {
        curr_taxonomy_names = paste(" g__", paste("s__",curr_taxa,sep=""))
      }

      mean_stat$Taxa[i] = curr_taxonomy_names
    }
  }

  sort_index = sort(mean_stat[,function_index+1], index.return=TRUE, decreasing=TRUE)$ix
  df = data.frame(mean_stat[sort_index,function_index+1], std_stat[sort_index,function_index+1])
  names(df) =  c("mean","std")
  df$Taxa = factor(mean_stat$Taxa[sort_index], levels=mean_stat$Taxa[sort_index])
  df = df[c("Taxa","mean","std")]

  ## remove taxa that have a very low contribution on selected function
  taxa_with_contribution = which(abs(df$mean) > min_contribution)
  df = df[taxa_with_contribution,]
  
  # if needed, clean the taxa names a bit
  if (!add_phyla_names_to_species_label) {
    df$Taxa = as.character(df$Taxa)
    df$Taxa[grepl("g__", df$Taxa)] = gsub("s__","", gsub("g__","", gsub("p__.*g__","g__", df$Taxa[grepl("g__", df$Taxa)])))
    
    df$Taxa = factor(df$Taxa, levels=unique(df$Taxa))
  }
    
  ############################################
  ## STEPS PLOT
  ############################################

  stepsPlot = ggplot(df, aes(x=mean, y=Taxa))  +
    geom_vline(xintercept = original_stat[function_index,2], colour="blue") +
    geom_errorbarh(aes(xmax = mean+std, xmin = mean-std, height = .4), size=1.5) +
    geom_point(colour="#990000", size=2) +
    ggtitle(paste(original_stat$KO[function_index], class_label_1, "vs", class_label_2)) +
    xlab(paste(function_name,": Mean contribution to",input_score,"when assessing with",input_permutation)) 
  
  
  ############################################
  ## RING PLOT
  ############################################

  if (plot_type == "ring" || plot_type == "bar") {
    # extract only the positive contributions
    pos = df[df['mean'] > 0,]
    neg = df[df['mean'] < 0,]
    pos_sum = abs(sum(pos$mean))
    neg_sum = abs(sum(neg$mean))

    if (neg_sum > pos_sum) { # this is a control-enriched function, so swipe the pos/neg data.frames
      tmp = neg
      neg = pos
      pos = tmp
    }

    pos_sum = abs(sum(pos$mean))
    pos$mean = round(100 * abs(pos$mean) / pos_sum)
    pos = pos[with(pos, order(mean, decreasing=TRUE)), ]
    pos$ymax = cumsum(pos$mean)
    pos$ymin = c(0, head(pos$ymax, n=-1))

    neg_sum = abs(sum(neg$mean))
    neg$mean = round(100 * abs(neg$mean) / neg_sum)
    neg = neg[with(neg, order(mean, decreasing=TRUE)), ]
    neg$ymax = cumsum(neg$mean)
    neg$ymin = c(0, head(neg$ymax, n=-1))


    neg_pos_ratio = neg_sum / pos_sum
    if (neg_pos_ratio < 0.1){
      neg_pos_ratio = 0.1
    }

    pos$pos_xmin = 3
    pos$pos_xmax = 4
    neg$neg_xmin = 2
    neg$neg_xmax = neg$neg_xmin + neg_pos_ratio

    # create the text of the percentages
    df_to_plot = data.frame(c(pos$mean, neg$mean), c(as.character(pos$Taxa), as.character(neg$Taxa)),
                          c(pos$ymax, neg$ymax), c(pos$ymin, neg$ymin),
                          c(pos$pos_xmax, neg$neg_xmax), c(pos$pos_xmin, neg$neg_xmin),
                           c((pos$pos_xmin+pos$pos_xmax)/2, (neg$neg_xmin+neg$neg_xmax)/2),
                           c(pos$mean/2 + c(0, cumsum(pos$mean)[-length(pos$mean)]),
                             neg$mean/2 + c(0, cumsum(neg$mean)[-length(neg$mean)])))

    names(df_to_plot) = c("val","Taxa","ymax","ymin","xmax","xmin","x_text","y_text")
    df_to_plot$Taxa = factor(df_to_plot$Taxa, levels=df_to_plot$Taxa)

    RingPlot = ggplot(df_to_plot, aes(fill=Taxa, ymax=ymax, ymin=ymin, xmax=xmax, xmin=xmin)) +
      geom_rect(colour="grey30") +
      coord_polar(theta="y") +
      geom_text(aes(x = x_text, y = y_text, label = paste(val,"%",sep="")), size=7) +
      xlim(c(0, 4)) +
      theme_bw() +
      theme(panel.border = element_blank()) +
      theme(panel.grid=element_blank()) +
      theme(axis.text=element_blank()) +
      theme(axis.ticks=element_blank()) +
      theme(axis.title=element_blank()) +
      ggtitle(paste(function_name,": Taxa contribution to",input_score,"when assessing with",input_permutation)) +
      theme(plot.title = element_text(face="bold", size=20))

    BarPlot = ggplot(df_to_plot, aes(fill=Taxa, xmax=ymax, xmin=ymin, ymax=xmax, ymin=xmin)) +
      geom_rect(colour="grey30") +
      geom_text(aes(x = y_text, y = x_text, label = paste(val,"%",sep="")), size=7) +
      ylim(c(0, 4)) +
      theme_bw() +
      theme(panel.border = element_blank()) +
      theme(panel.grid=element_blank()) +
      theme(axis.text=element_blank()) +
      theme(axis.ticks=element_blank()) +
      theme(axis.title=element_blank()) +
      ggtitle(paste(function_name,": Taxa contribution to",input_score,"when assessing with",input_permutation)) +
      theme(plot.title = element_text(face="bold", size=20))
  }

  #################################################################################################
  # RETURN
  #################################################################################################

  if (plot_type == "steps") {
    return(stepsPlot)
  } else if (plot_type == "ring") {
    return(RingPlot)
  } else if (plot_type == "bar") {
    return(BarPlot)
  } else {
    print("Error: unknown plot type")
    return(NULL)
  }

}

