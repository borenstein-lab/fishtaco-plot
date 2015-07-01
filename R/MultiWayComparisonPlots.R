MultiWayComparisonPlots <- function(input_dirs=NULL, input_names=NULL, input_prefix=NULL, input_score=NULL,
                                    input_matrix=NULL, input_original=NULL, input_permutation=NULL,
                                    input_suffix=NULL, plot_type="bars", input_function_meta=NULL,
                                    min_contribution=0, input_taxa_taxonomy=NULL, min_cont_as_separate=0.2,
                                    input_function_filter_list=NULL, show_only_pos=FALSE) {


  # example input
#     input_dirs=c("/Volumes/ohadm/OhadM/METAFIT/HMP_DATA/TONGUE_DORSUM_vs_BUCCAL_MUCOSA/Output_METAPHLAN",
#                  "/Volumes/ohadm/OhadM/METAFIT/HMP_DATA/TONGUE_DORSUM_vs_SUPRAGINGIVAL_PLAQUE/Output_METAPHLAN",
#                  "/Volumes/ohadm/OhadM/METAFIT/HMP_DATA/TONGUE_DORSUM_vs_ANTERIOR_NARES/Output_METAPHLAN",
#                  "/Volumes/ohadm/OhadM/METAFIT/HMP_DATA/TONGUE_DORSUM_vs_RETROAURICULAR_CREASE/Output_METAPHLAN",
#                  "/Volumes/ohadm/OhadM/METAFIT/HMP_DATA/TONGUE_DORSUM_vs_STOOL/Output_METAPHLAN")
#
#     input_names=c("TONGUEvsMUCOSA", "TONGUEvsPLAQUE", "TONGUEvsNARES", "TONGUEvsCREASE", "TONGUEvsSTOOL")
#
#     input_prefix="pathway_with_t2f"
#     input_score="wilcoxon"
#     input_matrix="taxa_contributions"
#     input_original="original_value"
#     input_permutation="permute_all_but_i" #   "permute_all_but_i" permuted_shapley_orderings
#     input_suffix=".tab"
#     plot_type="polar_bars"
#     input_function_meta="/Volumes/ohadm/OhadM/MUSiCC/Matrices/PATHWAYvsNAME_BACTERIAL_KEGG_2013_07_15.lst"
#     min_contribution=0
#     min_cont_as_separate=0.1
#     input_taxa_taxonomy = "/Volumes/ohadm/OhadM/METAFIT/HMP_DATA/METAPHLAN_taxa_vs_TAXONOMY.tab"
#     input_function_filter_list=c("ko02020", "ko00053", "ko00903", "ko00642", "ko00363")
#     show_only_pos=TRUE


  ##########################################
  # read input
  ##########################################

  number_of_comparisons = length(input_dirs)
  functions = input_function_filter_list

  #######
  # read the metanames of the function, for future use
  #######
  if (!is.null(input_function_meta)) {
    meta = read.table(input_function_meta, sep = "\t", header=FALSE, stringsAsFactors=FALSE,  quote="")
    meta = meta[meta$V1 %in% functions,]
    rownames(meta) = meta[,1]
    meta[,1] = NULL
  }

  #######
  # first read files and create the dataframe with just the list of all taxa in all comparisons
  #######
  all_taxa_list = vector(mode="character")
  for (i in 1:number_of_comparisons){
    curr_df = read.table(paste(input_dirs[i],"/",input_prefix,"_STAT_",input_matrix,"_SCORE_",
                          input_score,"_ASSESSMENT_",input_permutation,input_suffix, sep=""),
                    sep = "\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
    all_taxa_list = c(all_taxa_list, curr_df$Taxa)
  }
  all_taxa_list = unique(all_taxa_list)

  #######
  # now create the dataframe to hold all the comparisons
  #######
  df_comp = data.frame(all_taxa_list)
  row.names(df_comp) = df_comp[,1]
  df_comp[,1] = NULL

  for (i in 1:number_of_comparisons){
    #print(i)
    df = read.table(paste(input_dirs[i],"/",input_prefix,"_STAT_",input_matrix,"_SCORE_",
                          input_score,"_ASSESSMENT_",input_permutation,input_suffix, sep=""),
                    sep = "\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)

    # select the functions of interest
    #curr_functions_in_df = names(df)[names(df) %in% input_function_filter_list]
    #missing_functions_in_df = input_function_filter_list[!(input_function_filter_list %in% curr_functions_in_df)]
    #df = df[, c("Taxa",curr_functions_in_df)]
    df = df[, c("Taxa",input_function_filter_list)]
    row.names(df) = df[,1]
    df[,1] = NULL
    # add missing functions with zeros
    #if (length(missing_functions_in_df) > 0) {
    #  for (j in 1:length(missing_functions_in_df)){
    #    df[missing_functions_in_df[j]] = 0
    #  }
    #}

    # merge with other studies
    old_names = names(df_comp)
    #print(old_names)
    df_comp = merge(df_comp, df, by="row.names", all=TRUE)
    row.names(df_comp) = df_comp[,1]
    df_comp[,1] = NULL
    #names(df_comp) = c(old_names, paste(c(curr_functions_in_df,missing_functions_in_df),i,sep="_"))
    names(df_comp) = c(old_names, paste(input_function_filter_list,i,sep="_"))
  }

  df_comp[is.na(df_comp)] = 0

  #######
  # read in the taxonomy file
  #######
  if (!is.null(input_taxa_taxonomy)) {
    taxonomy = read.table(input_taxa_taxonomy, sep = "\t", header=FALSE, stringsAsFactors=FALSE)
    rownames(taxonomy) = taxonomy[,1]
    taxonomy[,1] = NULL
    names(taxonomy) = c('kingdom','phylum','class','order','family','genus','species')
  }

  #######
  # read in the function DA for the function at hand
  #######
  if (!is.null(input_original)) {
    function_da = data.frame()
    function_da_names = vector(mode="character")
    for (i in 1:number_of_comparisons){
      #print(i)
      curr_da = read.table(paste(input_dirs[i],"/",input_prefix,"_STAT_",input_original,"_SCORE_",
                                     input_score,"_ASSESSMENT_",input_permutation,input_suffix, sep=""),
                               sep = "\t", header=TRUE, stringsAsFactors=FALSE)

      function_da = rbind(function_da, curr_da[curr_da$KO %in% functions,])
      function_da_names = c(function_da_names, paste(curr_da$KO[curr_da$KO %in% functions],i,sep="_"))
    }
    #print(function_da_names)
    row.names(function_da) = function_da_names
    function_da[,1] = NULL
  }

  #######
  # remove taxa that have a very low contribution on all functions
  #######
  max_val = apply(abs(df_comp[, -1]), 1, max)
  taxa_with_contribution = which(max_val > min_contribution)
  df_comp = df_comp[taxa_with_contribution,]

  #######
  # if we are only using a single function, use the input names as labels
  #######
  if (!is.null(input_names) && (length(functions) == 1)) {
    names(df_comp) = input_names
    row.names(function_da) = input_names
  }

  #######
  # bring back the "Taxa" column
  #######
  functions = names(df_comp)
  df_comp$Taxa = row.names(df_comp)
  df = df_comp[, c("Taxa", functions)]

  #######
  # set some variables
  #######
  num_of_functions = length(functions)
  num_of_taxa = length(df$Taxa)

  ##########################################
  # draw plot
  ##########################################

  #######
  # separate each function to pos and neg
  #######
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

  if (!is.null(input_taxa_taxonomy)) {
    # find the strong taxa that contribute the most and leave them intact, but try to find their taxonomy name
    cont_fraction_pos = sweep(df_pos[,-1], 2, colSums(df_pos[,-1]), FUN='/')
    cont_fraction_pos[apply(cont_fraction_pos, 2, is.nan)] = 0

    cont_fraction_neg = sweep(df_neg[,-1], 2, colSums(df_neg[,-1]), FUN='/')
    cont_fraction_neg[apply(cont_fraction_neg, 2, is.nan)] = 0

    strong_contributors = sort(unique(c(which(apply(cont_fraction_pos, 1, max) > min_cont_as_separate),
                                        which(apply(cont_fraction_neg, 1, max) > min_cont_as_separate))))

    # always add "unknown" to strong contributors if we have it
    if (sum(row.names(cont_fraction_pos) == "Unknown") > 0) {
      strong_contributors = sort(unique(c(strong_contributors, which(row.names(cont_fraction_pos) == "Unknown"))))
    }


    df_pos_name_of_strongest = NULL
    # get the sorted list of sum of contributions, if we want to sort by the strongest contributor
    sorted_list_of_contributors = sort(rowSums(df_pos[,-1]), decreasing = TRUE,index.return = TRUE)$ix
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
        if (!is.na(taxonomy[curr_taxa,"species"])){
          #print(paste(taxonomy[curr_taxa,"genus"],taxonomy[curr_taxa,"species"]))
          if (taxonomy[curr_taxa,"species"] != " s__") {
            curr_taxonomy_names = paste(taxonomy[curr_taxa,"genus"], taxonomy[curr_taxa,"species"])
          } else if (taxonomy[curr_taxa,"genus"] != " g__") {
            curr_taxonomy_names = paste(taxonomy[curr_taxa,"genus"],paste("s__",curr_taxa,sep=""))
          } else if (taxonomy[curr_taxa,"family"] != " f__") {
            curr_taxonomy_names = paste(taxonomy[curr_taxa,"family"],paste("s__",curr_taxa,sep=""))
          } else if (taxonomy[curr_taxa,"order"] != " o__") {
            curr_taxonomy_names = paste(taxonomy[curr_taxa,"order"],paste("s__",curr_taxa,sep=""))
          } else if (taxonomy[curr_taxa,"order"] != " c__") {
            curr_taxonomy_names = paste(taxonomy[curr_taxa,"class"],paste("s__",curr_taxa,sep=""))
          } else {
            curr_taxonomy_names = paste(taxonomy[curr_taxa,"kingdom"],taxonomy[curr_taxa,"phylum"],curr_taxa)
          }
        } else if (curr_taxa == "Unknown"){
          curr_taxonomy_names = curr_taxa
        } else {
          curr_taxonomy_names = paste(" g__", paste("s__",curr_taxa,sep=""))
        }

        df_pos$Taxa[strong_contributors[i]] = curr_taxonomy_names
        df_neg$Taxa[strong_contributors[i]] = curr_taxonomy_names

        if (strong_contributors[i] == sorted_list_of_contributors[1]){

          df_pos_name_of_strongest = curr_taxonomy_names
          #print(df_pos_name_of_strongest)
        }

      }
    }

    #print(df_pos$Taxa)
    #print(taxa_da$Taxa)

    # for all the other taxa (that contribute less), merge them by pyhlum level
    for (i in 1:length(df_pos$Taxa)) {
      if (i %in% strong_contributors){ # skip stron contributors as they already have taxonomy assigned
        next
      }
      curr_taxa = df_pos$Taxa[i]
      #print(curr_taxa)
      if ( !is.na(taxonomy[curr_taxa,"phylum"])){

        ### option 1: add POS/NEG to the Phyla level
        #if (max(df_pos[i,-1]) * min(df_neg[i,-1]) < 0){
        #  df_pos$Taxa[i] = paste(taxonomy[curr_taxa,"phylum"])
        #  df_neg$Taxa[i] = paste(taxonomy[curr_taxa,"phylum"])
        #} else if (max(df_pos[i,-1]) > 0) {
        #  df_pos$Taxa[i] = paste(taxonomy[curr_taxa,"phylum"], "POS")
        #  df_neg$Taxa[i] = paste(taxonomy[curr_taxa,"phylum"], "POS")
        #} else {
        #  df_pos$Taxa[i] = paste(taxonomy[curr_taxa,"phylum"], "NEG")
        #  df_neg$Taxa[i] = paste(taxonomy[curr_taxa,"phylum"], "NEG")
        #}
        ### option 2: DON'T add POS/NEG to the Phyla level
        df_pos$Taxa[i] = paste(taxonomy[curr_taxa,"phylum"])
        df_neg$Taxa[i] = paste(taxonomy[curr_taxa,"phylum"])

      }
      else {
        df_pos$Taxa[i] = "Other taxa"
        df_neg$Taxa[i] = "Other taxa"
      }
    }

    #print(df_pos$Taxa)
    #print(taxa_da$Taxa)

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
    taxa_counts$Full[grepl("(1)", taxa_counts$Full, fixed=TRUE)] = taxa_counts$Taxa[grepl("(1)", taxa_counts$Full, fixed=TRUE)]

    # find the full name of the strongest contributor
    #print(taxa_counts$Taxa)
    #print(df_pos$Taxa)
    #print(strong_contributors)
    #print(df_pos_name_of_strongest)
    #print(which(taxa_counts$Taxa %in% df_pos_name_of_strongest))
    if (!is.null(df_pos_name_of_strongest)) {
      final_name_of_strong = taxa_counts$Full[which(taxa_counts$Taxa %in% df_pos_name_of_strongest)]
      #print(final_name_of_strong)
    }

    original_names = names(df_pos)
    df_pos = aggregate(df_pos[,-1], by=list(df_pos$Taxa), sum)
    names(df_pos) = original_names
    df_pos = merge(taxa_counts, df_pos, by="Taxa")
    df_pos$Taxa = NULL
    names(df_pos) = original_names

    df_neg = aggregate(df_neg[,-1], by=list(df_neg$Taxa), sum)
    names(df_neg) = original_names
    df_neg = merge(taxa_counts, df_neg, by="Taxa")
    df_neg$Taxa = NULL
    names(df_neg) = original_names

  }

  # if we want to show only the positive values, remove all rows with only zero values in the pos
  if (show_only_pos) {
    df_pos = df_pos[rowSums(df_pos[,-1]) > 0, ]
  }

  # set the color pallete
  num_of_final_taxa = length(df_pos$Taxa)
  is_other = sum(taxa_counts$Taxa == "Other taxa")
  is_unknown = sum(taxa_counts$Taxa == "Unknown")
  if (is_other && is_unknown){
    cPalette = c(rainbow(num_of_final_taxa-2), "#828282", "#000000")
  } else if (is_other) {
    cPalette = c(rainbow(num_of_final_taxa-1), "#828282")
  } else if (is_unknown) {
    cPalette = c(rainbow(num_of_final_taxa-1), "#000000")
  } else {
    cPalette = rainbow(num_of_final_taxa)
  }

  if (plot_type == "percentage_bars") {
    df_pos[,-1] = sweep(df_pos[,-1], 2, colSums(df_pos[,-1]), FUN='/')
    df_neg[,-1] = -1 * sweep(df_neg[,-1], 2, colSums(df_neg[,-1]), FUN='/')
  }

  df_pos.m = melt(df_pos)
  df_neg.m = melt(df_neg)

  BarPlot = ggplot() +
    geom_bar(data=df_pos.m, aes(x=variable, fill=Taxa, y=value), stat = "identity", colour="black") +
    coord_flip()

  if (!show_only_pos) {
    BarPlot = BarPlot + geom_bar(data=df_neg.m, aes(x=variable, fill=Taxa, y=value), stat = "identity", colour="black")
  }

  BarPlot = BarPlot +
    geom_abline(intercept = 0, slope=0, size = 2)  +
    scale_fill_manual(values=cPalette) +
    ggtitle(paste(input_function_filter_list,":","    #Taxa:",num_of_taxa,"#Func:",num_of_functions,"Taxa contribution to",input_score,"when assessing with",input_permutation)) +
    ylab(paste(input_score, "test statistic")) +
    xlab("Function") +
    theme_bw() +
    theme(plot.title = element_text(face="bold", size=20)) +
    theme(axis.title = element_text(face="bold", size=20)) +
    theme(axis.text.x = element_text(face="bold", size=15)) +
    theme(panel.grid=element_blank())

  if (num_of_functions < 50) {
    BarPlot = BarPlot +  theme(axis.text.y = element_text(face="bold", size=12))
  } else {
    BarPlot = BarPlot +  theme(axis.text.y = element_text(face="bold"))
  }

  if (plot_type == "percentage_bars") {
    BarPlot = BarPlot +
      ylab(paste(input_score, "test statistic fraction")) +
      scale_y_continuous(breaks=seq(-1,1,0.1))
  }

  if (!is.null(input_original) && (plot_type != "percentage_bars")) { # we have DA information for functions, so add line and markers
    # add markers
    #     if (!is.null(input_function_meta)) {
    #       function_da = merge(meta, function_da, by="row.names")[, 2:3]
    #       rownames(function_da) = function_da$V2
    #       function_da$V2 = NULL
    #     }
    #print(functions)
    #print(function_da)
    function_da$y_vals = function_da[functions, input_score]
    #print(function_da)
    function_da$x_val = 1:length(function_da$y_vals)
    function_da[input_score] = NULL
    rownames(function_da) = functions

    BarPlot = BarPlot + geom_point(data=function_da, aes(x_val, y_vals), colour="black", fill="red", size=4, shape=23)
  }


  #################################################################################################
  # Polar bar plot
  #################################################################################################

  if (plot_type == "polar_bars"){
    polar_df = df_pos.m
    polar_df$family = "None"
    #print(polar_df)
    # for each comparison, set the $family value to be the appropriate comparison
    for (i in 1:number_of_comparisons){
      #print(paste0("_",i))
      if (!is.null(input_names)) {
        polar_df$family[grepl(paste0("_",i), polar_df$variable)] = input_names[i]
      } else {
        polar_df$family[grepl(paste0("_",i), polar_df$variable)] = i
      }
    }

    if (!is.null(input_function_meta)) {
      polar_df$itemName = "None"
      for (i in 1:dim(meta)[1]){
        polar_df$itemName[grepl(row.names(meta)[i], polar_df$variable)] = meta$V2[i]
      }
    }
    names(polar_df) = c("score","item", "value","family", "itemName")
    #print(polar_df)

    PolarBarPlot = polarHistogram(polar_df, familyLabels=TRUE, pallete=cPalette)

  }



  if (plot_type == "bars" || plot_type == "percentage_bars") {
    return(BarPlot)
  } else if (plot_type == "polar_bars") {
    return(PolarBarPlot)
  } else {
    print("Error: unknown plot type")
    return(NULL)
  }

}




