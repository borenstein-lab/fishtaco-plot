PlotShapleyRandomOrderings <- function() {

  # plot shaply ordering info
    input_dir="/Volumes/ohadm/OhadM/METAFIT/HMP_DATA/TONGUE_DORSUM_vs_BUCCAL_MUCOSA/Output"
    input_prefix="pathway_with_t2f"
    input_score="wilcoxon"
    input_matrix="shapley_orderings"
    input_permutation="permuted_shapley_orderings"
    input_suffix=".tab"

    # read input
    df = read.table(paste(input_dir,"/",input_prefix,"_STAT_",input_matrix,"_SCORE_",
                        input_score,"_ASSESSMENT_",input_permutation,input_suffix, sep=""),
                  sep = "\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)

    rownames(df) = df[,1]
    df[,1] = NULL

    num_of_taxa = length(names(df))

    counts_per_taxa = matrix(0, num_of_taxa, num_of_taxa)
    # for each taxa, count the number of times it appears in each position
    for (i in 1:num_of_taxa){
      positions_for_taxa = df == i
      counts_per_taxa[i,] = colSums(positions_for_taxa)
    }

    # show average of positions for all taxa
    print("Average times a taxa appears in position:")
    print(colMeans(counts_per_taxa))

}






