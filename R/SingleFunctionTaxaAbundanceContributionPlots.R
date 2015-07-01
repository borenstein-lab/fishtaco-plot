#######################################################################
# CONTRIBUTION TO FUNCTION MEASURED ABUNDANCE
#
# plots one of the following plots for a single function:
# - Steps plot with mean and std contributions for all taxa
# - Ring plot where the contributions are shown as a ring
# - Bar plot where the contributions are shown as a bar
#######################################################################

SingleFunctionTaxaAbundanceContributionPlots <- function(input_function=NULL,
                                                         input_genomic_content=NULL,
                                                         input_function_abundance=NULL,
                                                         input_taxa_abundance=NULL,
                                                         input_sample_label=NULL) {



  # # sample data
#   input_function="ko00020"
#   input_genomic_content="/Volumes/ohadm/OhadM/METAFIT/HMP_DATA/METAPHLAN_taxa_vs_PATHWAY.tab"
#   input_function_abundance="/Volumes/ohadm/OhadM/METAFIT/HMP_DATA/WGS_PATHWAY_vs_SAMPLE_MUSiCC.tab"
#   input_taxa_abundance="/Volumes/ohadm/OhadM/METAFIT/HMP_DATA/TONGUE_DORSUM_vs_BUCCAL_MUCOSA/METAPHLAN_taxa_vs_SAMPLE.tab"
#   input_sample_label="/Volumes/ohadm/OhadM/METAFIT/HMP_DATA/TONGUE_DORSUM_vs_BUCCAL_MUCOSA/SAMPLE_vs_CLASS.tab"

  # read files
  taxa_abundance = read.table(input_taxa_abundance, sep = "\t", header=TRUE, stringsAsFactors=FALSE)
  samples_with_taxa_abun = names(taxa_abundance)[-1]
  row.names(taxa_abundance) = taxa_abundance[, 1]
  taxa_abundance[, 1] = NULL

  # first, cut out the measured function abundance for only this function
  # and select only intersected samples
  function_abundance = read.table(input_function_abundance, sep = "\t", header=TRUE, stringsAsFactors=FALSE)
  samples_with_function_abun = names(function_abundance)[-1]
  samples_with_both_info = intersect(samples_with_taxa_abun, samples_with_function_abun)
  num_of_samples = length(samples_with_both_info)
  #print(paste("number of samples with both info: ", num_of_samples))
  row.names(function_abundance) = function_abundance[, 1]
  function_abundance[, 1] = NULL
  function_abundance = melt(t(function_abundance[input_function, samples_with_both_info]))
  function_abundance = function_abundance[,c(1,3)]
  names(function_abundance) = c("SampleName", "Measured")

  # filter taxa abundance for intersecting samples
  # sort samples by their function abundace value
  sample_label = read.table(input_sample_label, sep = "\t", header=TRUE, stringsAsFactors=FALSE)
  sample_label = sample_label[sample_label$Sample %in% samples_with_both_info,]
  names(sample_label) = c("Sample", "label")
  sample_cases = sample_label$Sample[sample_label$label == 1]
  sample_control = sample_label$Sample[sample_label$label == 0]
  number_of_cases = length(sample_cases)
  index_cases_in_function_abun = which(function_abundance$SampleName %in% sample_cases)
  index_control_in_function_abun = which(function_abundance$SampleName %in% sample_control)
  sort_index_cases = index_cases_in_function_abun[sort(function_abundance$Measured[function_abundance$SampleName %in% sample_cases], decreasing=TRUE, index.return=TRUE)$ix]
  sort_index_control = index_control_in_function_abun[sort(function_abundance$Measured[function_abundance$SampleName %in% sample_control], decreasing=TRUE, index.return=TRUE)$ix]
  sort_index = c(sort_index_cases, sort_index_control)
  taxa_abundance = taxa_abundance[, samples_with_both_info[sort_index]]

  # first, cut out the genomic content for only this function
  genomic_content = read.table(input_genomic_content, sep = "\t", header=TRUE, stringsAsFactors=FALSE)
  function_ind = which(names(genomic_content) == input_function)
  genomic_content = genomic_content[, c(1,function_ind)]
  row.names(genomic_content) = genomic_content[, 1]
  genomic_content[, 1] = NULL

  # now merge taxa abundance with genomic content
  merge_content_abun = merge(genomic_content, taxa_abundance, by="row.names")
  for (i in 1:num_of_samples) {
    merge_content_abun[, (i+2)] = merge_content_abun[, 2] * merge_content_abun[, (i+2)]
  }
  merge_content_abun = merge_content_abun[, c("Row.names", samples_with_both_info[sort_index])]
  merge_content_abun.m = melt(merge_content_abun)
  names(merge_content_abun.m) = c("Taxa","SampleName","Contribution")
  merge_content_abun.m$Taxa = factor(merge_content_abun.m$Taxa)

  # plot resulting figure
  BarPlot = ggplot() +
    geom_bar(data=merge_content_abun.m, aes(x=SampleName,y=Contribution, fill=Taxa), stat="identity", colour="black") +
    geom_point(data=function_abundance, aes(x=SampleName, y=Measured), colour="black", fill="white", size=4, shape=23) +
    ylab("Contribution to function abundance") + xlab("Samples") +
    guides(fill = guide_legend(ncol = 3)) +
    geom_vline(xintercept = number_of_cases + 0.5, size = 2, linetype = "longdash") +
    theme_bw() + theme(axis.text.x = element_blank(),
                       axis.title = element_text(size=20), axis.text.y = element_text(size=15))

  return(BarPlot)
}
