
#' get de novo data from the 2012 Iossifov et al autism exome study in Neuron
#'
#' Supplementary table 1 (where the non-coding SNVs have been excluded) and
#' supplementary table 2 from:
#' Iossifov et al. (2012) Neuron 74:285-299
#' doi: 10.1016/j.neuron.2012.04.009
#'
#' @export
#'
#' @return data frame of de novos, with standardised genome coordinates and VEP
#      consequences for each variant
iossifov_neuron_de_novos <- function() {
    url = "http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S0896627312003406/1-s2.0-S0896627312003406-mmc2.xlsx/272195/FULL/S0896627312003406/26c5ba3b72a2410ef43fec52a40f35e6/mmc2.xlsx"
    snvs = gdata::read.xls(url, sheet="SNV.v4.1-normlized", stringsAsFactors=FALSE)
    
    url = "http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S0896627312003406/1-s2.0-S0896627312003406-mmc4.xlsx/272195/FULL/S0896627312003406/6caa42b35609c2ed5910b5381ddd5335/mmc4.xlsx"
    indels = gdata::read.xls(url, sheet="ID.v4.1-normlized", stringsAsFactors=FALSE)
    
    # trim out the low quality de novos (as defined by a flag in the table)
    snvs = snvs[snvs$SNVFilter == 1, ]
    indels = indels[indels$IndelFilter == 1, ]
    
    # merge the SNV and indel de novo calls
    snvs = subset(snvs, select = c("quadId", "location", "variant", "effectGenes", "effectType", "inChild"))
    indels = subset(indels, select = c("quadId", "location", "variant", "effectGenes", "effectType", "inChild"))
    variants = rbind(snvs, indels)
    
    # get the coordinates and VEP consequence
    variants = fix_coordinates_with_allele(variants, "location", "variant")
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    variants$person_id = variants$quadId
    variants$study_code = "iossifov_neuron_2012"
    variants$publication_doi = "10.1016/j.neuron.2012.04.009"
    variants$study_phenotype = "autism"
    
    variants$sex = substr(variants$inChild, 4, 4)
    variants$sex[variants$sex == "M"] = "male"
    variants$sex[variants$sex == "F"] = "female"
    
    variants = subset(variants, select=c("person_id", "sex", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence",
        "study_code", "publication_doi", "study_phenotype"))
    
    return(variants)
}
