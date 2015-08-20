
#' get de novo data from the Sanders et al autism exome study
#'
#' Supplementary table 2 (where the excel sheets for the probands and
#' siblings have been combined) from:
#' Sanders et al. (2012) Nature 485:237-241
#' doi: 10.1038/nature10945
#'
#' @export
#' 
#' @return data frame of de novos, with standardised genome coordinates and VEP
#      consequences for each variant
sanders_de_novos <- function() {
    
    url = "http://www.nature.com/nature/journal/v485/n7397/extref/nature10945-s2.xls"
    samples = gdata::read.xls(url, stringsAsFactors=FALSE)
    samples = samples[!(samples$Role %in% c("Mother", "Father")), ]
    samples$sex = tolower(samples$Gender)
    samples$person_id = samples$Sample
    gender = samples[, c("person_id", "sex")]
    
    url = "http://www.nature.com/nature/journal/v485/n7397/extref/nature10945-s3.xls"
    sanders_probands = gdata::read.xls(url, sheet="Probands", stringsAsFactors=FALSE)
    sanders_siblings = gdata::read.xls(url, sheet="Siblings", stringsAsFactors=FALSE)
    variants = rbind(sanders_probands, sanders_siblings)
    
    variants$chrom = gsub("chr", "", variants$Chr.1)
    variants$start_pos = gsub(" ", "", variants$Pos..hg19.)
    variants$end_pos = variants$start_pos
    variants$ref_allele = variants$Ref
    variants$alt_allele = variants$Alt
    
    # get the correct ref and alt alleles for indels
    indels = grep(":", variants$alt_allele)
    temp_distance = nchar(sapply(strsplit(variants$alt_allele[indels], ":"), "[[", 2))
    variants$alt_allele[indels] = apply(variants[indels, ], 1, get_sequence_in_region)
    variants$end_pos[indels] = as.numeric(variants$end_pos[indels]) + temp_distance
    variants$ref_allele[indels] = apply(variants[indels, ], 1, get_sequence_in_region)
    
    variants = fix_het_alleles(variants)
    
    # get the HGNC symbol and VEP consequence for each variant
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    variants$person_id = variants$Child_ID
    variants$study_code = "sanders_nature_2012"
    variants$publication_doi = "10.1038/nature10945"
    variants$study_phenotype = "autism"
    
    variants = merge(variants, gender, by="person_id", all.x=TRUE)
    
    variants = subset(variants, select=c("person_id", "sex", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence",
        "study_code", "publication_doi", "study_phenotype"))
    
    return(variants)
}
