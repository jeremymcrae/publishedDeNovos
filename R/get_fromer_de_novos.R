
#' get de novo data from Fromer et al. schizophrenia exome study
#'
#' De novo mutation data sourced from Supplementary table 1:
#' Fromer et al. (2014) Nature 506:179-184
#' doi: 10.1038/nature12929
#'
#' @export
#' 
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
fromer_de_novos <- function() {
    
    url = "http://www.nature.com/nature/journal/v506/n7487/extref/nature12929-s2.xlsx"
    variants = gdata::read.xls(url, stringsAsFactors=FALSE)
    
    variants$temp = strsplit(variants$Locus, ":")
    variants$chrom = sapply(variants$temp, "[", 1)
    variants$chrom = gsub("chr", "", variants$chrom)
    variants$start_pos = sapply(variants$temp, "[", 2)
    variants$end_pos = variants$start_pos
    
    # fix indel ranges
    indels = grepl("\\.\\.", variants$start_pos)
    variants$start_pos[indels] = sapply(strsplit(variants$start_pos[indels], "\\.\\."), "[", 1)
    variants$end_pos[indels] = sapply(strsplit(variants$end_pos[indels], "\\.\\."), "[", 2)
    
    variants$ref_allele = variants$Reference.allele
    variants$alt_allele = variants$Alternate.allele
    
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    variants$person_id = variants$Proband.ID
    variants$study_code = "fromer_nature_2014"
    variants$publication_doi = "10.1038/nature12929"
    variants$study_phenotype = "schizophrenia"
    
    # get the sex for each variant
    variants$sex = variants$Proband.gender
    
    variants = subset(variants, select=c("person_id", "sex", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence",
        "study_code", "publication_doi", "study_phenotype"))
    
    return(variants)
}
