
#' get de novo data from the 2014 De Rubeis et al. autism exome study in Nature
#'
#' De novo mutation data sourced from Supplementary table 3:
#' De Rubeis et al. (2013) Nature 515:209-215
#' doi: 10.1038/nature13772
#'
#' @export
#' 
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
de_rubeis_de_novos <- function() {
    url = "http://www.nature.com/nature/journal/v515/n7526/extref/nature13772-s4.xlsx"
    
    variants = gdata::read.xls(url, sheet="De Novo", stringsAsFactors=FALSE)
    
    # exclude the final row, which contains a footnote
    variants = variants[1:(nrow(variants) - 1), ]
    
    # rename columns to match the other de novo datasets, and strip whitespace
    variants$start_pos = gsub("[ \t]", "", variants$Pos)
    variants$chrom = gsub("[ \t]", "", variants$Chr)
    variants$person_id = gsub("[ \t]", "", variants$Child_ID)
    variants$ref_allele = gsub("[ \t]", "", variants$Ref)
    variants$alt_allele = gsub("[ \t]", "", variants$Alt)
    
    # get the end position
    variants$end_pos = as.character(as.numeric(variants$start_pos) + nchar(variants$ref_allele) - 1)
    
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    # figure out the sex for each variant
    variants$sex = variants$Child_Sex
    variants$sex[variants$sex == "2"] = "female"
    variants$sex[variants$sex == "1"] = "male"
    
    variants$study_code = "derubeis_nature_2014"
    variants$publication_doi = "10.1038/nature13772"
    variants$study_phenotype = "autism"
    
    variants = subset(variants, select=c("person_id", "sex", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence",
        "study_code", "publication_doi", "study_phenotype"))
    
    return(variants)
}
