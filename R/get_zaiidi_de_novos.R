
#' get de novo data from Zaidi et al. congenital heart disease exome study
#'
#' De novo mutation data sourced from Supplementary table 4:
#' Zaidi et al. (2013) Nature 498:220-223
#' doi: 10.1038/nature12141
#'
#' @export
#' 
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
zaidi_de_novos <- function() {
    
    url = "http://www.nature.com/nature/journal/v498/n7453/extref/nature12141-s1.pdf"
    
    # obtain the supplementary material
    path = tempfile()
    download.file(url, path)
    
    # extract the supplementary tables from the pdf
    test = tm::readPDF(control=list(text = "-layout"))(elem = list(uri = path), language = "en", id = "id1")
    unlink(path)
    
    table_s4 = test$content[314:912]
    table_s4 = table_s4[table_s4 != ""]
    
    # clean up table S4
    table_s4 = gsub("^[ \t\f]+", "", table_s4) # trim leading whitespace
    table_s4 = gsub(" bp beyond exon ", "_bp_beyond_exon_", table_s4)
    table_s4 = gsub(" bp up of exon ", "_bp_up_of_exon_", table_s4)
    
    # drop the section breaks from the table, as well as the lines that are part
    # of the line breaks
    table_s4 = table_s4[!grepl("[Mm]utations", table_s4)]
    table_s4 = table_s4[!grepl("N A T U R E", table_s4)]
    table_s4 = table_s4[!grepl("SUPPLEMENTARY INFORMATION", table_s4)]
    
    # cull it down to the first few entries in each line
    split_strings = strsplit(table_s4, "[ \t]+")
    variants = data.frame(t(sapply(split_strings, "[", 1:11)))
    names(variants) = c("person_id", "category", "hgnc", "type", "aa_change", "dbSNP", "transcript", "protein", "chrom", "start_pos", "alleles")
    
    # drop out the controls
    variants = variants[variants$category != "Control", ]
    
    variants = fix_zaiidi_coordinates(variants, "alleles")
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    variants$study_code = "zaidi_nature_2013"
    variants$publication_doi = "10.1038/nature12141"
    variants$study_phenotype = "congenital_heart_disease"
    
    # make up sex for Zaiidi samples, in proportion to the sexes in their cases.
    # That way when we exclude samples who are likely diagnostic, we will remove
    # samples in proportion to the study sex ratio.
    person_ids = unique(variants$person_id)
    sex = runif(length(person_ids))
    male_fraction = 220/(220+142)
    sex[sex < male_fraction] = "male"
    sex[sex != "male"] = "female"
    variants$sex = sex[match(variants$person_id, person_ids)]
    
    variants = subset(variants, select=c("person_id", "sex", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence",
        "study_code", "publication_doi", "study_phenotype"))
    
    return(variants)
}
