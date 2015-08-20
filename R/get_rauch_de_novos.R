
#' get de novo data for Rauch et al. intellectual disability exome study
#'
#' De novo mutation data sourced from supplementary tables 2 and 3 from
#' Rauch et al. (2012) Lancet 380:1674-1682
#' doi: 10.1016/S0140-6736(12)61480-9
#'
#' @export
#'
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
rauch_de_novos <- function() {
    
    url = "http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S0140673612614809/1-s2.0-S0140673612614809-mmc1.pdf/271074/FULL/S0140673612614809/55b26043f4a279334b3a5ec00b9faf4b/mmc1.pdf"
    
    # obtain the supplementary material
    path = tempfile()
    download.file(url, path)
    
    # extract the supplementary tables from the pdf
    test = tm::readPDF(control=list(text = "-layout"))(elem = list(uri = path), language = "en", id = "id1")
    unlink(path)
    
    # get the lines corresponding to each table
    table_s1_text = test$content[333:529]
    table_s2_text = test$content[769:847]
    table_s3_text = test$content[879:890]
    
    # clean up table S1
    table_s1_text = gsub("^[ \t\f]+", "", table_s1_text) # trim leading whitespace
    split_strings = strsplit(table_s1_text, "[ \t]{2,}")
    split_strings = split_strings[sapply(split_strings, length) > 1]
    split_strings = lapply(split_strings, function(x) x[1:3])
    table_s1 = data.frame(t(data.frame(split_strings)))
    row.names(table_s1) = 1:nrow(table_s1)
    
    # only select the lines with a single-letter sex code
    use_lines = apply(table_s1, 1, function(x) x[2] %in% c("f", "m") | x[3] %in% c("f", "m"))
    table_s1 = table_s1[use_lines, ]
    table_s1$sex = apply(table_s1, 1, function(x) x[x == "f" | x == "m"])
    names(table_s1) = c("person_id", "column2", "column3", "sex")
    table_s1 = table_s1[, c("person_id", "sex")]
    table_s1$sex = gsub("m", "male", table_s1$sex)
    table_s1$sex = gsub("f", "female", table_s1$sex)
    table_s1$person_id = gsub("-", "", table_s1$person_id)
    table_s1$person_id = as.character(table_s1$person_id)
    table_s1$person_id[table_s1$person_id == "TUTLN112014"] = "TUTLN"
    
    # clean up table S2
    table_s2_text = gsub("^[ \t\f]+", "", table_s2_text) # trim leading whitespace
    split_strings = strsplit(table_s2_text, "[ \t]+")
    split_strings = iconv(unlist(split_strings), "latin1", "ASCII", sub="")
    table_s2 = as.data.frame(matrix(split_strings, ncol=12, byrow=TRUE))
    names(table_s2) = c("person_id", "hgnc", "type", "hgvs_genomic")
    
    # clean up table S3
    table_s3_text = gsub("^[ \t\f]+", "", table_s3_text) # trim leading whitespace
    split_strings = strsplit(table_s3_text, "[ \t]+")
    split_strings[7][[1]] = c(split_strings[[7]], "") # one row lacks an entry
    split_strings = iconv(unlist(split_strings), "latin1", "ASCII", sub="")
    table_s3 = as.data.frame(matrix(split_strings, ncol=9, byrow=TRUE))
    names(table_s3) = c("person_id", "hgnc", "hgvs_genomic")
    table_s3$type = "synonymous"
    
    # standardise and merge the two tables
    table_s2 = subset(table_s2, select=c("person_id", "hgnc", "type", "hgvs_genomic"))
    table_s3 = subset(table_s3, select=c("person_id", "hgnc", "type", "hgvs_genomic"))
    variants = rbind(table_s2, table_s3)
    
    variants = fix_coordinates_with_hgvs_genomic(variants, "hgvs_genomic")
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    # include the sex codes for each de novo
    variants = merge(variants, table_s1, by="person_id", all.x=TRUE)
    
    # define the study details
    variants$study_code = "rauch_lancet_2012"
    variants$publication_doi = "10.1016/S0140-6736(12)61480-9"
    variants$study_phenotype = "intellectual_disability"
    
    variants = subset(variants, select=c("person_id", "sex", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence",
        "study_code", "publication_doi", "study_phenotype"))
    
    return(variants)
}
