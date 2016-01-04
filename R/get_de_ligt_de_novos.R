
#' get de novo data for De Ligt et al. intellectual disability exome study
#'
#' De novo mutation data sourced from supplementary table 3 from
#' De Ligt et al. (2012) N Engl J Med 367:1921-1929
#' doi: 10.1056/NEJMoa1206524
#'
#' @export
#'
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
deligt_de_novos <- function() {
    
    url = "http://www.nejm.org/doi/suppl/10.1056/NEJMoa1206524/suppl_file/nejmoa1206524_appendix.pdf"
    pdf = get_deligt_text(url)
    gender = get_deligt_sex(pdf)
    variants = get_deligt_denovo_table(pdf)
    
    variants = fix_coordinates_with_hgvs_genomic(variants, "hgvs_genomic")
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    variants = merge(variants, gender, by="person_id")
    
    variants$study_code = "deligt_nejm_2012"
    variants$publication_doi = "10.1056/NEJMoa1206524"
    variants$study_phenotype = "intellectual_disability"
    
    variants = subset(variants, select=c("person_id", "sex", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence",
        "study_code", "publication_doi", "study_phenotype"))
    
    return(list(variants=variants, gender=gender))
}

#' load the supplementary data for the de ligt de novos
#'
#' @param url url of the de ligt supplementary pdf.
#'
#' @return list of text in supplementary file
get_deligt_text <- function(url) {
    
    temp = "temp.html"
    cookie = "cookie.txt"
    system(paste("wget --cookies=on --keep-session-cookies --save-cookies=", cookie, " ", url, " -O", temp, sep=""))
    
    # obtain the supplementary material
    path = "corrupted_pdf.pdf"
    system(paste("wget --referer=", url, " --cookies=on --load-cookies=", cookie, " --keep-session-cookies --save-cookies=", cookie, " ", url, " -O ", path, sep=""))
    
    # repair the pdf with ghostscript
    repaired = "repaired.pdf"
    system(paste("gs -o ", repaired, " -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress ", path, sep=""))
    
    test = tm::readPDF(control=list(text = "-layout"))(elem = list(uri=repaired), language = "en", id = "id1")
    
    # delete the pdf and temporary files
    unlink(temp)
    unlink(cookie)
    unlink(path)
    unlink(repaired)
    
    return(test)
}

#' load the supplementary table for the de ligt sex information
#'
#' @param pdf list of text lines in supplementary pdf.
#'
#' @return dataframe of person_ids and sex codes
get_deligt_sex <- function(pdf) {
    start_line = 174
    end_line = 1160
    text = pdf$content[start_line:end_line]
    
    # trim leading whitespace
    text = gsub("^[ \t\f]+", "", text)
    gender = data.frame(text=text[grepl("^Trio", text)], sex=NA)
    gender$text = as.character(gender$text)
    
    # standardise the sex codes, based on the table contents
    gender$sex[grepl("[Tt]his (boy|male)+", gender$text)] = "male"
    gender$sex[grepl("[Tt]his (girl|female)+", gender$text)] = "female"
    
    # figure out the individual IDs
    gender$person_id = sapply(strsplit(gender$text, "( |-)+"), "[", 2)
    gender$sex[gender$person_id == "69"] = "female"
    
    return(gender)
}

#' load the supplementary table for the de ligt de novo variants
#'
#' @param pdf list of text lines in supplementary pdf.
#'
#' @return list of text in supplementary file
get_deligt_denovo_table <- function(pdf) {
    # clean up table S3
    start_line = 1709
    end_line = 1795
    text = pdf$content[start_line:end_line]
    
    # trim leading whitespace, then remove blank lines, commented out lines,
    # and short lines
    text = gsub("^[ \t\f]+", "", text)
    text = text[text != ""]
    text = text[text != "#"]
    text = text[lapply(text, nchar) > 5]
    split_strings = strsplit(text, "[ \t]+")
    
    # cull it down to the first few entries in each line
    variants = data.frame(t(sapply(split_strings, "[", 1:6)))
    names(variants) = c("person_id", "hgnc", "hgvs_genomic", "transcript", "hgvs_cdna", "hgvs_protein")
    
    # fix the hgvs genomic string
    variants$hgvs_genomic = gsub("\\(GRCh37\\)g", ":g", variants$hgvs_genomic)
    variants$hgvs_genomic = gsub("\\(GRCh37\\)", "", variants$hgvs_genomic)
    variants$hgvs_genomic = gsub("\\(GRCH37\\)", "", variants$hgvs_genomic)
    variants$hgvs_genomic = gsub("Chr", "chr", variants$hgvs_genomic)
    variants$hgvs_genomic = gsub("-", "_", variants$hgvs_genomic)
    
    # convert a position with NCBI36 genome assembly coordinates
    variants$hgvs_genomic[27] = "chr19:g.53958839G>C"
    
    return(variants)
}
