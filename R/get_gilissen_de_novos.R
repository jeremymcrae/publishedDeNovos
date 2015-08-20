
#' get de novo data for Gilissen et al. intellectual disability exome study
#'
#' De novo mutation data sourced from supplementary table 8 from:
#' Gilissen et al. (2014) Nature 511:344-347
#' doi: 10.1038/nature13394
#'
#' NOTE: the number of males and females will form part of the De Ligt counts,
#' NOTE: since the Gilissen samples were in the De Ligt study, they just didn't
#' receive a diagnosis in the original study
#'
#' @param deligt dataframe of de novos variants from the De Ligt N Engl J Med
#'        367:1921-1929 study, so we can get the sample IDs and sex for the
#'        current dataset.
#' @param deligt_sex dataframe of sample IDs and sex codes for all of the
#'        children in the de ligt study, not just those who had de novo variants.
#'
#' @export
#'
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
gilissen_de_novos <- function(deligt, deligt_sex) {
    
    url = "http://www.nature.com/nature/journal/v511/n7509/extref/nature13394-s1.pdf"
    
    # obtain the supplementary material
    path = tempfile()
    download.file(url, path)
    
    # extract the supplementary tables from the pdf
    test = tm::readPDF(control=list(text = "-layout"))(elem = list(uri = path), language = "en", id = "id1")
    unlink(path)
    
    # get the lines corresponding to each table
    table_s8 = test$content[1434:1563]
    
    # clean up table S2
    table_s8 = gsub("^[ \t]+", "", table_s8) # trim leading whitespace
    
    # fix the lines that come after the page breaks
    table_s8 = table_s8[table_s8 != ""]
    table_s8 = table_s8[table_s8 != "#"]
    table_s8 = table_s8[!grepl("NATURE.COM", table_s8)]
    table_s8 = table_s8[!grepl("SUPPLEMENTARY INFORMATION", table_s8)]
    
    # one line in the table is spread across three lines, rather than fixing
    # this in code, just re-insert the correct line
    table_s8[50] = "25  SATB2     Chr2(GRCh37):g.200213667_200213668insGTTGCCTTACAA    NM_001172517.1:c.929_930insTTGTAAGGCAAC  p.Gln310delinsHis_CysLys_AlaThr  Insertion  no  K  F  D  C  B  G  M  Known"
    
    # trim the too short lines (which only contain page numbers)
    table_s8 = table_s8[lapply(table_s8, nchar) > 5]
    
    # split the table, and drop the erroneous lines from the bad line
    split_strings = strsplit(table_s8, "[ \t]+")
    split_strings[50] = NULL
    split_strings[49] = NULL
    
    # cull it down to the first few entries in each line
    variants = data.frame(t(sapply(split_strings, "[", 1:6)))
    names(variants) = c("person_id", "hgnc", "hgvs_genomic", "hgvs_transcript", "hgvs_protein", "type")
    
    # fix the hgvs genomic string
    variants$hgvs_genomic = gsub("\\(GRCh37\\)g", ":g", variants$hgvs_genomic)
    variants$hgvs_genomic = gsub("\\(GRCh37\\)", "", variants$hgvs_genomic)
    variants$hgvs_genomic = gsub("\\(GRCH37\\)", "", variants$hgvs_genomic)
    variants$hgvs_genomic = gsub("Chr", "chr", variants$hgvs_genomic)
    variants$hgvs_genomic = gsub("-", "_", variants$hgvs_genomic)
    
    variants = fix_coordinates_with_hgvs_genomic(variants, "hgvs_genomic")
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    variants$study_code = "gilissen_nature_2014"
    variants$publication_doi = "10.1038/nature13394"
    variants$study_phenotype = "intellectual_disability"
    
    variants = merge(variants, deligt_sex, by="person_id", all.x=TRUE)
    
    variants = subset(variants, select=c("person_id", "sex", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence",
        "study_code", "publication_doi", "study_phenotype"))
    
    # remove the variants that were present in the De Ligt study
    temp = rbind(deligt, variants)
    temp$dups = duplicated(temp[, 2:7])
    variants = variants[!(temp$dups[temp$study_code == "gilissen_nature_2014"]), ]
    
    return(variants)
}
