
#' get de novo data for the Epi4K epilepsy exome study
#'
#' De novo mutation data from the most recent EPI4K publication:
#' Supplementary table 1:
#' American Journal of Human Genetics (2014) 95:360-370
#' doi: 10.1016/j.ajhg.2014.08.013
#'
#' This incorporates the de novo mutation data from supplementary table 2 of:
#' Allen et al. (2013) Nature 501:217-221
#' doi: 10.1038/nature12439
#'
#' @export
#' 
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
epi4k_de_novos <- function() {
    
    system("wget http://catalog.coriell.org/0/Excel/7235051.xls")
    clinical_data = read.table("7235051.xls", sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
    clinical_data$sex = NA
    clinical_data$sex[clinical_data$Gender == "Female"] = "female"
    clinical_data$sex[clinical_data$Gender == "Male"] = "male"
    clinical_data = clinical_data[, c("Ref", "sex")]
    unlink("7235051.xls")
    
    # obtain the supplementary material
    url = "http://www.nature.com/nature/journal/v501/n7466/extref/nature12439-s1.pdf"
    path = tempfile()
    download.file(url, path)
    
    # extract the supplementary tables from the pdf
    test = tm::readPDF(control=list(text = "-layout"))(elem = list(uri = path), language = "en", id = "id1")
    unlink(path)
    
    table_s1 = test$content[11:175]
    table_s1 = gsub("^[ \t\f]+", "", table_s1)
    table_s1 = table_s1[!grepl("N A T U R E", table_s1)]
    table_s1 = table_s1[!grepl("SUPPLEMENTARY INFORMATION", table_s1)]
    table_s1 = table_s1[!grepl("Trio", table_s1)]
    table_s1 = table_s1[!grepl("Coriell ID", table_s1)]
    table_s1 = table_s1[!grepl("Proband", table_s1)]
    table_s1 = table_s1[nchar(table_s1) > 1]
    table_s1 = strsplit(table_s1, "[ ]+")
    table_s1 = data.frame(t(data.frame(table_s1)))
    row.names(table_s1) = 1:nrow(table_s1)
    part1 = table_s1[, 1:4]
    names(part1) = c("TRIO.ID", "proband", "father", "mother")
    part2 = table_s1[, 5:8]
    names(part2) = c("TRIO.ID", "proband", "father", "mother")
    
    samples = rbind(part1, part2)
    samples = merge(samples, clinical_data, by.x="proband", by.y="Ref", all.y=TRUE)
    gender = samples[, c("proband", "TRIO.ID", "sex")]
    
    # now get the supplementary material from the AJHG paper, for the probands
    # who do not have sex info in the Nature paper supplementary material
    url = "http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S0002929714003838/1-s2.0-S0002929714003838-mmc1.pdf/276895/FULL/S0002929714003838/4a001cfb9c56ad9275203cfbe75e9e9e/mmc1.pdf"
    
    path = tempfile()
    download.file(url, path)
    
    # extract the supplementary tables from the pdf
    test = tm::readPDF(control=list(text = "-layout"))(elem = list(uri = path), language = "en", id = "id1")
    unlink(path)
    
    # clean up the text that forms part of supplementary table S2.
    table_s2 = test$content[200:320]
    table_s2 = gsub("^[ \t\f]+", "", table_s2)
    table_s2 = table_s2[!grepl("^Trio", table_s2)]
    table_s2 = strsplit(table_s2, "[ ]+")
    table_s2 = table_s2[sapply(table_s2, length) > 14]
    table_s2 = sapply(table_s2, function(x) x[c(1, 3)])
    
    table_s2 = data.frame(t(table_s2))
    names(table_s2) = c("TRIO.ID", "sex")
    table_s2$TRIO.ID = gsub("(LGS|IS)", "", table_s2$TRIO.ID)
    table_s2$sex = sapply(strsplit(as.character(table_s2$sex), "/"), "[", 1)
    table_s2$sex[table_s2$sex == "F"] = "female"
    table_s2$sex[table_s2$sex == "M"] = "male"
    table_s2$proband = NA
    table_s2 = table_s2[, c("proband", "TRIO.ID", "sex")]
    
    # One of the probands from the supplementary table is mislabelled, this
    # would give the wrong gender to the child. Checking against the excel
    # table from the American Journal of Human Genetics (2014) 95:360-370
    # identifies the correct proband, offset by one, which has the matching de
    # novo in the correct gene.
    table_s2$TRIO.ID[table_s2$TRIO.ID == "ci"] = "cj"
    
    gender = rbind(gender, table_s2)
    
    variants = gdata::read.xls("http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S0002929714003838/1-s2.0-S0002929714003838-mmc2.xlsx/276895/FULL/S0002929714003838/bf21945d72e3297fc44969dc0296f4f1/mmc2.xlsx", stringsAsFactors=FALSE)
    
    # the excel table contains three final tail rows which do not contain
    # tabular data, so we remove these
    variants = variants[1:(nrow(variants) - 3), ]
    
    variants = fix_coordinates_with_allele(variants,
        "hg19.coordinates..chr.position.", "Ref.Alt.alleles")
    
    variants$TRIO.ID = gsub("\\*", "", variants$TRIO.ID)
    variants$Child.ID = gsub("\\*", "", variants$Child.ID)
    
    # get the hgnc symbol, and clean any anomalies
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    variants$person_id = variants$Child.ID
    variants$study_code = "epi4k_ajhg_2014"
    variants$publication_doi = "10.1016/j.ajhg.2014.08.013"
    variants$study_phenotype = "epilepsy"
    
    # get a set of IDs that match the Coriell IDs
    matches = regexpr("ND[0-9]*", toupper(variants$Child.ID))
    start = as.vector(matches)
    stop = start + attributes(matches)$match.length - 1
    variants$altered_id = substr(toupper(variants$Child.ID), start, stop)
    variants$altered_id[variants$altered_id == ""] = toupper(variants$Child.ID)[variants$altered_id == ""]
    
    # get the sex code, using the different types of ID
    variants = merge(variants, gender[, c("TRIO.ID", "sex")], by="TRIO.ID", all.x=TRUE)
    variants = merge(variants, gender[, c("proband", "sex")], by.x="altered_id", by.y="proband", all.x=TRUE)
    
    # collapse the sex codes into a single value
    variants$sex = apply(variants[, c("sex.x", "sex.y")], 1, function(x) unique(na.omit(x)))
    variants$sex[sapply(variants$sex, function(x) identical(x, character(0)))] = NA
    variants$sex = unlist(variants$sex)
    
    variants = subset(variants, select=c("person_id", "sex", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence",
        "study_code", "publication_doi", "study_phenotype"))
    
    return(variants)
}
