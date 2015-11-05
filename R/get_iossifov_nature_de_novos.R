
#' get de novo data from the 2014 Iossifov et al. autism exome study in Nature
#'
#' De novo mutation data sourced from Supplementary table 2:
#' Iossifov et al. (2014) Nature 498:216-221
#' doi: 10.1038/nature13908
#'
#' @export
#'
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
iossifov_nature_de_novos <- function() {
    tmpdir = tempdir()
    path = tempfile(tmpdir=tmpdir)
    
    # obtain the dataframe of de novo variants
    download.file("http://www.nature.com/nature/journal/v515/n7526/extref/nature13908-s2.zip", path)
    unzip(path, files=c("nature13908-s2/Supplementary Table 1.xlsx",
        "nature13908-s2/Supplementary Table 2.xlsx"), exdir=tmpdir)
    variants = gdata::read.xls(file.path(tmpdir, "nature13908-s2", "Supplementary Table 2.xlsx"), stringsAsFactors=FALSE)
    families = gdata::read.xls(file.path(tmpdir, "nature13908-s2", "Supplementary Table 1.xlsx"), stringsAsFactors=FALSE)
    unlink(path)
    
    variants = fix_coordinates(variants, "location", "vcfVariant")
    
    # exclude the de novos from the unaffected sibs
    variants = variants[!grepl("^s", variants$inChild), ]
    
    # get the sex of the probands
    variants$sex = substr(variants$inChild, 2, 2)
    variants$sex[variants$sex == "M"] = "male"
    variants$sex[variants$sex == "F"] = "female"
    
    # NOTE: the variant with the most severe consequence might not necessarily
    # NOTE: be within the listed gene. I haven't accounted for this yet.
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    variants$person_id = variants$familyId
    variants$study_code = "iossifov_nature_2014"
    variants$publication_doi = "10.1038/nature13908"
    variants$study_phenotype = "autism"
    
    # identify the family IDs of the probands with low IQ, as we shall
    # reclassify these probands as having intellectual disability, since these
    # probands are more likely to be enriched for de novos for intellectual
    # diability than for autism.
    low_iq = families$familyId[families$probandVIQ < 70 & families$probandNVIQ < 70]
    low_iq = low_iq[!is.na(low_iq)]
    variants$study_phenotype[variants$person_id %in% low_iq] = "intellectual_disability"
    
    variants = subset(variants, select=c("person_id", "sex", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence",
        "study_code", "publication_doi", "study_phenotype"))
    
    return(variants)
}
