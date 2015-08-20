
#' get de novo data from the O'Roak et al autism exome study
#'
#' Supplementary table 3 from:
#' O'Roak et al. (2012) Nature 485:246-250
#' doi: 10.1038/nature10989
#'
#' @export
#' 
#' @return data frame of de novos, with standardised genome coordinates and VEP
#      consequences for each variant
oroak_de_novos <- function() {
    url = "http://www.nature.com/nature/journal/v485/n7397/extref/nature10989-s2.xls"
    variants = gdata::read.xls(url, sheet="Supplementary Table 3", stringsAsFactors=FALSE)
    samples = gdata::read.xls(url, sheet="Supplementary Table 1", stringsAsFactors=FALSE)
    samples$person_id = samples$child
    gender = samples[, c("person_id", "sex")]
    
    # the excel table contains two final tail rows which do not contain
    # tabular data, so we remove these
    variants = variants[1:(nrow(variants) - 2), ]
    
    # standardise the chrom, position and allele column names
    variants$chrom = gsub(" ", "", variants$Chromosome)
    variants$start_pos = gsub(" ", "", variants$Position..hg19.)
    variants$end_pos = variants$start_pos
    variants$ref_allele = gsub(" ", "", variants$Ref)
    variants$alt_allele = gsub(" ", "", variants$Allele)
    
    # sort out the "complex" allele events
    variants$alt_allele[61] = "S"
    variants$alt_allele[78] = "1D, -G"
    variants$alt_allele[116] = "1I, +C"
    variants$alt_allele[133] = "R"
    variants$alt_allele[185] = "1D, -G"
    
    # fix the alleles and positions for insertions and deletions
    deletions = grep("D, *-", variants$alt_allele)
    insertions = grep("I, *\\+", variants$alt_allele)
    
    # find the reference sequence at the site. Deletions use this as the
    # alternate allele, whereas the insertions use this as the reference allele
    alt_dels = apply(variants[deletions, ], 1, get_sequence_in_region)
    ref_ins = apply(variants[insertions, ], 1, get_sequence_in_region)
    
    # find the sequence at the site + the distance of the deletion
    variants$ref_allele[deletions] = paste(alt_dels, sapply(strsplit(variants$alt_allele[deletions], "D, *-"), "[", 2), sep="")
    variants$alt_allele[deletions] = alt_dels
    variants$ref_allele[insertions] = ref_ins
    variants$alt_allele[insertions] = paste(ref_ins, sapply(strsplit(variants$alt_allele[insertions], "I, *\\+"), "[", 2), sep="")
    
    # get the end coordinate, including those for the insertions and deletions
    variants$end_pos[deletions] = as.numeric(variants$end_pos[deletions]) + nchar(variants$ref_allele[deletions]) - 1
    
    variants = fix_het_alleles(variants)
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = sapply(vep, "[", 1)
    variants$hgnc = sapply(vep, "[", 2)
    
    variants$person_id = variants$Person
    variants$study_code = "oroak_nature_2012"
    variants$publication_doi = "10.1038/nature10989"
    variants$study_phenotype = "autism"
    
    variants = merge(variants, gender, by="person_id", all.x=TRUE)
    
    variants = subset(variants, select=c("person_id", "sex", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence",
        "study_code", "publication_doi", "study_phenotype"))
    
    return(variants)
}
