
#' get de novo data for Lelieveld et al. intellectual disability exome study
#'
#' De novo mutation data sourced from supplementary table 1 from:
#' Lelieveld et al. (2016) Nature Neuroscience 19:1194-1196
#' doi: 10.1038/nn.4352
#'
#' Note that the paper says that the variants were aligned to hg19, but their
#' table of de novo variants is definitely for hg18 (GRCh37).
#'
#' @export
#'
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
lelieveld_de_novos <- function() {
    
    url = "http://www.nature.com/neuro/journal/v19/n9/extref/nn.4352-S3.xlsx"
    variants = gdata::read.xls(url, stringsAsFactors=FALSE)
    variants = as.data.frame(variants)
    
    variants['person_id'] = as.character(unlist(variants['Sample.ID']))
    variants['hgnc'] = as.character(unlist(variants['Gene.name']))
    variants['chrom'] = as.character(gsub('chr', '', unlist(variants['Chromosome'])))
    variants['start_pos'] = variants['Start.position']
    variants['end_pos'] = variants['End.position']
    variants['ref_allele'] = as.character(unlist(variants['Reference.Allele']))
    variants['alt_allele'] = as.character(unlist(variants['Variant.Allele']))
    
    # get the reference allele for insertions
    ins = variants['ref_allele'] == ''
    seq = as.character(apply(variants[ins, ], 1, get_sequence_in_region, verbose=TRUE))
    variants['ref_allele'][ins] = seq
    variants['alt_allele'][ins] = paste(seq, as.character(unlist(variants['alt_allele'][ins])), sep='')
    
    dels = variants['alt_allele'] == ''
    variants['alt_allele'][dels] = '-'
    
    vep = apply(variants, 1, get_vep_consequence, verbose=TRUE)
    variants$consequence = unlist(sapply(vep, "[", 1))
    variants$hgnc = unlist(sapply(vep, "[", 2))
    
    variants$study_code = "lelieveld_nature_neuroscience_2016"
    variants$publication_doi = "10.1038/nn.4352"
    variants$study_phenotype = "intellectual_disability"
    
    # make up sex for the samples, in proportion to the sexes in their cases.
    # That way when we exclude samples who are likely diagnostic, we will remove
    # samples in proportion to the study sex ratio.
    person_ids = unique(variants$person_id)
    sex = runif(length(person_ids))
    male_fraction = 461/(461+359)
    sex[sex < male_fraction] = "male"
    sex[sex != "male"] = "female"
    variants$sex = sex[match(variants$person_id, person_ids)]
    
    variants = subset(variants, select=c("person_id", "sex", "chrom", "start_pos",
        "end_pos", "ref_allele", "alt_allele", "hgnc", "consequence",
        "study_code", "publication_doi", "study_phenotype"))
    
    return(variants)
}
