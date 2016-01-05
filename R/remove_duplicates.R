# consequence list, as sorted at http://www.ensembl.org/info/genome/variation/predicted_data.html
consequences = c("transcript_ablation", "splice_donor_variant",
    "splice_acceptor_variant", "stop_gained", "frameshift_variant",
    "start_lost", "initiator_codon_variant", "stop_lost",
    "transcript_amplification", "inframe_insertion", "inframe_deletion",
    "missense_variant", "protein_altering_variant", "splice_region_variant",
    "incomplete_terminal_codon_variant", "stop_retained_variant",
    "synonymous_variant", "coding_sequence_variant", "mature_miRNA_variant",
    "5_prime_UTR_variant", "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant", "intron_variant",
    "NMD_transcript_variant", "non_coding_transcript_variant",
    "upstream_gene_variant", "downstream_gene_variant", "TFBS_ablation",
    "TFBS_amplification", "TF_binding_site_variant",
    "regulatory_region_ablation", "regulatory_region_amplification",
    "regulatory_region_variant", "feature_elongation", "feature_truncation",
    "intergenic_variant")
severity = data.frame(consequence=consequences, rank=seq(1:length(consequences)),
    stringsAsFactors=FALSE)

#' remove duplicate variants from the same individual
#'
#' The autism de novos come from multiple studies. Some probands are in multiple
#' studies, so we have to remove duplicate entries, while accounting for
#' potential differences in sample IDs, and possibly different coordinates (for
#' indels).
#'
#' @param de_novos dataframe of variants
#'
#' @return dataframe of variants, without duplicates
#' @export
#'
#' @examples
#' variants = read.table(header=TRUE, text="
#'     person_id    chrom start_pos ref_allele alt_allele consequence        study_phenotype
#'     sample_01    1     10000     AAA        A          frameshift_variant autism
#'     sample_01.p1 1     10001     AAA        A          frameshift_variant normal_iq_autism
#'     sample_02    1     10001     AAA        A          frameshift_variant normal_iq_autism",
#'     stringsAsFactors=FALSE)
#' remove_duplicates(variants)
remove_duplicates <- function(de_novos) {
    
    # remove de novos that have been duplicated between studies. These are easily
    # spotted as individuals who have IDs that are nearly identical between
    # different studies eg 14323.p1 vs 14323.
    de_novos$temp_id = sapply(strsplit(de_novos$person_id, "\\."), "[", 1)
    
    for (pos in 1:nrow(de_novos)) {
        row = de_novos[pos, ]
        # Sometimes the positions differs slightly between studies, so we check
        # for near matches
        matches = find_matching_sites(row, de_novos)
        
        if (sum(matches) > 1) {
            # make sure the probands from the Iossifox Nature 2014 paper that
            # have been classified as normal_iq_autism, transfer their
            # classifications to the earlier de novo entries in case we select
            # that row.
            if (any(de_novos$study_phenotype[matches] == "normal_iq_autism")) {
                de_novos$study_phenotype[matches] = "normal_iq_autism"
            }
            
            most_severe = unique(de_novos$consequence[matches])
            if (length(most_severe) > 1) {
                sev = min(severity$rank[severity$consequence %in% most_severe])
                most_severe = severity$consequence[severity$rank == sev]
            }
            
            # pick one of the rows, and swap it in for all the possible matches.
            # This will ensure we can exclude the duplicates. Some entries will
            # get checked twice (due to the duplicate row), but there aren't
            # enough to worry about doing this in a single pass.
            severe_row = matches & as.vector(unlist(de_novos$consequence)) == unlist(most_severe)
            de_novos[matches, ] = de_novos[rep(which(severe_row)[1], sum(matches)), ]
        }
    }
    
    # remove duplicates, then trim the temporary ID
    de_novos = de_novos[!duplicated(de_novos[, c("temp_id", "chrom", "start_pos")]), ]
    de_novos$temp_id = NULL
    
    return(de_novos)
}

#' this identifies variants which match the current row
#'
#' Some variants in the same individual changed coordinates between studies.
#' These are all indels, so are due to realignment calling the variant at
#' slightly different locations. We can identify the initial position, by
#' finding the variant that is closest to the validation coordinates.
#'
#' @param row dataframe of a single row for a candidate, with columns for
#'        person_id, chrom, start_pos, ref_allele
#' @param de_novos dataframe of all variants
#'
#' @return boolean vector for sites that match the current row (the current
#'         row should match at minimum).
#' @export
#'
#' @examples
#' sample = read.table(header=TRUE, text="
#'     person_id    chrom start_pos ref_allele alt_allele consequence study_phenotype
#'     sample_01    1     10000     AAA        A          frameshift_variant autism",
#'     stringsAsFactors=FALSE)
#'
#' variants = read.table(header=TRUE, text="
#'     person_id    chrom start_pos ref_allele alt_allele consequence study_phenotype
#'     sample_01    1     10000     AAA        A          frameshift_variant autism
#'     sample_01.p1 1     10001     AAA        A          frameshift_variant normal_iq_autism
#'     sample_02    1     10001     AAA        A          frameshift_variant normal_iq_autism",
#'     stringsAsFactors=FALSE)
#' find_matching_sites(sample, variants)
find_matching_sites <- function(row, de_novos) {
    # identify the candidate sites for the proband on the same chromosome as the
    # given variant, the figure out how far they are from the given variant. We
    # expect the correct position to be within the distance of the ref allele.
    person_matches = de_novos$temp_id == row$temp_id & de_novos$chrom == row$chrom
    rows.delta = abs(as.numeric(as.character(de_novos$start_pos)) - as.numeric(as.character(row$start_pos)))
    distance_matches = rows.delta < nchar(de_novos$ref_allele) | rows.delta < nchar(de_novos$alt_allele)
    
    return(person_matches & distance_matches)
}
