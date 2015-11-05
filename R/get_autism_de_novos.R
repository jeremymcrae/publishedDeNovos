
#' get de novo data from autism exome studies
#'
#' I think this data was obtained from studies published using the Simon's
#' Simplex Collection data (http://sfari.org/resources/simons-simplex-collection)
#'
#' Supplementary table 2 (where the excel sheets for the probands and
#' siblings have been combined) from:
#' Sanders et al. (2012) Nature 485:237-241
#' doi: 10.1038/nature10945
#'
#' Supplementary table 3 from:
#' O'Roak et al. (2012) Nature 485:246-250
#' doi: 10.1038/nature10989
#'
#' Supplementary table 1 (where the non-coding SNVs have been excluded) and
#' supplementary table 2 from:
#' Iossifov et al. (2012) Neuron 74:285-299
#' doi: 10.1016/j.neuron.2012.04.009
#'
#' @export
#'
#' @return data frame of de novos, including gene symbol, functional consequence
#'     (VEP format), chromosome, nucleotide position and SNV or INDEL type
autism_de_novos <- function() {
    
    # load the publication specific de novo datasets
    sanders = sanders_de_novos()
    oroak = oroak_de_novos()
    iossifov_neuron = iossifov_neuron_de_novos()
    iossifov_nature = iossifov_nature_de_novos()
    derubeis_nature = de_rubeis_de_novos()
    
    # exclude de novos identified in previous studies, but make sure the probands
    # from the Iossifox Nature 2014 paper that have been classified as
    # intellectual_disability, transfer their classifications to the earlier de
    # novo entries.
    key_1 = paste(iossifov_neuron$person_id, iossifov_neuron$start_pos)
    key_2 = paste(iossifov_nature$person_id, iossifov_nature$start_pos)
    intersection = key_1[key_1 %in% key_2]
    iossifov_neuron[["study_phenotype"]][match(intersection, key_1)] = iossifov_nature[["study_phenotype"]][match(intersection, key_2)]
    iossifov_nature = iossifov_nature[!(key_2 %in% key_1), ]
    
    autism_de_novos = rbind(iossifov_nature, sanders, oroak, iossifov_neuron,
        derubeis_nature)
    
    # remove de novos that have been dupicated between studies. These are easily
    # spotted as individuals who have IDs that are nearly identical between
    # different studies eg 14323.p1 vs 14323
    temp = autism_de_novos
    temp$person_id = sapply(strsplit(temp$person_id, "\\."), "[", 1)
    to_remove =  temp[duplicated(temp[, c("person_id", "chrom", "start_pos")]), ]
    autism_de_novos = autism_de_novos[!row.names(autism_de_novos) %in% row.names(to_remove), ]
    
    return(autism_de_novos)
}
