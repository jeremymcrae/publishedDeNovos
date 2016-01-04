
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
    
    autism_de_novos = rbind(iossifov_nature, sanders, oroak, iossifov_neuron,
        derubeis_nature)
    
    autism_de_novos = remove_duplicates(autism_de_novos)
    
    return(autism_de_novos)
}
