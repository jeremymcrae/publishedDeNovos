% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ensembl_rest.R
\name{get_sequence_in_region}
\alias{get_sequence_in_region}
\title{find genomic sequence within a region}
\usage{
get_sequence_in_region(variant, build = "grch37", verbose = FALSE)
}
\arguments{
\item{variant}{data frame or list for a variant, containing columns named
"chrom", "start_pos", and "end_pos" for a single variant}

\item{build}{genome build to find consequences on}

\item{verbose}{flag indicating whether to print variants as they are checked}
}
\value{
a character string containing the HGNC symbol.
}
\description{
find genomic sequence within a region
}
\examples{
get_gene_id_for_variant(data.frame(chrom=c("1"), start_pos=c("1000000"),
    end_pos=c("1000000")))
get_gene_id_for_variant(list(chrom="1", start_pos="1000000",
    end_pos="1000000"))
}

