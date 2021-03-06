% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fix_alleles.R
\name{fix_ins_alleles}
\alias{fix_ins_alleles}
\title{fix insertion allele codes}
\usage{
fix_ins_alleles(variants, allele_column)
}
\arguments{
\item{variants}{data frame of variants}

\item{allele_column}{name of column containing allele information}
}
\value{
a data frame with chrom, start_pos, end_pos and allele columns.
}
\description{
fix the allele column for variants with allele columns that are structured
like "ins(1)"
}

