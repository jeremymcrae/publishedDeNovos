
library(publishedDeNovos)

# load all the de novo datasets from the published articles
deligt_data = deligt_de_novos()
deligt = deligt_data$variants
deligt_sex = deligt_data$gender
gilissen = gilissen_de_novos(deligt, deligt_sex)
rauch = rauch_de_novos()
epi4k = epi4k_de_novos()
autism = autism_de_novos()
fromer = fromer_de_novos()
zaidi = zaidi_de_novos()

# join the de novo datasets togther, and standardise the variant type column
variants = rbind(rauch, deligt, gilissen, epi4k, autism, fromer, zaidi)
variants$type = "snv"
variants$type[nchar(variants$ref_allele) != 1 | nchar(variants$alt_allele) != 1] = "indel"

# make sure all the columns are character
variants[] = lapply(variants, as.character)

save(variants, file="data/variants.rda", compress="xz")

# also save a copy, so that non-R code can also access the data
gz = gzfile("data-raw/variants.txt.gz", "w")
write.table(variants, gz, sep="\t", quote=FALSE, row.names=FALSE)
