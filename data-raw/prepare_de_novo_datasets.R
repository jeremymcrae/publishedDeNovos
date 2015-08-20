
library(publishedDeNovos)

# load all the de novo datasets fro the published articles
rauch = rauch_de_novos()
deligt_data = deligt_de_novos()
deligt = deligt_data$variants
deligt_sex = deligt_data$gender
gilissen = gilissen_de_novos(deligt, deligt_sex)
epi4k = epi4k_de_novos()
autism = autism_de_novos()
fromer = fromer_de_novos()
zaidi = zaidi_de_novos()

# join the de novo datasets togther, and standardise the variant type column
published_de_novos = rbind(rauch, deligt, gilissen, epi4k, autism, fromer, zaidi)
published_de_novos$type = "snv"
published_de_novos$type[nchar(published_de_novos$ref_allele) != 1 | nchar(published_de_novos$alt_allele) != 1] = "indel"

# make sure all the columns are character
published_de_novos[] = lapply(published_de_novos, as.character)

save(published_de_novos, file="data/test_published_de_novos.rda", compress="xz")
