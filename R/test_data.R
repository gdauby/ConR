# Script to read dictionaries and generate the R/sysdata.rda

# loading packages
#library(stringr)
library(readr)
path = "C://Users//renato//Documents//raflima//Pos Doc//Manuscritos//Artigo Extincao na MA//data analysis"
dic_files <- list.files(path = path,
                        pattern = "csv",
                        full.names = TRUE)
dic_files = dic_files[grepl("teste", dic_files)]
encoding <- "ISO-8859-15" # substituir aqui pelo encoding correto
dic <- lapply(dic_files, read_csv, locale = locale(encoding = encoding))
lapply(dic, head)

# transforma em data.frame
dic <- lapply(dic, as.data.frame)

# criterion Ataxonomists:
example_criterionA <- dic[[1]]

# só checando como estao os arquivos
head(example_criterionA)

#Saving the sysdata
save(example_criterionA,
     file = "data/example_criterionA.rda",
     compress = "xz")


