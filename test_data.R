# Script to read dictionaries and generate the R/sysdata.rda

# loading packages
library(stringr)
library(readr)

#### Reading and pre-editing files ####
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

#### criterion A ####
example_criterionA <- dic[[1]]

# só checando como estao os arquivos
head(example_criterionA)

#Random data for different patterns of decline
set.seed(50)
ys = seq(0, 31, 1)
##linear
#lm(c(10000, 8400, 6700) ~ c(0,8,16), data=x)
error = runif(length(ys), -150, 150)
pop.lin = round(10000 - (200 * ys) + error, 0)
#plot(ys, pop.lin, xlim=c(0,18), ylim = c(4000,10500))
#curve(10000 - 200*x, lwd=2, col=1, add=TRUE)
##exponential
#points(c(10100, 8500, 7500, 7000, 6800) ~ c(0,4,8,12,16), pch=19)
#nls(c(10100, 8500, 7500, 7000, 6800) ~ a*exp(b*c(0,4,8,12,16)), start=list(a=1000,b=-0.1))
error = runif(length(ys), -150, 150)
pop.exp = round(9750 * exp(-0.0265*ys) + error, 0)
#points(pop.exp ~ ys, col = 2)
#curve(9750 *exp(-0.0265*x), lwd=2, col=2, add=TRUE)
##accelerating
#points(c(10100, 9700, 9200, 8300, 6900) ~ c(0,4,8,12,16), pch=19)
#nls(c(10100, 9700, 9200, 8300, 6900) ~ a/(b+c(0,4,8,12,16)^c), start=list(a=1000000,b=1000,c=2.5))
error = runif(length(ys), -150, 150)
pop.acc = round(1.852e+07/ (1.854e+03 + ys^2.42) + error, 0)
#points(pop.acc ~ ys, col = 4)
#curve(1.852e+07/ (1.854e+03 + x^2.42), add=TRUE, lwd=2, col=4)
##piece-wise
#points(c(10000, 9950, 9850, 9850, 9100, 8100, 7200, 7100, 7000, 6950) ~ c(0,2,4,6,8,10,12,14,16,18), pch=19, col=2)
anos = c(0,2,4,6,8,10,12,14,16,18)
md <- stats::lm(c(10000, 9950, 9850, 9850, 9100, 8100, 7200, 7100, 7000, 6950) ~ 
                  anos) 
piece2 <- segmented::segmented(md, seg.Z = ~anos, psi = c(5, 7), 
                               control = segmented::seg.control(display = FALSE), it.max = 100, n.boot = 50)
error = runif(length(ys), -50, 40)
pop.piece = round(predict(piece2, data.frame(anos = ys)) + error, 0)
#plot(piece2, add=TRUE, lwd= 2, col=5)
#points(pop.piece ~ ys, col = 5, pch=19)

new.pops = rbind.data.frame(c("species 3", pop.lin),
                            c("species 4", pop.exp),
                            c("species 5", pop.acc),
                            c("species 6", pop.piece),
                            deparse.level = 0, 
                            stringsAsFactors = FALSE)
names(new.pops) = names(example_criterionA)
new.pops[2] = round((as.numeric(new.pops[,2]) + 20000)/3,0)
#new.pops[17] = round((as.numeric(new.pops[,17]) + 6700)/2,0)
example_criterionA = rbind.data.frame(example_criterionA,
                                      new.pops,
                                      deparse.level = 0, 
                                      stringsAsFactors = FALSE)

## Saving the sysdata ##
save(example_criterionA,
     file = "data/example_criterionA.rda",
     compress = "xz")

load("data/example_criterionA.rda")

#### criterion C ####
example_criterionC_subpops <- dic[[2]][,1:8]

# Overall pop sizes
example_criterionC <- rowsum(example_criterionC_subpops[,c(2:8)],
       example_criterionC$species)

# checking
example_criterionC_subpops
example_criterionC

save(example_criterionC_subpops,
     file = "data/example_criterionC_subpops.rda",
     compress = "xz")

save(example_criterionC,
     file = "data/example_criterionC.rda",
     compress = "xz")

#### FLUCTUATIONS ####
example_fluctuation <- dic[[3]]

example_fluctuation <- example_fluctuation[!apply(example_fluctuation[,-1], 1, function(x) all(is.na(x))), ]

save(example_fluctuation,
     file = "data/example_fluctuation.rda",
     compress = "xz")



