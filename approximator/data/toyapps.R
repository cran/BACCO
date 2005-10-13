require(emulator, quietly=TRUE, save=FALSE)
require(approximator, quietly=TRUE, save=FALSE)

"D1.toy" <-
  matrix(
         c(0.525, 0.075, 0.275, 0.675, 0.625, 0.775, 0.025, 
           0.975, 0.825, 0.475, 0.425, 0.225, 0.725, 0.575, 0.125, 0.875, 
           0.325, 0.375, 0.175, 0.925, 0.375, 0.175, 0.325, 0.825, 0.575, 
           0.025, 0.475, 0.125, 0.675, 0.975, 0.075, 0.875, 0.925, 0.775, 
           0.725, 0.225, 0.525, 0.275, 0.425, 0.625, 0.375, 0.775, 0.925, 
           0.125, 0.625, 0.475, 0.725, 0.975, 0.025, 0.075, 0.825, 0.225, 
           0.275, 0.325, 0.575, 0.675, 0.425, 0.525, 0.875, 0.175),
         ncol=3
         )
colnames(D1.toy) <- letters[1:3]
rownames(D1.toy) <- paste("point",1:nrow(D1.toy),sep=".")

level1 <- 1:20
level2 <- c(3, 4, 5, 8, 13, 15, 18, 19)
level3 <- c(3,       8, 13, 15, 18, 19)
level4 <- c(3,       8,     15, 18, 19)

"subsets.toy" <-
  list(
       level1 = level1,
       level2 = level2,
       level3 = level3,
       level4 = level4
       )
names(subsets.toy[[1]]) <- paste("level.1_run",1:length(level1),sep=".")
names(subsets.toy[[2]]) <- paste("level.2_run",1:length(level2),sep=".")
names(subsets.toy[[3]]) <- paste("level.3_run",1:length(level3),sep=".")
names(subsets.toy[[4]]) <- paste("level.4_run",1:length(level4),sep=".")

rm(level1)
rm(level2)
rm(level3)
rm(level4)




jj  <- generate.toy.observations(export.truth=TRUE)
hpa.toy <- jj$hpa
betas.toy <- jj$betas
rm(jj)


set.seed(0)
"z.toy" <- 
generate.toy.observations(D1=D1.toy, subsets=subsets.toy, basis.fun=basis.toy)
