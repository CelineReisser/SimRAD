#################################################
## SimRAD code for simulation genome digestion ##
#################################################
# Celine MO Reisser, May 2017

#source("https://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
#biocLite("ShortRead")
#install.packages("SimRAD",dependencies=T)



library(SimRAD)


#simulate the genome
simseq <- sim.DNAseq(size=100000000, GCfreq=0.40)


#######################################################
## Single digest followed by shearing: classical RADseq

#PstI#
cs_5p1 <- "G"
cs_3p1 <- "AATTC"

#digestion of the "simseq" genome:
simseq.dig <- insilico.digest(simseq, cs_5p1, cs_3p1, verbose=TRUE)


# doing the same code, but doing 10 iteration to obtain a min max and average
x<-c()
for (i in 1:10) {
  simseq <- sim.DNAseq(size=100000000, GCfreq=0.40)
  simseq.dig <- insilico.digest(simseq, cs_5p1, cs_3p1, verbose=F)
  x[i]=length(simseq.dig)-1
}

min<-paste("The minimum number of loci is",min(x),"for a genome of 100Mb and 40% GC.",sep=" ")
max<-paste("The maximum number of loci is",max(x),"for a genome of 100Mb and 40% GC.",sep=" ")
ave<-paste("The average number of loci is", mean(x),"for a genome of 100Mb and 40% GC.",sep=" ")
min
max
ave


########################################################################
## ddRAD double restriction enzyme digestion followed by size selection:

#Define the restriction enzyme 1 recognition pattern:
#PstI#
cs_5p1 <- "G"
cs_3p1 <- "AATTC"

#Define the restriction enzyme 2 recognition pattern:
#MspI :  C'CGG
cs_5p2 <- "C"
cs_3p2 <- "CGG"

simseq <- sim.DNAseq(size=100000000, GCfreq=0.40)

#Simulation of the digestion just like before, but by adding the new recognition site
simseq.dig <- insilico.digest(simseq, cs_5p1, cs_3p1, cs_5p2, cs_3p2, verbose=F)

#selecting fragment with ends corresponding to each enzyme: E1--E2 and E2--E1 (versus E1--E1 and E2--E2)
simseq.sel <- adapt.select(simseq.dig, type="AB+BA", cs_5p1, cs_3p1, cs_5p2, cs_3p2)

#selecting the fraction of fragments between 450 and 530bp:
nar.simseq <- size.select(simseq.sel,  min.size = 450, max.size = 530, graph=T, verbose=T)


#Now do that again, also doing 10 iterations to get min max and average:
x<-c()
for (i in 1:10) {
  simseq <- sim.DNAseq(size=100000000, GCfreq=0.40)
  simseq.dig <- insilico.digest(simseq, cs_5p1, cs_3p1, cs_5p2, cs_3p2, verbose=F)
  simseq.sel <- adapt.select(simseq.dig, type="AB+BA", cs_5p1, cs_3p1, cs_5p2, cs_3p2)
  nar.simseq <- size.select(simseq.sel,  min.size = 450, max.size = 530, graph=F, verbose=F)
  x[i]=length(nar.simseq)
}

min<-paste("The minimum number of loci is", min(x),"for a genome of 100Mb and 40 GC.",sep=" ")
max<-paste("The maximum number of loci is", max(x),"for a genome of 100Mb and 40% GC.",sep=" ")
ave<-paste("The average number of loci is", mean(x),"for a genome of 100Mb and 40% GC.",sep=" ")
min
max
ave



