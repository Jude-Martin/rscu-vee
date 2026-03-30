#file:cb_config.R

# Define ORF names
region<-c("orfns","orfst")

#ORF locations (start, stop, number of sections, number of this section) last 2 only if dealing with gaps. 

# VEEV
# orfns <-c(44,7540,1,1)
# VEEV
# orfst <-c(7576,11340,1,1)

#CHIKV

orfns <-c(26,7450,1,1)
orfst <-c(7516,11262,1,1)



region_start <- c(orfns[1],orfst[1])   
region_stop <- c(orfns[2],orfst[2])

#directory with input csv files

#input_dir<-"./alignments/HFtoHF/sorted/pysamstats"
#output_dir<-"./alignments/HFtoHF/sorted/pysamstats/output_files"

input_dir<-"./alignments/WT/sorted/pysamstats"
output_dir<-"./alignments/WT/sorted/pysamstats/output_files"

