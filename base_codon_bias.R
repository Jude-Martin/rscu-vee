## ---------------------------
##
## Script name: base_codon_bias.R 
##
## Purpose of script: Population Codon Bias caluculation
##
## Author: J Martin
##
## 
## Email: 
##
## ---------------------------
##
## Notes:
##   For future  could add a pre-read or selection of ORF  
##
## ---------------------------


library (data.table)
library(BiocGenerics)
library(S4Vectors)
library(XVector)
library(Biostrings)
library(IRanges)
library(readr)



# define variables


#positions of open reading frames, note these correspond to config  tsv files, check if using new input 
#  Using cd_config.R file now instead to define variables
#orfns <-c(44,7540,1,1)
#orfst <-c(7576,11340,1,1)
#orftf <-c(9978,10058,1,1)


#debug  check wd here 

source("./r_script/cb_config.R")

#print(input_dir)
#getwd()

# function for assigning frequencies to nucleotides in ORF
init_freq<-function(count){                                                                                                                                                                                                
  
  for (f in 1:3)
  {
    
    
    assign(paste0('f', f,'a'), as.numeric(rep(0,times=count)),envir = .GlobalEnv)
    assign(paste0('f', f,'c'), as.numeric(rep(0,times=count)),envir = .GlobalEnv)
    assign(paste0('f', f,'g'), as.numeric(rep(0,times=count)),envir = .GlobalEnv)
    assign(paste0('f', f,'t'), as.numeric(rep(0,times=count)),envir = .GlobalEnv)
    
    assign(paste0('p', f,'a'), as.numeric(rep(0,times=count)),envir = .GlobalEnv)
    assign(paste0('p', f,'c'), as.numeric(rep(0,times=count)),envir = .GlobalEnv)
    assign(paste0('p', f,'g'), as.numeric(rep(0,times=count)),envir = .GlobalEnv)
    assign(paste0('p', f,'t'), as.numeric(rep(0,times=count)),envir = .GlobalEnv)
    
    
  }
  
}

cleanup_freq<-function()
{  for (f in 1:3)
{
  rm(list=paste0("f", f,"a"),envir = .GlobalEnv)
  rm(list=paste0("f", f,"c"),envir = .GlobalEnv)
  rm(list=paste0("f", f,"g"),envir = .GlobalEnv)
  rm(list=paste0("f", f,"t"),envir = .GlobalEnv)
  
  
  
  rm(list=paste0("p", f,"a"),envir = .GlobalEnv)
  rm(list=paste0("p", f,"c"),envir = .GlobalEnv)
  rm(list=paste0("p", f,"g"),envir = .GlobalEnv)
  rm(list=paste0("p", f,"t"),envir = .GlobalEnv)
  
  
  
}  
}     



init_results<-function()
{  
  #make a dataframe for results
  df_results_by_codon<<-data.frame(col1,"AA"=col2)
  
  #add columns with 0 values
  df_results_by_codon$syn_count<<-0
  df_results_by_codon$codon_count<<-0
  df_results_by_codon$syn_sum<<-0
  #df_results_by_codon$RSCU<-0
  
  #renaming first column may need to do this if issues with setting up dataframe
  names(df_results_by_codon)[1]<<-"codon"
  
  #use names to make a lookuptable
  #this is a dataframe object from conversion to vector by unlist  then table to tabulate counts as.data.frame converts array to df
  codon_synonym_count <<-as.data.frame(table(unlist(df_results_by_codon$AA)))
  cods<<-codon_synonym_count$Var1
  names(cods)<<-codon_synonym_count$Freq
  for(r_length in 1:length(df_results_by_codon$AA))
    
  {
    
    #populate the synonymous count
    val<<-df_results_by_codon$AA[r_length]
    lookup<<-names(cods)[val]
    
    df_results_by_codon$syn_count[r_length]<<-lookup
  } 
  
} 


#Define genes and features    

#Will need to define a gapped region or region made up of multiple sections




#create dataframe to define regions

df_regions<-data.frame(region,region_start,region_stop)

#length of regions by count of amino acid residues
df_regions$length<-((region_stop-region_start+1)/3)

#put an error-check here to catch if ORF length is not divisible by 3


#defining amino acids in aatable

aa<-c("alanine","arganine","asparagine","aspartate","cysteine","glutamine","glutamate","glycine","histidine","isoleucine","leucine","lysine","methionine","phenylalanine","proline","serine","threonine","tryptophan","tyrosine","valine","stop")

aashort<- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","*")

codons<-c("GCT,GCC,GCA,GCG","CGT,CGC,CGA,CGG,AGA,AGG","AAT,AAC","GAT,GAC","TGT,TGC","CAA,CAG","GAA,GAG","GGT,GGC,GGA,GGG","CAT,CAC","ATT,ATC,ATA","CTT,CTC,CTA,CTG,TTA,TTG","AAA,AAG","ATG","TTT,TTC","CCT,CCC,CCA,CCG","TCA,TCC,TCT,TCG,AGT,AGC","ACA,ACT,ACG,ACC","TGG","TAC,TAT","GTA,GTC,GTT,GTG","TGA,TAG,TAA")

aatable <- data.frame(aa,aashort,codons)

all_codon_st<-paste(codons, collapse=',' )
col1<-strsplit(all_codon_st,",")
all_codon_st<- gsub('[ ,]','',all_codon_st)
coldna<-DNAString(all_codon_st)
col2<-translate(coldna)   

# files to load from a specified directory 
# now in config file   input_dir<-"./test" 
files <- list.files(path=input_dir, pattern="*.csv", full.names=TRUE, recursive=FALSE)

# parse into temp df record subs
# use fread from data.table to make a temp dataframe with required data

#looping through files in directory
#Loop 1

for (filenumber in 1:length(files))
  
{
  # reading columns from file, change this to be more informative
  
  # to check progress uncomment print(filenumber) 
  # There is an issue here ...
  temp <- fread(files[filenumber],select=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22),fill=TRUE,sep = ",")
  
  #Previously was having issues reading in vivo files from VEEV dataset?  Check that  this may be due to zero values? 
  # get name for file being worked on
  
  workingfile<-files[filenumber]
  
  
  #loop here to go through regions
  
  #Do translation of base sequence This is working
  # Base is ref sequence values at each position from file,gsub cuts out commas to make a dna string to translate
  
  seq <- toString(temp$ref)
  seq<- gsub('[ ,]','',seq)
  
  
  
  #Have to make entries sequentially into this table and then save/write it 
  
  
  ############ CALL FUNCTION TO MAKE RESULTS DataFrame
  #### This was causing an issue  not initializing, or not being referenced correctly, now ok
  init_results()
  
  
  
  #translate
  
  #Choose which orf to use
  #loop through orfs,genes here
  #Do something here to define/handle gapped regions or concat all
  
  #Loop 3
  
  
  for (region_gene in 1:length(df_regions$region)) 
  {
    #Can this next line be expanded to sections?
    
    orfcount<-df_regions$length[region_gene]
    orf<-df_regions$region[region_gene]
    thisorf<-as.character(orf)
    position <-(1:orfcount)
    orfstart <-df_regions$region_start[region_gene]
    orfend <- df_regions$region_stop[region_gene]
    
    
    
    
    
    # Is there something here to fix? Can put in check that sequence is divisible by 3 and check number of stop positions
    ##translate the orf
    
    transl <-toString(temp$ref[df_regions$region_start[region_gene]:df_regions$region_stop[region_gene]])
    transl <- gsub('[ ,]','',transl)
    transl <-DNAString(transl)
    translation <-translate(transl)
    
    print (translation)
    
    #consensus translated sequence
    consen <- base::strsplit(as.character(translation),"")
    
    print (consen)
    
    
    #Used to set up dataframe with  frequency and probability of a,c,g,t at each nucleotide position
    init_freq(orfcount)
    
    for (f in 1:3)
    {
      
      assign(paste0('f', f,'a'), as.numeric(rep(0,times=orfcount)))
      assign(paste0('f', f,'c'), as.numeric(rep(0,times=orfcount)))
      assign(paste0('f', f,'g'), as.numeric(rep(0,times=orfcount)))
      assign(paste0('f', f,'t'), as.numeric(rep(0,times=orfcount)))
      
      assign(paste0('p', f,'a'), as.numeric(rep(0,times=orfcount)))
      assign(paste0('p', f,'c'), as.numeric(rep(0,times=orfcount)))
      assign(paste0('p', f,'g'), as.numeric(rep(0,times=orfcount)))
      assign(paste0('p', f,'t'), as.numeric(rep(0,times=orfcount)))
    }
    
    
    #generate dataframe orfdf with numeric positions in orf, this is translated sequence, consensus sequence, frequencies of       nucleotides at position, probabilities of nucleotides at position   

    
    orfdf<-data.frame(position,consen,f1a,f1c,f1g,f1t,p1a,p1c,p1g,p1t,f2a,f2c,f2g,f2t,p2a,p2c,p2g,p2t,f3a,f3c,f3g,f3t,p3a,p3c,p3g,p3t,stringsAsFactors = FALSE)
    names(orfdf)[2]<-"consen"
    

    
     #orfdf<-data.frame(position,consen,f1a,f1c,f1g,f1t,p1a,p1c,p1g,p1t,f2a,f2c,f2g,f2t,p2a,p2c,p2g,p2t,f3a,f3c,f3g,f3t,p3a,p3c,p3g,p3t,stringsAsFactors = FALSE)
     #orfdf<-data.frame(position,consen,f1a,f1c,f1g,f1t,p1a,p1c,p1g,p1t,f2a,f2c,f2g,f2t,p2a,p2c,p2g,p2t,f3a,f3c,f3g,f3t,p3a,p3c,p3g,p3t,stringsAsFactors = FALSE)

#     Error in names(orfdf) <- `*vtmp*` : 
#       'names' attribute [2] must be the same length as the vector [1]
#     Execution halted
     
          
####Error in data.frame(position, consen, f1a, f1c, f1g, f1t, p1a, p1c, p1g,  : 
####                      arguments imply differing number of rows: 2475, 0
####                  Execution halted

    
    
    
    ## generate codon probabilities from matrix of base values This works
    #examine position 1 to 3 calculate
    
    #In pysamstats output 23  headings are:  chrom,pos,ref,reads_all,reads_pp,matches,matches_pp,mismatches,mismatches_pp,deletions,deletions_pp,insertions,insertions_pp,A,A_pp,C,C_pp,T,T_pp,G,G_pp,N,N_pp
    
    for (p in seq(orfstart,orfend,by=3))
    {
      
      p1 <- p
      p2 <- p+1
      p3 <- p+2
      
      #mc 3x4 matrix of counts  
      
      mc <- matrix(nrow=3,ncol=4, dimnames = list  (c("pos1","pos2","pos3"),c("A","C","G","T")))
      #mp 3x4 matrix of probabilities
      
      mp <- matrix(nrow=3,ncol=4,dimnames = list  (c("pos1","pos2","pos3"),c("p(A)","p(C)","p(G)","p(T)")))
      
      for (c_pos in p1:p3)
      {
        total <- temp$reads_all[c_pos]
        #from reference get the consensus/expected at position
        cbase <- temp$ref [c_pos]
        
        # debug section
        # print (pos1total)
        # print (cbase)
        
        acount <-temp$A[c_pos]
        ccount <-temp$C[c_pos]
        gcount <-temp$G[c_pos]
        tcount <-temp$T[c_pos]
        
        
        #if (cbase == 'A'){acount=temp$PM[c_pos]  + acount}
        #if (cbase == 'C'){ccount=temp$PM[c_pos]  + ccount}
        #if (cbase == 'G'){gcount=temp$PM[c_pos] + gcount}
        #if (cbase == 'T'){tcount=temp$PM[c_pos]  + tcount}
        
        #vector and matrix entry of counts
        vc <- c(acount,ccount,gcount,tcount)
        
        #assign counts to matrix
        mc [(c_pos%%p1)+1,] <-vc
        
        pa <- as.double(acount/total)
        pa[is.na(pa)]<-0
        
        pc <- as.double(ccount/total)
        pc[is.na(pc)]<-0
        
        pg <- as.double(gcount/total)
        pg[is.na(pg)]<-0
        
        pt <- as.double(tcount/total)
        pt[is.na(pt)]<-0
        
        # vector and matrix entry of probabilities
        
        vp <-c(pa,pc,pg,pt)
        
        # assign probabilities to matrix  
        mp [(c_pos%%p1)+1,] <-vp
        
        
      }
      
      # orfdf[dfnum,"mp"] <-mp
      
      dfnum <-(((p-orfstart)/3)+1)
      
      orfdf[dfnum,"p1a"] <- mp[1,1]
      orfdf[dfnum,"p1c"] <- mp[1,2]
      orfdf[dfnum,"p1g"] <- mp[1,3]
      orfdf[dfnum,"p1t"] <- mp[1,4]
      
      orfdf[dfnum,"f1a"] <- mc[1,1]
      orfdf[dfnum,"f1c"] <- mc[1,2]
      orfdf[dfnum,"f1g"] <- mc[1,3]
      orfdf[dfnum,"f1t"] <- mc[1,4]
      
      orfdf[dfnum,"p2a"] <- mp[2,1]
      orfdf[dfnum,"p2c"] <- mp[2,2]
      orfdf[dfnum,"p2g"] <- mp[2,3]
      orfdf[dfnum,"p2t"] <- mp[2,4]
      
      orfdf[dfnum,"f2a"] <- mc[2,1]
      orfdf[dfnum,"f2c"] <- mc[2,2]
      orfdf[dfnum,"f2g"] <- mc[2,3]
      orfdf[dfnum,"f2t"] <- mc[2,4]
      
      orfdf[dfnum,"p3a"] <- mp[3,1]
      orfdf[dfnum,"p3c"] <- mp[3,2]
      orfdf[dfnum,"p3g"] <- mp[3,3]
      orfdf[dfnum,"p3t"] <- mp[3,4]
      
      orfdf[dfnum,"f3a"] <- mc[3,1]
      orfdf[dfnum,"f3c"] <- mc[3,2]
      orfdf[dfnum,"f3g"] <- mc[3,3]
      orfdf[dfnum,"f3t"] <- mc[3,4]
      
      
      
      
    }
    
    # checked orfdf files, seem to be generating them ok, but may have issue with summation of values and rounding errors changed   sum/count section to integer values, represents ok in tables
    
    ##This works now
    #May need to swap in "Stop" for asterix otherwise things go screwy in orfdf
    
    #orfdf[orfdf == "*"] <- "Stop"
    
    #Can write orfdf to file for error checking here if needed
    
    syn_codons<-as.character(rep("NA",times=orfcount),stringsAsFactors = FALSE)
    
    #output handling starts here
    
    
    df_output<-data.frame(orfdf[,2],syn_codons,stringsAsFactors = FALSE)
    #Next part loops through orfdf
    #cthing will be codon
    
    for (row_numbers in 1:nrow(orfdf))
    {
      cthing<-(orfdf[row_numbers,2])
      
      syn_list<- as.character(aatable$codons[aatable$aashort==cthing])
      
      df_output[row_numbers,2]<-syn_list
    }
    
    
    
    #This works
    
    #Use probabilities as weights in RSCU calculation
    #RSCU=(number of codons coding aa)*(count of particular codon)/sum(total count of all synonymous codons)
    #Neep to Display as some kind of heatmap or bar diagram/histogram
    
    #Adding columns to output df
    #C columns for codon V columns for count of each codon
    
    df_output$C1<-"NA"
    df_output$C2<-"NA"
    df_output$C3<-"NA"
    df_output$C4<-"NA"
    df_output$C5<-"NA"
    df_output$C6<-"NA"
    
    
    df_output$V1<-as.double(0)
    df_output$V2<-as.double(0)
    df_output$V3<-as.double(0)
    df_output$V4<-as.double(0)
    df_output$V5<-as.double(0)
    df_output$V6<-as.double(0)
    
    # identify probability of each variant codon
    #Issue must be occuring somewhere in this section? corrected
    
    #working through df_output row by row
    for (output_row in 1:nrow(df_output))
    {
      #split codon list into individual codons
      
      tocut<-df_output$syn_codons[output_row]
      
      #making a list of synonymous codon variants
      
      ncodons<-as.list(strsplit(tocut,",")[[1]])
      
      #for each synonymous codon at position check orfdf for probabilities at each position
      
      for (num_codons in 1:length(ncodons))
      {
        
        toadd<-(ncodons[[num_codons]])
        
        df_output[output_row,paste0("C",num_codons)]<-toadd
        
        #split each codon into bases
        
        bases_resolve<-as.list(strsplit(toadd,"")[[1]])
        p1lc<-tolower(bases_resolve[[1]])
        p2lc<-tolower(bases_resolve[[2]])
        p3lc<-tolower(bases_resolve[[3]])
        
        prob1<-orfdf[output_row,paste0("p1",p1lc)]
        prob2<-orfdf[output_row,paste0("p2",p2lc)]
        prob3<-orfdf[output_row,paste0("p3",p3lc)] 
        pc<-as.double(prob1*prob2*prob3)
        df_output[output_row,paste0("V",num_codons)]<-pc
      }
      
    }
    
    dfoutfile<- gsub(".tsv", "dfout",workingfile)
    dfoutfile<-gsub ("input", "output_files",dfoutfile)
    dfoutfile<-gsub(".csv","",dfoutfile)
    dfoutfile<- paste0(dfoutfile,orf,".tsv")
    
    write_tsv(df_output, dfoutfile)
    #Up to about here to create df_output dataframe. Can we split in two here.
    
    
    
    for(r_length2 in 1:length(df_results_by_codon$AA))  
      
    {  
      
      
      #designate a codon to query df_results dataframe with
      codon_q<-(df_results_by_codon$codon[r_length2])
      sum_codon_values<-as.double(0)
      
      #loop through the reading frame
      
      for (j in 1:length(df_output$C1))
      {
        
        # loop through the columns   
        for (sum_rows in 1:6)
        {
          #    print(k)
          q_k<-paste0("df_output$C",sum_rows,"[j]")
          
          
          
          if (eval(parse(text=q_k))==codon_q) 
            
            # add to sum values
          {
            q_v<-paste0("sum_codon_values+df_output$V",sum_rows,"[j]")
            sum_codon_values<-(eval(parse(text=q_v)))
          }
          
          
          else
          {sum_codon_values<-sum_codon_values+0} 
          #TRYING as.integer
          
          # total of the values to results table
          
        }
        
        
      }
      #print(sum_codon_values)
      #here to by orf
      
      
      
      df_results_by_codon$codon_count[r_length2]<-sum_codon_values
    }
    
    
    
    
    #  Need to make function/loop  to use  for total of all syn codons
    
    
    
    
    
    
    #Fill in synonymous codon sum
    for (z in 1:length(df_results_by_codon$codon) )
    {
      index <- df_results_by_codon$AA[z]  
      sum<-sum(df_results_by_codon[which(df_results_by_codon$AA==index),"codon_count"])
      # (orfdf[, 2], syn_codons, stringsAsFactors = FALSE) 
      df_results_by_codon$syn_sum[z]<-sum
      
      #Then make a new column with the resulting RSCU
      
      RSCU_vector<-df_results_by_codon$codon_count/df_results_by_codon$syn_sum
      
      #df_results_by_codon$syn_count<-as.double(as.character()df_results_by_codon$syn_count)
      
      df_results_by_codon$syn_count<-as.integer(df_results_by_codon$syn_count)
      
      #TRYING as.integer as this should be a count of discrete number of reads
      
      RSCU_vector<-RSCU_vector*df_results_by_codon$syn_count
      
      #  this thing paste0(thisorf,"RSCU")
      
      current_RSCU<-paste0(thisorf,"RSCU")
      
      eval(parse(text = paste0("df_results_by_codon$",current_RSCU,"<-RSCU_vector")))
      
      
      
      
      
      #could write to results by codon dataframe with with columns for each region?
      #Have to make results by codon outside region loop
    } 
  }
  
  #Save files
  savefile<- gsub(".csv", "",workingfile)
  savefile<-gsub("input","output_files",savefile)
  write.csv(df_results_by_codon, paste0(savefile,"_RSCU.csv"), row.names=FALSE, quote=FALSE)     
  
  
}
