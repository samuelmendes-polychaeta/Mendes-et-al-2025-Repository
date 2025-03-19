#' Author: Samuel Lucas da S. Delgado Mendes, Paulo Cesar de Paiva,
#' Rodolfo L. Nascimento
#' Subject: Functional dispersion framework reveals that
#' depth and upwelling occurrence drive the assembly of 
#' soft-bottom polychaete communities
#' Journal: Marine Ecology Progress Series

# Packages ---------------------------------------------------------------------

{
lib <- .libPaths()[1] 
required.packages <- c("vegan", "picante",
                       "writexl", "openxlsx",
                       "readxl",  "here","tidyverse", 
                       "ade4", "hclust") 
i1 <- !(required.packages %in% row.names(installed.packages())) 
if(any(i1)) { 
  install.packages(required.packages[i1], dependencies = TRUE, lib = lib) 
} 
lapply(required.packages, require, character.only = TRUE)

} # required packages


# Import -----------------------------------------------------------------------

{
  traits<-read.xlsx(here("Data", "data.xlsx"), sheet = "R_trait_df") #Raw fuzzy coded traits
  
  abund<-read.xlsx(here("Data", "data.xlsx"), sheet = "abund") #Polychaete genera incidence per community
  
  abund[is.na(abund)] <- 0 #empty space in L matrix must be zero
  
  env<-read.xlsx(here("Data", "data.xlsx"), sheet = "env_variables")#enviromental variables per community
}


# Fuzzy scores calculation -----------------------------------------------------

{
  #Fuzzy scores calculation
  
  Genus = traits[1]
  
  jaw<-prep.fuzzy(traits[2:3], col.blocks = 2) #1
  branchiae<-prep.fuzzy(traits[4:7], col.blocks = 4) #2
  parapodial_rami<-prep.fuzzy(traits[8:10], col.blocks = 3) #3
  parapodial_apdg<-prep.fuzzy(traits[11:15], col.blocks = 5) #4
  parapodial_lb<-prep.fuzzy(traits[16:17], col.blocks = 2) #5
  chaetae<-prep.fuzzy(traits[18:29], col.blocks = 12) #6
  pharinx_morphology<-prep.fuzzy(traits[30:33], col.blocks = 4)#7
  palps<-prep.fuzzy(traits[34:37], col.blocks = 4) #8
  head_ap<-prep.fuzzy(traits[38:40], col.blocks = 3) #9
  sensory_organs<-prep.fuzzy(traits[41:49], col.blocks = 9) #10
  body_ap<-prep.fuzzy(traits[50:51], col.blocks = 2) #11
  body_surface<-prep.fuzzy(traits[52:54], col.blocks = 3) #12
  
  larval_dev<-prep.fuzzy(traits[55:57], col.blocks = 3) #13
  fate_ova<-prep.fuzzy(traits[58:63], col.blocks = 6) #14
  asex_rep<-prep.fuzzy(traits[64:65], col.blocks = 2) #15
  rep_mod<-prep.fuzzy(traits[66:68], col.blocks = 3) #16
  
  habitat<-prep.fuzzy(traits[69:71], col.blocks = 3) #17
  feeding<-prep.fuzzy(traits[72:76], col.blocks = 5) #18
  motility<-prep.fuzzy(traits[77:80], col.blocks = 4) #19
  
  size<-prep.fuzzy(traits[81:83], col.blocks = 3) #20
  seg_number<-prep.fuzzy(traits[84:85], col.blocks = 2) #21
  
  #db.trait dataframe 
  
  db.trait = cbind.data.frame(Genus, jaw, branchiae, parapodial_rami, 
                              parapodial_apdg, parapodial_lb, chaetae, 
                              pharinx_morphology, palps, head_ap, 
                              sensory_organs, body_ap, body_surface, 
                              larval_dev, fate_ova, asex_rep, rep_mod, 
                              habitat, feeding, motility, size, seg_number)
  
  
  rownames(db.trait) <- db.trait$Genus
  db.trait<-db.trait[,-1]
  str(db.trait)
  
} #preparing db.trait object with polychaete traits fuzzy scores



# Functional trait dispersion --------------------------------------------------

{
  #listing traits
  ktab_list<- ktab.list.df(list(jaw, branchiae, parapodial_rami, 
                                parapodial_apdg, parapodial_lb, chaetae, 
                                pharinx_morphology, palps, head_ap, 
                                sensory_organs, body_ap, body_surface, 
                                larval_dev, fate_ova, asex_rep, rep_mod,
                                habitat, feeding, motility, size, seg_number)) 
  #naming traits
  ktab_list$blo <- c(ncol(jaw), ncol(branchiae), ncol(parapodial_rami), 
                     ncol(parapodial_apdg), ncol(parapodial_lb), ncol(chaetae), 
                     ncol(pharinx_morphology), ncol(palps), ncol(head_ap), 
                     ncol(sensory_organs), ncol(body_ap), 
                     ncol(body_surface), ncol(larval_dev), ncol(fate_ova), 
                     ncol(asex_rep), ncol(rep_mod), 
                     ncol(habitat), ncol(feeding), ncol(motility), 
                     ncol(size), ncol(seg_number)) 
  
  #specifying that the traits are fuzzy coded
  dist_matrix <- dist.ktab(ktab_list, type = c("F", "F", "F", "F", "F", 
                                               "F", "F", "F", "F", "F", 
                                               "F", "F", "F", "F", "F", 
                                               "F", "F", "F", "F", "F", "F")) 
  
  dist_obj <- as.dist(as.matrix(dist_matrix)) #creating a distance object
  attr(dist_obj, "Labels") <- rownames(db.trait) #attributing taxa names 
  
  #Check for coherence among polychaete taxa functional
  #traits gower distance matrix
  plot(hclust(dist_obj)) 
  
} #1. Gower distance matrix calculation

{
  #SES MPD calculation, requiring the incidence matrix and the gower distance matrix
  ses_mpd<- ses.mpd(as.matrix(abund[,-c(1)]), as.matrix(dist_obj), 
                    runs = 999, null.model = "taxa.labels", abundance.weighted = F)
  
  #SES MNTD calculation, requiring the incidence matrix and the gower distance matrix
  ses_mntd<- ses.mntd(as.matrix(abund[,-c(1)]), as.matrix(dist_obj), 
                      runs = 999, null.model = "taxa.labels", abundance.weighted = F) 
  
  #Grouping the results in a dataframe
  FD_indexes <- data.frame(MPD = ses_mpd$mpd.obs.z, MNTD = ses_mntd$mntd.obs.z) 
  
  
} #2. SES MPD and MNTD  

# RLQ and fourthcorner analysis ------------------------------------------------

{
  #Data preparation
  {
    {
      R = data.frame(env[,-c(1:4,16:17)])
      rownames(R) <- env$comm 
      str(R) 
      
      L = data.frame(abund[-1]) 
      rownames(L) <- abund$comm
      str(L)
      
      # Combining all fuzzy objects into a single data frame
      Q <- data.frame(db.trait)
      
      # Defining col.blocks (number of trait modalities, i.e. columns, for each fuzzy trait)
      
      col.blocks<- c(1,1,
                     2,2,2,2,
                     3,3,3,
                     4,4,4,4,4,
                     5,5, 
                     6,6,6,6,6,6,6,6,6,6,6,6,
                     7,7,7,7,
                     8,8,8,8,
                     9,9,9,
                     10,10,10,10,10,10,10,10,10,
                     11,11,
                     12,12,12,
                     13,13,13,
                     14,14,14,14,14,14,
                     15,15,
                     16,16,16,
                     17,17,17,
                     18,18,18,18,18,
                     19,19,19,19,
                     20,20,20,
                     21,21)
      
      str(col.blocks)
      
      attr(Q, "col.blocks") <- col.blocks
      
      row.w <- rep(1, nrow(Q)) #equal weights for fuzzy coded traits' FCA
      
      attr(Q, "row.w") <- row.w
      
    }#for RLQ and fourthcorner (preparation of R, L and Q objects)
  }
  
  #1. ordinations
  {
    {
      L_dudi <- dudi.coa(L, scannf = FALSE, nf = 2) #COA for L matrix
      R_dudi <- dudi.pca(R, row.w = L_dudi$lw, scale = T, scannf = FALSE, nf = 2) #PCA for R matrix
      Q_dudi <- dudi.fca(Q, scannf = FALSE, nf = 2) #FCA for Q matrix
      
      #Correcting line weights after the ordinations for the RLQ analysis
      R_dudi$lw <- L_dudi$lw  
      Q_dudi$lw <- L_dudi$cw  
    }
  }
  
  #2. RLQ analysis and significance
  {
    RLQ <- rlq(dudiR = R_dudi, dudiL = L_dudi, dudiQ = Q_dudi, scannf = FALSE, nf = 2) #RLQ analysis
    plot(RLQ) #plotting the RLQ 
    summary(RLQ) 
    
    rand_test_abiotic = randtest(RLQ, nrepet = 49999, modeltype = 6) #Significance test for the multivariate pattern
    
    adjusted_p_value <- p.adjust(rand_test_abiotic$pvalue, method = "fdr") #Adjusting p values according to False Discovery Rate (FDR) method
    
    p_value_SRLQ = fourthcorner2(R, L, Q,
                                 modeltype = 6, p.adjust.method.G = "fdr", nrepet = 49999)#Fourthcorner statistic 
    
  }
  
  #3. Fourthcorner analysis
  {
    fc_adj = fourthcorner(tabR = R, tabL = L, tabQ = Q, 
                          modeltype = 6, p.adjust.method.G = "fdr", 
                          p.adjust.method.D = "fdr", nrepet = 49999) #Fourthcorner analysis (trait x enviroment correlations and their p values FDR corrected)
    
    plot(fc_adj)
    plot(fc_adj, x.rlq = RLQ, stat = "D2", type = "biplot")
    
    summary(fc_adj)
    
    fc_Q = fourthcorner.rlq(RLQ, modeltype = 6,
                            typetest = "Q.axes", nrepet = 49999, p.adjust.method.G = "fdr",
                            p.adjust.method.D = "fdr") #Significance of each trait correlation with RLQ axes (FDR corrected)
    
    
    fc_R = fourthcorner.rlq(RLQ, modeltype = 6,
                            typetest = "R.axes", nrepet = 49999, p.adjust.method.G = "fdr",
                            p.adjust.method.D = "fdr") #Significance of each enviromental variable correlation with the RLQ axes (FDR corrected)
  }
}



