# runs a simulation of a word frequency sequence effect experiment (no feedback) at some level of L, F, SF and dp,
# returning five data frames: 
# 1) the main data frame containing the results in terms of the mean ∆P("old") for targets and lures across trial quarters 
#       (and SEs); 
# 2) the quartile means (and SEs) for all four test item types (LF lures, etc.);
# 3) the four overall means and SEs; 
# 4) sample criteria across trials four two sample "subjects, one starting low (0.25 below expected
#       mean of the intensity distribution) and one starting high (0.25 above expected intensity mean); 
# 5) the distribution of intensities  across test items.

sequence<- function(SF1, dp, L, A, N_subjects){
  library(dplyr) 
  HF = 5; LF = 1
  Frequency<- 0    
  N_el = 100; N_items = 80
  N_seqSegment = 4                       # cutting quartet means into quarters (N_seqSegment = 4), deciles (N_seqSegment =10)?
  N_trialPres = (N_items/2)/N_seqSegment
  N_study = N_items*2
  N_types = 4; N_test = N_items*N_types      
  indices<- seq(1:N_items)
  LFind<- indices[ indices %% 4 == 1| indices %% 4 == 2]
  HFind<- indices[ indices %% 4 == 3| indices %% 4 == 0]  
  quartetOrder<- vector(, N_test/N_seqSegment); testOrder<- vector(, N_test)
  
  testItems<- matrix(0, nrow=N_test, ncol=N_el)
  studyItems<- matrix(0, nrow=N_study, ncol=N_el) 
  
  N_bins = 21
  Bin_size = 0.075
  Bin_start = -.525
  Bin_vector<- vector(,N_bins)
  Critter<- matrix(0, nrow=N_test,ncol=2)
  Distribution<- matrix(0, nrow=N_bins, ncol=N_types)
  
  crit<- 0; crit1<- 0; crit2<- 0; crit3<- 0    
  cStart<-0.5                           # mean for Gaussian controlling setting of initial criterion position
  cVar<-0.2                             # variance for above   
  
  Means<- vector(, N_types); meanError<- vector(, N_types)
  quartMeans<- matrix(0, nrow=N_seqSegment, ncol=N_types); quartError<- matrix(0, nrow=N_seqSegment, ncol=N_types) 
  dpOldMeans<- matrix(0, nrow=N_seqSegment, ncol=N_types/2); dpOldError<- matrix(0, nrow=N_seqSegment, ncol=N_types/2) 
  Summary<- matrix(0,nrow=N_subjects, ncol=N_types)
  summaryQ<- array(0, dim=c(N_seqSegment, N_types, N_subjects))
  dPoldQ<- array(0, dim=c(N_seqSegment, N_types/2, N_subjects))
  
  ######### functions   ##########################
  # standard error    
  se<- function(x){  
    sd(x)/sqrt(length(x))
  }
  qMean<- function(x){
    .Internal(mean(x))
  }
  # calculate 95% CI based on t-distribution using N_subjects-1 df
  CI<-function(x){
    y<- abs(qt(.025, N_subjects-1))
    return(x*y)
  }
  # Hintzman's activation function; cube of similarity, which is dot product of two vectors, divided by 
  # the number of relevant features (i.e, those features that are non-zero for at least one of the two vectors)
  getSim<-function(VectA, VectB) {
    VectC<-abs(VectA)+abs(VectB)
    nR<- length(VectC[VectC != 0])
    return((VectA %*% VectB)/nR)
  }  
  #constructs a single echo, normalizes it, computes discrepancy with studied item, and then
  # applies probabilistic weighting of elements for encoding. Is applied row-wise to memory matrix within makeEcho().
  constructE<- function(z){
    x<- getSim(Probe,z)^3 
    Echo<- Echo+(x*z)
    return(Echo)
  }
  
  #### set up distribution vector #######
  Bin_vector[1]<- Bin_start
  for (k in 2:N_bins) {
    Bin_vector[k]<- Bin_vector[k-1] + Bin_size 
  }
  #### start N_subject simulations  
  crit2<- 1
  crit3<- 0
  for (N_reps in 1: N_subjects) {
    crit<- rnorm(1, 0.5, 0.2)
    crit1<- crit                                                # crit1 used to record initial crit setting
    SF<- SF1
    
    Qholo<- matrix(sample(c(-1,1), replace=T, N_items*N_el), nrow=N_items, byrow=T)
    Qholn<- matrix(sample(c(-1,1), replace=T, N_items*N_el), nrow=N_items, byrow=T)
    Qhnlo<- matrix(sample(c(-1,1), replace=T, N_items*N_el), nrow=N_items, byrow=T)
    Qhnln<- matrix(sample(c(-1,1), replace=T, N_items*N_el), nrow=N_items, byrow=T)    
    Items<- array(0, dim=c(N_items,N_el,N_types) )
    
    Memory<- matrix(0,nrow=0, ncol=N_el)    
    Echo<- vector(, N_el)
    Probe<- vector(, N_el)
    
    HOLO<- 0; HOLN<- 0; HNLO<- 0; HNLN<- 0        
    LLcnt<- 0; HLcnt<- 0; LTcnt<- 0; HTcnt<- 0                                              # counts hits/fas
    HFoldLFold<- vector(,N_seqSegment);  HFoldLFnew<- vector(,N_seqSegment)                 # counts p("old") across different sequence conditons
    HFnewLFold<- vector(,N_seqSegment); HFnewLFnew<- vector(,N_seqSegment)  
    
    ######### test trial quartets consisting of 2 lf (old/new) and 2 hf (old/new) trials, 
    ######### then combine into Item array, with item types on separate pages. 
    Items[1:(N_items/2), ,1]<-            Qholo[LFind, ]
    Items[((N_items/2)+1):N_items, ,1]<-  Qhnlo[LFind, ]
    Items[1:(N_items/2), ,2]<-            Qholo[HFind, ]
    Items[((N_items/2)+1):N_items, ,2]<-  Qholn[HFind, ]
    Items[1:(N_items/2), ,3]<-            Qhnln[LFind, ]
    Items[((N_items/2)+1):N_items, ,3]<-  Qholn[LFind, ]
    Items[1:(N_items/2), ,4]<-            Qhnln[HFind, ]
    Items[((N_items/2)+1):N_items, ,4]<-  Qhnlo[HFind, ]
    
    ########## test array ##################
    # Make an array of individual test items by randomly combining trial quartets.
    # This ensures random ordering of trials (outside of quartets), while ensuring 
    # that there are an equal number of test items in each test condition 
    # weave together quartets from quartet arrays, and form testOrder vector (coding for individual trial conditions)  from quartetOrder vector     
    quartetOrder<- c(sample(rep(seq(1:N_types), each = N_items/N_types)), sample(rep(seq(1:N_types), each = N_items/N_types)),
                     sample(rep(seq(1:N_types), each = N_items/N_types)),sample(rep(seq(1:N_types), each = N_items/N_types)))
    j<-1; k<- 1; m<- 1; n<- 1;  v<- 1
    for (i in 1:(N_test/N_seqSegment) )  {
      switch(quartetOrder[i], 
             "1" = {testItems[j:(j+3), ]<- Qhnln[k:(k+3), ]
                    testOrder[j:(j+1)]<- 3                                 
                    testOrder[(j+2):(j+3)]<-  4                
                    k<- k+4 },
             "2" = {testItems[j:(j+3), ]<- Qhnlo[m:(m+3), ]
                    testOrder[j:(j+1)] <- 1    
                    testOrder[(j+2):(j+3)]<- 4
                    m<- m+4 },
             "3" = {testItems[j:(j+3), ]<- Qholn[n:(n+3), ]
                    testOrder[j:(j+1)]<- 3    
                    testOrder[(j+2):(j+3)]<- 2
                    n<- n+4 },
             "4" = {testItems[j:(j+3), ] = Qholo[v:(v+3), ]
                    testOrder[j:(j+1)] = 1    
                    testOrder[(j+2):(j+3)] = 2
                    v<- v+4 }
      )
      j<- j+4
    }
    #test QUARTET order: 1 = HnLn, 2 = HnLo, 3 = HoLn, 4 = HoLo
    # test ITEM order: 1= Low old, 2 = High old, 3 = Low new, 4 = High new 
    ####### pre-study: how accessible are memories from before the study period #######  
    j<- 0
    for (i in 1:N_types) {
      if (i == 1 || i == 3) Frequency<- LF
      if (i == 2 || i == 4) Frequency<- HF
      for (k in 1:Frequency) {        
        pV<- matrix(sample(c(0,1), replace=T, prob=c(1-A,A), N_items*N_el), nrow=N_items, ncol=N_el)
        Memory<- rbind(Memory, pV*Items[, , i])
      }
    }    
    ##### study #####       
    studyItems[1:N_items, ]<- Items[ , ,1]
    studyItems[(N_items+1):N_study, ]<- Items[ , ,2]   
    studyItems[sample(nrow(studyItems)),]
    for (j in 1:N_study) {
      Echo<- vector(, N_el)
      Probe<- studyItems[j, ]   
#       for (i in 1:nrow(Memory)){Echo<- Echo + Memory[i,]*getSim(Memory[i,], Probe )^3 }
      Echo<- apply(Memory, 1, constructE )
      Echo<- as.vector(t(rowSums(Echo)))
      Echo<- Echo/max(abs(Echo))
      pE<- sample(c(0,1),N_el,replace=T, prob=c(1-L,L))
      Echo<- pE*(Probe-Echo)      
      Memory<- rbind(Memory,Echo) 
    }   
    #### test ####   
    for (j in 1:N_test) {
      Echo<- vector(,N_el)
      Probe<- testItems[j,]  
      Echo<- apply(Memory, 1, constructE )
      Echo<- as.vector(t(rowSums(Echo)))
      Echo<- Echo/max(abs(Echo)) 
      Intensity<- 0.0
      Intensity<- getSim(Probe, Echo)             # clarity intensity
      Evidence<- Intensity - crit       
      
      if ( (crit2 == 1 && crit1 <= 0.25) || (crit2 == 2 && crit1  >= 0.75)) {
        Critter[j,crit2]<- crit    
        crit3<- 1
      }      
      # construction of presentation quarters: evaluate performance for 10 presentations of HF quartets: 
      # 2 HF trials per HF quartet, so "w" counts off 20-trial blocks (using standard values for N_items, etc.)
      w<- 0
      switch(testOrder[j],
             "1"={},
             "2"={
               if (j > 2 && testOrder[j-2] == 1) {
                 HOLO<- HOLO+1
                 w<- ceiling(HOLO/N_trialPres)                                    # N_trialPres = 10
               }else{
                 HOLN<- HOLN+1
                 w<- ceiling(HOLN/N_trialPres)
               }
             },
             "3"={},
             "4"={
               if (j > 2 && testOrder[j-2] == 1) {
                 HNLO<- HNLO+1
                 w<- ceiling(HNLO/N_trialPres)
               }else{
                 HNLN<- HNLN+1
                 w<- ceiling(HNLN/N_trialPres)
               }
             }
      )      
      # scoring a P('old') response
      if (Evidence > 0.0) {    
        switch(testOrder[j],
               "1"={
                 LTcnt<- LTcnt + 1.0 },
               "2"={
                 HTcnt<- HTcnt + 1.0
                 if (j > 2 && testOrder[j-2] == 1) {
                   HFoldLFold[w]<- HFoldLFold[w]+1.0
                 } else {
                   HFoldLFnew[w]<- HFoldLFnew[w]+1.0
                 }
               },
               "3"={
                 LLcnt<- LLcnt + 1.0 },
               "4"={
                 HLcnt<- HLcnt + 1.0
                 if (j > 2 && testOrder[j-2] == 1)  {
                   HFnewLFold[w]<- HFnewLFold[w]+1.0
                 }else{
                   HFnewLFnew[w]<- HFnewLFnew[w]+1.0
                 }
               }
        )
      }                                   # end big if (Evidence > 0)
      crit<- crit + SF*tanh(Evidence)^3
      SF<- SF - dp*SF
      for (k in 1:N_bins) {               
        if (Intensity >= Bin_vector[k]-Bin_size/2 && Intensity <= Bin_vector[k]+Bin_size/2 ) {
          Distribution[k,testOrder[j] ]<- Distribution[k,testOrder[j]] + 1.0 }
      }
      
    }                                                 # end of test loop
    # adjust crit counters to stock collecting after 2 good subjects
    if (crit3 == 1) {
      crit2<- crit2+1 
      crit3 <- 0
    }
    Summary[N_reps, ]<- c( LLcnt/N_items, HLcnt/N_items, HTcnt/N_items, LTcnt/N_items)
    for (i in 1:N_seqSegment ){
      summaryQ[ i, , N_reps]<- c( HFnewLFnew[i]/N_trialPres, HFnewLFold[i]/N_trialPres, 
                                  HFoldLFnew[i]/N_trialPres, HFoldLFold[i]/N_trialPres)
      dPoldQ[i, , N_reps]<- c( HFnewLFnew[i]/N_trialPres - HFnewLFold[i]/N_trialPres, 
                               HFoldLFnew[i]/N_trialPres - HFoldLFold[i]/N_trialPres)
    } 
  }                  # end of simulation (subject ) loop
  
  # averaging across subjects
  Means<- colMeans(Summary)
  meanError<- apply(Summary, 2, se)     
  quartMeans<- apply(summaryQ, c(1,2), qMean)  # each row of quartMeans holds means for one test quarter across test items, with columns holding means for a test item across test quaters
  quartError<- apply(summaryQ, c(1,2), se)
  dpOldMeans<- apply(dPoldQ, c(1,2), qMean)  # each row of quartMeans holds means for one test quarter across test items (Lures, Targets), with columns holding means for a test item across test quaters
  dpOldError<- apply(dPoldQ, c(1,2), se)
  #turn these into data frames  
  dPold<- data.frame(Quartiles=factor(8, levels = c(1,2,3,4)), "Test.Item"=factor(8, levels=c("Lure", "Target")), 
                     "Mean.deltaP" = numeric(8),SE = numeric(8))
  dPold[,1]<- rep(c(1:4), times=2)
  dPold[1:4,2]<- "Lure"
  dPold[5:8,2]<- "Target"
  dPold[1:4,3]<- dpOldMeans[,1]
  dPold[5:8,3]<- dpOldMeans[,2]
  dPold[1:4,4]<- dpOldError[,1]
  dPold[5:8,4]<- dpOldError[,2] 
  dPold<- dPold %>% mutate("CI"=CI(SE))
  
  qMeans<- data.frame(Quartiles=factor(16, levels = c(1,2,3,4)), 
                      "Test.Item"=factor(16,levels=c("HFnew(LFnew)", "HFnew(LFold)","HFold(LFnew)","HFold(LFold)" )),
                      "Mean.P.old" = numeric(16),SE = numeric(16))
  qMeans[,1]<- c(rep(seq(1:4), times=4))
  qMeans[1:4,2]<- "HFnew(LFnew)"
  qMeans[5:8,2]<- "HFnew(LFold)"
  qMeans[9:12,2]<- "HFold(LFnew)"
  qMeans[13:16,2]<- "HFold(LFold)"
  qMeans[1:4,3]<- quartMeans[,1]
  qMeans[5:8,3]<- quartMeans[,2]
  qMeans[9:12,3]<- quartMeans[,3]
  qMeans[13:16,3]<- quartMeans[,4]
  qMeans[1:4,4]<- quartError[,1]
  qMeans[5:8,4]<- quartError[,2]
  qMeans[9:12,4]<- quartError[,3]
  qMeans[13:16,4]<- quartError[,4]
  qMeans<- qMeans %>% mutate("CI" = CI(SE))
  
  allMeans<- data.frame("Test.Item"=factor(4, levels=c("LF Lure", "HF Lure","HF Target","LF Target" )), 
                        "Mean.P.old" = numeric(4),SE = numeric(4))
  allMeans[,1]<- as.factor(c("LF Lure", "HF Lure","HF Target","LF Target" ))
  allMeans[,2]<- as.vector( t(t(Means)) )                 # t(t()) because first transpose turns vector into matrix without changing shape
  allMeans[,3]<- as.vector( t(t(meanError)))
  allMeans$Test.Item<- factor(allMeans$Test.Item, levels=c("LF Lure", "HF Lure","HF Target","LF Target" ))
  allMeans<- allMeans %>% mutate("CI"=CI(SE))
  
  Crit<- data.frame(Trial = numeric(N_test*2), Type=factor(N_test*2, levels=c("Low", "High")), Criterion = numeric(N_test*2))
  Crit[,1]<- rep(seq(1:N_test), times = 2)   
  Crit[,2]<- as.factor(rep(c("Low", "High"), each = N_test ))
  Crit[1:N_test,3]<- Critter[,1]
  Crit[(N_test+1):(N_test*2),3]<- Critter[,2]
  
  Dist<- data.frame(bins=numeric(N_bins*N_types), 
                    "Test.Item"=factor(N_bins*N_types, levels=c("LF Target", "HF Target", "LF Lure","HF Lure" ) ),
                    Count = numeric(N_bins*N_types) )
  Dist[,1]<- rep(Bin_vector, times=N_types)
  Dist[,2]<- as.factor(rep(c("LF Target", "HF Target","LF Lure", "HF Lure" ), each = N_bins))
  Dist[1:N_bins,3]<- Distribution[,1]
  Dist[(N_bins+1):(N_bins*2),3]<- Distribution[,2]
  Dist[((N_bins*2)+1):(N_bins*3),3]<- Distribution[,3]
  Dist[((N_bins*3)+1):(N_bins*4),3]<- Distribution[,4]
  
  outs<- list(dPold, qMeans, allMeans, Crit, Dist)
  return(outs)
}                    # end of function

# set up multiple plots with ggplot, From Cookbook for R
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
