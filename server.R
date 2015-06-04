library(ggplot2)
source("helpers.R")  # sequence() loads dplyr, multiplot() loads grid

shinyServer(
  function(input, output) {
    seqOut<- reactive({
      sequence(input$SF, as.numeric(input$dp),input$L, input$A, as.numeric(input$N_subjects))
      })  # returns an object named "outs", a list of 5 data frames
    
    output$all<- renderPlot({
        if (input$runSim == 0) return()
        putOut<- isolate(seqOut())
        plotOutput(putOut)
        })

    plotOutput<- function(outs){
      # building up plot themes: first, panel and plot properties
      themePro<- theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),   
                       panel.border=element_rect(colour="black", fill=NA), 
                       panel.background=element_rect(fill=NA),                         
                       strip.background=element_rect(fill=NA)) 
      #now text formatting, axes & titles
      themePro<- themePro + theme(axis.text=element_text(vjust=-0.7, size=14, colour="black"),
                                  axis.title.x=element_text(vjust=-0.1,size=16), axis.title.y=element_text(vjust=1.5, size=16), 
                                  plot.title=element_text(size=18), strip.text.x=element_text(size=16),strip.text.y=element_text(size=16)) 
      #now legend formatting                     
      themePro<- themePro + theme(legend.key=element_rect(fill=NA),
                                  legend.title=element_text(size=16), 
                                  legend.text=element_text(size=14), strip.text.y=element_text(angle=0), 
                                  legend.position="top")
      
      dp<- ggplot(outs[[1]], aes(Quartiles, Mean.deltaP, group= Test.Item, linetype=Test.Item, shape=Test.Item))+
        geom_line(size=1.0)+geom_point(size=4)+ylim(min(outs[[1]]$Mean.delta-outs[[1]]$CI),
        max(outs[[1]]$Mean.deltaP+outs[[1]]$CI))+scale_shape_manual(name="HF Test item ", values=c(0,1))+
        scale_linetype_manual(name="HF Test item ", values=c(1,3))+geom_hline(aes(yintercept=0)) + 
        geom_errorbar(aes(ymin=Mean.deltaP-CI, ymax=Mean.deltaP+CI), width=0.15, lty=1)+
        labs(x="Trial Quartile", y=expression("Mean "~Delta~italic(P)~"('old')"))+themePro
      
      means<- ggplot(outs[[3]], aes(Test.Item, Mean.P.old))+geom_bar(stat="identity", fill=c("grey55", "black","black","grey55"))+
                        geom_errorbar(aes(ymin=Mean.P.old-SE, ymax=Mean.P.old+SE),width=0.15, 
                        colour=c("grey55", "black","black","grey55"))+ylim(0.0,1.0)+
                        labs(x="Test Items", y="Mean P('old')")+themePro+theme(legend.position="none")

      crit<- ggplot(outs[[4]], aes(Trial, Criterion, group=Type,linetype=Type, colour=Type))+geom_line(size=1.0)+
        geom_hline(yintercept=0.5, size=1.0)+scale_linetype_manual(name="Initial criterion position",values=c(1,1))+
        scale_colour_manual(name="Initial criterion position", values=c("blue","red"))+ylim(min(outs[[4]]$Criterion-.05),
        max(outs[[4]]$Criterion+0.05))+themePro
      
      outs[[5]]$Test.Item <- as.factor(outs[[5]]$Test.Item)
      outs[[5]]$Test.Item <- factor(outs[[5]]$Test.Item, levels=c("LF Lure", "HF Lure","HF Target","LF Target" ))
      dists<- ggplot(outs[[5]], aes(bins, Count, group=Test.Item,linetype=Test.Item, 
                        colour=Test.Item))+geom_line(size=0.75)+geom_vline(xintercept=0.5, size=1.0)+
                        scale_linetype_manual(name="Test item",values=c(2,2,1,1))+xlab("Test item intensity")+
                        scale_colour_manual(name="Test item", values=c("blue","red","red", "blue"))+xlim(0.225,1.05)+themePro
      
     multiplot(dp, means, crit, dists, cols=2, layout=matrix(c(1,2,3,4), nrow=2, byrow=T))
}
   
})
