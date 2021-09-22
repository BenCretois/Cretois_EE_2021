myINLAplot <- function(ModelList,
                    ModelNames = NULL,
                    VarNames = NULL, VarOrder = NULL,
                    Intercept = TRUE, Size = 1,
                    tips = 0.2){
  
  require(dplyr); require(ggplot2); require(INLA); require(MCMCglmm)
  
  Graphlist<-list()
  
  if(!class(ModelList) == "list"){
    ModelList <- list(ModelList)
  }
  
  for(i in 1:length(ModelList)){
    model<-ModelList[[i]]
    
      Graph<-as.data.frame(summary(model)$fixed)
      colnames(Graph)[which(colnames(Graph)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
      colnames(Graph)[which(colnames(Graph)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
      colnames(Graph)[which(colnames(Graph)%in%c("mean"))]<-c("Estimate")
    
    Graph$Model<-i
    Graph$Factor<-rownames(Graph)
    
    Graphlist[[i]]<-Graph
    
  }
  
  Graph <- bind_rows(Graphlist)
  
  Graph$Model <- as.factor(Graph$Model)
  
  if(!is.null(ModelNames)){
    levels(Graph$Model)<-ModelNames
  }
  
  position <- ifelse(length(unique(Graph$Model))  ==  1, "none", "right")
  
  if(is.null(VarOrder)) VarOrder <- rev(unique(Graph$Factor))
  if(is.null(VarNames)) VarNames <- VarOrder
  
  Graph$Factor <- factor(Graph$Factor, levels = VarOrder)
  levels(Graph$Factor) <- VarNames
  
  Graph %<>% as.data.frame %>% filter(!is.na(Factor))
  
  if(Intercept == FALSE){
    
    VarNames <- VarNames[!str_detect(VarNames, "ntercept")]
    
    Graph <- Graph %>% filter(Factor %in% VarNames)
    
  }
  
  min<-min(Graph$Lower,na.rm = T)
  max<-max(Graph$Upper,na.rm = T)
  
  return(Graph)
  #ggplot(Graph,
  #       aes(x = as.factor(Factor),
  #           y = Estimate,
  #           group = Model,
  #           colour = Model)) +
  #  geom_errorbar(position = position_dodge(w = 0.5),
  #                aes(ymin = Lower, ymax = Upper), size = 0.3,
  #                width = tips) +
  #  geom_hline(aes(yintercept = 0),lty = 2) + labs(x = NULL) + coord_flip() +
  #  theme(legend.position = position) +
  #    geom_point(position = position_dodge(w = 0.5), size = 1) +
  #  scale_color_manual(values = c('#2EB872', '#f0a73a', '#364F6B')) +
  #  labs(color = 'RSF')
}
