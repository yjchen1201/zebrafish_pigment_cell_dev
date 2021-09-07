### Generate bar plots representing frequency and distribution of iridophore-associated DM/ARs with a particular TF motif ###
# ALX1_ALX3_ALX4_GBX2_TFEC_SOX10_MOTIF_occurrence_in_Iri_DMR_DAR_DMAR plot
data <- read.table("MOTIF_occurrence_in_24vsIri_DMAR_table.txt", sep = "\t", header = T, stringsAsFactors=F)
data<-as.data.frame(data)
data$MOTIF<-as.factor(data$MOTIF)
data$MOTIF <- factor(data$MOTIF, levels = c("ALX1","ALX3","ALX4","GBX2","TFEC","SOX10"))
data$DMAR_type <- factor(data$DMAR_type, levels = c("solo_closeDAR","solo_openDAR","solo_hyperDMR","solo_hypoDMR","hypo_closing_DMAR","hypo_opening_DMAR"))

mypalette <- brewer.pal(n = 12, name = "Paired")
#mypalette <-c("#740001","#ae0001","#e35d6a","#ffdfba","#ffffba","#baffc9","#bae1ff","#428bca","#d896ff")
p1 <- ggplot(data, aes(x=MOTIF, y=DMAR_wMOTIF,fill=DMAR_type,label = DMAR_wMOTIF)) +
  geom_bar(position="stack",stat="identity",colour = "black")+
  ggtitle("Motifs occurrence")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Motif", y = "# of DM/ARs with motif")+scale_fill_manual(values = mypalette[c(1,2,3,4,9,10)])+scale_color_manual(values=mypalette[c(1,2,3,4,7,8)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ 
  theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))+geom_text(size = 3, fontface = "bold",position = position_stack(vjust = 0.5))+coord_flip()
p1

p2<- ggplot(data, aes(x=MOTIF, y=MOTIF_pct, label=MOTIF_pct, fill=DMAR_type)) +
  geom_bar(position="stack",stat="identity",colour = "black")+
  ggtitle("Motifs occurrence")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Motif", y = "pct")+scale_fill_manual(values = mypalette[c(1,2,3,4,9,10)])+scale_color_manual(values=mypalette[c(1,2,3,4,7,8)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ 
  theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))+geom_text(size = 3, fontface = "bold",position = position_stack(vjust = 0.5))+coord_flip()
p2

pdf("ALX1_ALX3_ALX4_GBX2_TFEC_SOX10_MOTIF_occurrence_in_Iri_DMR_DAR_DMAR.pdf", width = 10, height = 8)
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
dev.off()

# DMR_DAR_DMAR_occurence plot
data <- read.table("total_DMAR_distribution.txt", sep = "\t", header = T, stringsAsFactors=F)
data<-as.data.frame(data)
data$DMAR_type <- factor(data$DMAR_type, levels = c("solo_closeDAR","solo_openDAR","solo_hyperDMR","solo_hypoDMR","hypo_closing_DMAR","hypo_opening_DMAR"))
data$DMAR <-"DMAR"
p3 <- ggplot(data, aes(x=DMAR, y=DMAR_number, fill=DMAR_type,label=DMAR_number)) +
  geom_bar(position="stack",stat="identity",colour = "black")+
  ggtitle("DM/AR number")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(y = "pct")+scale_fill_manual(values = mypalette[c(1,2,3,4,9,10)])+scale_color_manual(values=mypalette[c(1,2,3,4,7,8)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ 
  theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))+geom_text(size = 3, fontface = "bold",position = position_stack(vjust = 0.5))+coord_flip()
p4 <- ggplot(data, aes(x=DMAR, y=pct, fill=DMAR_type,label=pct)) +
  geom_bar(position="stack",stat="identity",colour = "black")+
  ggtitle("DM/AR number")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs( y = "pct")+scale_fill_manual(values = mypalette[c(1,2,3,4,9,10)])+scale_color_manual(values=mypalette[c(1,2,3,4,7,8)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ 
  theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))+geom_text(size = 3, fontface = "bold",position = position_stack(vjust = 0.5))+coord_flip()

pdf("DMR_DAR_DMAR_occurence.pdf", width = 12, height = 4)
grid.draw(rbind(ggplotGrob(p3), ggplotGrob(p4), size = "last"))
dev.off()