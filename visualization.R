packages <- c("circlize", "ggplot2", "grid", "plyr")
options(repos = c(CRAN = "https://cloud.r-project.org"))
install_pac=setdiff(packages, rownames(installed.packages()))
if (length(install_pac)>0){
	install.packages(install_pac)
}
library(circlize)
library(ggplot2)
library(grid)
library(plyr)

args=commandArgs((trailingOnly=TRUE))
outdir=args[1]
virus_len=as.numeric(args[2])
hg_version=args[3]


dir.create(paste(outdir,"plot",sep='/'))
file_list=list.files(pattern="\\.AVID.final.txt",path=outdir,recursive = TRUE)
file_list=paste(outdir,file_list,sep='/')
merge_dat=lapply(file_list,function(f) {dat=read.table(f,sep='\t',header=T,stringsAsFactors = F) 
                                        dat$sample=gsub(".virusclip2.final.txt","",basename(f))
                                        return (dat)})
merge_dat <- do.call(rbind, merge_dat)
merge_dat_filter=merge_dat[(merge_dat$repeat_mark=="No")&(merge_dat$support_reads_softclip>0)&(merge_dat$chr_human!="chrM"),]
merge_dat_filter$Virus=rep("Virus",nrow(merge_dat_filter))

dd=read.cytoband(species = hg_version)
cytoband.df=dd$df
cytoband.df=rbind(cytoband.df,c("Virus",1,virus_len*100000,0,0))
cytoband.df$V2=as.numeric(cytoband.df$V2)
cytoband.df$V3=as.numeric(cytoband.df$V3)

#chr_col=data.frame(col=c("#f6cb52","#b56576","#50b99b","#802754","#487bea","#8447ff","#76A18AFF","#7cd5f3","#ffb700","#027da9","#B57CB2FF","#fb8500","#9BE965FF","#523C49FF","#ff7bc4","#64A89AFF","#3D3A1AFF","#833446FF","#0B5D2DFF","#d12e64","#f9f110","#ABAAB3FF","#ff9b85","#8ecae6","#de1616"))
#chr_col=data.frame(col=c("orange3","rosybrown4","palegreen4","palevioletred4","royalblue3","purple1","honeydew4","lightblue","orange","steelblue4","thistle4","orange1","lightgreen","mediumorchid4","lightpink2","lightcyan4","snow4","plum4","olivedrab4","maroon4","lightgoldenrod1","lightgray","mistyrose2","lightblue3","red"))
#chr_col=data.frame(col=c("#FF000080", "#00FF0080", "#0000FF80","#ABC4E0FF","#843CB9FF","purple1","honeydew4","lightblue","orange","steelblue4","thistle4","orange1","lightgreen","mediumorchid4","lightpink2","lightcyan4","snow4","plum4","olivedrab4","maroon4","lightgoldenrod1","lightgray","mistyrose2","lightblue3","red"))
#rownames(chr_col)=c(paste("chr",1:22,sep=''),"chrX","chrY","Virus")


circos_human=merge_dat_filter[,c("chr_human","breakpoint_human")]
circos_virus=merge_dat_filter[,c("Virus","breakpoint_virus")]

circos_human$breakpoint_human=as.numeric(circos_human$breakpoint_human)
circos_virus$breakpoint_virus=as.numeric(circos_virus$breakpoint_virus)*100000

circos_human$end=circos_human$breakpoint_human
circos_virus$end=circos_virus$breakpoint_virus

circos_human$value=rep(1,nrow(circos_human))
circos_virus$value=rep(1,nrow(circos_virus))

colnames(circos_human)=c("chr","start","end","value")
colnames(circos_virus)=c("chr","start","end","value")




pdf(paste(outdir,"/plot/circos_plot.pdf",sep='/'),6,6)
p=circos.initializeWithIdeogram(cytoband.df,plotType = c("labels"),chromosome.index=c(paste("chr",1:22,sep=''),"chrX","chrY","Virus"),sort.chr = F)+
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    circos.rect(xlim[1], 0, xlim[2], 1, col =rand_color(1))
    breaks = seq(0, 1e10, by = 1e8)
  }, track.height = 0.15, bg.border = NA)+circos.genomicLink(circos_human, circos_virus,lwd=0.001)
print(p)
dev.off()



hist_dat=circos_virus
hist_dat$breakpoint=circos_virus$start/100000
hist_dat=as.data.frame(table(hist_dat$breakpoint),stringsAsFactors = F)
hist_dat$Freq=hist_dat$Freq/sum(hist_dat$Freq)*100
test=data.frame(Var1=as.character(1:virus_len),Freq=rep(0,virus_len))
hist_dat=rbind(hist_dat,test[!(test$Var1%in%hist_dat$Var1),])
hist_dat$Var1=factor(hist_dat$Var1,levels = as.character(1:virus_len))

pdf(paste(outdir,"/plot/virus_breakpoint_histogram.pdf",sep='/'),5,2)
p=ggplot(hist_dat, aes(x=Var1,y=Freq)) +
  geom_bar(stat ="identity", fill="blue",position="identity",width = 15)+
  theme_bw()+theme(panel.background = element_blank(),panel.grid = element_blank(),panel.border = element_rect(linewidth=1))+scale_x_discrete(breaks=seq(0,virus_len,1000))+xlab("Virus genome")+ylab("% Virus integration")
print(p)
dev.off()

















































