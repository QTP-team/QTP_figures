library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggthemes)
library(ggThemeAssist)
library(lemon)
library(tomtom)
library("gapminder")
library(ggnewscale)
library(ggforce)


##################################################
############   Figure E6
{

df=read.table("0.data/dset.txt",
              sep="\t",check.names = F,header = 1)#,fileEncoding = "UTF16")

set1=c("Tibetanhorse","Tibetanass","N1","N5","N4","N3","Tibetancattle","Yak","N2","Tibetansheep",
       "Tibetanantelope")
abb1=c("H","A","OddToed","Primary","EvenToed","Bovine","C","Y","Ovisae","S","AN")
new_abb=data.frame(row.names=abb1)
new_abb$abb=c("TH","TA","N1","N5","N4","N3","TC","Yak","N2","TS","TAN")
set2=c("Tibetanhorse","Tibetanass","N1","N4","N3","Tibetancattle","Yak","N2","Tibetansheep",
       "Tibetanantelope")
df$Source1=new_abb[df$Source,'abb']
df$Target1=new_abb[df$Target,'abb']
df$vsgroup1=paste(df$Target1,df$Source1,sep  = " vs. ")
df$PathwayName
p2=ggplot(df,aes(x=vsgroup1,y=PathwayName))+
  geom_point(aes(x=vsgroup1,y=PathwayName,size=-log10(FDR),fill=UpDown),
             shape = 21,col="grey", stroke = 1)+
  scale_fill_manual(values  = c('#013369','#D50A0A'))+
  theme_bw()+
  theme(plot.margin = unit(c(0,-.5,0,.5), "lines"))+
  facet_grid(Level2~VSclass,
                     space="free",
                     scales = "free")
p2=p2+  theme(axis.text.x = element_text(angle = 45, vjust = 1,size = 12, hjust = 1))
ggsave("1.result/sel5.CAZy.KEGG.pdf",plot = p2,width = 10,height = 15,limitsize = FALSE)

inter_df=df[df$VSclass=='Inter-group',]

for (i in rownames(inter_df)){
  mtx=inter_df[i,'matrix']
  tmp=as.numeric(unlist(strsplit(mtx,split=";")))
  inter_df[i,'Source_rel']=tmp[1]/tmp[2]
  inter_df[i,'Source_num']=tmp[1]
  inter_df[i,'Target_num']=tmp[3]
  inter_df[i,'Target_rel']=tmp[3]/tmp[4]
  if(inter_df[i,'p-value']>0){
    inter_df[i,'Source_alpha']=0.8
    inter_df[i,'Target_alpha']=0.2
    inter_df[i,'Source_shape']=19
    inter_df[i,'Target_shape']=13
    
  }else{
    inter_df[i,'Source_alpha']=0.2
    inter_df[i,'Target_alpha']=0.8
    inter_df[i,'Source_shape']=13
    inter_df[i,'Target_shape']=19
  }
}
inter_df$Target_shape=as.numeric(inter_df$Target_shape)
inter_df$`-log10(FDR)`=-log10(inter_df$FDR)


CAZy_class=c("Multiple polysaccharides","Pectin","Chitin",
             "Cellulose and other beta-glucans, hemicelluloses",
             "Polysialic_acid","Plant cell wall","Starch and storage glycans",
             "Xylan","Host O- and N-linked glycans","Glycosaminoglycans","Mannan",
             "Fungal cell wall polysaccharides")
CAZy_class_abb=c("M1","Pectin","Chitin","M5","M7","M8","M2","M3","M4","M9","Mannan","M6")
CAZy_abb_df1=data.frame(row.names = CAZy_class_abb)
CAZy_abb_df1$abb=CAZy_class
n=inter_df[inter_df$Level2=='CAZy','Pathway']
inter_df[inter_df$Level2=='CAZy','PathwayName']=CAZy_abb_df1[n,'abb']
inter_df=inter_df[inter_df$FDR<0.05,]
inter_df$type="T2"
inter_df[inter_df$vsgroup1 %in% c("N4 vs. N1","N2 vs. N3"),'type']="T1"


df_test=dcast(inter_df,Pathway+PathwayName+Level2~vsgroup1,value.var='FDR')

rownames(df_test)=df_test$Pathway
pth_order=c("M5","M8","M4","Chitin","M7","M1","M9","Pectin","M6","M2",
            "Mannan","M3","map00010","map00020","map00640","map00040",
            "map00051","map00520","map00562","map00660","map00620","map00030",
            "map00650","map00052","map00500","map00053","map00630","map00561",
            "map00061","map00140","map01040","map00600","map00072","map00590",
            "map00071","map00100","map00564","map00592","map00062","map00120",
            "map00121","map00565","map00260","map00400","map00340","map00280",
            "map00310","map00290","map00380","map00300","map00270","map00250",
            "map00360","map00220","map00330","map00780","map00790","map00770",
            "map00730","map00760","map00130","map00785","map00830","map00750",
            "map00740","map00860","map00670","")
pth_order=c("M4","Chitin","M7","M1","M5","M8","M9","Pectin","M6","M2",
  "Mannan","M3","map00040","map00010","map00020","map00640",
  "map00051","map00520","map00562","map00660","map00620","map00030",
  "map00650","map00052","map00500","map00053","map00630","map00561",
  "map00061","map00140","map01040","map00600","map00072","map00590",
  "map00071","map00100","map00564","map00592","map00062","map00120",
  "map00121","map00565","map00260","map00400","map00340","map00280",
  "map00310","map00290","map00380","map00300","map00270","map00250",
  "map00360","map00220","map00330","map00780","map00790","map00770",
  "map00730","map00760","map00130","map00785","map00830","map00750",
  "map00740","map00860","map00670")
pth_order=c("M4","Chitin","M7","M1","M5","M8","M9","Pectin","M6","M2",
            "Mannan","M3","map00040","map00010","map00020","map00640",
            "map00051","map00520","map00562","map00660","map00620","map00030",
            "map00650","map00052","map00500","map00053","map00630","map00561",
            "map00061","map00140","map01040","map00600","map00072","map00590",
            "map00071","map00100","map00564","map00592","map00062","map00120",
            "map00121","map00565","map00260","map00400","map00340","map00280",
            "map00310","map00290","map00380","map00300","map00270","map00250",
            "map00360","map00220","map00330","map00780","map00790","map00770",
            "map00730","map00760","map00130","map00785","map00830","map00750",
            "map00740","map00860","map00670",
            "map00405","map00521","map00943","map00311","map00333","map00950",
            "map00332","map00997","map00941","map00945","map00261","map00998",
            "map00401","map00960","map00944","map00940","map00524","map00525",
            "map00966","map00930","map00633","map00626","map00621","map00362",
            "map00622","map00643","map00625","map00361","map00623","map00364",
            "map00980","map00642","map00791","map00983","map00624","map00982",
            "map00984","map00540","map00511","map00571","map00531","map00604",
            "map00513","map00541","map00603","map00572","map00550","map00510",
            "map00514","map00515","map00532","map00410","map00480","map00471",
            "map00473","map00440","map00450","map00430","map00472","map00460",
            "map00906","map00908","map01054","map01055","map01051","map00523",
            "map00903","map00909","map01053","map00281","map00900","map01052",
            "map00253","map01057","map00680","map00720","map00920","map00190",
            "map00910","map00240","map00230","map03070","map02060","map02010")
pth_name_order=df_test[pth_order,'PathwayName']

level2_order=c("CAZy","Carbohydrate metabolism","Lipid metabolism",
               "Amino acid metabolism","Metabolism of cofactors and vitamins")
level2_order=c("CAZy","Carbohydrate metabolism",
               "Amino acid metabolism","Metabolism of cofactors and vitamins",
               "Xenobiotics biodegradation and metabolism")

inter_df$Level2=factor(inter_df$Level2,levels = level2_order)
inter_df$Pathway=factor(inter_df$Pathway,levels = rev(pth_order))


sel_tar_enrich=inter_df[inter_df$Target_shape==19,]
sel_sou_enrich=inter_df[inter_df$Source_shape==19,]
dim(sel_tar_enrich)
dim(sel_sou_enrich)
dim(inter_df)
p3=ggplot(inter_df,aes(x=vsgroup1,y=Pathway))+
  geom_point(aes(x=vsgroup1,y=Pathway,size=Source_num,col=-log10(FDR)),
             shape=21,position = position_nudge(x = 0.1))+
  geom_point(data=sel_sou_enrich,aes(x=vsgroup1,y=Pathway,size=Source_num,col=-log10(FDR)),
             shape=19,position = position_nudge(x = 0.1))+
  scale_color_gradient2(high="green4",mid = "lightgreen")+
  new_scale_color() +

  geom_point(aes(x=vsgroup1,y=Pathway,size=Target_num,col=-log10(FDR)),
             shape=21,position = position_nudge(x = -0.1))+
  geom_point(data=sel_tar_enrich,aes(x=vsgroup1,y=Pathway,size=Target_num,col=-log10(FDR)),
             shape=19,position = position_nudge(x = -0.1))+
  
  scale_color_gradient2( high="blue4",mid="royalblue1")+
  scale_y_discrete(position = "right",breaks=pth_order,labels=pth_name_order) +
  facet_grid(Level2~type,space="free",scales = "free",switch = "y")+
  theme_bw()+ 
  theme(legend.position="top")+
  
  theme(strip.placement = "outside")
p3
ggsave("Figure6.v0905.0.pdf",plot = p3,width = 9,height = 11)
if (FALSE){
  p3=ggplot(inter_df,aes(x=vsgroup1,y=Pathway))+
    geom_point(aes(x=vsgroup1,y=Pathway,size=Source_num,col=-log10(FDR)),
               shape=21,position = position_nudge(x = 0.1))+
    geom_point(data=sel_sou_enrich,aes(x=vsgroup1,y=Pathway,size=Source_num,col=-log10(FDR)),
               shape=19,position = position_nudge(x = 0.1))+
    scale_color_gradient2(high="orangered4",mid = "orange1")+
    new_scale_color() +
    geom_point(aes(x=vsgroup1,y=Pathway,size=Target_num,col=-log10(FDR)),
               shape=21,position = position_nudge(x = -0.1))+
    geom_point(data=sel_tar_enrich,aes(x=vsgroup1,y=Pathway,size=Target_num,col=-log10(FDR)),
               shape=19,position = position_nudge(x = -0.1))+
    scale_color_gradient2( high="mediumpurple4",mid="mediumorchid1")+
    scale_y_discrete(position = "right",breaks=pth_order,labels=pth_name_order) +
    facet_grid(Level2~type,space="free",scales = "free",switch = "y")+
    theme_bw()+ 
    theme(legend.position="top")+
    
    theme(strip.placement = "outside")
  p3
  ggsave("1.result/Figure6.v0905.1.pdf",plot = p3,width = 9,height = 11)
}

}
