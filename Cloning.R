source("./findMed_LINCS.R")

Cloning<-function(){
    
    strGEOID="Cloning"
    cls="SCNT2IVF"
    chip="EmbryoGENE"
    from="SCNT"
    to="IVF"
    suffix_disease="SCNT2IVF"
    findMed_LINCS(strGEOID=strGEOID,cls=cls,chip=chip,from=from,to=to,suffix_disease=suffix_disease)
    
    cls="SCNT2IVF1"    
    to="IVF1"
    suffix_disease="SCNT2IVF1"
    findMed_LINCS(strGEOID=strGEOID,cls=cls,chip=chip,from=from,to=to,suffix_disease=suffix_disease)
    merge_results_cloning()
    }
######
merge_results_cloning<-function(){
    results_SCNT2IVF_all=read.csv("./Results/Cloning_SCNT2IVF.txt",sep="\t",header = TRUE);
    results_SCNT2IVF1=read.csv("./Results/Cloning_SCNT2IVF1.txt",sep="\t",header = TRUE);
    
    merged = merge(results_SCNT2IVF_all,results_SCNT2IVF1,by="instance_id")
    ii=order(as.numeric(merged[,"rank.x"]))
    merged=merged[ii,]
    nn=which(merged[,"rank.x"]<30 & merged[,"rank.y"]<300)
    
    candidates=merged[nn,];
    kk=grep("BRD",candidates[,"cmap.name.x"])
    if(length(kk)!=0){candidates=candidates[-kk,]}

    ii=order(candidates[,"rank.x"]+candidates[,"rank.y"])
    candidates=candidates[ii,]

    print("Best candidates for imroving SCNT efficancy are:")
    print(candidates[,"cmap.name.x"]);
    
    print("Best candidate based on IVF samples for embryo1 and both embryos is:")
    print(candidates[1,"cmap.name.x"]);
    
    write.table(candidates,"candidates.txt");
    xx=1
}    

######
Cloning()
