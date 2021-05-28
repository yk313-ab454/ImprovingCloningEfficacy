library("stringr")
##
findMed_LINCS<-function(strGEOID,cls=strGEOID,chip,from,to,suffix_disease){

    lincsNamestbl=read.table("./data/lincs/LINCS_compounds.txt",sep="\t",row.names=2)
    lincsNames=lincsNamestbl[,1]
    names(lincsNames)=rownames(lincsNamestbl)
    ##
    for(i in 0:19){
      gmt_up_file= paste0("./data/lincs/lincs_upsig_symb_",i,"of20.gmt")
      gmt_down_file=paste0("./data/lincs/lincs_downsig_symb_",i,"of20.gmt")
      suffix=paste0("_",suffix_disease,"_",i)
      findMed(strGEOID=strGEOID,cls=cls,chip=chip,gmt_up_file=gmt_up_file,gmt_down_file=gmt_down_file,from=from,to=to,suffix=suffix,lincsNames)
      resFile1=paste0("./Results/",strGEOID,suffix,".txt")
    }
    mergResults(strGEOID,suffix_disease);
}
##
mergResults<-function(strGEOID,suffix_disease){
  tbl_res=c();
  for(i in 0:19){
    suffix=paste0("_",suffix_disease,"_",i)
    resFile=paste0("./Results/",strGEOID,suffix,".txt")    
    tbl= read.table(resFile,header=TRUE,sep="\t")
    file.remove(resFile)
    tbl_res=rbind(tbl_res,tbl)
  }
  ii=order(as.numeric(tbl_res[,"score"]));
  tbl_res=tbl_res[ii,];
  tbl_res[,"rank"]=1:nrow(tbl_res)
  outFile=paste0("./Results/",strGEOID,"_",suffix_disease,".txt");
  write.table(tbl_res,outFile,sep="\t",row.names=FALSE,quote=FALSE)
  
}
##
findMed<-function(strGEOID,cls=strGEOID,chip,gmt_up_file,gmt_down_file,from,to,suffix="",lincsNames){
  
  strGct_DB= paste0("./data/",strGEOID,".gct");
  strCls=paste0("./data/cls_",cls,".cls");
  strPlatform=paste0("./data/",chip,".chip");
  
  dir.create(paste0("./output/",strGEOID), showWarnings = FALSE)

  for(strDirection in c("up","down")){
    strFolderName = paste0("output/" , strGEOID) ;
    strFolderName_dir = paste0( strFolderName , "/" , strDirection);
    strGMT= if(strDirection=="up") {gmt_up_file}else{gmt_down_file} 

    strCMD = paste0("java -cp ./gsea/gsea2-2.2.4.jar -Xmx5000m xtools.gsea.Gsea -res " , strGct_DB ," -cls " , strCls , paste0("#",from)  , paste0("_versus_",to) ,  " -gmx " , strGMT , " -collapse true -mode Max_probe -norm meandiv -nperm 0 -permute gene_set -rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis -metric Signal2Noise -sort real -order descending -chip " , strPlatform , " -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 0 -rnd_seed timestamp -save_rnd_lists false -zip_report false -out " , strFolderName_dir , " -gui false -set_min 1")
    system(strCMD)
  }
  mergScores(strGEOID,from,to,suffix,lincsNames);
}
##
mergScores<-function(strGEOID,from,to,suffix,lincsNames){
  
  dfUp <- file.info(list.files(paste0("./output/",strGEOID,"/up/"), full.names = T));
  upFolder=rownames(dfUp)[which.max(dfUp$mtime)][1]
  
  dfdown <- file.info(list.files(paste0("./output/",strGEOID, "/down/"), full.names = T));
  downFolder=rownames(dfdown)[which.max(dfdown$mtime)][1]
  
  file_up_h=  dir(upFolder,pattern = paste0("gsea_report_for_",to,"_.*xls"),full.names = TRUE);
  file_up_d=dir(upFolder,pattern = paste0("gsea_report_for_",from,"_.*xls"),full.names = TRUE);
  file_down_h= dir(downFolder,pattern = paste0("gsea_report_for_",to,"_.*xls"),full.names = TRUE);
  file_down_d=  dir(downFolder,pattern = paste0("gsea_report_for_",from,"_.*xls"),full.names = TRUE);
  
  tbl_up_h=read.table(file_up_h,header=TRUE,sep="\t")
  tbl_up_d=read.table(file_up_d,header=TRUE,sep="\t")
  tbl_up= rbind(tbl_up_h,tbl_up_d);
  tbl_up=tbl_up[order(tbl_up[,1]),];
  
  tbl_down_h=read.table(file_down_h,header=TRUE,sep="\t")
  tbl_down_d=read.table(file_down_d,header=TRUE,sep="\t")
  tbl_down=rbind(tbl_down_h,tbl_down_d);
  tbl_down=tbl_down[order(tbl_down[,1]),];
  
  score_Up=tbl_up[,"ES"];
  score_down=tbl_down[,"ES"];
  score=(as.numeric(score_Up) - as.numeric(score_down))/2;
  sp= str_split( tbl_up[,"NAME"],"_|-");
  instance_ID = tbl_up[,"NAME"];
  id=paste0("BRD-", unlist(lapply(sp, function(x)return(x[5]))))
  cmap_name = lincsNames[id]
  
  tbl_results = cbind(instance_ID,cmap_name,score_Up,score_down,score);
  tbl_results = tbl_results[order(score),];
  rank=1:nrow(tbl_results);
  cID=rep(NA,nrow(tbl_results));
  tbl_results=cbind(rank,tbl_results,cID);
  colnames(tbl_results)= c("rank","instance_id","cmap name","score_up","scoreDown","score","cID");
  write.table(tbl_results,paste0("./Results/",strGEOID,suffix,".txt"),sep="\t",row.names = FALSE,quote = FALSE);
}
###
