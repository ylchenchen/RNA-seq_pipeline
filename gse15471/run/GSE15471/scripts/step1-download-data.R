
if( file.exists(f) ){
  load(f)
}else{
  ##输入要下载的芯片数据集；若太慢，请手动
  eSet <- geoChina( gse_number )
  class(eSet)
  length(eSet)
  eSet = eSet[[1]]
  probeM <- exprs(eSet)  #注为方便区分，探针表达矩阵用probeM表示
  #(2)获取临床信息（为了告诉R你的分组信息是什么）
  pd <- pData(eSet)
  #(3)先获取芯片平台编号，再获取芯片注释信息 
  gpl_number <- eSet@annotation;gpl_number
  ids <- AnnoProbe::idmap(gpl_number)
  save(gse_number,pd,probeM,ids,
       gpl_number,file = f )
}               