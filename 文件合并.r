rm(list = ls())

library(tidyverse)

acuce.kegg = NULL

files = dir('./KEGG mapping download/')

for (i in files) {
  f = data.table::fread(paste0('./KEGG mapping download/',i),header = FALSE) %>% 
    mutate(GeneID = '')
  colnames(f)[1:6] = c('Query',	'KO',	'Definition',	'Score',	'Second', 'best')
  
  f = dplyr::select(f, c('GeneID','Query',	'KO',	'Definition',	'Score',	'Second', 'best'))
  
  for (j in 1:nrow(f)) {
    f$GeneID[j] = stringr::str_split(f$Query[j],'\\.')[[1]][]
  }
  
  acuce.kegg = rbind(acuce.kegg,f)
}

write.table(acuce.kegg,file = 'Acuce.protein.KEGG.mapping.txt', sep = '\t',row.names = FALSE,quote = FALSE)

keggid = unique(acuce.kegg$KO) %>% 
  as.data.frame() 
colnames(keggid) = 'keggID'
  
keggid %>% 
  dplyr::filter(keggID != '') %>% 
  mutate(id = 1:nrow(.)) %>% 
  dplyr::select('id','keggID') %>% 
  write.table(file = 'KEGG_ID.txt',sep = '\t',row.names = FALSE, quote = FALSE)


# 合并爬虫结果文件
dir2 = dir('./爬虫结果(带Class)/')

kegg_res = NULL

for (i in dir2) {
  f2 = data.table::fread(paste0('./爬虫结果(带Class)/',i),header = TRUE,sep = '\t')
  colnames(f2)[1] = 'KO'
  
  kegg_res = rbind(kegg_res,f2)
}

acuce.kegg.2 = acuce.kegg %>% 
  dplyr::select(GeneID,'KO') %>% 
  dplyr::filter(KO != '') %>% 
  dplyr::mutate(temp = paste0(GeneID,KO)) %>% 
  group_by(temp) %>% 
  dplyr::distinct(GeneID,KO) %>% 
  ungroup() %>% 
  dplyr::select(GeneID,KO)


merge(acuce.kegg.2, kegg_res, by = 'KO',all = TRUE) %>% 
  dplyr::select(GeneID,KO,KEGG_pathway_ID,KEGG_pathway_Description,
                `KEGG First Class`,`kegg_class_second`,`PMID Link`,
                Species, kingdom, Organisms) %>% 
  dplyr::mutate(temp = paste0(GeneID,KEGG_pathway_ID,KEGG_pathway_Description)) %>% 
  group_by(temp) %>% 
  dplyr::distinct(GeneID,KO,KEGG_pathway_ID,KEGG_pathway_Description,
                  `KEGG First Class`,`kegg_class_second`,`PMID Link`,
                  Species, kingdom, Organisms) %>% 
  ungroup() %>% 
  dplyr::select(-temp) %>% 
  dplyr::filter(kingdom %in% c('Plants','')) %>% 
  write.table('Acuce.KEGG.all.results.20210626.txt',row.names = FALSE,quote = FALSE,sep = '\t')





















