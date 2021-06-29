import os
import time
from numpy import fabs, inf
import pandas as pd
import requests
from bs4 import BeautifulSoup
from requests.api import head

os.chdir('D:/坚果云/科研相关/数据/纯系/基因组/genome/Acuce蛋白序列处理/月亮谷蛋白KEGG mapping结果/')

kegg_species = pd.read_table('D:/坚果云/科研相关/数据/纯系/基因组/genome/Acuce蛋白序列处理/月亮谷蛋白KEGG mapping结果/KEGG species code.txt',header=0)
koid = pd.read_table('D:/坚果云/科研相关/数据/纯系/基因组/genome/Acuce蛋白序列处理/月亮谷蛋白KEGG mapping结果/KEGG_ID.txt',header=0)

for i in range(3236,koid.shape[0]):
    kegg_res = pd.DataFrame(columns=['keggID','KEGG_pathway_ID','KEGG_pathway_Description'])

    url = 'https://www.kegg.jp/entry/ko:' + koid['keggID'][i]

    res = requests.get(url)
    soup = BeautifulSoup(res.text,'html.parser')
    pathway = soup.find_all('table',style='border:0;border-collapse:collapse;')
    info = soup.find_all('div',style='width:555px;overflow-x:auto;overflow-y:hidden')

    for j in pathway:
        item = j.get_text()
        if str.startswith(item,'k'):
            detail = str.split(item,'\xa0')  
            kegg_res.loc[kegg_res.shape[0]] = {'keggID': koid['keggID'][i], 'KEGG_pathway_ID':detail[0], 'KEGG_pathway_Description':detail[2]}
    

    if True:
        keggid = [koid['keggID'][i]]
        pmid = []
        species = []

        for m in info:
            term = m.get_text()

            if str.startswith(term,'PMID'):
                pmid.append('https://pubmed.ncbi.nlm.nih.gov/' + term.split(':')[1] + '/')
            if str.startswith(term,'['):
                species.append(term.replace('[','').split(':')[0])
        
        max_num = max(len(pmid),len(species))

        keggid = keggid*max_num 

        if max_num == 0:
            keggid.append(koid['keggID'][i])
            pmid.append('NaN')
            species.append('NaN')
        
        if len(pmid) == 0:
            pmid.append('NaN')
        if len(species) == 0:
            species.append('NaN')
        
        kegg_info = pd.concat([pd.DataFrame(keggid),pd.DataFrame(pmid),pd.DataFrame(species)],axis=1)

        #print(kegg_info)

        if False:
            if kegg_info.shape[0] == 0:
                temp_id = koid['keggID'][i]
                kegg_info = pd.DataFrame({'keggID':temp_id,'PMID Link':'NaN','Species':'NaN'},index = 'a')
                #print(kegg_info)

                #kegg_info = kegg_info[kegg_info.EncodedPixels.notnull()]

        kegg_info.columns = ['keggID','PMID Link','Species']

        kegg_info = kegg_info.fillna('None')
        kegg_info = pd.merge(kegg_info,kegg_species,on='Species',how='left')
        #print(kegg_info)


    if kegg_res.shape[0] != 0:

        # 抓取pathway class
        kegg_class_df = pd.DataFrame()
        for k in kegg_res['KEGG_pathway_ID']:

            #print(k)

            url_2 = 'https://www.genome.jp/entry/' + k

            res_2 = requests.get(url_2)
            soup_2 = BeautifulSoup(res_2.text,'html.parser')
            class_info = soup_2.find_all('div',style='width:555px;overflow-x:auto;overflow-y:hidden')

            for j in class_info:
                kegg_class = j.get_text().replace('\n','')
                if str.endswith(kegg_class,'BRITE hierarchy'):
                    if kegg_class == 'BRITE hierarchy':
                        kegg_class_first = 'None'
                        kegg_class_second = 'None'
                    else:
                        kegg_class = kegg_class.replace('BRITE hierarchy','').split(';')
                        kegg_class_first = kegg_class[0]
                        kegg_class_second = kegg_class[1].replace(' ','',1)

                    kegg_class_temp = pd.DataFrame()
                    kegg_dict = {'KEGG_pathway_ID':[k],'KEGG First Class':[kegg_class_first],'kegg_class_second':[kegg_class_second]}
                    kegg_class_temp = kegg_class_temp.append(pd.DataFrame(kegg_dict))
                    kegg_class_df = pd.concat([kegg_class_df, kegg_class_temp],axis=0)

        kegg_res = pd.merge(kegg_res,kegg_class_df,on='KEGG_pathway_ID',how='left')

        kegg_res = pd.merge(kegg_res,kegg_info,on='keggID',how='left')

        filenames = 'D:/坚果云/科研相关/数据/纯系/基因组/genome/Acuce蛋白序列处理/月亮谷蛋白KEGG mapping结果/爬虫结果(带Class)/第' + str(koid['id'][i]) + '(' + str(koid.shape[0] + 1 ) + ')个_'+ koid['keggID'][i] + '.txt'
        kegg_res.to_csv(filenames,sep='\t',index=False,header=True)

    files = os.listdir('D:/坚果云/科研相关/数据/纯系/基因组/genome/Acuce蛋白序列处理/月亮谷蛋白KEGG mapping结果/爬虫结果(带Class)/')

    print('已完成第' + str(koid['id'][i]) + '/' + str(koid.shape[0] + 1 ) + '个！')
    print('已保存' + str(len(files)) + '个！')
    print('---------------------------------------------------------------------')

    time.sleep(0)