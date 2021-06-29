# coding: utf-8
import os
from numpy import heaviside, isin, unicode_
import pandas as pd
import sys
import subprocess
import requests
from bs4 import BeautifulSoup
import re
import time

pep = pd.read_table('D:/坚果云/科研相关/数据/纯系/基因组/genome/Python处理Acuce蛋白序列/pep_seq.txt',sep=' ',header=None,names=['pepID','geneID','sequence'])

path = str.replace(os.getcwd(),'\\','/') + '/'

# 提取蛋白序列
file_name = path + sys.argv[1]

gene_name = pd.read_table(file_name,header=0,names=['gene'])

df = pep.loc[pep['geneID'].isin(gene_name['gene'])]
df.reset_index(drop=True,inplace=True) 

res_name = path + str.split(sys.argv[1],'.')[0] + '_seq.fasta'

res = open(res_name,'w')
for i in range(0,df.shape[0]):
    res.write(df['pepID'][i] + '\n')
    res.write(df['sequence'][i] + '\n')
res.close()

pri = 'echo -----------------------------------------------\n'
subprocess.call(pri,shell=True)
pri = 'echo Protein sequence extraction has been completed！\n'
subprocess.call(pri,shell=True)
pri = 'echo -----------------------------------------------\n'
subprocess.call(pri,shell=True)

# 蛋白序列比对
blast_out = str.split(res_name,'_')[0] + '_swissprot_blast.txt'
blast_out = '/' +str.split(blast_out,':')[0].lower() + str.split(blast_out,':')[1]

res_name = '/' +str.split(res_name,':')[0].lower() + str.split(res_name,':')[1]

blast_shell = 'blastp -query ' + res_name + ' -db /e/Blast+/database/swissprot/blastp/swissprot -num_threads 12 -evalue 1e-5 -outfmt 6 -out ' +  blast_out

bash = open('run_blastp.sh','w')
bash.write(blast_shell)
bash.close()

if str(sys.argv[4]) == 'T':
    subprocess.call('bash run_blastp.sh',shell=True)
else:
    pri = 'echo 已经完成Blast，无需再次Blast!\n'
    subprocess.call(pri,shell=True)
    pri = 'echo -----------------------------------------------\n'
    subprocess.call(pri,shell=True)


os.remove('./run_blastp.sh')

# 蛋白序列比对结果筛选
blast_final_file_name = sys.argv[1].split('.')[0] + '_swissprot_blast.txt'
if str(sys.argv[4]) == 'T':
    blast_final_file_name = sys.argv[1].split('.')[0] + '_swissprot_blast.txt'
    col_names = ['Query id','Subject id','identity','alignment length','mismatches','gap openings','q.Start','q.End','s.Start','s.End','E value','score']

    pep_blast = pd.read_table(blast_final_file_name,sep='\t',header=None,names=col_names)
    pep_blast = pep_blast.loc[pep_blast['identity'] >= float(sys.argv[2])]
    pep_blast.reset_index(drop=True,inplace=True) 

    pep_blast.to_csv(blast_final_file_name,header=True,index=False)

    pri = 'echo Protein sequence alignment has been completed！\n'
    subprocess.call(pri,shell=True)
    pri = 'echo -----------------------------------------------\n'
    subprocess.call(pri,shell=True)

# 蛋白信息爬虫
pep_file = pd.read_table(blast_final_file_name,sep=',',header=0)

bash = open('run_temp.sh','w')
bash.write('echo Crawler started......')
bash.close()

subprocess.call('bash run_temp.sh',shell=True)
os.remove('./run_temp.sh')

if str(sys.argv[3]) == '0':
    os.makedirs('spider_res')

for i in range(int(sys.argv[3]),pep_file.shape[0]):
    if (pep_file.shape[0] - int(sys.argv[2])) > 50:
        time.sleep(5)
    spider_res = pd.DataFrame(columns=('Subject id','pep id','organism','gene name','protein','status','pes seq original','pep seq'))

    pep_id = str.split(pep_file['Subject id'][i],'|')[1]

    spider_link = 'https://www.uniprot.org/uniprot/' + pep_id
    res_spider = requests.get(spider_link)
    soup = BeautifulSoup(res_spider.text,'html.parser')

    organism = soup.find('div',id = 'content-organism', class_ = 'entry-overview-content').get_text()
    gene_name = soup.find('div',id = 'content-gene', class_ = 'entry-overview-content').get_text()
    protein = soup.find('div',id = 'content-protein', class_ = 'entry-overview-content').get_text()
    status = soup.find('div', id = 'content-status',class_ = 'entry-overview-content').find('span', class_ = 'context-help tooltipped-click').get_text()
    res_str = re.findall('<p>(.*?)</p>',status)
    status = status.replace(res_str[0],'').replace('\n','.').replace('                                    <p></p>','').replace('-','').replace('leveli','level')
    status = str.replace(status,'.','')

    # 蛋白序列爬取
    fasta_link = 'https://www.uniprot.org/uniprot/' + pep_id + '.fasta'
    res_fasta = requests.get(fasta_link)
    soup_fasta = BeautifulSoup(res_fasta.text,'html.parser').get_text()

    res = []
    res.append(soup_fasta)

    sequence_original = res[0].replace('\n','')
    p1 = re.compile(r'[>](.*)[=]', re.S)
    start_index = len(re.findall(p1, sequence_original)[0]) + 3
    end_index = int(len(sequence_original))
    seq_final = sequence_original[start_index:end_index]
    
    spider_res.loc[spider_res.shape[0]] = [pep_file['Subject id'][i],pep_id,organism,gene_name,protein,status,sequence_original,seq_final]
    
    fmt = format((i + 1)/(pep_file.shape[0])*100,'.2f')
    info = 'echo Note:' + str(i+1) + ' out of ' + str(pep_file.shape[0]) + ' has been completed[' + str(fmt) + '%]!' +'\n'
    bash = open('run_temp.sh','w')
    bash.write(info)
    bash.close()

    subprocess.call('bash run_temp.sh',shell=True)
    os.remove('./run_temp.sh')

    spider_name = './spider_res/' + sys.argv[1].split('.')[0] + '_'+ str(i) +'_blast_spider.txt'
    spider_res.to_csv(spider_name,sep='\t',header=True,index=False)

    if i == (pep_file.shape[0]-1):
        files_spider = os.listdir('./spider_res/')

        res_spider_final = pd.DataFrame(columns=('Subject id','pep id','organism','gene name','protein','status','pes seq original','pep seq'))

        for file in files_spider:  
            file_path = './spider_res/' + file
            file_temp = pd.read_table(file_path,sep='\t',header=0)
            
            res_spider_final = pd.concat([res_spider_final,file_temp],axis=0)
    
        spider_res = pd.merge(pep_blast,res_spider_final,on='Subject id',how='left')

        spider_name = sys.argv[1].split('.')[0] + '_blast_spider.txt'
        spider_res.to_csv(spider_name,sep='\t',header=True,index=False)
        spider_name = sys.argv[1].split('.')[0] + '_blast_spider.xlsx'
        spider_res.to_csv(spider_name,header=True,index=False)

        pri = 'echo Uniprot crawler has been completed！\n'
        subprocess.call(pri,shell=True)
        pri = 'echo -----------------------------------------------\n'
        subprocess.call(pri,shell=True)


