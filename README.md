# “Acuce”蛋白序列处理相关代码

`pep_process.py`功能：输入“Acuce”的基因ID、相似度、爬虫起始点及是否Blast等参数，返回`UbiProt`的爬虫结果。相似度参数是用于筛选Blastp的结果。

```shell
python /d/坚果云/GitHub/protein-sequence-processing/pep_process.py degs.txt 50 543 F
```

`KEGG爬虫.py`用于爬取KEGG相关信息。

`文件合并.r`用于合并爬虫后的结果。

