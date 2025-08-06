
## ./Rfam/

本目录下存储与Rfam数据库数据相关的scripts，以及一些不占用太大空间的表格数据

## 文件说明

Rfam.pdb 来源 https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.seed.gz

rfam_seed_headers.csv 处理来自Rfam中所有家族seed文件的header

rfam_pdb_labeled.csv 合并Rfam.pdb和rfam_seed_headers.csv部分内容得到的，方便后续进行家族特异性的分析，其中CL和TP是比较好的官方分类，当然rfam_acc(Rfam accession)也是，就是类分的有点太细

## TODO

./structure_download_label.py
- download的逻辑似乎和utils/pdb_utils.py有重合，但考虑到start和end的存在，需要至少重构一版数据
- 重命名rfam_pdb_labeled.csv内的所有列名，并写文档说明各列含义