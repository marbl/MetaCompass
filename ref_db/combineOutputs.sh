#!/bin/bash

database_dir=${1}

cd $database_dir/RefSeq_V2_db
mkdir marker_index

mkdir marker_index/ArgS_COG0018
mkdir marker_index/CysS_COG0215
mkdir marker_index/Ffh_COG0541
mkdir marker_index/FtsY_COG0552
mkdir marker_index/Gtp1_COG0012
mkdir marker_index/HisS_COG0124
mkdir marker_index/LeuS_COG0495
mkdir marker_index/PheS_COG0016
mkdir marker_index/RplA_COG0081
mkdir marker_index/RplB_COG0090
mkdir marker_index/RplC_COG0087
mkdir marker_index/RplD_COG0088
mkdir marker_index/RplE_COG0094
mkdir marker_index/RplF_COG0097
mkdir marker_index/RplK_COG0080
mkdir marker_index/RplM_COG0102
mkdir marker_index/RplN_COG0093
mkdir marker_index/RplO_COG0200
mkdir marker_index/RplP_COG0197
mkdir marker_index/RplR_COG0256
mkdir marker_index/RplV_COG0091
mkdir marker_index/RpoA_COG0202
mkdir marker_index/RpoB_COG0085
mkdir marker_index/RpsB_COG0052
mkdir marker_index/RpsC_COG0092
mkdir marker_index/RpsD_COG0522
mkdir marker_index/RpsE_COG0098
mkdir marker_index/RpsG_COG0049
mkdir marker_index/RpsH_COG0096
mkdir marker_index/RpsI_COG0103
mkdir marker_index/RpsK_COG0100
mkdir marker_index/RpsL_COG0048
mkdir marker_index/RpsM_COG0099
mkdir marker_index/RpsO_COG0184
mkdir marker_index/RpsQ_COG0186
mkdir marker_index/RpsS_COG0185
mkdir marker_index/SecY_COG0201
mkdir marker_index/SerS_COG0172
mkdir marker_index/TsaD_COG0533
mkdir marker_index/ValS_COG0525

for gene in $(ls marker_index); do  
   echo $gene
   head -n 1 marker_indexaa/${gene}/${gene}_genome.tsv > marker_index/${gene}/${gene}_genome.tsv
   for file in marker_index[a-z][a-z]/${gene}/${gene}.fna; do 
      cat $file >> marker_index/${gene}/${gene}.fna
      cat ${file/fna/faa} >> marker_index/${gene}/${gene}.faa
      tail -n+2 ${file/.fna/_genome.tsv} >> marker_index/${gene}/${gene}_genome.tsv
   done
done