#!/usr/bin/env bash


if [[ ! -d "classifier_gen" ]]; then
    mkdir classifier_gen
fi

cd classifier_gen

if [ ! -f "Silva_132_release.zip" ]; then 
    wget https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_132_release.zip
fi

if [[ ! -d "SILVA_132_QIIME_release/" ]]; then
    unzip Silva_132_release.zip
fi


qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path $PWD/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna \
--output-path silva_132_99_16s.qza

qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path $PWD/SILVA_132_QIIME_release/taxonomy/16S_only/99/taxonomy_7_levels.txt \
--output-path ref-taxonomy.qza 

echo "Classifier over the whole 16s sequence"
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads silva_132_99_16s.qza \
--i-reference-taxonomy ref-taxonomy.qza \
--o-classifier 16s-whole-seq-classifier.qza


echo "Classifier trained on 515f 806r region"
qiime feature-classifier extract-reads \
--i-sequences silva_132_99_16s.qza \
--p-f-primer GTGCCAGCMGCCGCGGTAA \
--p-r-primer GGACTACHVGGGTWTCTAAT \
--p-trunc-len 120 \
--p-min-length 100 \
--p-max-length 400 \
--o-reads 515-806-ref-seqs.qza

qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads 515-806-ref-seqs.qza \
--i-reference-taxonomy ref-taxonomy.qza \
--o-classifier 515-806-classifier.qza
