scripts_dir=$(pwd)
cd data/

================================================================================
# GTEx cis-eQTL data
wget https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar
# Extract lung and blood data
tar -xf GTEx_Analysis_v8_eQTL.tar --wildcards \
    "*/Lung.v8.signif_variant_gene_pairs.txt.gz" \
    "*/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz"
# Delete source file
rm GTEx_Analysis_v8_eQTL.tar
================================================================================


================================================================================

================================================================================



cd "${scripts_dir}"
