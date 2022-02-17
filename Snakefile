configfile: "config.yaml"

rule all:
    input: 
        #afc_unconstrained = config["results"]["unconstrained_afc"],
        unconstrained_for_afc = config["results"]["unconstrained_for_afc"]
        #annotations_snps = config["results"]["annotations_snps"],
        #annotations_genes = config["results"]["annotations_genes"],
        #population = config["results"]["population"],
        #genes = config["results"]["druggable_genes"],
        #drug_disease = config["results"]["drug_disease"],
        #asthma_drugs_ppi = config["results"]["asthma_drugs_ppi"]

rule comorbidity_analysis:
    input: config["data"]["gwas"]
    #output:
    #snps = config["analysis"]["comorbid"]
    #shell: "Rscript scripts/get_gwas_snps.R -i {input} -o {output.snps}"
           
rule annotate_snps:
    input:
        eqtls_constrained = expand("{eqtls}{tissue}/significant_eqtls.txt",
                      eqtls=config["results"]["constrained_eqtls"],
                      tissue=config["tissues"]["list"]),
        eqtls_unconstrained = expand("{eqtls}{tissue}/significant_eqtls.txt",
                      eqtls=config["results"]["unconstrained_eqtls"],
                      tissue=config["tissues"]["list"])
        #config["analysis"]["curated_snps"]
    output:
        config["results"]["annotations_snps"]

rule annotate_genes:
    input:
        eqtls_constrained = expand("{eqtls}{tissue}/significant_eqtls.txt",
                      eqtls=config["results"]["constrained_eqtls"],
                      tissue=config["tissues"]["list"]),
        eqtls_unconstrained = expand("{eqtls}{tissue}/significant_eqtls.txt",
                      eqtls=config["results"]["unconstrained_eqtls"],
                      tissue=config["tissues"]["list"])
        #gencode = config["data"]["gencode"]
    output:
        config["results"]["annotations_genes"], config["analysis"]["genes_by_drugs"]
    shell: "python scripts/get_gene_type.py -o {output}"

           
rule population_analysis:
    input:
        annotations_snps = config["results"]["annotations_snps"]
    output:
        config["results"]["population"]

        
rule analyse_drug_interactions:
    input:
        genes = config["results"]["annotations_genes"],
        drugbank = config["data"]["drugbank"],
        out_dir = config["results"]["drugbank"]
    output:
        genes = config["results"]["druggable_genes"],
        drugs = config["results"]["drugs"],
        asthma_drugs = config["results"]["asthma_drugs"],
        asthma_drug_interactions = config["results"]["asthma_drug_interactions"],
        bootstrap_genes = config["results"]["bootstrap_genes"],
        bootstrap_drugs = config["results"]["bootstrap_drugs"]

rule map_drugs_to_disease:
    input: config["results"]["drugs"]
    output: config["results"]["drug_disease"]
    shell: "python scripts/map_drugs_to_disease.py --drugs {input} -o {output}"

rule get_protein_interactions:
    input:
        asthma_drugs = config["results"]["asthma_drugs"], #results/drugbank/asthma_drugs.txt'
        genes_by_drugs = config["analysis"]["genes_by_drugs"], #analysis/data_tidying/genes_by_drugs.txt'
        out = config["results"]["string_dir"]
    output:
        asthma_drugs_ppi = config["results"]["asthma_drugs_ppi"]#results/stringdb/string_asthma_drugs.txt
    shell:
        "python scripts/string_analysis.py -d {input.asthma_drugs} -g {input.genes_by_drugs} -o {input.out}"
        
rule map_eqtls:
    input: config["analysis"]["curated_snps"]
    output:
        eqtls_constrained = expand("{eqtls}{tissue}/significant_eqtls.txt",
                      eqtls=config["results"]["constrained_eqtls"],
                      tissue=config["tissues"]["list"]),
        eqtls_unconstrained = expand("{eqtls}{tissue}/significant_eqtls.txt",
                      eqtls=config["results"]["unconstrained_eqtls"],
                      tissue=config["tissues"]["list"])

rule get_gwas_variants:
    input: config["data"]["gwas"]
    output:
        snps = config["analysis"]["gwas_snps"]
    shell: "Rscript scripts/get_gwas_snps.R -i {input} -o {output.snps}"

rule calc_afc_unconstrained:
    input: 
        unconstrained_for_afc = config["results"]["unconstrained_for_afc"]
    output:
        afc_unconstrained = config["results"]["unconstrained_afc"]
    shell: """Rscript scripts/collate_unconstrained_eqtls_for_afc.R \
        -i {input.unconstrained_for_afc} \
        -o {output.afc_unconstrained}
    """

rule prep_for_afc_unconstrained:
    input: 
        eqtls_unconstrained = expand("{eqtls}{tissue}/significant_eqtls.txt",
                      eqtls=config["results"]["unconstrained_eqtls"],
                      tissue=config["tissues"]["list"])
    output:
        unconstrained_for_afc = config["results"]["unconstrained_for_afc"]
    shell: """Rscript scripts/collate_unconstrained_eqtls_for_afc.R \
        -i {input.eqtls_unconstrained} \
        -o {output.unconstrained_for_afc}
    """

           
rule curate_all_variants:
    input:
        gwas = config["analysis"]["gwas_snps"],
        combined = config["analysis"]["combined_snps"]
    output: config["analysis"]["curated_snps"]
    shell: "Rscript scripts/curate_all_variants.R -i {input} -o {output}"
     	

rule download_data:
    output:
    	gtex = expand("{gtex_data}{tissue}.v8.signif_variant_gene_pairs.txt.gz", \
	    gtex_data=config["data"]["GTEx"], tissue=config["tissues"]["list"]), 
	gwas = config["data"]["gwas"],
        expression = config["data"]["expression"],
        gencode = config["data"]["gencode"]
    shell: "bash scripts/download_data.sh"
