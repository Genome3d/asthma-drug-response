#! /usr/bin/env python

import pandas as pd
import stringdb
import sys
import argparse

def query_string(asthma_drugs_fp, egenes_fp, output_fp):
    asthma_drugs = pd.read_csv(asthma_drugs_fp, sep='\t')
    egenes = pd.read_csv(egenes_fp, sep='\t')
    df = []
    for _, drug in asthma_drugs.iterrows():
        print(drug['drug'])
        try:
            genes = drug['gene_target'].split(';')
            string_ids = stringdb.get_string_ids(genes)
            res = stringdb.get_interaction_partners(
                string_ids['queryItem'], required_score=700)
            egene_overlap = res[res.preferredName_B.isin(egenes.gene)]
            if not egene_overlap.empty:
                egene_overlap['drug'] = drug['drug']
                df.append(egene_overlap)
        except:
            pass
    df = pd.concat(df)
    print(df)
    df.to_csv(output_fp, sep='\t', index=False)

def bootstrap_query_string(asthma_drugs_fp, egenes_fp, gene_ref_fp, output_fp):
    asthma_drugs = pd.read_csv(asthma_drugs_fp, sep='\t')
    egenes = pd.read_csv(egenes_fp, sep='\t')
    gene_ref = pd.read_csv(gene_ref_fp, sep='\t', header=None,
                           names=['chrom', 'start', 'end', 'gene', 'gencode_id'])
    df = []
    sim = 0
    num_sim = 1000
    print('Bootstrapping...')
    for i in range(num_sim):
        print(f'\tsimulation {i}')
        sim_df = []
        for _, drug in asthma_drugs.iterrows():
            #print(drug['drug'])
            try:
                num_genes = len(drug['gene_target'].split(';'))
                genes = gene_ref['gene'].sample(n=num_genes, random_state=sim).tolist()
                sim += 1
                string_ids = stringdb.get_string_ids(genes)
                res = stringdb.get_interaction_partners(
                    string_ids['queryItem'], required_score=700)
                egene_overlap = res[res.preferredName_B.isin(egenes.gene)]
                if not egene_overlap.empty:
                    egene_overlap['drug'] = drug['drug']
                    sim_df.append(egene_overlap)
            except:
                pass
        sim_df = pd.concat(sim_df)
        df.append((i, sim_df['preferredName_A'].nunique(), sim_df['preferredName_B'].nunique()))
    df = pd.DataFrame(df, columns=['sim', 'preferredName_A', 'preferredName_B'])
    print(df)
    df.to_csv(output_fp, sep='\t', index=False)


def query_string(eqtls_fp, output_fp):
    gene_df = pd.read_csv(eqtls_fp, sep='\t')
    df = []

    try:
        genes = gene_df['gene'].drop_duplicates().to_list()
        string_ids = stringdb.get_string_ids(genes)
        res = stringdb.get_interaction_partners(
            string_ids['queryItem'], required_score=700)
        print(res)
        res = gene_df.merge(res, how = 'inner', left_on = 'gene', right_on = 'preferredName_A')
        print(res)
        res.to_csv(output_fp, sep='\t', index=False)
    except:
        pass
    #df = pd.concat(df)
    #print(df)
    #df.to_csv(output_fp, sep='\t', index=False)

if __name__=='__main__':
    pd.options.mode.chained_assignment = None 

    parser = argparse.ArgumentParser(description='')
    parser.add_argument(
        '-d', '--drugs', required=True,
        help='File containing asthma drugs from Drugbank.')
    parser.add_argument(
        '-g', '--genes', required=True,
        help='File containing eGenes and the category of drugs of their eQTLs.')
    parser.add_argument(
        '-o', '--out', required=True,
        help='Directory to write results.')
    args = parser.parse_args()
    
    output_fp = os.path.join(args.out,  'string_asthma_drugs.txt')
    asthma_drugs_fp = '../results/drugbank/asthma_drugs.txt'
    egenes_fp = '../analysis/data_tidying/genes_by_drugs.txt'
    eqtls_fp = '../analysis/data_tidying/eqtls_by_drugs.txt'
    gene_ref_fp = '../data/gene_reference.bed'
    
    query_string(asthma_drugs_fp, egenes_fp, output_fp)
    output_fp = os.path.join(args.out,  'string_asthma_drugs_bootstrap.txt')
    bootstrap_query_string(asthma_drugs_fp, egenes_fp, gene_ref_fp, output_fp)

    output_fp = os.path.join(args.out,  'string_genes_by_drugs.txt')
    query_string(eqtls_fp, output_fp)
