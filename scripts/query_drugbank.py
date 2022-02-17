#! /usr/bin/env python

import pandas as pd
import xml.etree.ElementTree as et
import os
import sys
import re
import istarmap
import multiprocessing as mp
from tqdm import tqdm
from itertools import repeat
import numpy as np
import argparse


def get_pubmed_ids(drug):
    pubmed_ids = []
    for article in [ids for ids in drug.find(
            'drugbank:general-references', ns).find(
                'drugbank:articles', ns).findall(
                    'drugbank:article', ns)]:
        pubmed_id = article.find('drugbank:pubmed-id', ns).text
        pubmed_ids.append(pubmed_id if pubmed_id else '')
    return '; '.join(pubmed_ids) if len(pubmed_ids) > 0 else ''
    # return pubmed_ids


def get_description(drug):
    description = drug.find('drugbank:description', ns).text
    return description.replace(
        '\n', '').replace('\r', '') if description else ''


def get_indication(drug):
    indication = drug.find('drugbank:indication', ns).text
    return indication.replace(
        '\n', '').replace('\r', '') if indication else ''


def get_pharmacodynamics(drug):
    pharmocodynamics = drug.find('drugbank:pharmacodynamics', ns).text
    return pharmocodynamics.replace(
        '\n', '').replace('\r', '') if pharmocodynamics else ''


def get_mechanism_of_action(drug):
    mechanism = drug.find('drugbank:mechanism-of-action', ns).text
    return mechanism.replace(
        '\n', '').replace('\r', '') if mechanism else ''


def get_products(drug):
    products = []
    for product in [product for product in drug.find(
            'drugbank:products', ns).findall(
                'drugbank:product', ns)]:
        product_name = product.find('drugbank:name', ns).text
        products.append(product_name if product_name else '')
    return '; p'.join(products) if len(products) > 0 else ''


def get_drug_interactions(drug):
    interactions = []
    descriptions = []
    for interaction in [interaction for interaction in drug.find(
            'drugbank:drug-interactions', ns).findall(
                'drugbank:drug-interaction', ns)]:
        interactions.append(interaction.find('drugbank:name', ns).text)
        descriptions.append(interaction.find('drugbank:description', ns).text)
    return interactions, descriptions


def get_food_interactions(drug):
    interactions = [interaction.text for interaction in drug.find(
        'drugbank:food-interactions', ns).findall(
        'drugbank:food-interaction', ns)]
    return '; '.join(interactions) if len(interactions) > 0 else ''

def get_identifiers(drug):
    ttd_id = ''
    kegg_drug_id = ''
    pubchem_sid = ''
    pharmgkb_id = ''
    ids = {}
    identifier_root = drug.find(
        'drugbank:external-identifiers', ns)
    if identifier_root == None:
        return ids
    identifiers = identifier_root.findall(
        'drugbank:external-identifier', ns)
    for identifier in identifiers:
        #drug.find(
        #'drugbank:external-identifiers', ns).findall(
        #'drugbank:external-identifier', ns):
        resource = identifier.find('drugbank:resource', ns).text
        ids[resource] = identifier.find('drugbank:identifier', ns).text
        '''
        if resource == 'Therapeutic Targets Database':
            ttd_id = identifier.find('drugbank:identifier', ns).text
        if resource == 'KEGG Drug':
            kegg_drug_id = identifier.find('drugbank:identifier', ns).text
        if resource == 'PubChem Substance':
            pubchem_id = identifier.find('drugbank:identifier', ns).text
        if resource == 'PharmGKB':
            pharmgkb_id = identifier.find('drugbank:identifier', ns).text
        '''

    #return ttd_id, kegg_drug_id, pubchem_id, pharmgkb_id
    return ids


def get_atc_code(drug):
    try:
        return drug.find('drugbank:atc-codes', ns).find(
            'drugbank:atc-code', ns).get('code')
    except:
        return ''

def get_inchikey(drug, ns):
    try:
        inchikeys = []
        salts = drug.find('drugbank:salts', ns)
        for salt in salts.findall('drugbank:salt', ns):
            inchikey = salt.find('drugbank:inchikey', ns).text
            if inchikey != '':
                inchikeys.append(inchikey)
        if len(inchikeys) > 0:
            return ', '.join(inchikeys)
        else:
            return ''
    except:
        return ''

    
def get_action_on_target(target):
    return [action.text for action in target.find('drugbank:actions', ns).findall(
        'drugbank:action', ns)]


def get_target_polypeptide(target):
    name = ''
    gene = ''
    polypeptide = target.find('drugbank:polypeptide', ns)
    if polypeptide:
        name = polypeptide.find('drugbank:name', ns).text
        gene = polypeptide.find('drugbank:gene-name', ns).text
    return(name, gene)


def get_targets(drug):
    polypeptides = []
    genes = []
    actions = []
    for target in [target for target in drug.find(
            'drugbank:targets', ns).findall(
                'drugbank:target', ns)]:
        t_actions = ', '.join(get_action_on_target(target))
        actions.append(t_actions)
        polypeptide, gene = get_target_polypeptide(target)

        polypeptides.append(polypeptide if polypeptide else '')
        genes.append(gene if gene else '')
    if len(actions) > 0:
        actions = '; '.join(actions)
    else:
        actions = ''
    if len(polypeptides) > 0:
        polypeptides = '; '.join(polypeptides)
    else:
        polypeptides = ''
    if len(genes) > 0:
        genes = '; '.join(genes)
    else:
        genes = ''
    return(polypeptides, genes, actions)


def get_xroot(drug_bank_fp):
    print('Parsing Drugbank database...')
    xtree = et.parse(drug_bank_fp)
    return xtree.getroot()


def get_drug_by_name(name, xroot, ns):
    for drug in xroot:
        if name == drug.find('drugbank:name', ns).text:
            return drug


def get_drugs(gene_list, xroot, ns, output_dir):
    drug_df = []
    genes = {}
    interactions_df = []
    ttd = 'Therapeutic Targets Database'
    kegg_drug = 'KEGG Drug'
    pubchem_s = 'PubChem Substance'
    pubchem_c = 'PubChem Compound'
    pharmgkb = 'PharmGKB'
    chebi = 'ChEBI'
    chembl = 'ChEMBL'

    for drug in xroot:
        name = drug.find('drugbank:name', ns).text
        if name is None:
            continue
        target_polypeptides, target_genes, target_actions = get_targets(
            drug)
        gene_inter = set(gene_list).intersection(
            set([gene.strip() for gene in target_genes.split(';')]))

        if len(gene_inter) < 1:
            continue
        print(f'\t\t{name}')
        for gene in gene_inter:
            if gene not in genes.keys():
                genes[gene] = []
            genes[gene].append(name)
        #ttd_id, kegg_drug_id, pubchem_id, pharmgkb_id = get_identifiers(drug)
        ids = get_identifiers(drug)
        atc_code = get_atc_code(drug)
        inchikey = get_inchikey(drug, ns)
        description = get_description(drug)
        pubmed_ids = get_pubmed_ids(drug)
        indication = get_indication(drug)
        pharmacodynamics = get_pharmacodynamics(drug)
        mechanism_of_action = get_mechanism_of_action(drug)
        products = get_products(drug)
        food_interactions = get_food_interactions(drug)
        drug_interactions, drug_interaction_descriptions = get_drug_interactions(
            drug)
        num_drug_interactions = len(drug_interactions)
        drug_df.append(
            [name,
             ids[ttd] if ttd in ids.keys() else '',
             ids[kegg_drug] if kegg_drug in ids.keys() else '',
             ids[pubchem_s] if pubchem_s in ids.keys() else '',
             ids[pubchem_c] if pubchem_c in ids.keys() else '',
             ids[pharmgkb] if pharmgkb in ids.keys() else '',
             ids[chebi] if chebi in ids.keys() else '',
             ids[chembl] if chembl in ids.keys() else '',
             #ttd_id, kegg_drug_id, pubchem_id, pharmgkb_id,
             atc_code, inchikey, description, pubmed_ids, indication,
             pharmacodynamics, mechanism_of_action, products, num_drug_interactions,
             food_interactions, target_polypeptides, target_genes, target_actions])
        for i in range(len(drug_interactions)):
            i_name = get_drug_by_name(drug_interactions[i], xroot, ns)
            if i_name == None:
                continue
            i_ids = get_identifiers(i_name)
            interactions_df.append(
                [name, drug_interactions[i],
                 i_ids[ttd] if ttd in i_ids.keys() else '',
                 i_ids[kegg_drug] if kegg_drug in i_ids.keys() else '',
                 i_ids[pubchem_s] if pubchem_s in i_ids.keys() else '',
                 i_ids[pubchem_c] if pubchem_c in i_ids.keys() else '',
                 i_ids[pharmgkb] if pharmgkb in i_ids.keys() else '',
                 i_ids[chebi] if chebi in i_ids.keys() else '',
                 i_ids[chembl] if chembl in i_ids.keys() else '',
                 drug_interaction_descriptions[i]])
        
    drug_cols = ['drug',
                 'ttd_id', 'kegg_drug_id', 'pubchem_sid', 'pubchem_cid', 'pharmgkb_id', 'chebi_id', 'chembl_id',
                 'atc_code', 'inchikey', 'description', 'pubmed_ids',
                 'indication', 'pharmacodynamics', 'mechanism_of_action', 'products', '#drug_interactions',
                 'food_interactions', 'target_polypeptides', 'target_genes', 'target_actions']

    drug_df = pd.DataFrame(drug_df, columns=drug_cols)
    drug_df = drug_df.replace({r'\\n': ''}, regex=True)
    drug_df.to_csv(os.path.join(output_dir, 'drugs.txt'), sep='\t', index=False)
    #print(drug_df)

    interactions_cols = ['drug', 'drug_interactions', 'ttd_id', 'kegg_drug_id',
                         'pubchem_sid', 'pubchem_cid', 'pharmgkb_id', 'chebi_id',
                         'chembl_id', 'drug_interaction_descriptions']
    interactions_df = pd.DataFrame(interactions_df, columns=interactions_cols)
    interactions_df.to_csv(os.path.join(output_dir, 'drug_interactions.txt'), sep='\t', index=False)
    #print(interactions_df)

    genes_out_df = []
    for gene in genes:
        for drug in genes[gene]:
            genes_out_df.append([gene, drug])
    genes_out_df = pd.DataFrame(genes_out_df, columns=['gene', 'drug'])
    #print(genes_out_df)
    genes_out_df.to_csv(os.path.join(output_dir, 'genes.txt'), sep='\t', index=False)


def get_asthma_drugs(xroot, ns, output_dir):
    drug_df = []
    interactions_df = []
    regex = re.compile('asthma')
    for drug in xroot:
        indication = get_indication(drug)
        if re.search('asthma', indication, re.IGNORECASE):
            name = drug.find('drugbank:name', ns).text
            drug_interactions, drug_interaction_descriptions = get_drug_interactions(
                drug)
            print(f'\t\t{name}')
            target_polypeptides, target_genes, target_actions = get_targets(
                drug)
            drug_df.append((name, target_genes, indication, len(drug_interactions)))
            for i, inter in enumerate(drug_interactions):
                inter_drug = get_drug_by_name(inter, xroot, ns)
                atc_code = get_atc_code(inter_drug)
                inchikey = get_inchikey(inter_drug, ns)

                interactions_df.append((name, inter, atc_code, inchikey, drug_interaction_descriptions[i]))
    drug_cols = ['drug', 'gene_target', 'indication', 'num_interactions']
    inter_cols = ['drug', 'drug_interaction', 'atc_code', 'inchikey', 'interaction_desc']
    drug_df = pd.DataFrame(drug_df, columns=drug_cols)
    interactions_df = pd.DataFrame(interactions_df, columns=inter_cols)
    drug_df.to_csv(os.path.join(output_dir, 'asthma_drugs.txt'), sep='\t', index=False)
    interactions_df.to_csv(os.path.join(output_dir, 'asthma_drug_interactions.txt'), sep='\t', index=False)


def bootstrap_genes(gene_list, xroot, ns, gene_ref_fp, output_dir):
    gene_ref = pd.read_csv(gene_ref_fp, sep='\t', header=None,
                           names=['chrom', 'start', 'end', 'gene', 'gencode_id'])
    num_sim = 10000
    sim_dir = os.path.join(output_dir, 'simulations')
    manager = mp.Manager()
    num_proc = 2
    desc = 'Simulating set...'
    bar_format = '{desc}: {percentage:3.0f}% |{bar}| {n_fmt}/{total_fmt} {unit}'
    print(desc)
    for i in range(6999, 10000):
        print(i)
        sim(i,
            gene_ref,
            sim_dir,
            xroot,
            ns)
    '''
    args = zip(range(num_sim),
               repeat(gene_ref),
               repeat(sim_dir),
               repeat(xroot),
               repeat(ns))
    with mp.Pool(num_proc) as pool:
        for _ in tqdm(
                pool.istarmap(sim, args),
                total=num_sim,
                desc=desc,
                bar_format=bar_format,
                unit='sim',
                ncols=80
        ):
            pass
    '''
        
def sim(i, gene_ref, sim_dir, xroot, ns):#    for i in range(num_sim):
    #print(f'\t{i}')
    i_dir = os.path.join(sim_dir, str(i))
    os.makedirs(i_dir,exist_ok=True)
    if os.path.isfile(os.path.join(i_dir, 'genes.txt')):
        return
    rand_gene_list = gene_ref['gene'].sample(n=len(gene_list), random_state=i).tolist()
    get_drugs(rand_gene_list, xroot, ns, i_dir)
        

def collate_sims(output_dir):
    eqtls_fp = '../analysis/data_tidying/eqtls_by_drugs.txt'
    eqtls = pd.read_csv(eqtls_fp, sep='\t')
    asthma_inter_fp = '../results/drugbank/asthma_drug_interactions.txt'
    asthma_inter = pd.read_csv(asthma_inter_fp, sep='\t')
    asthma_drugs_fp = '../results/drugbank/asthma_drug_interactions.txt'
    asthma_drugs = pd.read_csv(asthma_drugs_fp, sep='\t')
    sim_dir = os.path.join(output_dir, 'simulations')
    df = []
    print('Collating simulations...')
    for i in range(1000):
        i_dir = os.path.join(sim_dir, str(i))
        gene_df = pd.read_csv(os.path.join(i_dir, 'genes.txt'), sep='\t')
        gene_df = gene_df.rename(columns={'drug': 'inter_drug'})
        gene_df = gene_df.merge(asthma_inter[['drug', 'drug_interaction']], how='inner',
                                left_on='inter_drug', right_on='drug_interaction')
        df.append((i, gene_df['gene'].nunique(), gene_df['drug'].nunique(), gene_df['inter_drug'].nunique()))
        #df.append((i, gene_df['gene'].nunique(), gene_df['drug'].nunique()))
    df = pd.DataFrame(df, columns=['sim', 'genes', 'drugs', 'inter_drugs'])
    df.to_csv(os.path.join(output_dir, 'bootstrap_summary.txt'), sep='\t', index=False)
    print('Collation completed.')


def bootstrap_drugs(xroot, ns, output_dir):
    asthma_drugs_fp = '../results/drugbank/asthma_drug_interactions.txt'
    study_drugs_fp = '../results/drugbank/genes.txt'
    drugs_df =[]
    for drug in xroot:
        drugs_df.append((drug.find('drugbank:name', ns).text))
    drugs_df = pd.DataFrame(drugs_df, columns = ['drug']).drop_duplicates()
    print(f'Total number of drugs in database: {len(drugs_df)}')
    asthma_drugs = pd.read_csv(asthma_drugs_fp, sep='\t')
    study_drugs = pd.read_csv(study_drugs_fp, sep='\t')

    num_drugs = len(study_drugs['drug'].drop_duplicates())
    print(f'Number of drugs in study: {num_drugs}')
    num_sim = 10000
    df = []
    for i in tqdm(range(num_sim),
                  desc='Boostrapping drugs',
                  unit='simulations',
                  ncols=80):
        rand_drug_list = drugs_df['drug'].sample(n=num_drugs, random_state=i).tolist()
        overlap = asthma_drugs[asthma_drugs['drug_interaction'].isin(rand_drug_list)]
        df.append((i, overlap['drug'].nunique(), overlap['drug_interaction'].nunique()))
    df = pd.DataFrame(df, columns=['sim', 'asthma_drug', 'inter_drug'])
    print(df)
    df.to_csv(os.path.join(output_dir, 'bootstrap_summary_drugs.txt'), sep='\t', index=False)
    sys.exit()


def parse_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-g', '--genes', required=True,
                        help='File containing eGenes.')
    parser.add_argument('-d', '--drugbank', required=True,
                        help='Drugbank interactions in XML format.')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='Directory to write results.')
    
        
if __name__ == '__main__':
    output_dir = args.output_dir #'../results/drugbank'
    gene_fp = args.genes #'../analysis/data_tidying/genes_by_drugs.txt'
    gene_ref_fp = '../data/gene_reference.bed'
    drug_bank_fp = args.drugbank #'../data/full database.xml'
    #drug_bank_fp = '../../data/test_cyclosporine.xml'

    gene_df = pd.read_csv(gene_fp, sep='\t')
   # global gene_list
    gene_list = gene_df['gene'].tolist()
    #gene_list = ['PPP3R2', 'CAMLG']

    xtree = et.parse(drug_bank_fp)
    xroot = get_xroot(drug_bank_fp)
   # global ns
    ns = {
        'drugbank': 'http://www.drugbank.ca',
        'description': '{http://www.drugbank.ca}description',
    }
    
    get_drugs(gene_list, xroot, ns, output_dir)
    get_asthma_drugs(xroot, ns, output_dir)
    bootstrap_genes(gene_list, xroot, ns, gene_ref_fp, output_dir)
    collate_sims(output_dir)
    bootstrap_drugs(xroot, ns, output_dir)

