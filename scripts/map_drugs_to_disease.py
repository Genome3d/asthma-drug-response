#!/usr/bin/env python
import pandas as pd
import sys
import re
import argparse

def get_ttd_id(drug_df, drug_col, ttd_df):
    ids = {}
    for idx, row in drug_df.iterrows():
        chebi = ttd_df[ttd_df['value'] == 'ChEBI:{}'.format(row['chebi_id'])]
        chembl = ttd_df[ttd_df['value'] == row['chembl_id']]
        name = ttd_df[ttd_df['value'].str.contains(re.escape(row[drug_col]), case=False, regex=True)]
        if row[drug_col] not in ids.keys():
            ids[row[drug_col]] = []
        if not chebi.empty:
            ids[row[drug_col]] += chebi['ttd_id'].tolist()
        elif not chembl.empty:
            ids[row[drug_col]] += chembl['ttd_id'].tolist()
        elif not name.empty:
            ids[row[drug_col]] += name['ttd_id'].tolist()
    df = []
    for drug in ids:
        for ttd_id in ids[drug]:
            df.append([drug, ttd_id])
    df = pd.DataFrame(df, columns=[drug_col, 'ttd_id'])
    return df
    

def get_icd11_codes(ttd_ids, icd11_df):
    ids = {}
    df = icd11_df.loc[(icd11_df['ttd_id'].isin(ttd_ids['ttd_id'])) &
                      (icd11_df['variable'] == 'INDICATI')]
    replacements = {
        'value': {r'[]': '|'}}
    df['value'] = df['value'].apply(lambda x: x.replace('[', '|').replace(']', '|'))
    df[['disease', 'icd11', 'status']] = df['value'].str.split('|', expand=True)
    df['icd11'] = df['icd11'].apply(lambda x: x.replace('ICD-11:', '').strip())
    df['disease'] = df['disease'].str.strip()
    df = df.drop(['variable', 'value'], axis=1)
    return df.merge(ttd_ids, how='left', on=['ttd_id']).drop_duplicates()


def read_drugs(drug_fp):
    return pd.read_csv(drug_fp, sep='\t')


def read_ttd():
    ttd_fp = '../data/ttd/P1-03-TTD_crossmatching.txt'
    ttd_df = pd.read_csv(ttd_fp, skiprows=27, header=None, sep='\t')
    ttd_df.columns = ['ttd_id', 'variable', 'value']
    return ttd_df


def read_icd11():
    icd11_fp = '../data/ttd/P1-05-Drug_disease.txt'
    icd11_df = pd.read_csv(icd11_fp, skiprows=21, header=None, sep='\t')
    icd11_df.columns = ['ttd_id', 'variable', 'value']
    return icd11_df


def map_disease(drug_fp, drug_col, output_fp, ttd_df, icd11_df):
    print(f'Processing {drug_fp}')
    drug_df = read_drugs(drug_fp)
    print(drug_df)

    ttd_ids = get_ttd_id(drug_df, drug_col, ttd_df)
    print(ttd_ids)
    
    disease_df = get_icd11_codes(ttd_ids, icd11_df)
    print(disease_df)
    disease_df.to_csv(output_fp, sep='\t', index=False)


def map_promiscuous(drug_fp, drug_col):
    drug_info_fp = '../data/promiscuous/drug_information.csv'
    icd_drugs_fp = '../data/promiscuous/icd_drugs.csv'
    atc_class_fp = '../data/promiscuous/atc_classification.csv'
    output_fp = '../results/drugbank/drug_disease_promiscuous.txt'
    drug_df = read_drugs(drug_fp)

    atc_class_df = pd.read_csv(atc_class_fp, sep=',')
    atc_res = atc_class_df.merge(drug_df, how='inner', left_on=['level5'], right_on=['atc_code'])
    df = atc_res[['pid', drug_col]]
    drug_info_df = pd.read_csv(drug_info_fp, sep='|', usecols=['pid', 'name', 'inchikey'])
    name_res = drug_info_df.merge(drug_df, how='inner', left_on=['name'], right_on=[drug_col])
    df = df.append(name_res[['pid', drug_col]], ignore_index=True)
    inchikey_res = drug_info_df.merge(drug_df[drug_df['inchikey'].notna()], how='inner', on=['inchikey'])
    df = df.append(inchikey_res[['pid', drug_col]], ignore_index=True)

    icd_df = pd.read_csv(icd_drugs_fp, sep=',')
    df = df.merge(icd_df[['pid', 'indication', 'icd10']], how='inner', on=['pid'])
    df = df[df['icd10'].notna()]
    df['uid'] = df.index

    codes = pd.DataFrame(df['icd10'].str.split(',').tolist(), index=df['uid']).stack()
    codes = codes.reset_index()[[0, 'uid']]
    codes.columns = ['icd', 'uid']
    codes = pd.DataFrame(codes['icd'].str.split('-').tolist(), index=codes['uid']).stack()
    codes = codes.reset_index()[[0, 'uid']]#.drop_duplicates()
    codes.columns = ['icd', 'uid']
    codes['icd'] = codes['icd'].str.strip().replace({r'\.': ''}, regex=True)

    codes = codes.merge(df, on=['uid'], how='left')
    df = codes[['drug', 'icd', 'indication']].drop_duplicates()
    df.to_csv(output_fp, sep='\t', index=False)

    
if __name__=='__main__':
    pd.set_option('mode.chained_assignment', None)
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-d', '--drugs', required=True,
        help='File of drugs that target eGenes.')
    parser.add_argument('-o', '--out', required=True,
        help='Filepath to write results.')
    args = parser.parse_args()
    '''
    ttd_df = read_ttd()
    icd11_df = read_icd11()
    '''
    
    # Drugs
    drug_fp = args.drugs #'../results/drugbank/drugs.txt'
    output_fp = args.out #'../results/drugbank/drug_disease.txt'
    drug_col = 'drug'
    '''
    ttd_df = read_ttd()
    icd11_df = read_icd11()
    '''
    
    # Drugs
    drug_fp = '../results/drugbank/drugs.txt'
    drug_col = 'drug'
    output_fp = '../results/drugbank/drug_disease.txt'
    map_promiscuous(drug_fp, drug_col)

    '''
    #map_disease(drug_fp, drug_col, output_fp, ttd_df, icd11_df)
    
    # Drug interactions
    drug_fp = '../results/drugbank/drug_interactions.txt'
    drug_col = 'drug_interactions'
    output_fp = '../results/drugbank/drug_interaction_disease.txt'
    map_disease(drug_fp, drug_col, output_fp, ttd_df, icd11_df)
    '''
