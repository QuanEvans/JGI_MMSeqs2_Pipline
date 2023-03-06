import os
import argparse
import pandas as pd



def main(query, target, createdb, tmp, min_seq_id=0.75, threads=1, default_params = '-c 0.8 --cov-mode 1'):

    if createdb == 'y':
        query_without_fasta = query.split('.fasta')[0]
        create_db_cmd = f'mmseqs createdb {query} {query_without_fasta}'
        os.system(create_db_cmd)
        target_without_fasta = target.split('.fasta')[0]
        create_db_cmd = f'mmseqs createdb {target} {target_without_fasta}'
        os.system(create_db_cmd)
        query = query_without_fasta
        target = target_without_fasta

    search_cmd = f'mmseqs linsearch {query} {target} {query}_vs_target {tmp} \
        --min-seq-id {min_seq_id} {default_params} --threads {threads}'
    os.system(search_cmd)
    convertail_cmd = f'mmseqs convertalis {query} {target} {query}_vs_target {query}_vs_target.tsv --format-output "query"'
    os.system(convertail_cmd)
    get_filtered_index(f'{query}_vs_target.tsv', f'{query}.lookup', query)
    createsubdb_cmd = f'mmseqs createsubdb {query}_vs_target.index {query} {query}_vs_target_filtered --subdb-mode 1'
    os.system(createsubdb_cmd)
    conver2fasta_cmd = f'mmseqs convert2fasta {query}_vs_target_filtered out.fasta'
    os.system(conver2fasta_cmd)
    post_del(query,'y')

def get_filtered_index(query2filter,query_lookup,query_name):
    """
    Create a file with the index of the query sequences that are not in the query2filter file.
    """
    colunm_names = ['idx', 'query', 'id']
    query_df = pd.read_csv(query_lookup, sep='\t', names=colunm_names)
    query2filter_df = pd.read_csv(query2filter, sep='\t', names=['query_id'])
    
    all_query = set(query_df['query'].tolist())
    query2filter = set(query2filter_df['query_id'].tolist())
    filtered_query = all_query - query2filter
    filtered_query_idx = query_df[query_df['query'].isin(filtered_query)]['idx'].tolist()
    
    out_name = f'{query_name}_vs_target.index'
    with open(out_name, 'w') as f:
        for idx in filtered_query_idx:
            f.write(f'{idx}\n')


def post_del(query,createdb):
    """
    delete all input files and rename output files to original name
    """
    os.system(f'rm {query}*')
    os.system(f'mv out.fasta {query}.fasta')
    
    if createdb == 'y':
        os.system(f'mmseqs createdb {query}.fasta {query}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compare two files')
    parser.add_argument('query', type=str, help='query cluster DB')
    parser.add_argument('target', type=str, help='target cluster DB')
    parser.add_argument('createdb', type=str, help='whether to create mmseqs DBs (y/n)')
    parser.add_argument('--id', type=float, default=0.75, help='minimum sequence identity (default: 0.75)')
    parser.add_argument('--threads', type=int, default=32, help='number of threads (default: 1)')
    parser.add_argument('--tmp', type=str, default='./tmp', help='temporary directory (default: ./tmp)')
    parser.add_argument('--params', type=str, default='-c 0.8 --cov-mode 1', help='mmseqs search parameters (default: -c 0.8 --cov-mode 1)')
                        
    args = parser.parse_args()
    main(args.query, args.target, args.createdb, args.tmp, args.id, args.threads, args.params)