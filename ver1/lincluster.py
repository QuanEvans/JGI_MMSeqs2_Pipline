import os
import argparse


def main(query, createdb, tmp, min_seq_id=0.75, threads=1, default_params = '-c 0.8 --cov-mode 1'):
    default_params = '-c 0.8 --cov-mode 1'

    if createdb == 'y':
        query_without_fasta = query.split('.fasta')[0]
        create_db_cmd = f'mmseqs createdb {query} {query_without_fasta}'
        os.system(create_db_cmd)
        query = query_without_fasta

    lincluster_cmd = f'mmseqs linclust {query} {query}_clu {tmp} --min-seq-id {min_seq_id} {default_params} --threads {threads}'
    os.system(lincluster_cmd)
    createsubdb_cmd = f'mmseqs createsubdb {query}_clu {query} {query}_clu_rep --subdb-mode 1'
    os.system(createsubdb_cmd)
    conver2fasta_cmd = f'mmseqs convert2fasta {query}_clu_rep out.fasta'
    os.system(conver2fasta_cmd)
    post_del(query,'y')


def post_del(query,createdb):
    """
    delete all input files and rename output files to original name
    """
    os.system(f'rm {query}*')
    os.system(f'mv out.fasta {query}.fasta')
    
    if createdb == 'y':
        os.system(f'mmseqs createdb {query}.fasta {query}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--query', help='query fasta file', required=True)
    parser.add_argument('-c', '--createdb', help='create database for query', required=False, default='y')
    parser.add_argument('-id', '--min_seq_id', help='minimum sequence identity', required=False, default=0.75)
    parser.add_argument('-t', '--threads', help='number of threads', required=False, default=1)
    parser.add_argument('--tmp', help='temporary directory', required=False, default='./tmp')
    parser.add_argument('-d', '--default_params', help='default parameters', required=False, default='-c 0.8 --cov-mode 1')
    args = parser.parse_args()
    main(args.query, args.createdb, args.tmp, args.min_seq_id, args.threads, args.default_params)