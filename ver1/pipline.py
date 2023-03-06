import os
import psutil
import pandas as pd


class Pipeline:

    def __init__(self, 
                 working_dir,
                 input_fasta,
                 available_cores=16, 
                 available_memory=1024 * 1024 * 1024 * 64, 
                 min_seq_id=0.75):
        
        self.working_dir = self.get_wd(working_dir)
        self.input_fasta = input_fasta
        self.min_seq_id = min_seq_id

        self.available_cores = available_cores
        self.available_memory = available_memory

        # get chunk size, usually 1/4 of available memory
        self.chunk_size = self.available_memory // 4
        self.num_chunks = None

        # initialize class path variables
        self.path_var_init()

    def path_var_init(self):
        """
        initialize class variables
        """
        # initialize class variables
        self.main_dir = f'{self.working_dir}/main'
        self.chunks_dir = f'{self.working_dir}/chunks'
        self.results_dir = f'{self.working_dir}/results'
        self.logs_dir = f'{self.working_dir}/logs'
        self.mainDB = f'{self.main_dir}/mainDB' # main database
        # splited mainDB index file, unsed for creating chunks of mainDB
        self.chunks_index = f'{self.working_dir}/chunks_index'

    def get_wd(self, wd):
        """
        get working directory
        """
        # check if working directory is absolute path
        if wd[0] == '/':
            return wd
        # if not, get absolute path
        else:
            full_path = os.path.abspath(wd)
        return full_path
    
    def file_setup(self):
        """
        create required directories
        """
        # create required directories
        os.system(f'cp {self.input_fasta} {self.main_dir}/{self.input_fasta}')
        self.input_fasta = f'{self.main_dir}/{self.input_fasta}'
        # create main database from input fasta file using mmseqs
        os.system(f'mmseqs createdb {self.input_fasta} {self.mainDB}')
        self.split_index()
        self.split_file() # split main database into chunks


    def init_setup_dir(self, delete=True):
        """
        remove all files and directories created by the pipeline
        initialize directories
        """
        if delete:
            # remove all files and directories created by the pipelines
            os.system(f'rm -rf {self.working_dir}/main')
            os.system(f'rm -rf {self.working_dir}/chunks')
            os.system(f'rm -rf {self.working_dir}/results')
            os.system(f'rm -rf {self.working_dir}/logs')
            os.system(f'rm -rf {self.working_dir}/chunks_index')

        # create required directories
        os.system(f'mkdir -p {self.working_dir}/main')
        os.system(f'mkdir -p {self.working_dir}/chunks')
        os.system(f'mkdir -p {self.working_dir}/results')
        os.system(f'mkdir -p {self.working_dir}/logs')
        os.system(f'mkdir -p {self.working_dir}/chunks_index')
    
    def split_index(self, chunk_size=None):
        if chunk_size is None:
            chunk_size = self.chunk_size
        # get file size
        file_size = os.path.getsize(f'{self.mainDB}')
        # get number of chunks
        if self.num_chunks is None:
            num_chunks = (file_size // chunk_size) + 1
        # count lines in file
        idx_file = self.mainDB + '.index'
        with open(idx_file, 'r') as f:
            num_lines = sum(1 for line in f)
        # get number of lines per chunk
        lines_per_chunk = (num_lines // num_chunks) + 1

        # split file using split command
        os.system(f'split -l {lines_per_chunk} {idx_file} {self.chunks_index}/chunk_')

    def split_file(self):
        chunck_names = os.listdir(self.chunks_index)

        for chunck in chunck_names:
            # create a directory for each chunk
            os.makedirs(f'{self.chunks_dir}/{chunck}_dir')
            # create the chunk database
            os.system(f'mmseqs createsubdb {self.chunks_index}/{chunck} {self.mainDB} {self.chunks_dir}/{chunck}_dir/{chunck}_DB --subdb-mode 1')
            # create the chunk fasta file
            os.system(f'mmseqs convert2fasta {self.chunks_dir}/{chunck}_dir/{chunck}_DB {self.chunks_dir}/{chunck}_dir/out.fasta')
            # post processing
            self.post_del(f'{self.chunks_dir}/{chunck}_dir/{chunck}_DB',f'{self.chunks_dir}/{chunck}_dir/out.fasta', createdb='y')

    def post_del(self,query,out_fasta,createdb):
        """
        delete all input files and rename output files to original name
        """
        os.system(f'rm {query}*')
        os.system(f'mv {out_fasta} {query}.fasta')

        if createdb == 'y':
            os.system(f'mmseqs createdb {query}.fasta {query}')

    def lincluster(self):
        """
        run linclust on all chunks
        """
        # get all chunks
        chunks = os.listdir(self.chunks_dir)
        for chunk in chunks:
            db_name = chunk.replace('_dir', '')
            query = f'{self.chunks_dir}/{chunk}/{db_name}_DB'
            tmp_dir = f'{self.chunks_dir}/{chunk}/tmp'
            os.system(f'python lincluster.py -q {query} -c n -id {self.min_seq_id} -t {self.available_cores} \
                      -t {self.available_cores} --tmp {tmp_dir}')
            # remove tmp directory
            os.system(f'rm -rf {tmp_dir}')
            os.system(f'echo {chunk} finished >> {self.chunks_dir}/{chunk}/finished.txt')
    
    def compare2db(self):
        chunks = os.listdir(self.chunks_dir)
        chunks_name = [chunk.replace('_dir', '') for chunk in chunks]
        num_chunks = len(chunks)
        for i in range(1,num_chunks):
            target = f'{self.chunks_dir}/{chunks[i-1]}/{chunks_name[i-1]}_DB'
            # targert need to be indexed
            os.system(f'mmseqs createlinindex {target} {self.chunks_dir}/{chunks[i-1]}/tmp')
            # remove tmp directory
            os.system(f'rm -rf {self.chunks_dir}/{chunks[i-1]}/tmp')
            for j in range(i,num_chunks):
                query = f'{self.chunks_dir}/{chunks[j]}/{chunks_name[j]}_DB'
                tmp_dir = f'{self.chunks_dir}/{chunks[j]}/tmp'
                cmd =f'python compare2db.py {query} {target} \
                          n --id {self.min_seq_id} --tmp {tmp_dir} --threads {self.available_cores}'
                os.system(f'echo {cmd} >> {self.logs_dir}/cmd.txt')
                os.system(cmd)
                # remove tmp directory
                os.system(f'rm -rf {tmp_dir}')

    def run(self):
        """
        run the pipeline
        """
        self.init_setup_dir()
        self.file_setup()
        self.lincluster()
        self.compare2db()

if __name__ == '__main__':
    wd = './wd'
    input_fasta = 'DB.fasta'
    available_cores = 30
    min_seq_id = 0.75
    # 100 GB total memory
    #total_memory = 100 * 1024 * 5
    instance = Pipeline(wd, input_fasta, available_cores=available_cores, min_seq_id=min_seq_id)
    instance.chunk_size = 1024 * 1024 * 100
    instance.run()