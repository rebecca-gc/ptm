import sys
import os
import rfc
from joblib import Parallel, delayed

sys.path.append(os.path.abspath("data_preprocess"))
sys.path.append(os.path.abspath("Source"))

import data_pipeline
import ican


dirs = ['glycosylation', 's_nitrosylation', 'acetylation', 'methylation']
#dirs = ['s_nitrosylation']

def ican_parallel(dir):
    seqs = os.path.join('data', dir, 'smiles.smi')
    output = os.path.join('data', dir)
    sys.argv = ['ican.py', f'--output_path={output}', seqs]
    ican.main()
    rfc.main(output)
    return 42

def main():
    #data_pipeline.main()
    hello = Parallel(n_jobs=-1)(delayed(ican_parallel)(dir) for dir in dirs)
    print(hello)

if __name__ == '__main__':
    main()
