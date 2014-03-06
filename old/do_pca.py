# pypy script <in_fasta> <out_fasta>

import sys
from libs.bio_file_parsers import fasta_parser
import cPickle
import numpy as np

def main(args):
    
    print 'Loading pickle...'
    
    #~ import cPickle
    #~ matrix = cPickle.load(open('results/kmer_obj.pickle', 'r'))
    
    kmer_array = np.load('results/kmer_obj.pickle')
    
    print 'Loading matrix...'
    #~ kmer_array = np.array(matrix)
    kmer_array = np.transpose(kmer_array)
    #~ del matrix
    
    #~ kmer_array = kmer_array[:1000]
    print kmer_array.shape
    
    print 'PCA...'
    import time
    #~ time.sleep(200)
    
    # Do PCA using scikit
    import pylab as pl
    from mpl_toolkits.mplot3d import Axes3D
    from sklearn import decomposition
    from sklearn.preprocessing import scale
    
    fig = pl.figure(1, figsize=(4, 3))
    pl.clf()
    ax = Axes3D(fig)
    
    kmer_array = scale(kmer_array, axis=0)
    
    pl.cla()
    pca = decomposition.PCA(n_components=400)
    pca.fit(kmer_array)
    
    print pca.explained_variance_ratio_
    cum = np.cumsum(pca.explained_variance_ratio_)
    print cum
    
    X = pca.transform(kmer_array)
    print X.shape
    
    print 'Plotting...'
    ax.scatter(X[:, 0], X[:, 1], X[:, 2], cmap=pl.cm.spectral)
    pl.show()
    
    #~ fig = pl.figure()
    #~ pl.plot(range(1, len(cum) + 1), cum)
    #~ pl.show()

def read_phenotype(pheno_file):
    """ Reads pheno file into df.
    """
    pheno_df = pd.read_csv(pheno_file,
                           sep=',',
                           index_col=None,
                           header=0,
                           na_values=['NA', 0])
    # Delete the first column
    del pheno_df['row']
    # Use wlid as index then delete it
    pheno_df.index = pheno_df['w1id']
    del pheno_df['w1id']

    return pheno_df
    

if __name__ == '__main__':
    
    args = {}
    #~ args['in_fasta'] = sys.argv[1]
    
    main(args)
