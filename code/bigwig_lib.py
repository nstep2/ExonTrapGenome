

bigwig_path = exp_output_path.bigwig_path+'phastcons_17way/hg38.phastCons17way.bw'



import numpy as np
import pyBigWig

bw = pyBigWig.open(bigwig_path)

score_array = bw.stats("chr1", 100000, 100050,nBins=51)
score_array = np.array(score_array)

class bigwig_cons():
    def __init__(self, file_path):
        self.bw = pyBigWig.open(file_path)
    def get_base_scores(self, chrom, start, end):
        end += 1
        bases_len = end-start+1
        score_array = self.bw.stats(chrom, start, end,nBins=bases_len)
        score_array = np.array(score_array)

        return score_array[1:]
    def check_if_none(self, score_array):
        if None in score_array:
            return True
        else:
            return False






