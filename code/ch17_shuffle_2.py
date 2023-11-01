
from experiment_paths.experiment_paths import *




import os



import pyfaidx
input_fasta = exp_output_path.chr17_shuffle_input_folder + 'chr17_shuffle.fa'
genome_fasta = pyfaidx.Fasta(input_fasta)


import tensorflow.compat.v1 as tf

tf.config.threading.set_intra_op_parallelism_threads(6)
tf.config.threading.set_inter_op_parallelism_threads(6)



from tensorflow.compat.v1.keras.models import load_model
tf.disable_v2_behavior()



import time
import numpy as np






def one_hot_encode(seq):  #take from SpliceAI github

    map = np.asarray([[0, 0, 0, 0],
                      [1, 0, 0, 0],
                      [0, 1, 0, 0],
                      [0, 0, 1, 0],
                      [0, 0, 0, 1]])

    seq = seq.upper().replace('A', '\x01').replace('C', '\x02')
    seq = seq.replace('G', '\x03').replace('T', '\x04').replace('N', '\x00')

    return map[np.fromstring(seq, np.int8) % 5]





input_sequence = 'ACAGATATTATATATTATATGAGAGCACCCAAATACTAACATATATATCCATTCTATCTATTTATACATTTTTCAGGACGAAGAAAAGACACTGCGCTGTGTGCGAAGAACGAGATATGATATTAGATATGACCCGCGCCGCATGACGAAGAAGAGAGAACCTCATCTATCGTACGGCTAGCTGTCTAGCTCGGTAAGTTATAGATGTAGGgGGGGATGATGAGCGTGTTGATGGTAGTGCGCATGGCGCTGTACTAGATGTGGTCGTAGTGT'

context = 10000
model_paths = [exp_output_path.spliceAI_model_path+'spliceai{:}.h5'.format(x) for x in range(1,6)]

models = [load_model(path) for path in model_paths]
x = one_hot_encode('N'*(context//2) + 'ACAGATATTATATATTATATGAGAGCACCCAAATACTAACATATATATCCATTCTATCTATTTATACATTTTTCAGGACGAAGAAAAGACACTGCGCTGTGTGCGAAGAACGAGATATGATATTAGATATGACCCGCGCCGCATGACGAAGAAGAGAGAACCTCATCTATCGTACGGCTAGCTGTCTAGCTCGGTAAGTTATAGATGTAGGgGGGGATGATGAGCGTGTTGATGGTAGTGCGCATGGCGCTGTACTAGATGTGGTCGTAGTGT' + 'N'*(context//2))[None, :]
y = np.mean([models[m].predict(x) for m in range(5)], axis=0)





x_stack = [one_hot_encode('N'*(context//2) + input_sequence + 'N'*(context//2))[None, :] for ii, input_sequence in enumerate([input_sequence,input_sequence,input_sequence,input_sequence,input_sequence])]
x_stack=np.vstack(x_stack)


start=time.time()
y = np.mean([models[m].predict(x_stack, batch_size=40 ) for m in range(5)], axis=0)
end=time.time()-start
print(end)


def spliceAI_score_seq(seq_stack):
    x_stack = list()

    x_stack = [one_hot_encode('N'*(context//2) + input_sequence + 'N'*(context//2))[None, :] for ii, input_sequence in enumerate(seq_stack)]
    x_stack=np.vstack(x_stack)

        
    y_stack = np.mean([models[m].predict(x_stack, batch_size=1 ) for m in range(5)], axis=0)
    
    acceptor_prob_list = list()
    donor_prob_list    = list()
    for ii, y in enumerate(y_stack):
        acceptor_prob = y[:, 1]
        donor_prob = y[:, 2]
        acceptor_prob_list.append(acceptor_prob)
        donor_prob_list.append(donor_prob)

        
    return acceptor_prob_list, donor_prob_list










chr_17_len=len(genome_fasta['chr17_shuffle'])
import datetime

chrom_stack = list()

strand_flag='+'
file_path = exp_output_path.chr17_shuffle_folder + 'chr17_shuffle_%s_spliceAI_computations.txt' % (strand_flag)
f = open(file_path,'wt')


ii_batch = list()
context = 11000
number_fragments = int(chr_17_len/1000)
for ii in range(number_fragments+1):
    ii_batch.append(ii)
    lower = int(max(ii*1000-context/2+500, 0))
    upper = int(min((ii+1)*1000+context/2-500, chr_17_len))
    lower_pad = -1*int(min((ii*1000-context/2+500), 0))*'N'
    upper_pad = int(max((context/2-500+(ii+1)*1000-chr_17_len), 0))*'N'
    if strand_flag == '+':
        seq = lower_pad+str(genome_fasta['chr17_shuffle'][lower:upper])+upper_pad
    else:
        seq = upper_pad+str(genome_fasta['chr17_shuffle'][lower:upper].reverse.complement)+lower_pad
    
    chrom_stack.append(str(seq))

    if ii % 10000 == 0:
        #print('iteration:',ii)
        print('%d fragments and at time:  %s' % (ii,datetime.datetime.now().strftime('%H:%M:%S')))

        
    if ii % 200 == 0  and ii != 0 or ii == len(genome_fasta['chr17_shuffle'])-1:
        acceptor_prob_list, donor_prob_list = spliceAI_score_seq(chrom_stack)
        
        for kk, val in enumerate(acceptor_prob_list):
            acceptor_prob=acceptor_prob_list[kk]
            donor_prob=donor_prob_list[kk]
            
            acceptor_prob = acceptor_prob[5000:6001]
            donor_prob = donor_prob[5000:6001]
            threshold_index_3ss = acceptor_prob > 0.2
            threshold_index_5ss = donor_prob > 0.2

            if ii % 1000 == 0 and kk == 0:
                print('%d fragments and at time:  %s' % (ii,datetime.datetime.now().strftime('%H:%M:%S')))

            if strand_flag == '+':
                for jj in range(0,1000):        
                    out = '%d\t%s\t%s\n' % (jj+ii_batch[kk]*1000,str(acceptor_prob[jj]),str(donor_prob[jj]))
                    f.write(out)

            else:                
                for jj in range(999,-1,-1):        
                    out = '%d\t%s\t%s\n' % (jj+ii_batch[kk]*1000,str(acceptor_prob[jj]),str(donor_prob[jj]))
                    f.write(out)
        
        #reset
        ii_batch=list()
        chrom_stack=list()
    
print("%d/%d iterations completed" % (ii,number_fragments))

f.close()


















