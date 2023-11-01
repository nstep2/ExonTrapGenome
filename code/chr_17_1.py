'''
hdf5_format.py
.decode('utf-8')) replace with )#.decode('utf-8')) or remove

'''



from experiment_paths.experiment_paths import *


strand_flag = '+'


import tensorflow.compat.v1 as tf

tf.config.threading.set_intra_op_parallelism_threads(6)
tf.config.threading.set_inter_op_parallelism_threads(6)

from tensorflow.compat.v1.keras.models import load_model
tf.disable_v2_behavior()

training = tf.placeholder(tf.bool, name='training')

import time
import numpy as np


def one_hot_encode(seq):  

    map = np.asarray([[0, 0, 0, 0],
                      [1, 0, 0, 0],
                      [0, 1, 0, 0],
                      [0, 0, 1, 0],
                      [0, 0, 0, 1]])

    seq = seq.upper().replace('A', '\x01').replace('C', '\x02')
    seq = seq.replace('G', '\x03').replace('T', '\x04').replace('N', '\x00')

    return map[np.fromstring(seq, np.int8) % 5]








context = 10000

model_paths = [exp_output_path.spliceAI_model_path+'spliceai{:}.h5'.format(x)for x in range(1,6)]
models = [load_model(path) for path in model_paths]







import random

num_to_base={0:'A',1:'C',2:'G',3:'T'}

def make_random_sequence(length):
    seq=''
    for ii in range(length):
        b = random.randint(0, 256)%4
        seq = '{:}{:}'.format(seq, num_to_base[b])
    
    return seq





def spliceAI_score_seq(input_sequence):
    x = one_hot_encode('N'*(context//2) + input_sequence + 'N'*(context//2))[None, :]
    y = np.mean([models[m].predict(x) for m in range(5)], axis=0)
    acceptor_prob = y[0, :, 1]
    donor_prob = y[0, :, 2]
    
    acceptor_prob = y[0, :, 1]
    donor_prob = y[0, :, 2]
    
    return acceptor_prob, donor_prob




def spliceAI_score_seq_vstack(input_vstack):
    batch_size=1
    
    test_vstack = list()
    for ii, input_sequence in enumerate(input_vstack):
        x = one_hot_encode(input_sequence)[None, :]
        test_vstack.append(x)
    
    test_vstack=np.vstack(test_vstack)
    y = np.mean([models[m].predict(test_vstack, batch_size=batch_size) for m in range(5)], axis=0)
    
    
    acceptor_prob_list = list()
    donor_prob_list = list()
    for ii, seq_result in enumerate(y):
        acceptor_prob_list.append(seq_result[:][:,1])
        donor_prob_list.append(seq_result[:][:,2])
    
    
    return acceptor_prob_list, donor_prob_list, y





nn = 124
xx=1000
seq = genome_fasta['chr17'][xx*nn:xx*(nn+1)]
acceptor_prob, donor_prob = spliceAI_score_seq(str(seq))

import time
start = time.time()
for ii in range(10):
    acceptor_prob, donor_prob = spliceAI_score_seq(str(seq))
end = time.time()-start
print(end)








file_path = exp_output_path.SpliceAI_ET_exons + 'chr17_%s_spliceAI_computations_tmp.txt' % (strand_flag)

f = open(file_path,'wt')







import datetime

strand_flag = '+'

chr_17_len = len(genome_fasta['chr17'])

context = 11000
number_fragments = int(chr_17_len/1000)

query_vstack = list()
ii_vstack = list()

for ii in range(number_fragments+1):
    
    lower = int(max(ii*1000-context/2+500, 0))
    upper = int(min((ii+1)*1000+context/2-500, chr_17_len))
    lower_pad = -1*int(min((ii*1000-context/2+500), 0))*'N'
    upper_pad = int(max((context/2-500+(ii+1)*1000-chr_17_len), 0))*'N'
    if strand_flag == '+':
        seq = lower_pad+str(genome_fasta['chr17'][lower:upper])+upper_pad
    else:
        seq = upper_pad+str(genome_fasta['chr17'][lower:upper].reverse.complement)+lower_pad
    
    
    query_vstack.append(str(seq))
    ii_vstack.append(ii)
    
    if ii % 1 == 0  or ii == number_fragments:
        acceptor_prob_list, donor_prob_list, y = spliceAI_score_seq_vstack(query_vstack)
        
        for kk, val in enumerate(acceptor_prob_list):
            acceptor_prob, donor_prob = acceptor_prob_list[kk], donor_prob_list[kk]
            cur_ii = ii_vstack[kk]
            
            
            if cur_ii % 10000 == 0:
                print('%d fragments and at time:  %s' % (cur_ii,datetime.datetime.now().strftime('%H:%M:%S')))

            if strand_flag == '+':
                for jj in range(0,1000):        
                    out = '%d\t%s\t%s\n' % (jj+cur_ii*1000,str(acceptor_prob[jj]),str(donor_prob[jj]))
                    f.write(out)

            else:                
                for jj in range(999,-1,-1):        
                    out = '%d\t%s\t%s\n' % (jj+cur_ii*1000,str(acceptor_prob[jj]),str(donor_prob[jj]))
                    f.write(out)

        


        query_vstack=list()
        ii_vstack=list()

print("%d/%d iterations completed" % (ii,number_fragments))

 



f.close()









