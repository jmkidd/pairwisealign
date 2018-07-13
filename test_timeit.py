import sys
import nwalign as nw
import random
from align import aligner
import pairwisealign_1
import pairwisealign_c1
import pairwisealign_c2
import pairwisealign_c3
import pairwisealign_c4
import pairwisealign_c5


# setup matrix for using alignment...
MY_MATRIX = {}

MY_MATRIX['A'] = {'A':1, 'C':-1,'T':-1,'G':-1,'N':-1 }
MY_MATRIX['C'] = {'A':-1, 'C':1,'T':-1,'G':-1,'N':-1 }
MY_MATRIX['T'] = {'A':-1, 'C':-1,'T':1,'G':-1,'N':-1 }
MY_MATRIX['G'] = {'A':-1, 'C':-1,'T':-1,'G':1,'N':-1 }
MY_MATRIX['N'] = {'A':-1, 'C':-1,'T':-1,'G':-1,'N':1 }


#####################################################################
def random_seqs(seqLen,numPairs):
    data = []
    for i in range(numPairs):
        s1 = ''
        s2 = ''
        for j in range(seqLen):
            c = random.choice(['A','C','T','G','-'])
            if c != '-':
                s1 += c
            c = random.choice(['A','C','T','G','-'])
            if c != '-':
                s2 += c
        data.append([s1,s2])
    return data
#####################################################################
        
        



#####################################################################
def run_align_nw(data):
    for s in data:
        nw.global_align(s[0],s[1])    
#####################################################################    
def run_align_aligner(data):
    for s in data:
        aligner(s[0],s[1],method='global',gap_open=-1,gap_extend=-1,matrix=MY_MATRIX)        
#####################################################################           


#####################################################################    
def run_align_pw_python(data):
    for s in data:
        pairwisealign_1.pw_align_python(s[0],s[1])
#####################################################################           
def run_align_pw_c1(data):
    for s in data:
        pairwisealign_c1.pw_align(s[0],s[1])
#####################################################################           
def run_align_pw_c2(data):
    for s in data:
        pairwisealign_c2.pw_align(s[0],s[1])
#####################################################################           
def run_align_pw_c3(data):
    for s in data:
        pairwisealign_c3.pw_align(s[0],s[1])
#####################################################################           
#####################################################################           
def run_align_pw_c4(data):
    for s in data:
        pairwisealign_c4.pw_align(s[0],s[1])
#####################################################################           
def run_align_pw_c5(data):
    for s in data:
        pairwisealign_c5.pw_align(s[0],s[1])
#####################################################################           





