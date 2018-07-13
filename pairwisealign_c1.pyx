#pairwisealign_c1
# attempt at cython usage...

import numpy as np
cimport numpy as np

DTYPE = np.int8
ctypedef np.int8_t DTYPE_t



################################################
def pw_align(seqA,seqB,gapScore=-1, matchScore=1, mismatchScore=-1):
    cdef DTYPE_t MAX_LEN = 50

    cdef DTYPE_t UP = 1
    cdef DTYPE_t LEFT = 2
    cdef DTYPE_t DIAG = 3
    cdef DTYPE_t NONE = 4
    
    cdef DTYPE_t lenA = len(seqA)  
    cdef DTYPE_t lenB = len(seqB)
    
    cdef DTYPE_t gap = gapScore
    cdef DTYPE_t match = matchScore
    cdef DTYPE_t mismatch = mismatchScore

        
    # A is rows, B is columns
    scores = np.zeros((lenA+1,lenB+1),dtype = DTYPE)
    pointers = np.zeros((lenA+1,lenB+1),dtype = DTYPE)
    
    # initalize
    pointers[0, 0] = NONE
    scores[0,0] = 0 

  
    pointers[0, 1:] = LEFT
    pointers[1:,0] = UP

    # setup all gaps row and colum
    scores[0, 1:] = gap * np.arange(1, lenB + 1, dtype=DTYPE)
    scores[1:, 0] = gap * np.arange(1, lenA + 1, dtype=DTYPE)

    # do the loop through....
    for a in range(1,lenA+1):
        cA = seqA[a-1]  # get the char
        for b in range(1,lenB+1,):
            cB = seqB[b-1] # get the char
            
            if cA == cB:
                diag_score = scores[ a -1, b -1] + match
            else:
                diag_score = scores[ a -1, b -1] + mismatch
            
            up_score = scores[a-1,b] + gap
            left_score = scores[a,b-1] + gap
            
            
            # figure out scores..
            if diag_score >= up_score:
                if diag_score >= left_score:
                    scores[a,b] = diag_score
                    pointers[a,b] = DIAG
                else: # left is max
                    scores[a,b] = left_score
                    pointers[a,b] = LEFT
            else: # up score is >
                if up_score >= left_score:
                    scores[a,b] = up_score
                    pointers[a,b] = UP
                else:
                    scores[a,b] = left_score
                    pointers[a,b] = LEFT
                
    align_a = np.zeros(MAX_LEN,dtype = DTYPE)
    align_b = np.zeros(MAX_LEN,dtype = DTYPE)

    # follow the seqs.
    p = pointers[a,b]
    align_counter = 0
    while p != NONE:
        if p == DIAG:
            a -= 1
            b -= 1
            align_a[align_counter] = ord(seqA[a])
            align_b[align_counter] = ord(seqB[b])
        elif p == LEFT:
            b -= 1
            align_a[align_counter] = ord('-')
            align_b[align_counter] = ord(seqB[b])
        elif p == UP:
            a -= 1
            align_a[align_counter] = ord(seqA[a])
            align_b[align_counter] = ord('-')
        align_counter += 1
        p = pointers[a,b]
        
    align_a = align_a[0:align_counter]
    align_b = align_b[0:align_counter]

    align_a = align_a[::-1] # reverse
    align_b = align_b[::-1] 
    
    align_a_list= [chr(i) for i in align_a]
    align_b_list= [chr(i) for i in align_b]
    
    return(align_a_list,align_b_list)
            
################################################






