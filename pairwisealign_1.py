#pairwisealign_1

import numpy as np

# naive, uses all python

################################################
def pw_align_python(seqA,seqB,gap=-1,match=1,mismatch=-1):
    MAX_LEN = 50

    UP = 1
    LEFT = 2
    DIAG = 3
    NONE = 4
    
    lenA = len(seqA)  
    lenB = len(seqB)
        
    # A is rows, B is columns
    scores = np.zeros((lenA+1,lenB+1),dtype = np.int8)
    pointers = np.zeros((lenA+1,lenB+1),dtype = np.int8)
    
    # initalize
    pointers[0, 0] = NONE
    scores[0,0] = 0 

  
    pointers[0, 1:] = LEFT
    pointers[1:,0] = UP

    # setup all gaps row and colum
    scores[0, 1:] = gap * np.arange(1, lenB + 1, dtype=np.int)
    scores[1:, 0] = gap * np.arange(1, lenA + 1, dtype=np.int)

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
                
    align_a = np.zeros(MAX_LEN,dtype = np.int8)
    align_b = np.zeros(MAX_LEN,dtype = np.int8)

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






