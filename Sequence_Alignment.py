#variables defined
crt_bst=list()#current best score: a list of list
cnt_seq=0
i=0# the ith base in sequence A
j=0# the jth base in sequence B
gap_s=0
gap_e=0
v_s=0
v_e=0
v_se=0
tmp_seq=''
tmp_l=list()
tmp_s=list()
tmp_e=list()
tmp_se=list()
dir_min=0
pnt_dit=list()
###########read data from the web page
########### set the path in the computer
import os
os.getcwd()
#change path
#os.chdir('')

# read in fasta file
seq_fasta= open('Sequence.fasta')

seq=list()
cnt_seq=0
for line in seq_fasta:
    line=line.rstrip()
    if line.startswith('>'):
        if cnt_seq==0:
            cnt_seq=cnt_seq+1
            continue
        seq.append(tmp_seq[:])
        tmp_seq=''
        continue
    tmp_seq=tmp_seq+line

seq.append(tmp_seq[:])
tmp_seq=''
#print(len(seq[0]),len(seq[1]))

########alignment forward###############
#first build a list for current best value of penalty score
#initiation of first row in current best score
i=0
while i<=len(seq[0]):
    tmp_l.append([i*5,['e']])
    i=i+1
crt_bst.append(tmp_l[:])
#calculate the folling rows
j=1
while j<=len(seq[1]):
    tmp_l.clear()#store current value in line j
    #add gap to the last row +5
    if (j == len(seq[1])):
        gap_e = 5
    else:
        gap_e = 9
    i = 0
    while i<=len(seq[0]):
        if (i==0):# for column zero, go down south and penalize 5
            tmp_l.append([j*5,['s']])
            i=i+1
            continue # colum zero over
        #other colunmns
        #add gap to the last column
        if (i==len(seq[0])):
            gap_s=5
        else:
            gap_s=9
        #compare s,e,se for colunm i and row j
        v_s=crt_bst[j - 1][i][0]+ gap_s
        v_e=tmp_l[i - 1][0] + gap_e
        tmp_s.append([v_s,'s'])
        tmp_e.append([v_e,'e'])

        if seq[0][i-1]==seq[1][j-1]:
            v_se=crt_bst[j-1][i-1][0]
            tmp_se.append([v_se,'se'])
        else:
            v_se=crt_bst[j-1][i-1][0]+6
            tmp_se.append([v_se,'se'])

        #max of 3 different direction, if 3 scores are the same
        dir_min=min(tmp_e[0][0],tmp_s[0][0],tmp_se[0][0])
        for v in (tmp_s[0], tmp_e[0],tmp_se[0]):
            if (v[0]==dir_min):
                pnt_dit.append(v[1][:])
        tmp_l.append([dir_min,pnt_dit[:]])
        tmp_e.clear()
        tmp_s.clear()
        tmp_se.clear()
        pnt_dit.clear()
        i=i+1
    crt_bst.append(tmp_l[:])
    j=j+1

#The largest score
#print(crt_bst[-1][-1][0])
#trace back
j=len(seq[1])
i=len(seq[0])
alg_seqA=list()
alg_seqB=list()
#SE is a matched or mismatched pair,
#S is sequence A with a gap
#E is sequence B with a gap
while (j>0 and i>0):
    if('se'in crt_bst[j][i][1]):
        alg_seqA.append(seq[0][i-1])
        alg_seqB.append(seq[1][j-1])
        i = i - 1
        j = j - 1
    else:#insert gaps in the short sequence
        if(len(seq[0])>=len(seq[1])):
            if('e' in crt_bst[j][i][1]):
                alg_seqA.append(seq[0][i-1])
                alg_seqB.append('-')
                i = i - 1
                j = j
            else:
                alg_seqA.append('-')
                alg_seqB.append(seq[1][j-1])
                j = j - 1
                i = i
        else:
            if ('s' in crt_bst[j][i][1]):
                alg_seqA.append('-')
                alg_seqB.append(seq[1][j-1])
                j = j - 1
                i = i
            else:
                alg_seqA.append(seq[0][i-1])
                alg_seqB.append('-')
                i = i - 1
                j = j
#trace back stops at the first row or first column
if(j==0 and i>0):
    while(i>0):
        alg_seqA.append(seq[0][i-1])
        alg_seqB.append('-')
        i=i-1
elif(j>0 and i==0):
    while(j>0):
        alg_seqA.append('-')
        alg_seqB.append(seq[1][j-1])
        j=j-1

#number of initial gaps,3'
n_ig=0
n_ig=n_ig+min(alg_seqA.index('A'),alg_seqA.index('T'),alg_seqA.index('C'),alg_seqA.index('G'))
n_ig=n_ig+min(alg_seqB.index('A'),alg_seqB.index('T'),alg_seqB.index('C'),alg_seqB.index('G'))

alg_seqA.reverse()
alg_seqB.reverse()
#initial gap,5'
n_ig=n_ig+min(alg_seqA.index('A'),alg_seqA.index('T'),alg_seqA.index('C'),alg_seqA.index('G'))
n_ig=n_ig+min(alg_seqB.index('A'),alg_seqB.index('T'),alg_seqB.index('C'),alg_seqB.index('G'))

#check score
k=0
check_scr=0
while k<len(alg_seqA):
    if(alg_seqA[k]!=alg_seqB[k]):
        check_scr=check_scr+6
        if(alg_seqA[k]=='-' or alg_seqB[k]=='-'):
            check_scr=check_scr+3
    k=k+1
#The largest score
#print('Optimal Score: ',crt_bst[-1][-1][0])
check_scr=check_scr-4*n_ig
#print('Check Score: ', check_scr)
import math
a=math.floor(len(alg_seqA)/10)
k=0
print('\nOne of the optimal alignment\n')
while(k<a):
    print('Paired sequence basepair:', 10*k+1, 'to', 10*(k+1))
    print('Sequence A',alg_seqA[(10 * k):(10 * (k+1))])
    print('Sequence B',alg_seqB[(10 * k):(10 * (k+1))])
    print('')
    k=k+1
print('Paired sequence basepair:', 10*a+1, 'to', len(alg_seqA))
print('Sequence A',alg_seqA[(10 * a):])
print('Sequence B',alg_seqB[(10 * a):])
print('\n')
#The largest score
print('Optimal Score: ',crt_bst[-1][-1][0])
