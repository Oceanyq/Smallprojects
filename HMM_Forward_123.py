#Set working path
import os
import math
cwd=os.getcwd()
#print(cwd)

#read in data from txt
seq_txt=open("HW2_seq.txt")

seq="";
for line in seq_txt:
    line = line.rstrip()
    seq=seq+line
#print(len(seq))

#emission probability A0,A1,T0,T1----
A=[0.3,0.2]
T=[0.3,0.2]
G=[0.2,0.3]
C=[0.2,0.3]
#transition probability 00,01,10,11
trans_p = [[[0.98, 0.02], [0.05, 0.95]], [[0.8, 0.2], [0.5, 0.5]], [[0.51, 0.49], [0.49, 0.51]]]

##calculate of the current total probability of observing the sequence
p_x=list()
p_y=list()
log_p_x=list()
log_p_y=list()
k = 0
while (k < len(trans_p)):
    p_x.clear()
    p_y.clear()
    log_p_x.clear()
    log_p_y.clear()
    #initiation--initial probability of state p(y=0)=1; position 0(sequence position 1)
    if seq[0]=="A":
        p_x.append(1 * A[0])
    elif seq[0]=="T":
        p_x.append(1 * T[0])
    elif seq[0]=="G":
        p_x.append(1 * G[0])
    else:
        p_x.append(1 * C[0])

    p_y.append((1,0))
    l_x0=math.log(p_x[0])
    log_p_x.append(l_x0)
    log_p_y.append((0,None))
    # position 1(sequence position 2): transition probability and emission probability
    y0=0
    y1=0
    #position 1 is AT-rich[0]|GC-rich[0]
    if seq[1]=="A":
        y0 = p_x[0] * trans_p[k][0][0] * A[0]
        y1 = p_x[0] * trans_p[k][0][1] * A[1]
    elif seq[1]=="T":
        y0 = p_x[0] * trans_p[k][0][0] * T[0]
        y1 = p_x[0] * trans_p[k][0][1] * T[1]
    elif seq[1]=="G":
        y0 = p_x[0] * trans_p[k][0][0] * G[0]
        y1 = p_x[0] * trans_p[k][0][1] * G[1]
    else:
        y0 = p_x[0] * trans_p[k][0][0] * C[0]
        y1 = p_x[0] * trans_p[k][0][1] * C[1]
    p_x.append(y0 + y1)
    p_y.append((y0,y1))
    l_x1=math.log(p_x[1])
    l_y1_0=math.log(p_y[1][0])
    l_y1_1=math.log(p_y[1][1])
    log_p_x.append(l_x1)
    log_p_y.append((l_y1_0,l_y1_1))

    #from position 2 to the last position
    i=2
    while(i<len(seq)):
        m_L_yi_1=max(log_p_y[i-1][0],log_p_y[i-1][1])
        m_yi_1=math.exp(m_L_yi_1)
        if seq[i] == "A":
            y00 = math.exp(-m_L_yi_1 + log_p_y[i - 1][0] + math.log(trans_p[k][0][0] * A[0]))
            y01 = math.exp(-m_L_yi_1 + log_p_y[i - 1][0] + math.log(trans_p[k][0][1] * A[1]))
            y10 = math.exp(-m_L_yi_1 + log_p_y[i - 1][1] + math.log(trans_p[k][1][0] * A[0]))
            y11 = math.exp(-m_L_yi_1 + log_p_y[i - 1][1] + math.log(trans_p[k][1][1] * A[1]))
        elif seq[i] == "T":
            y00 = math.exp(-m_L_yi_1 + log_p_y[i - 1][0] + math.log(trans_p[k][0][0] * T[0]))
            y01 = math.exp(-m_L_yi_1 + log_p_y[i - 1][0] + math.log(trans_p[k][0][1] * T[1]))
            y10 = math.exp(-m_L_yi_1 + log_p_y[i - 1][1] + math.log(trans_p[k][1][0] * T[0]))
            y11 = math.exp(-m_L_yi_1 + log_p_y[i - 1][1] + math.log(trans_p[k][1][1] * T[1]))
        elif seq[i] == "G":
            y00 = math.exp(-m_L_yi_1 + log_p_y[i - 1][0] + math.log(trans_p[k][0][0] * G[0]))
            y01 = math.exp(-m_L_yi_1 + log_p_y[i - 1][0] + math.log(trans_p[k][0][1] * G[1]))
            y10 = math.exp(-m_L_yi_1 + log_p_y[i - 1][1] + math.log(trans_p[k][1][0] * G[0]))
            y11 = math.exp(-m_L_yi_1 + log_p_y[i - 1][1] + math.log(trans_p[k][1][1] * G[1]))
        else:
            y00 = math.exp(-m_L_yi_1 + log_p_y[i - 1][0] + math.log(trans_p[k][0][0] * C[0]))
            y01 = math.exp(-m_L_yi_1 + log_p_y[i - 1][0] + math.log(trans_p[k][0][1] * C[1]))
            y10 = math.exp(-m_L_yi_1 + log_p_y[i - 1][1] + math.log(trans_p[k][1][0] * C[0]))
            y11 = math.exp(-m_L_yi_1 + log_p_y[i - 1][1] + math.log(trans_p[k][1][1] * C[1]))

        l_yi_0=m_L_yi_1+math.log(y00+y10)
        l_yi_1=m_L_yi_1+math.log(y01+y11)
        l_xi=m_L_yi_1+math.log(y00+y01+y10+y11)
        log_p_x.append(l_xi)
        log_p_y.append((l_yi_0, l_yi_1))
        i=i+1

    print("Question ",k+1,"\nLog-likelihood value is ",log_p_x[1199])
    #print(math.exp(log_p_x[300]))
    #print(math.exp(log_p_x[400]))
    #print(math.exp(log_p_x[500]))
    k=k+1


#Question  1 :
#Log-likelyhood value is -1660.5708078732996
#3.673776266433154e-181
#2.8408481591006164e-241
#8.559407805705273e-302
#Question  2 :
#Log-likelyhood value is -1669.383184158238
#3.822118624687534e-183
#2.636358914558528e-243
#2.3122392091840993e-304
#Question  3 :
#Log-likelyhood value is -1663.6983892282558
#4.959952742429479e-182
#3.1068836008111504e-242
#1.9738631375102624e-302