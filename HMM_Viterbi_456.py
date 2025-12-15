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
trans_p=[[[0.98,0.02],[0.05,0.95]],[[0.8,0.2],[0.5,0.5]],[[0.51,0.49],[0.49,0.51]]]
##calculate of the current most possible path 0,1
mp_x=list() #current most possible path generate xi
mp_y=list() #all possiblities calculate path at current position 0,1
log_mp_x=list() #log transformation
log_mp_y=list() #log transformation
v_y=list()

k = 0
while (k < len(trans_p)):
    mp_x.clear()
    mp_y.clear() #all possiblities calculate path at current position 0,1
    log_mp_x.clear() #log transformation
    log_mp_y.clear() #log transformation
    v_y.clear()
    #initiation--initial most possible state p(y=0)=1; position 0(sequence position 1)
    if seq[0]=="A":
        mp_x.append(1 * A[0])
    elif seq[0]=="T":
        mp_x.append(1 * T[0])
    elif seq[0]=="G":
        mp_x.append(1 * G[0])
    else:
        mp_x.append(1 * C[0])

    l_mp_x0 = math.log(mp_x[0])
    log_mp_x.append(l_mp_x0)
    mp_y.append((1,0))
    log_mp_y.append((0,None))


    # position 1 (sequence position 2): transition probability and emission probability
    # most probable path
    y0=0
    y1=0
    #position 1 is AT-rich[0]|GC-rich[0]
    if seq[1]=="A":
        y0 = mp_x[0] * trans_p[k][0][0] * A[0]
        y1 = mp_x[0] * trans_p[k][0][1] * A[1]
    elif seq[1]=="T":
        y0 = mp_x[0] * trans_p[k][0][0] * T[0]
        y1 = mp_x[0] * trans_p[k][0][1] * T[1]
    elif seq[1]=="G":
        y0 = mp_x[0] * trans_p[k][0][0] * G[0]
        y1 = mp_x[0] * trans_p[k][0][1] * G[1]
    else:
        y0 = mp_x[0] * trans_p[k][0][0] * C[0]
        y1 = mp_x[0] * trans_p[k][0][1] * C[1]
    mp_x.append(max(y0,y1))
    mp_y.append((y0,y1))
    l_mp_x1=math.log(mp_x[1])
    log_mp_x.append(l_mp_x1)
    l_mp_y1_0=math.log(mp_y[1][0])
    l_mp_y1_1=math.log(mp_y[1][1])
    log_mp_y.append((l_mp_y1_0,l_mp_y1_1))

    #from position 2 to the last position
    i=2
    while(i<len(seq)):
        l_mp_yi_1=log_mp_x[i-1]
        mp_path=log_mp_y[i-1].index(max(log_mp_y[i-1][0],log_mp_y[i-1][1]))
        if mp_path == 0:
            if seq[i] == "A":
                y0 = math.exp(-l_mp_yi_1 + log_mp_y[i - 1][0] + math.log(trans_p[k][0][0] * A[0]))
                y1 = math.exp(-l_mp_yi_1 + log_mp_y[i - 1][0] + math.log(trans_p[k][0][1] * A[1]))
            elif seq[i] == "T":
                y0 = math.exp(-l_mp_yi_1 + log_mp_y[i - 1][0] + math.log(trans_p[k][0][0] * T[0]))
                y1 = math.exp(-l_mp_yi_1 + log_mp_y[i - 1][0] + math.log(trans_p[k][0][1] * T[1]))
            elif seq[i] == "G":
                y0 = math.exp(-l_mp_yi_1 + log_mp_y[i - 1][0] + math.log(trans_p[k][0][0] * G[0]))
                y1 = math.exp(-l_mp_yi_1 + log_mp_y[i - 1][0] + math.log(trans_p[k][0][1] * G[1]))
            else:
                y0 = math.exp(-l_mp_yi_1 + log_mp_y[i - 1][0] + math.log(trans_p[k][0][0] * C[0]))
                y1 = math.exp(-l_mp_yi_1 + log_mp_y[i - 1][0] + math.log(trans_p[k][0][1] * C[1]))
        elif mp_path == 1:
            if seq[i] == "A":
                y0 = math.exp(-l_mp_yi_1 + log_mp_y[i - 1][1] + math.log(trans_p[k][1][0] * A[0]))
                y1 = math.exp(-l_mp_yi_1 + log_mp_y[i - 1][1] + math.log(trans_p[k][1][1] * A[1]))
            elif seq[i] == "T":
                y0 = math.exp(-l_mp_yi_1 + log_mp_y[i - 1][1] + math.log(trans_p[k][1][0] * T[0]))
                y1 = math.exp(-l_mp_yi_1 + log_mp_y[i - 1][1] + math.log(trans_p[k][1][1] * T[1]))
            elif seq[i] == "G":
                y0 = math.exp(-l_mp_yi_1 + log_mp_y[i - 1][1] + math.log(trans_p[k][1][0] * G[0]))
                y1 = math.exp(-l_mp_yi_1 + log_mp_y[i - 1][1] + math.log(trans_p[k][1][1] * G[1]))
            else:
                y0 = math.exp(-l_mp_yi_1 + log_mp_y[i - 1][1] + math.log(trans_p[k][1][0] * C[0]))
                y1 = math.exp(-l_mp_yi_1 + log_mp_y[i - 1][1] + math.log(trans_p[k][1][1] * C[1]))
        l_yi_0 = l_mp_yi_1 + math.log(y0)
        l_yi_1 = l_mp_yi_1 + math.log(y1)
        l_xi=max(l_yi_0,l_yi_1)
        log_mp_x.append(l_xi)
        log_mp_y.append((l_yi_0, l_yi_1))
        i=i+1

    i=len(seq)-1
    while(i>0):
        v_y.append(log_mp_y[i].index(max(log_mp_y[i][0],log_mp_y[i][1])))
        i=i-1
    #print(len(v_y))
    v_y.append(0)
    #print(len(v_y))
   # print(v_y)
    v_y.reverse()
    print("Question ",k+4,"\nLogarithm of the probability of generating the sequence data from the most probable hidden path is:\n",log_mp_x[1199])
    print("The most probable hidden path is:\n",v_y)
    #print(math.exp(log_mp_x[300]))
    #print(math.exp(log_mp_x[400]))
    #print(math.exp(log_mp_x[500]))
    k=k+1