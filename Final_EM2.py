import math
pi_A=0.25
pi_C=0.25
pi_G=0.25
pi_T=0.25
s=0.5
i=1
AA=1715
AC=799
AG=536
AT=402

CA=828
CC=1359
CG=478
CT=333

GA=580
GC=451
GG=776
GT=243

TA=397
TC=343
TG=261
TT=499
ATCG=AA+AC+AG+AT+CA+CC+CG+CT+GA+GC+GG+GT+TA+TC+TG+TT
#print("total pairs", ATCG)
#print("Initialization\n", "Cycle 0" , "substitution rate", s)
#print("pi_A:",pi_A," pi_C:",pi_C, " pi_G:",pi_G, " pi_T:",pi_T)
while(i<50):
	## E-step
	KAA=(s*pi_A)/(math.exp(-s)+(1-math.exp(-s))*pi_A)
	KCC=(s*pi_C)/(math.exp(-s)+(1-math.exp(-s))*pi_C)
	KGG=(s*pi_G)/(math.exp(-s)+(1-math.exp(-s))*pi_G)
	KTT=(s*pi_T)/(math.exp(-s)+(1-math.exp(-s))*pi_T)
	K01=s/(1-math.exp(-s))
	K=(AA*KAA+CC*KCC+GG*KGG+TT*KTT+(AC+AG+AT+CA+CG+CT+GA+GC+GT+TA+TC+TG)*K01)/10000
	#i,j,w
	#E(c(w))=delta(i,j)*delta(i,w)(math.exp(-s)/(math.exp(-s)+(1-math.exp(-s))*pi_j))
	#+ delta(i,j)*(2*delta(i,w)+ (s/(1-math.exp(-s))-1)*pi_w)*(1-math.exp(-s))*pi_j/(math.exp(-s)+(1-math.exp(-s))*pi_j)
	#+ (1-delta(i,j))*(delta(i,w)*(1-delta(j,w))+(1-delta(i,w))*delta(j,w) +(s/(1-math.exp(-s))-1)*pi_w)*1
	##w= A
	cnt_A = AA * ((math.exp(-s) / (math.exp(-s) + (1 - math.exp(-s)) * pi_A)) + (2 + (s / (1 - math.exp(-s)) - 1) * pi_A) * (1 - math.exp(-s)) * pi_A / ( math.exp(-s) + (1 - math.exp(-s)) * pi_A)) + CC * ((s / (1 - math.exp(-s)) - 1) * pi_A) * (1 - math.exp(-s)) * pi_C / ( math.exp(-s) + (1 - math.exp(-s)) * pi_C) + GG * ((s / (1 - math.exp(-s)) - 1) * pi_A) * ( 1 - math.exp(-s)) * pi_G / (math.exp(-s) + (1 - math.exp(-s)) * pi_G) + TT * ((s / (1 - math.exp(-s)) - 1) * pi_A) * ( 1 - math.exp(-s)) * pi_T / (math.exp(-s) + (1 - math.exp(-s)) * pi_T) + (AC + AG + AT + CA + GA + TA) * ((1 + (s / (1 - math.exp(-s)) - 1) * pi_A) * 1) + (CG + CT + GC + GT + TC + TG) * ((s / (1 - math.exp(-s)) - 1) * pi_A * 1)

	cnt_C = CC * ((math.exp(-s) / (math.exp(-s) + (1 - math.exp(-s)) * pi_C)) + (2 + (s / (1 - math.exp(-s)) - 1) * pi_C) * (1 - math.exp(-s)) * pi_C / ( math.exp(-s) + (1 - math.exp(-s)) * pi_C)) + AA * ((s / (1 - math.exp(-s)) - 1) * pi_C) * (1 - math.exp(-s)) * pi_A / ( math.exp(-s) + (1 - math.exp(-s)) * pi_A) + GG * ((s / (1 - math.exp(-s)) - 1) * pi_C) * ( 1 - math.exp(-s)) * pi_G / (math.exp(-s) + (1 - math.exp(-s)) * pi_G) + TT * ((s / (1 - math.exp(-s)) - 1) * pi_C) * ( 1 - math.exp(-s)) * pi_T / (math.exp(-s) + (1 - math.exp(-s)) * pi_T) + (CA + CG + CT+ AC + GC + TC) * ((1 + (s / (1 - math.exp(-s)) - 1) * pi_C) * 1) + (AG + AT + GA + GT + TA + TG) * ((s / (1 - math.exp(-s)) - 1) * pi_C * 1)

	cnt_G = GG * ((math.exp(-s) / (math.exp(-s) + (1 - math.exp(-s)) * pi_G)) + (2 + (s / (1 - math.exp(-s)) - 1) * pi_G) * (1 - math.exp(-s)) * pi_G / ( math.exp(-s) + (1 - math.exp(-s)) * pi_G)) + AA * ((s / (1 - math.exp(-s)) - 1) * pi_G) * (1 - math.exp(-s)) * pi_A / ( math.exp(-s) + (1 - math.exp(-s)) * pi_A) + CC * ((s / (1 - math.exp(-s)) - 1) * pi_G) * ( 1 - math.exp(-s)) * pi_C / (math.exp(-s) + (1 - math.exp(-s)) * pi_C) + TT * ((s / (1 - math.exp(-s)) - 1) * pi_G) * ( 1 - math.exp(-s)) * pi_T / (math.exp(-s) + (1 - math.exp(-s)) * pi_T) + (GA + GC + GT+ AG + CG + TG) * ((1 + (s / (1 - math.exp(-s)) - 1) * pi_G) * 1) + (AC + AT + CA + CT + TA + TC) * ((s / (1 - math.exp(-s)) - 1) * pi_G * 1)

	cnt_T = TT * ((math.exp(-s) / (math.exp(-s) + (1 - math.exp(-s)) * pi_T)) + (2 + (s / (1 - math.exp(-s)) - 1) * pi_T) * (1 - math.exp(-s)) * pi_T / ( math.exp(-s) + (1 - math.exp(-s)) * pi_T)) + AA * ((s / (1 - math.exp(-s)) - 1) * pi_T) * (1 - math.exp(-s)) * pi_A / ( math.exp(-s) + (1 - math.exp(-s)) * pi_A) + CC * ((s / (1 - math.exp(-s)) - 1) * pi_T) * ( 1 - math.exp(-s)) * pi_C / (math.exp(-s) + (1 - math.exp(-s)) * pi_C) + GG * ((s / (1 - math.exp(-s)) - 1) * pi_T) * ( 1 - math.exp(-s)) * pi_G / (math.exp(-s) + (1 - math.exp(-s)) * pi_G) + (TA + TC + TG+ AT + CT + GT) * ((1 + (s / (1 - math.exp(-s)) - 1) * pi_T) * 1) + (AC + AG + CA + CG + GA + GC) * ((s / (1 - math.exp(-s)) - 1) * pi_T * 1)

	#print("Total events", cnt_A + cnt_C + cnt_G + cnt_T)
	pi_A=cnt_A/(cnt_A + cnt_C + cnt_G + cnt_T)
	pi_C=cnt_C/(cnt_A + cnt_C + cnt_G + cnt_T)
	pi_G=cnt_G/(cnt_A + cnt_C + cnt_G + cnt_T)
	pi_T=cnt_T/(cnt_A + cnt_C + cnt_G + cnt_T)
	s=K
	if (i < 11):
		# print("Total events", cnt_A + cnt_C + cnt_G + cnt_T)
		print("Cycle", i, "substitution rate", s)
		print("    pi_A:",pi_A,"\n    pi_C:",pi_C, "\n    pi_G:",pi_G, "\n    pi_T:",pi_T)
	# else:
	#	print("Cycle", i, "substitution rate", s)
	#	print("pi_A: ", pi_A, "pi_C: ", pi_C, "pi_G: ", pi_G, "pi_T: ", pi_T)
	i = i + 1

#print((AC+AG+AT+CA+CG+CT+GA+GC+GT+TA+TC+TG)/10000)