import math
import numpy
import statistics
import pandas as pd
import matplotlib.pyplot as plt
import time
import time
import random
start_time=time.time()

random.seed(0)
#MH algorithm
#probability of purines
theta_list=list()
ln_th=list()

#current data
AG=2
TC=5

delta_v=0.1

rep=0
while(rep<10):
	#initialization of Theta with the current available data 2 purines and 5 pyrimidines
	#theta_list.append(AG/(AG+TC))
	theta_0=numpy.random.uniform(0,1,1)[0]
	theta_list.append(theta_0)
	ln_th.append(math.log(theta_0))
	i=1
	while(i<150000):
		#randomly proposing new theta_astrerisk from a jumping function J(theta_astrerisk|theta)
		#calculate the jumping function
		H_pgc=min(delta_v,1-theta_list[i-1])
		L_pgc=min(delta_v,theta_list[i-1])
		uni_smp1=numpy.random.uniform(0,1,1)[0]
		theta_prop=(theta_list[i-1]-L_pgc)+uni_smp1*(H_pgc+L_pgc)

		#receptive probability=
		#(likelihood_thetalist*
		#prior probility of theta_proposed*
		#junmp_function of the current value given the proposed value)
		#/
		#(likelihood*
		# prior probability of theta[0]*
		#jumpfunction of the proposed value given the current value)
		j_pgc=1/(H_pgc+L_pgc)
		#jump function of the proposed value given the current value
		#(the same as previous)
		H_cgp=min(delta_v,1-theta_prop)
		L_cgp=min(delta_v,theta_prop)
		j_cgp=1/(H_cgp+L_cgp)
		# total number of possible combination
		#cb=math.comb(AG+TC,AG)
		#r=(cb*math.pow(theta_list[0],AG)*math.pow((1-theta_list[0]),TC)*1*j_cgp)/(cb*math.pow(theta_prop,AG)*math.pow((1-theta_prop),TC)*1*j_pgc)
		if (theta_list[i - 1] >= 0 and theta_list[i - 1] <= 0.3):
			prior_theta = 10 / 3 - (100 / 9) * theta_list[i - 1]
		else:
			prior_theta = 100 * theta_list[i - 1] / 49 - 30 / 49

		if (theta_prop >= 0 and theta_prop <= 0.3):
			prior_theta_prop = 10 / 3 - (100 / 9) * theta_prop
		else:
			prior_theta_prop = 100 * theta_prop / 49 - 30 / 49

		r = (math.pow(theta_prop, AG) * math.pow((1 - theta_prop), TC) * prior_theta_prop * j_cgp) / (
			 math.pow(theta_list[i - 1], AG) * math.pow((1 - theta_list[i - 1]), TC) * prior_theta * j_pgc)

		p_r=min(r,1)
		uni_smp2=numpy.random.uniform(0,1,1)[0]
		if (uni_smp2<p_r):
			theta_list.append(theta_prop)
			ln_th.append(math.log(theta_prop))
		else:
			cur_thetha=theta_list[i-1]
			theta_list.append(cur_thetha)
			ln_th.append(math.log(cur_thetha))
		i=i+1

	#print("Current algorithm runs for", time.time()-start_time)

	E_theta=statistics.mean(theta_list[15000:])
	var_theta=statistics.variance(theta_list[15000:])
	sd_theta=statistics.stdev(theta_list[15000:])
	print(rep,"Theta0",theta_0,"\nMean",E_theta,"Variance", var_theta, "Standard Deviation", sd_theta)

	e_th=sum(theta_list[15000:])/len(theta_list[15000:])
	sd_th=math.sqrt(sum((theta_list[15000:]-e_th)*(theta_list[15000:]-e_th))/len(theta_list[15000:]))

	print("Check Mean",e_th,"Standard Deviation", sd_th)

	#dt=numpy.array(theta_list[15000:])
	#plt.hist(dt,color="b")
	#plt.title("Histogram of theta values")
	#plt.xlabel("Theta")
	#plt.ylabel("Frequency")
	#plt.show()
	#ln_dt=numpy.array(ln_th)
	#plt.plot(ln_dt,color="b")
	#plt.show()
	rep=rep+1
