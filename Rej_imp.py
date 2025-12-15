import numpy as np
import statistics
import matplotlib.pyplot as plt
import math
import random
random.seed(0)
d = list()
w = list()
s = list()
sw = list()
ssw = list()
i = 0
k = 0
while (k < 1000000):
    # sample theta value from uniform distribution g(theta)
    smp = np.random.uniform(0, 1, 1)[0]

    # set c=4 since the p(theta) has 3.3 as the highest value
    # rate of accept
    if (smp >= 0 and smp <= 0.3):
        post_theta = 10 / 3 - (100 / 9) * smp
    else:
        post_theta = 100 * smp / 49 - 30 / 49

    c = 4
    r = post_theta / c
    s.append(smp)
    w.append(post_theta / 1)
    sw.append(smp * post_theta / 1)
    ssw.append(smp * smp * post_theta / 1)
    # whether recept or not
    gm = np.random.uniform(0, 1, 1)[0]
    if (gm < r):
        d.append(smp)
        k = k + 1

    i = i + 1

E_theta = statistics.mean(d[:])
var_theta = statistics.variance(d[:])
sd_theta = statistics.stdev(d[:])

print("Rejection sampling:\nMean", E_theta, "\nVariance", var_theta, "\nStandard Deviation",sd_theta)
e_th=sum(d[:])/len(d[:])
sd_th=math.sqrt(sum((d[:]-e_th)*(d[:]-e_th))/len(d[:]))
print("Check:\nMean", e_th, "\nStandard Deviation", sd_th)

dt = np.array(d[:])
#print(plt.hist(dt))
plt.hist(dt, color="b")
plt.title("Rejection sampling:histogram of theta values")
plt.xlabel("Theta")
plt.ylabel("Frequency")
plt.show()

t_w = sum(w[:])
E_w=statistics.mean(w[:])
var_w=statistics.variance(w[:])
sd_w=statistics.stdev(w[:])
print("Weigt of importance sampling:\nMean", E_w, "\nVariance", var_w, "\nStandard Deviation", sd_w)

E_theta_i = sum(sw[:]) / t_w
var_theta_i = sum(ssw[:]) / t_w - E_theta_i * E_theta_i
sd_th_i=math.sqrt(var_theta_i)
print("number of improtanct samples", len(w))
print("Importance sampling:\nMean", E_theta_i, "\nVariance", var_theta_i, "\nStandard Deviation",sd_th_i)

#print(plt.hist(df))
#plt.plot(s,w)
#plt.title("Importance sampling:histogram of theta values")
#plt.xlabel("Theta")
#plt.ylabel("Weight")
#plt.show()


