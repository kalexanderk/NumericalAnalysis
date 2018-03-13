import numpy as np
import math
import pylab

A = 0.9
L = 0.8
Q = 1.5
P = math.pi*1.
C = 1.
alpha = A*pow(P, 2)/pow(L, 2) + Q

def u(t, x):
    global L, P, alpha, C
    return C*math.exp(-t*alpha)*math.sin(P*x/L)

def uCalculate(n, m, h, k, u0, e):
    global A, Q
    R = A*k/(2.*pow(h, 2))

    matr = np.zeros([n + 1, n + 1])
    b = np.zeros(n + 1)
    matr[0, 0] = 1.
    matr[n, n] = 1.
    for i in range(1, n):
        matr[i, i - 1] = -R
        matr[i, i] = 1. + 2.*R + Q*k/2.
        matr[i, i + 1] = -R
    b[0] = 0.
    b[n] = 0.
    for i in range(1, n):
        b[i] = u0[i]*(1. - 2.*R - Q*k/2.) + u0[i - 1]*R + u0[i + 1]*R

    return np.linalg.solve(matr, b)

def KNMethod(n=10, m=50):
    h = L/n
    k = L/m

    x = np.zeros(n + 1)
    for i in range(n + 1):
        x[i] = i*h
    t = np.zeros(m + 1)
    for i in range(m + 1):
        t[i] = i*k

    u_ = [[0.]]*(m + 1)

    tmp = np.zeros(n + 1)
    for i in range(n + 1):
        tmp[i] = u(0, x[i])
    u_[0] = tmp

    for i in range(1, m + 1):
        u_[i] = uCalculate(n, m, h, k, u_[i - 1], i)

    u_r = [[0.]]*(m + 1)
    for i in range(m + 1):
        tmp = np.zeros(n + 1)
        for j in range(n + 1):
            tmp[j] = u(t[i], x[j])
        u_r[i] = tmp

    return x, t, np.array(u_), np.array(u_r)


x, t, u_, u_r = KNMethod(100, 500)

print "t            Residual"
x_h = np.zeros(len(u_[0]))
t_h = np.zeros(len(u_))
for i in range(len(u_)):
    tmp = abs(u_[i] - u_r[i])
    t_h[i] = np.amax(tmp)
for i in range(len(u_)):
    print t[i], t_h[i]
print "\n"
print "x            Residual"
u_x = np.transpose(u_)
u_rx = np.transpose(u_r)
for i in range(len(u_x)):
    tmp = abs(u_x[i] - u_rx[i])
    x_h[i] = np.amax(tmp)
for i in range(len(u_x)):
    print x[i], x_h[i]


pylab.figure(1)
pylab.plot(t, t_h)
pylab.figure(2)
pylab.plot(x, x_h)
pylab.show()
