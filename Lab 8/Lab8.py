import numpy as np
import matplotlib.pyplot as plt

A = 0.267; B = 0.871; C = -0.194; D = -1.093; E = -1.572
P = -3.; Q = lambda x: 1./x

def y(x):
    global A; global B; global C; global D; global E
    return A*pow(x, 2) + B*x + C + 1./(D*x + E)

def dy(x):
    global A; global B; global D; global E
    return 2.*A*x + B - D/pow((D*x + E), 2)

def d2y(x):
    global A; global D; global E
    return 2.*A + 2.*pow(D, 2)/pow((D*x + E), 3)

def f(x):
    global P
    return d2y(x) + P*dy(x) + Q(x)*y(x)

def fdMethod(n=30, x0=0.4, xn=0.7):
    h = float((xn - x0)/n)

    x = np.zeros(n + 1); x[0] = x0
    for i in range(1, n + 1):
        x[i] = x[i - 1] + h

    y_ = np.zeros(n + 1)
    for i in range(n + 1):
        y_[i] = y(x[i])

    matr = np.array([[0.]*(n + 1)]*(n + 1))
    matr[0, 0] = -2./pow(h, 2) + Q(x[0]) - 3./h + 1.5*P
    matr[0, 1] = 2./pow(h, 2)
    matr[n, n] = 2.
    for i in range(1, n):
        matr[i, i - 1] = 1./pow(h, 2) - P/(2*h)
        matr[i, i] = -2./pow(h, 2) + Q(x[i])
        matr[i, i + 1] = 1./pow(h, 2) + P/(2*h)

    b = np.array([0.]*(n + 1))
    b[0] = f(x[0]) - 2.*(1.5*y(x0) - dy(x0))/h + P*(1.5*y(x0) - dy(x0))
    b[n] = 2.*y(xn)
    for i in range(1, n):
        b[i] = f(x[i])

    return x, np.linalg.solve(matr, b), y_, h

N = 10
BEST = 4
n = (30, 40, 60, 80, 100, 150, 300, 400, 500, 3000)

x_ = [[0.]]*N
y_fd = [[0.]]*N
y_ = [[0.]]*N
h = np.zeros(N)
for i in range(N):
    x_[i], y_fd[i], y_[i], h[i] = fdMethod(n=n[i])

fdAcc = [[0]]*N
for i in range(N):
    fdAcc[i] = abs(y_[i] - y_fd[i])

hfdAcc = np.zeros(N)
for i in range(N):
    hfdAcc[i] = np.amax(fdAcc[i])


print "x            y                y_fd               Acc"
for i in range(n[BEST] + 1):
    print x_[BEST][i], "   ", y_[BEST][i], "  ", y_fd[BEST][i], "    ", fdAcc[BEST][i]

print "n        h                    e(h)"
for i in range(N):
    print n[i], "    ", h[i], "     ", hfdAcc[i]


# plt.plot(x_[BEST], y_[BEST], x_[BEST], y_fd[BEST])
# plt.title("y and y_fd for n = %i"%n[BEST])
# plt.grid(True)
# plt.xlabel('x')
# plt.ylabel('y')
# plt.show()
# plt.plot(x_[BEST], fdAcc[BEST])
# plt.title("Residual for n = %i"%n[BEST])
# plt.grid(True)
# plt.xlabel('x')
# plt.ylabel('fdAcc')
# plt.show()
#
# plt.plot(h, hfdAcc)
# plt.title("Residual max for all n")
# plt.grid(True)
# plt.xlabel('h')
# plt.ylabel('hfdAcc')
# plt.show()

# print " "
# aa=[0.]*(len(h))
# for i in range(len(h)):
#     aa[i] = hfdAcc[i]/(h[i]**2)
# aak = np.mean(aa)
# print aa
# print aak
