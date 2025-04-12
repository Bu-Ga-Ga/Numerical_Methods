import math as mth
import numpy as np
#отрезок [a,b]
a = 0
b = 1
eps = 0.001

y0 = [1.,1.]
x0 = 0
h = 1
ph = 0
n = 0

def func(x,y):
    global n
    n += 1
    #return x,x
    #return 2*x,2*x,2*x
    return y[1],y[0]




#Первый пункт
def k11(x,y,h):
    return np.multiply(func(x,y),h)

def k21(x,y,h):
    return np.multiply(func(x+h,np.add(y,k11(x,y,h))),h)

def yn1(x,y,h):
    return np.add(y,np.multiply(np.add(k21(x,y,h),k11(x,y,h)),0.5))
x1 = x0
y1 = y0
h1 = h
cnt = 0
while x1 < b:
    if x1 + h1 > b:
        h1 = b - x1
    yi1 = yn1(x1,y1,h1)
    yi21 = yn1(x1,y1,h1/2)
    yi22 = yn1(x1+h1/2,yi21,h1/2)
    ph = 0
    for i in range(len(yi1)):
        if ph < abs(yi22[i]-yi1[i]) / (1 - 2**(-2)):
            ph = abs(yi22[i]-yi1[i]) / (1 - 2**(-2))
    cnt += 1
    if ph > eps*h1/(b-a):
        h1 = h1/2
    else:
        x1 = x1 + h1
        y1 = yi1
        h1 = 2*h1
print (x1 , y1, n, cnt)
#print (yn(1,1,0.1))





#Второй пункт
def k12(x,y,h):
    return np.multiply(func(x,y),h)

def k22(x,y,h):
    return np.multiply(func(x+0.5*h,np.add(y,np.multiply(k11(x,y,h),0.5))),h)

def k32(x,y,h):
    # return h*func(x+h, y-k21(x,y,h)+2*k22(x,y,h))
    return np.multiply(func(x+h, np.subtract( np.add(y,np.multiply(k22(x,y,h),2)), k12(x,y,h))),h)

def yn2(x,y,h):
    # return y + 1/6*(k12(x,y,h) + 4* k22(x,y,h) + k32(x,y,h))
    return np.add(y, np.multiply(np.add(np.add(k12(x,y,h),np.multiply(k22(x,y,h),4)),k32(x,y,h)),1/6))

x2 = 0
y2 = [1.,1.]
h2 = 1
n = 0
ct = 0
while x2 < b:
    if x2 + h2 > b:
        h2= b - x2
    sol1 = yn1(x2,y2,h2)
    sol2 = yn2(x2,y2,h2)
    ph = 0
    for i in range(len(sol1)):
        if ph <  abs(sol2[i]-sol1[i]) / (1 - 2**(-3)):
            ph = abs(sol2[i]-sol1[i]) / (1 - 2**(-3))
    ct += 1
    if ph > eps*h2/(b-a):
        h2 = h2/2
    else:
        x2 = x2 + h2
        y2 = sol1
        h2 = 2*h2
print(x2,y2,n,ct)






#Третий пункт
def k13(x,y,h):
    return np.multiply(func(x,y),h)

def k23(x,y,h):
    return np.multiply(func(x+0.5*h,np.add(y,np.multiply(k13(x,y,h),1/2))),h)

def k33(x,y,h):
    return np.multiply(func(x+0.5*h,np.add(y,np.multiply(k23(x,y,h),0.5))),h)

def k43(x,y,h):
    return np.multiply(func(x+h,np.add(y,k33(x,y,h))),h)

def yn3(x,y,h):
    return np.add(y,np.multiply(np.add(np.add(k13(x,y,h),np.multiply(np.add(k23(x,y,h),k33(x,y,h)),2)),k43(x,y,h)),1/6))

def E(x,y,h):
    return 2/3*(k13(x,y,h)-k23(x,y,h)-k33(x,y,h)+k43(x,y,h))

x3 = x0
y3 = y0
h3 = h
ct = 0
n = 0
while x3 < b:
    if x3 + h3 > b:
        h3 = b - x3
    sol3 = yn3(x3,y3,h3)
    ct += 1
    ph = 0
    for i in range(len(sol3)):
        if E(x3,y3,h3)[i] > eps*h3/(b-a):
            ph += 1
    if ph == len(sol3):
        h3 = h3/2
    else:
        x3 = x3 + h3
        y3 = sol3
        h3 = 2*h3
print(x3,y3,n,ct)




