import numpy as np
import math

eps = 0.01
p = 5 




#11
a = 0.00001  
b = 1.00001
def solution(x):
    #return 1  
    return x + math.exp(-x)

def func(x):
    #return 1 + (1/x) * (math.cos(x/2) - 1)
    return math.exp(-x)

def K(x, s):
    #return math.sin(x * s)
    return 0.5 * x*math.exp(s)





"""
#12
a = 0.00001  
b = 0.50001
def solution(x):
    return 1  

def func(x):
    return 1 + (1/x) * (math.cos(x/2) - 1)

def K(x, s):
    return math.sin(x * s)
"""




"""
#15
a = 0.00001  
b = 2 * 3.1415
def solution(x):
    return math.cos(2*x)  

def func(x):
    return math.cos(2*x)

def K(x, s):
    return math.sin(x)*math.cos(s)
"""



"""
#17
a = 0.00001  
b = 1.00001
def solution(x):
    return 1 + (4/9)*x  

def func(x):
    return 1

def K(x, s):
    return x*(s**2)
"""



"""
#15
a = 0.00001  
b = 2 * 3.1415
def solution(x):
    return math.cos(2*x)  

def func(x):
    return math.cos(2*x)

def K(x, s):
    return math.sin(x)*math.cos(s)
"""



"""
#18
a = 0.00001  
b = 1.00001
def solution(x):
    return x 

def func(x):
    return (5/6) * x

def K(x, s):
    return 0.5 * x * s
"""


"""
#19
a =-1.00001  
b = 1.00001
def solution(x):
    return 1 

def func(x):
    return 1 - x * (math.exp(x) - math.exp(-x))

def K(x, s):
    return (x**2)*math.exp(x*s)
"""











def u_n(x, u):
    h = (b - a) / (len(u) - 1)
    res = 0
    for i in range(len(u)):
        res += 2 * h * K(x, a + i * h) * u[i]
    return res + func(x)


def bode_integral(alpha, betha, u1, u2, n):
    if n % 4 != 0:
        raise ValueError("Error")
    
    h = (betha - alpha) / n
    integral = 0
    
    for i in range(0, n, 4):
        x0 = alpha + i * h
        x1 = x0 + h
        x2 = x0 + 2 * h
        x3 = x0 + 3 * h
        x4 = x0 + 4 * h
        
        u1_0 = u_n(x0 + h / 2, u1)
        u1_1 = u_n(x1 + h / 2, u1)
        u1_2 = u_n(x2 + h / 2, u1)
        u1_3 = u_n(x3 + h / 2, u1)
        u1_4 = u_n(x4 + h / 2, u1)
        
        u2_0 = u_n(x0 + h / 2, u2)
        u2_1 = u_n(x1 + h / 2, u2)
        u2_2 = u_n(x2 + h / 2, u2)
        u2_3 = u_n(x3 + h / 2, u2)
        u2_4 = u_n(x4 + h / 2, u2)
        
        integral += (2 * h / 45) * (7 * (u1_0 - u2_0) ** 2 + 32 * (u1_1 - u2_1)**2 + 12 * (u1_2 - u2_2)**2 + 32 * (u1_3 - u2_3)**2 + 7 * (u1_4 - u2_4) ** 2)
    
    return integral

def L2_norm(u1, u2):
    sep = 10
    h = (b - a) / sep
    mod_eps = eps * h / (b - a)
    alpha = a
    betha = a
    res = 0
    while abs(b - betha) > eps:
        alpha = betha
        betha = min(betha + 2 * (b - a) / sep, b)
        n = 4
        I_h = bode_integral(alpha, betha, u1, u2, n)
        I_2h = bode_integral(alpha, betha, u1, u2, 2 * n)
        delt = abs(I_h - I_2h) / (1 - 2 ** (1 - p))
        while delt > mod_eps:
            I_h = I_2h
            n *= 2
            I_2h = bode_integral(alpha, betha, u1, u2, 2 * n)
            delt = abs(I_2h - I_h) / (1 - 2 ** (1 - p))
            mod_eps = eps * h / (b - a)
        res += I_h
    return res ** 0.5

def sle(n):
#    h = (b - a) / (n - 1)
    h = (b-a)/n
    A = np.eye(n)
    f = np.zeros((n, 1))
    
    for i in range(n):
        xi = a + i * h
        f[i] = func(xi)
        for j in range(0, n, 4):
            if j + 4 < n:
                xj_0 = a + j * h
                xj_1 = a + (j + 1) * h
                xj_2 = a + (j + 2) * h
                xj_3 = a + (j + 3) * h
                xj_4 = a + (j + 4) * h
                
                A[i][j] -= 2 * h / 45 * (7 * K(xi, xj_0) + 32 * K(xi, xj_1) + 12 * K(xi, xj_2) + 32 * K(xi, xj_3) + 7 * K(xi, xj_4))
    
    u1 = np.linalg.solve(A, f)
    
    n *= 2
 #   n -= 1
 #   h = (b - a) / (n - 1)
    h = (b - a)/n
    A = np.eye(n)
    f = np.zeros((n, 1))
    
    for i in range(n):
        xi = a + i * h
        f[i] = func(xi)
        for j in range(0, n, 4):
            if j + 4 < n:
                xj_0 = a + j * h
                xj_1 = a + (j + 1) * h
                xj_2 = a + (j + 2) * h
                xj_3 = a + (j + 3) * h
                xj_4 = a + (j + 4) * h
                
                A[i][j] -= 2 * h / 45 * (7 * K(xi, xj_0) + 32 * K(xi, xj_1) + 12 * K(xi, xj_2) + 32 * K(xi, xj_3) + 7 * K(xi, xj_4))
    
    u2 = np.linalg.solve(A, f)
    
    integ = L2_norm(u1, u2)

    while integ > eps:
        n *= 2
#        n -= 1
#        h = (b - a) / (n - 1)
        h = (b - a)/n
        A = np.eye(n)
        f = np.zeros((n, 1))
        
        for i in range(n):
            xi = a + i * h
            f[i] = func(xi)
            for j in range(0, n, 4):
                if j + 4 < n:
                    xj_0 = a + j * h
                    xj_1 = a + (j + 1) * h
                    xj_2 = a + (j + 2) * h
                    xj_3 = a + (j + 3) * h
                    xj_4 = a + (j + 4) * h

                    A[i][j] -= 2 * h / 45 * (7 * K(xi, xj_0) + 32 * K(xi, xj_1) + 12 * K(xi, xj_2) + 32 * K(xi, xj_3) + 7 * K(xi, xj_4))
        u2 = u1
        u1 = np.linalg.solve(A, f)
        
        integ = L2_norm(u1, u2)
        print(f'n = {(n + 1) // 2}:   = {integ}')
    
    return u1

def compare(u):
    res = 0
    h = (b - a) / (len(u) - 1)
    for i in range(1, len(u), 2):
        res += abs(u[i] - solution(a + i * h)) ** 2 * h
    otv = res ** 0.5
    print('u_n - u = ', otv)

compare(sle(5))
