import numpy as np
import control as con
import optim_tools
import matplotlib.pyplot as plt

# Solution directly from polyplacement char. polynomial to wanted polynomial (fastest) with given polynomical coeff a
def k_explizit_a(roots_p1, roots_pmin, pmin, a):
    n = len(roots_p1)
    
    k0 = np.zeros((n,))
    k1 = np.zeros((n,))

    a_tilde_p1 = np.matrix(np.poly(roots_p1)[:0:-1]).T # get the characteristical polynomial backwards
    a_tilde_pmin = np.matrix(np.poly(roots_pmin)[:0:-1]).T # get the characteristical polynomial backwards
    
    k1 = (a_tilde_p1 - a_tilde_pmin) / (1.0-1.0/pmin)
    k0 = a_tilde_p1 - a - k1

    return np.matrix(k0), np.matrix(k1)
#%timeit k_explizit_a([1,2,4], [2,4,8], 0.1, a)

# Do canonical form first (costy)
def k_explizit_Ab(roots_p1, roots_pmin, pmin, A, b):
    (A_R, _, _, _), _, _ = optim_tools.get_Steuerungsnormalform(A, b, b, 0)
    a = -A_R[-1][:].T
    n = len(roots_p1)
    
    k0 = np.zeros((n,))
    k1 = np.zeros((n,))

    a_tilde_p1 = np.matrix(np.poly(roots_p1)[:0:-1]).T # get the characteristical polynomial backwards
    a_tilde_pmin = np.matrix(np.poly(roots_pmin)[:0:-1]).T # get the characteristical polynomial backwards
    
    k1 = (a_tilde_p1 - a_tilde_pmin) / (1.0-1.0/pmin)
    k0 = a_tilde_p1 - a - k1

    return np.matrix(k0), np.matrix(k1)
#%timeit k_explizit_Ab([1,2,4], [2,4,8], 0.1, A, b)

# Use python control to place poles and then interpolate (faster with A and b)
def k_explizit_Ab2(roots_p1, roots_pmin, pmin, A, b):
    r0 = con.place(A, b, roots_p1) #k(p=1)
    r1 = con.place(A, b, roots_pmin) #k(p=pmin)

    # This seems to work as expected
    k1 = 1.0/(1.0-1.0/pmin) * (r0 - r1)
    k0 = r0 - k1
    return np.matrix(k0).T, np.matrix(k1).T

#%timeit k_explizit_Ab2([1,2,4], [2,4,8], 0.1, A, b)

#"""
# Aemo of plotting complex functions in python.
#
# Jim M | Feb 2011 | GPL
#"""
# Plotting functions ; see the example below
# and http://matplotlib.sourceforge.net/
#from matplotlib.pyplot import plot, legend

# Complex math (cmath) python functions ;
# see  see http://docs.python.org/library/cmath.html
#from cmath import sin, cos, exp, pi, log, polar, rect, phase, sqrt

# Note that python represents imaginary numbers like "3i" as "3j",
# where "j" meaning "sqrt(-1) must be preceded by a number,
# so "sqrt(-1)" alone would in python be "1j".
#
#     (3,4) complex rectangular form:     z = 3 + 4j
#     (x,y) complex rectangular form :    z = x + y * 1j
#     polar form :                        z = r * exp(1j * theta)
#     abs(z)   is length of complex number = r
#     phase(z) is angle of complex number = theta
#     z.real   is real part
#     z.imag   is imaginary part
#
# abs() is a python built-in; as are complex numbers themselves.
# But the other functions needed to be imported in their complex versions.
# The numeric constant pi can be imported from math or cmath.

# Remember that
# 1. lambda(x: ...) is an anyonymous function of x, e.g. lambda(x: 2*x+1)
# 2. map(f, [a, b, c, ...])  # returns [f(a), f(b), f(c), ...]
#
#
# == So here are a few utility functions for multiplying scalars and vectors.
#

# return real part of a vector
def real_vector(vector):
    return map(lambda x: x.real, vector)

# return imaginary part of a vector
def imag_vector(vector):
    return map(lambda x: x.imag, vector)

def plot_moving_poles(A, b, c, d, k_0, k_1, pmin=0.1):
    poles = []
    for p in np.arange(1, pmin, -0.001):
        sys_closed = con.ss(A-b*(k_0+1.0/p*k_1).T, b, c, d)
        pole = con.matlab.pole(sys_closed)
        poles.append(pole)
    
    # another approach to plot
    real_part = real_vector(poles)
    imag_part = imag_vector(poles)

    # Display a window with a plot of real, imag
    plt.plot(real_part, imag_part, 'b-')
    plt.plot(real_part[0], imag_part[0], 'b*')
    plt.plot(real_part[-1], imag_part[-1], 'rx')
    plt.show
    
def narf():
    # a scalar times a vector returns a vector
    def scale_vector(scale, vector):
        result = [0]*len(vector)
        for i in range(len(result)):
            result[i] = scale * vector[i]
        return result

    # dot product of two vectors = sum(x[0]*y[0] + ... + x[n-1]*y[n-1])
    def vector_dot(vector1, vector2):
        result = 0
        for i in range(len(vector1)):
            result += vector1[i] * vector2[i]
        return result



    from cmath import sin, cos, exp, pi, log, polar, rect, phase, sqrt
    # Generate numbers around the complex unit circle.
    # (These are the same numbers that show up in the Fourier Transform.)
    N = 128
    theta = scale_vector(2*pi/N, range(N))
    exp_theta = map(lambda x: exp(1j * x), theta)

    real_part = real_vector(exp_theta)
    imag_part = imag_vector(exp_theta)

    # Display a window with a plot of real, imag
    plt.plot(theta, real_part, '-', label='real')
    plt.plot(theta, imag_part, '--', label='imag')
    plt.legend(loc='lower left', title='exp(i theta)')
    plt.show()


    # TEST k_explizit functions

    # uboot
    A_x = np.matrix([[0., 1., 0.],
                  [0., 0., 1.],
                  [0., 0., -0.005]])
    a_x = -A_x[-1][:].T #!!!!

    b_x = np.matrix([[0], [0], [1.]])

    d_x = 0
    c_x = np.matrix([[1], [0], [0]])

    sys_x = con.ss(A_x, b_x, c_x.T, d_x)
    #mag, phase, omega = bode(sys1)
    #arr1, arr2 = control.step(sys)
    #plt.plot(arr2, arr1)

    #poles, zeros = control.matlab.pzmap(sys, True)
    plt.axis([-5,.1,-3,3])
    #plt.show

    roots_p1_x = [-1, -1+1j, -1-1j]
    roots_pmin_x = [-3, -3+2j, -3-2j]
    pmin_x = 0.1

    k0_x, k1_x = k_explizit(roots_p1_x, roots_pmin_x, pmin_x, a_x)

    k0_x2, k1_x2 = k_explizit_x(roots_p1_x, roots_pmin_x, pmin_x, A_x, b_x)
    k0_x3, k1_x3 = k_explizit2(roots_p1_x, roots_pmin_x, pmin_x, A_x, b_x)

    #print k0_x, k1_x
    #print k0_x2, k1_x2
    #print k0_x3, k1_x3

    assert np.allclose(k0_x, k0_x2)
    assert np.allclose(k1_x, k1_x2)
    assert np.allclose(k0_x2, k0_x3)
    assert np.allclose(k1_x2, k1_x3)


    poles = []
    for p in np.arange(1, 0.09, -0.01):
        #sys_cl = control.ss(A-b*(r0), b, c.T, d)
        #poles2, zeros2 = control.matlab.pzmap(sys_cl, True)
        #sys_cl = control.ss(A-b*(r1), b, c.T, d)
        #poles2, zeros2 = control.matlab.pzmap(sys_cl, True)
        #print p
        sys_cl = con.ss(A-b*(k0_x3+1.0/p*k1_x3).T, b, c, d)
        pole, zeros = con.matlab.pzmap(sys_cl, True)
        #print pole[0]
        poles.append(pole)
        #print poles2
    plt.show
    
    # another approach to plot
    real_part = real_vector(poles)
    imag_part = imag_vector(poles)

    # Display a window with a plot of real, imag
    plt.plot(real_part, imag_part, '-x')
    plt.show