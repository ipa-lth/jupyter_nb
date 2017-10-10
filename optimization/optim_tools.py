#!/usr/bin/env python
import cvxpy

'''
## Example bisection code (MATLAB)
## https://see.stanford.edu/materials/lsocoee364a/hw6sol.pdf (p.3)

## Example bisection code
## https://www.coursehero.com/file/p3eovvo/To-solve-the-problem-we-can-use-a-bisection-method-solving-an-LP-feasibility/

import numpy as np
import matplotlib.pyplot as plt
import cvxpy as cvx

k = 201
t = -3 + 6 * np.arange(k) / (k - 1)
y = np.exp(t)

Tpowers = np.vstack((np.ones(k), t, t**2)).T

a = cvx.Variable(3)
b = cvx.Variable(2)
gamma = cvx.Parameter(sign='positive')
lhs = cvx.abs(Tpowers * a - (y[:, np.newaxis] * Tpowers) * cvx.vstack(1, b))
rhs = gamma * Tpowers * cvx.vstack(1, b)
problem = cvx.Problem(cvx.Minimize(0), [lhs <= rhs])

l, u = 0, np.exp(3) # initial upper and lower bounds
bisection_tol = 1e-3 # bisection tolerance
while u - l >= bisection_tol:
    gamma.value = (l + u) / 2
    # solve the feasibility problem
    problem.solve()
    if problem.status == 'optimal':
        u = gamma.value
        a_opt = a.value
        b_opt = b.value
        objval_opt = gamma.value
    else:
        l = gamma.value

y_fit = (Tpowers * a_opt / (Tpowers * np.vstack((1, b_opt)))).A1

plt.figure()
plt.plot(t, y, 'b', t, y_fit, 'r+')
plt.xlabel('t')
plt.ylabel('y')
plt.show()

plt.figure()
plt.plot(t, y_fit - y)
plt.xlabel('t')
plt.ylabel('err')
plt.show()
'''

''' Bisection function
Tries to find the max value of parameter of the given feasibility problem by bisection


l: lower bound with infeasible solution
u: upper bound with optimal solution
bisection_tol: bisection tolerance
sol: solver for problem
problem: cvxpy Problem
variables: array of variables in problems
parameter:

Note: make sure solution with l is infeasible, with u is optimal
Example: [[o_Q, o_z0, o_z1], o_g] = bisect(
                                        0, 2,
                                        prob, g,
                                        [Q, z0, z1],
                                        solver=cvxpy.SCS,
                                        verbose=False)
'''
def bisect(l, u, problem, parameter, variables,
           bisection_tol=1e-3, solver=cvxpy.CVXOPT, verbose=False):

    # cross check bound
    if (u < l):
        #print "upperBound < lowerBound"
        raise ValueError("upperBound({}) < lowerBound({})".format(u, l))

    # check validity solution of l is optimal, solution of u is infeasible
    parameter.value = l
    problem.solve(solver=solver)
    lStatus = problem.status

    parameter.value = u
    problem.solve(solver=solver)
    uStatus = problem.status

    if not ('optimal' in lStatus and 'optimal' not in uStatus):
        #print "UpperBound({})={}, LowerBound({})={}".format(u, uStatus, l, lStatus)
        raise ValueError("UpperBound({})={}, LowerBound({})={}".format(u, uStatus, l, lStatus))

    variables_opt = [None] * len(variables)
    while u - l >= bisection_tol:
        parameter.value = (l + u) / 2.0
        ## solve the feasibility problem
        problem.solve(solver=solver)
        if verbose:
            print "Range: {}-{}; parameter {} -> {}".format(l, u, parameter.value, problem.status)

        if 'optimal' in problem.status:
            l = parameter.value
            # update Variables
            for i in range(len(variables)):
                variables_opt[i] = variables[i].value
            #a_opt = a.value
            #b_opt = b.value
            # update Parameters
            objval_opt = parameter.value
        else:
            u = parameter.value
    return [variables_opt, objval_opt]


''' Helper functions'''
### Seite 37 (im Text zwischen (4.12) und (4.13))
def _N(n):
    return np.diag([p for p in xrange(-n, 0, 1)])

### Seite 55; (4.64)
def _M(n):
    return np.diag([p for p in xrange(0, n, 1)])

### Seite 55; (4.64)
def _P(l, k, n):
    I = np.eye(n)
    Mn = _M(n)
    P = I
    if k == 0:
        pass # P=I
    else:
        for q in xrange(0, k, 1):
            P = P * ((l-q)*I + Mn)
    return np.matrix(P)
