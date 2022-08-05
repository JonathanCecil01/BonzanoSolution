# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import sympy as sp
import numpy as np
from sympy import Eq, degree

def bijection(expr, a, b):
    c = (a+b)/2
    fa = expr.evalf(subs = {x: a})
    fb = expr.evalf(subs = {x: b})
    fc = expr.evalf(subs = {x: c})
    #print(fc)
    while(fc>0.001 or fc<-0.001):
        if fa<0 and fb>0:
            if fc<0:
                a = c
            else:
                b = c
        elif fa>0 and fb<0:
            if fc>0:
                a = c
            else:
                b = c
        c = (a+b)/2
        fa = expr.evalf(subs = {x: a})
        fb = expr.evalf(subs = {x: b})
        fc = expr.evalf(subs = {x: c})  
    print("Bisection : ")
    #print(fc)
    print(c)
    
def regulaFalse(expr, a ,b):
    fa = expr.evalf(subs = {x: a})
    fb = expr.evalf(subs = {x: b})
    h = abs(fa*(b-a))/abs(fa+fb)
    x1 = a+h    
    fx = expr.evalf(subs = {x: x1})
    while(fx>0.0001 or fx<-0.0001):
        if fa<0 and fb>0:
            if fx<0:
                a = x1
            else:
                b = x1
        elif fa>0 and fb<0:
            if fx>0:
                a = x1
            else:
                b = x1
        h = abs(fa*(b-a))/abs(fa+fb)
        x1 = a+h    
        fx = expr.evalf(subs = {x: x1})
        fa = expr.evalf(subs = {x: a})
        fb = expr.evalf(subs = {x: b})
    print("Regula false : ")
    #print(fx)
    print(x1)
 
def NewtonRalph(expr, a):
    dx = sp.diff(expr, x)
    fx = expr.evalf(subs = {x: a})
    fdx =  dx.evalf(subs = {x: a})
    h = -1*fx/fdx
    while(fx>0.0001 or fx<-0.0001):
        fx = expr.evalf(subs = {x: a})
        fdx =  dx.evalf(subs = {x: a})
        h = -1*fx/fdx
        a = a+h
    print("Newton Raphson : ")
    #print(fx)
    print(a)

def fixedPoint(expr, a, b):
    if degree(expr)>0:
        deg = degree(expr)
        coeffi = expr.coeff(x**deg)
        expr = expr - coeffi*x**deg
        expr = ((expr/coeffi)**(1/deg))*-1
        dx = sp.diff(expr, x)
        print("Fixed point : ")
        if abs(dx.evalf(subs = {x:(b+a)/2}))>=1:
            print("Convergence not Saisfied")
            return
        fx = expr.evalf(subs = {x: a})-a
        #print(dx)
        while(fx>0.0001 or fx<-0.0001):
            fx = expr.evalf(subs = {x: a})-a
            a = expr.evalf(subs = {x: a})
        print(a)
        
   
x,y = sp.symbols('x y')
#expr = 2*x - sp.log(x)-7
expr = 3*x**1-sp.cos(x)-1
#expr = 10**x+x-4
bijection(expr,0.5,5)
regulaFalse(expr,0.5,5)
NewtonRalph(expr, 0.5)
fixedPoint(expr, 0, 1)

