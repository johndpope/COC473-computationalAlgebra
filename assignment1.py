from numpy.linalg import inv
import numpy as np
from math import *


tolerance = 0.0001


def testFuction1(x,derivative = 0):
	# root found on Wolfram Alpha = 277 or -277
	# returns function value at x if derivative = 0 (default)
	# returns derivate of f ate x if dirivative = 1 
	g = 9.806
	k = 0.00341
	functionValue = log(cosh(x*sqrt(g*k))) - 50
	derivativeValue = sqrt(g*k)*tanh(x*sqrt(g*k))
	if (derivative == 1):
		return derivativeValue
	return functionValue

def testFuction2(x,derivative = 0):
	# root found on Wolfram Alpha = [-14.1,-10.9,-7.85,-4.71,-1.55,0.597]
	functionValue = 4*cos(x)- exp(2*x)
	derivativeValue = -2*(exp(2*x)+2*sin(x))
	if (derivative == 1):
		return derivativeValue
	return functionValue

## Bissection Method
def bissectionMethod(function,a,b,tol):
	iterations = 0
	while(abs(a-b) > tol and iterations <100):
		rootPoint = (a+b)/2
		rootValue = function(rootPoint)
		if(rootValue > 0):
			b = rootPoint 
		else:
			a = rootPoint 
		iterations+= 1
	return rootPoint			
# print(bissectionMethod(testFuction2,-30.0,10.0,tolerance ))
# print(bissectionMethod(testFuction1,-1000.0,1000.0,tolerance ))

## Newton Method
## TO DO  - Metodo esta convergindo para valores que nao deveria (eg -50)
def newtonMetod(function,initialValue,tol):
	iterations = 1000
	rootPoint = initialValue
	for i in range (iterations):
		lastRoot = rootPoint  
		rootPoint  = rootPoint - (function(rootPoint)/function(rootPoint,1))
	if(abs(rootPoint - lastRoot) < tol):
		return rootPoint
	else:
		return "convergence not reached"

# print(newtonMetod(testFuction1,200,tolerance))
# print(newtonMetod(testFuction2,15,tolerance))
# print(newtonMetod(testFuction2,-10,tolerance))
# print(newtonMetod(testFuction2,-15,tolerance))


## Secant Method

def secantMethod(function,initialValue,tol):
	iterations = 1000
	rootPoint = initialValue
	delta = 0.001
	lastRoot = rootPoint
	rootPoint = lastRoot + delta
	fa = function(lastRoot)
	for i in range (iterations):
		fi = function(rootPoint	)
		rootPoint = rootPoint - (fi * (rootPoint-lastRoot))/(fi-fa)
		if(abs(rootPoint - lastRoot) < tol):
			return rootPoint
		else:
			fa = fi
	return "convergence not reached"

# print(newtonMetod(testFuction1,200,tolerance))
# print(newtonMetod(testFuction2,15,tolerance))
# print(newtonMetod(testFuction2,-10,tolerance))
# print(newtonMetod(testFuction2,-15,tolerance))



## Inverse Interpolation Method

def inverseInterpolationMethod(function,x1,x2,x3,tol):
	lastRoot = pow(10,36)
	iterations = 1000
	x = [x1,x2,x3]
	y = [0,0,0]
	for i in range (iterations):
		y[0],y[1],y[2] = function(x[0]),function(x[1]),function(x[2])
		rootPoint = (y[1]*y[2]*x[0])/((y[0]-y[1])*(y[0]-y[2])) + (y[0]*y[2]*x[1])/((y[1]-y[0])*(y[1]-y[2])) + (y[0]*y[1]*x[2])/((y[2]-y[0])*(y[2]-y[1]))
		if(abs(rootPoint-lastRoot) < tol):
			return rootPoint
		else:
			i = y.index(max(y))
			x[i] = rootPoint
			y[i] = function(x[i])
			lastRoot = rootPoint
	return "Convergence not reached"

# print(inverseInterpolationMethod(testFuction1,200,250,260,tolerance))

## Multi dimensional systems

# def multiDimensionalNewtonMethod(system):
# 	iterations = 1000
# 	systemDimension = system.shape
# 	lastX = np.random.random((systemDimension[0],systemDimension[1]))*10
# 	for i in range(iterations):



## Adjust of Nonlinear Functions

def jacobian(x):
	

a = np.array([(1.5,2,3), (4,5,6)])
b = np.random.random((2,3))*10

p = np.poly1d([2,3, 1, 2])
print(p)
q = p.deriv()
print (q)
print(q(1))