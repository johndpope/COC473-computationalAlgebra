from numpy.linalg import inv
import numpy as np
from math import *
from sympy import *
x,y,z = symbols('x y z')
init_printing(use_unicode=True)

tolerance = 0.0001
g = 9.806
k = 0.00341
testFunction1 = log(cosh(x*sqrt(g*k))) - 50
testFunction2 = 4*cos(x)- exp(2*x)


## Bissection Method
def bissectionMethod(function,a,b,tol):
	iterations = 0
	while(abs(a-b) > tol and iterations <100):
		rootPoint = (a+b)/2
		rootValue = function.subs(x,rootPoint)
		if(rootValue > 0):
			b = rootPoint 
		else:
			a = rootPoint 
		iterations+= 1
	return rootPoint			
#print(bissectionMethod(testFunction2,-30.0,10.0,tolerance))
#print(bissectionMethod(testFunction1,-1000.0,1000.0,tolerance))

## Newton Method
## TO DO  - Metodo esta convergindo para valores que nao deveria (eg -50)
def newtonMetod(function,initialValue,tol):
	iterations = 1000
	rootPoint = initialValue
	for i in range (iterations):
		lastRoot = rootPoint  
		rootPoint  = rootPoint - (function.subs(x,rootPoint)/diff(function,x).subs(x,rootPoint))
	if(abs(rootPoint - lastRoot) < tol):
		return rootPoint
	else:
		return "convergence not reached"

#print(newtonMetod(testFunction1,200,tolerance))
#print(newtonMetod(testFuction2,15,tolerance))
#print(newtonMetod(testFuction2,-10,tolerance))
#print(newtonMetod(testFuction2,-15,tolerance))


## Secant Method

def secantMethod(function,initialValue,tol):
	iterations = 1000
	rootPoint = initialValue
	delta = 0.001
	lastRoot = rootPoint
	rootPoint = lastRoot + delta
	fa = function(lastRoot)
	for i in range (iterations):
		fi = function.subs(x,rootPoint)
		rootPoint = rootPoint - (fi * (rootPoint-lastRoot))/(fi-fa)
		if(abs(rootPoint - lastRoot) < tol):
			return rootPoint
		else:
			fa = fi
	return "convergence not reached"

#print(newtonMetod(testFuction1,200,tolerance))
#print(newtonMetod(testFuction2,15,tolerance))
#print(newtonMetod(testFuction2,-10,tolerance))
#print(newtonMetod(testFuction2,-15,tolerance))



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
			y[i] = function.subs(x,x[i])
			lastRoot = rootPoint
	return "Convergence not reached"

#print(inverseInterpolationMethod(testFuction1,200,250,260,tolerance))

## Multi dimensional systems
## Adjust of Nonlinear Functions
symbolsList = [x,y]

def jacobian(functionArray):
	jacobian = []
	 
	# TODO: check if need to get symbolslist from functionArray
	# symbolsSet = []
	# for i in functionArray:
	# 	s = i.free_symbols
	# 	if len(s) > len(symbolsSet):
	# 		symbolsSet = s

	# symbolsList = list(symbolsSet)

	global symbolsList
	
	for i in range(functionArray.size):
		temp = []
		for j in range(len(symbolsList)):
			temp.append(diff(functionArray[i], symbolsList[j]))
		jacobian.append(temp)
	return jacobian

def changeValuesMatrix(matrix, valueArray):
	global symbolsList
	functionMatrix = matrix[:]
	for i in range(len(functionMatrix)):
		for j in range(len(functionMatrix[i])):
			for k in range(len(symbolsList)):
				functionMatrix[i][j] = functionMatrix[i][j].subs(symbolsList[k], valueArray[k])
	return functionMatrix

def changeValuesArray(functionArray, valueArray):
	global symbolsList
	for i in range(len(functionArray)):
		for k in range(len(symbolsList)):
				functionArray[i] = functionArray[i].subs(symbolsList[k], valueArray[k])
	return functionArray

functionArray1 = np.array([16*(x**4)+16*(y**4)+(z**4)-16, (x**2)+(y**2)+(z**2)-3, (x**3)-y+z-1])
functionArray2 = np.array([x+2*y-2,(x**2)+4*(y**2)-4])
def multiDimensionalNewtonMethod(functionArray, X0):
	iterations = 7
	jacob = jacobian(functionArray)
	lastX = X0

	for i in range(iterations):
		print("jacob", jacob)
		j = changeValuesMatrix(jacob, lastX)
		f = changeValuesArray(functionArray, lastX)

		j_np = np.array(j).astype(np.float64)
		f_np = np.array(f).astype(np.float64)
		print(lastX)
		deltaX = -np.dot(inv(j_np),f_np)
		lastX = lastX + deltaX

		tolk = np.linalg.norm(deltaX, ord=2) / np.linalg.norm(lastX, ord=2)
		if (tolk < tolerance):
			return lastX


	return "Convergence not reached"


print(multiDimensionalNewtonMethod(functionArray2,[2,3]))
# t = [-1415.66666667, -1415.66666667,  2834.33333333]

# print(16*(t[0]**4)+16*(t[1]**4)+(t[2]**4))
