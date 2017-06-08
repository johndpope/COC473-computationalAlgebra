from numpy.linalg import inv
import numpy as np
from math import *
from copy import deepcopy
from sympy import *

x,y,z = symbols('x y z')
init_printing(use_unicode=True)

tolerance = 0.0001
g = 9.806
k = 0.00341
testFunction1 = log(cosh(x*sqrt(g*k))) - 50
testFunction2 = 4*cos(x)- exp(2*x)


# Support Functions

def jacobian(funcArray):
	global symbolsList
	functionArray = deepcopy(funcArray)
	jacobian = []
	for i in range(functionArray.size):
		temp = []
		for j in range(len(symbolsList)):
			temp.append(diff(functionArray[i], symbolsList[j]))
		jacobian.append(temp)
	return jacobian


def changeValuesMatrix(matrix, valueArray):
	global symbolsList
	functionMatrix = deepcopy(matrix)
	for i in range(len(functionMatrix)):
		for j in range(len(functionMatrix[i])):
			for k in range(len(symbolsList)):
				functionMatrix[i][j] = functionMatrix[i][j].subs(symbolsList[k], valueArray[k])
	return functionMatrix


def changeValuesArray(array, valueArray):
	global symbolsList
	functionArray = deepcopy(array)
	for i in range(len(functionArray)):
		for k in range(len(symbolsList)):
				functionArray[i] = functionArray[i].subs(symbolsList[k], valueArray[k])
	return functionArray


def changeValuesArrayBroyden(array, valueArray):
	global symbolsList
	functionArray = deepcopy(array)
	for i in range(len(functionArray)):
		for j in range(len(functionArray[i])):
			for k in range(len(symbolsList)):
				functionArray[i][j] = functionArray[i][j].subs(symbolsList[k], valueArray[k])
	return functionArray




def bissectionMethod(function,a,b,tol):
	"""
	Implements: Bissection Method
	Arguments:
		function: function to apply the method (sympy function)
		a: begin of the analyzed interval (number)
		b: end of the analyzed interval (number)
		tol: error tolerance (number)
	Return: root of analyzed function (number) 
	"""
	print("## Running Bissection Method\n")
	iterations = 0
	while(abs(a-b) > tol and iterations < 1000):
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




#TODO  - Metodo esta convergindo para valores que nao deveria (eg -50)
def newtonMetod(function,initialValue,tol):
	"""
	Implements: Newthon Method
	Arguments:
		function: function to apply the method (sympy function)
		initialValue: starting point (number)
		tol: error tolerance (number)
	Return: root of analyzed function (number) 
	"""
	print("## Running Newton Method\n")
	iterations = 100
	rootPoint = initialValue
	for i in range (iterations):
		lastRoot = rootPoint  
		functionDerivative =  diff(function,x)
		rootPoint  = rootPoint - function.evalf(subs={x:rootPoint})/functionDerivative.evalf(subs={x:rootPoint})
	if(abs(rootPoint - lastRoot) < tol):
		return rootPoint
	else:
		return "convergence not reached"

print("Function 1 with initial point 200: ",newtonMetod(testFunction1,200,tolerance))
print("Function 1 with initial point -200: ",newtonMetod(testFunction1,-200,tolerance))
print("Function 2 with initial point -10: ",newtonMetod(testFunction2,-10,tolerance))
print("Function 2 with initial point -15: ",newtonMetod(testFunction2,-15,tolerance))
print("Function 2 with initial point -15: ",newtonMetod(testFunction2,-1,tolerance))


def secantMethod(function,initialValue,tol):
	"""
	Implements: Secant Method
	Arguments:
		function: function to apply the method (sympy function)
		initialValue: starting point (number)
		tol: error tolerance (number)
	Return: root of analyzed function (number) 
	"""
	print("## Running Secant Method\n")
	iterations = 1000
	rootPoint = initialValue
	delta = 0.001
	lastRoot = rootPoint
	rootPoint = lastRoot + delta
	fa = function(lastRoot)
	for i in range (iterations):
		fi = function.subs(x,rootPoint)
		# The line below is the main difference to Newton Method. We calculate the derivative using the definition
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




def inverseInterpolationMethod(function,x1,x2,x3,tol):
	"""
	Implements: Inverse Interpolation Method
	Arguments:
		function: function to apply the method (sympy function)
		x1,x2,x3: points to use in Langrange`s interpolation (number)
		tol: error tolerance (number)
	Return: root of analyzed function (number) 
	"""
	print("## Running Inverse Interpolation Method\n")
	# We start x* (root of last iteration) as 10^36
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




## 4 Multi Dimensional Newton Method
def multiDimensionalNewtonMethod(functionArray, X0, tol, symbols):
	print("## Running MultiDimensional Newton Method\n")
	iterations = 15000
	jacob = jacobian(functionArray)
	lastX = X0

	# define symbols that will be used
	global symbolsList
	symbolsList = symbols

	for i in range(iterations):
		j = changeValuesMatrix(jacob, lastX)
		f = changeValuesArray(functionArray, lastX)

		j_np = np.array(j).astype(np.float64)
		f_np = np.array(f).astype(np.float64)
		
		deltaX = -np.dot(inv(j_np),f_np)
		lastX = lastX + deltaX

		tolk = np.linalg.norm(deltaX, ord=2) / np.linalg.norm(lastX, ord=2)
		if (tolk < tol):
			return lastX

	return "Convergence not reached"

functionArray1 = np.array([16*(x**4)+16*(y**4)+(z**4)-16, (x**2)+(y**2)+(z**2)-3, (x**3)-y+z-1])
functionArray2 = np.array([x+2*y-2,(x**2)+4*(y**2)-4])
# print(multiDimensionalNewtonMethod(functionArray1, [1,1,1], tolerance, [x,y,z]))
# t = [ 0.79040954,  0.80688815,  1.31308198]

# print(16*(t[0]**4)+16*(t[1]**4)+(t[2]**4))




## 4 Multi Dimensional Broyden Method
def multiDimensionalBroydenMethod(functionArray, X0, B0, tol, symbols):
	print("## Running MultiDimensional Broyden Method\n")
	iterations = 15000
	X_list = [X0]
	B_list = [B0] 		# receive the jacobian of start

	# define symbols that will be used
	global symbolsList
	symbolsList = symbols

	for i in range(1, iterations+1):
		j_np = B_list[i-1]

		f_ant = changeValuesArrayBroyden(functionArray, X_list[i-1])

		f_np_ant = np.array(f_ant).astype(np.float64)
		
		deltaX = -np.dot(inv(j_np),f_np_ant)

		X_list.append(np.add(X_list[i-1], deltaX.transpose())[0])

		f = changeValuesArrayBroyden(functionArray, X_list[i])
		f_np = np.array(f).astype(np.float64)

		lastY = f_np - f_np_ant

		#print("### iteraction: ",i," B: ",j_np, " f_ant: ",f_ant, " f_np_ant: ",f_np_ant, " xlist -1: ", X_list[i-1], "deltaX: ",deltaX)

		tolk = np.linalg.norm(deltaX, ord=2) / np.linalg.norm(X_list[i], ord=2)
		if (tolk < tol):
			return X_list[i]
		else:
			deltaX_transp = deltaX.transpose()

			#print("### iteraction: ",i," deltaX_transp",deltaX_transp)

		b_dividendo = np.dot(lastY-np.dot(np.asmatrix(B_list[i-1]), deltaX), deltaX_transp)
		b_divisor = np.dot(deltaX_transp, deltaX)

		next_B = B_list[i-1] + np.divide(b_dividendo, b_divisor)

		#print("### nextB: ", next_B.tolist())

		B_list.append(next_B.tolist())

		print("### iteraction: ",i," ---- VALUE",X_list[i])

	return "Convergence not reached"

functionArray3 = np.array([[x+2*y-2], [(x**2)+4*(y**2)-4]])
b0_broyden = np.array([[1,2],[4,24]])
print(multiDimensionalBroydenMethod(functionArray3, [2,3], b0_broyden, tolerance, [x,y]))

