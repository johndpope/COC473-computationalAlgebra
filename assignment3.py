import numpy as np
from math import *
from copy import deepcopy
from sympy import *
y, yy, t = symbols('y yy t')

init_printing(use_unicode=True)


# Support Functions

def changeValuesInFuction(function, valueArray):
	global symbolsList
	func = deepcopy(function)
	for k in range(len(symbolsList)):
		func = func.subs(symbolsList[k], valueArray[k])
	return func



## 1 a) Euler Method
def integrateEulerMethod(function, delta, X0, maximumT):
	print("## Running Euler Method\n")

	X_list = [X0]
	t = 0
	index = 0

	while (t < maximumT):
		new_x = X_list[index] + delta * changeValuesInFuction(function, [t, X_list[index]])
		X_list.append(new_x)

		index += 1
		t = index * delta

	return X_list



## 1 b) Runge-Kutta second-order Method
def integrateRungeKutta2Method(function, delta, X0, maximumT):
	print("## Running Runge-Kutta second-order Method\n")

	X_list = [X0]
	t = 0
	index = 0

	while (t < maximumT):

		K1 = changeValuesInFuction(function, [t, X_list[index]])
		K2 = changeValuesInFuction(function, [t + delta, X_list[index] + delta * K1])

		new_x = X_list[index] + delta/2 * (K1 + K2)
		X_list.append(new_x)

		# print("#### t: ", t, " ## K1: ", K1: ", K1, " ## K2: ", K2)

		index += 1
		t = index * delta

	return X_list



## 1 c) Runge-Kutta fourth-order Method
def integrateRungeKutta4Method(function, delta, X0, maximumT):
	print("## Running Runge-Kutta fourth-order Method\n")

	X_list = [X0]
	t = 0
	index = 0

	while (t < maximumT):

		K1 = changeValuesInFuction(function, [t, X_list[index]])
		K2 = changeValuesInFuction(function, [t + delta/2, X_list[index] + delta/2 * K1])
		K3 = changeValuesInFuction(function, [t + delta/2, X_list[index] + delta/2 * K2])
		K4 = changeValuesInFuction(function, [t + delta, X_list[index] + delta * K3])

		new_x = X_list[index] + delta/6 * (K1 + 2*K2 + 2*K3 + K4)
		X_list.append(new_x)

		# print("#### t: ", t, " ## K1: ", K1, " ## K2: ", K2, " ## K3: ", K3, " ## K4: ", K4)

		index += 1
		t = index * delta

	return X_list





## 2 a) Second order Taylor approximation
def taylorSerieAproximation(function, delta, X0, XX0, maximumT):
	print("## Running Second order Taylor approximation\n")

	X_list = [X0]
	XX_list = [XX0]
	XXX_list = []
	t = 0
	index = 0

	while (t < maximumT):

		# x'' = x''(t) = f(t, x, x')
		current_xxx = changeValuesInFuction(function, [t, X_list[index], XX_list[index]])
		XXX_list.append(current_xxx)

		# x
		new_x = X_list[index] + XX_list[index] * delta + XXX_list[index]/2 * delta**2
		X_list.append(new_x)

		# x'
		new_xx = XX_list[index] + XXX_list[index] * delta
		XX_list.append(new_xx)

		index += 1
		t = index * delta

	return X_list



## 2 b) Runge Kutta Nystrom Method
def rungeKuttaNystrom(function, delta, X0, XX0, maximumT):
	print("## Running Runge Kutta Nystrom Method\n")

	X_list = [X0]
	XX_list = [XX0]
	t = 0
	index = 0

	while (t < maximumT):

		K1 = delta/2 * changeValuesInFuction(function, [t, X_list[index], XX_list[index]])
		Q  = delta/2 * (XX_list[index] + K1/2)
		K2 = delta/2 * changeValuesInFuction(function, [t + delta/2, X_list[index] + Q, XX_list[index] + K1])
		K3 = delta/2 * changeValuesInFuction(function, [t + delta/2, X_list[index] + Q, XX_list[index] + K2])
		L  = delta * (XX_list[index] + K3)
		K4 = delta/2 * changeValuesInFuction(function, [t + delta, X_list[index] + L, XX_list[index] + 2 * K3])

		new_x = X_list[index] + delta * (XX_list[index] + 1/3 * (K1 + K2 + K3))
		
		new_xx = XX_list[index] + 1/3 * (K1 + 2*K2 + 2*K3 + K4)

		X_list.append(new_x)
		XX_list.append(new_xx)

		# print("#### t: ", t, " ## K1: ", K1: ", K1, " ## Q: ", Q, " ## K2: ", K2, " ## K3: ", K3, " ## L: ", L, " ## K4: ", K4)

		index += 1
		t = index * delta

	return X_list





## 1 ODE First-Order
def integrateMethod(function, delta, X0, maximumT, order):
	# order 1 = Euler Method
	# order 2 = Runge-Kutta second-order Method
	# order 4 = Runge-Kutta fourth-order Method

	# define symbols that will be used
	global symbolsList
	symbolsList = [t, y]

	if (order == 1):
		integrateEulerMethod(function, delta, X0, maximumT)
	elif (order == 2):
		integrateRungeKutta2Method(function, delta, X0, maximumT)
	elif (order == 4):
		integrateRungeKutta4Method(function, delta, X0, maximumT)
	else:
		print("Order not defined. Please choose 1, 2 or 4.")



## 2 ODE Second-Order
def solveSecondOrder(function, delta, X0, XX0, maximumT, method):
	# method 0 = Second order Taylor approximation
	# method 1 = Runge Kutta Nystrom Method

	# define symbols that will be used
	global symbolsList
	symbolsList = [t, y, yy]

	if (method == 0):
		taylorSerieAproximation(function, delta, X0, XX0, maximumT)
	elif (method == 1):
		rungeKuttaNystrom(function, delta, X0, XX0, maximumT)
	else:
		print("Method not defined. Please choose 0 or 1.")



function = -2*t*y**2
y0 = 1

sampleFunction = t + y

integrateMethod(sampleFunction, 0.1, 0.0, 0.4, 1)

sampleFunction2 = 2 * t + y + yy

taylorSerieAproximation(sampleFunction2, 0.1, 0, 1, 2)
rungeKuttaNystrom(sampleFunction2, 0.1, 0, 1, 2)

