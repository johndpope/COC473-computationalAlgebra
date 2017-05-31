from numpy.linalg import inv
import numpy as np
from math import *
from copy import deepcopy
from sympy import *
x,y,z,t = symbols('x y z t')

init_printing(use_unicode=True)

# define symbols that will be used
symbolsList = [t, y]

# Support Functions

def changeValuesInFuction(function, valueArray):
	global symbolsList
	func = deepcopy(function)
	for k in range(len(symbolsList)):
		func = func.subs(symbolsList[k], valueArray[k])
	return func



## a) Euler Method
def integrateEulerMethod(function, delta, X0, maximumT):
	X_list = [X0]
	t = 0
	index = 0

	while (t < maximumT):
		new_x = X_list[index] + delta * changeValuesInFuction(function, [t, X_list[index]])
		X_list.append(new_x)

		print("#### index: ", index+1, " new_x: ", new_x)

		index += 1
		t = index * delta



## b) Runge-Kutta second-order Method
def integrateRungeKutta2Method(function, delta, X0, maximumT):
	X_list = [X0]
	t = 0
	index = 0

	while (t < maximumT):

		K1 = changeValuesInFuction(function, [t, X_list[index]])
		K2 = changeValuesInFuction(function, [t + delta, X_list[index] + delta * K1])

		new_x = X_list[index] + delta/2 * (K1 + K2)
		X_list.append(new_x)

		print("#### index: ", index+1, " new_x: ", new_x)
		print("#### K1: ", K1, " ## K2: ", K2)

		index += 1
		t = index * delta



## c) Runge-Kutta fourth-order Method
def integrateRungeKutta4Method(function, delta, X0, maximumT):
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

		print("#### index: ", index+1, " new_x: ", new_x)
		print("#### K1: ", K1, " ## K2: ", K2, " ## K3: ", K3, " ## K4: ", K4)

		index += 1
		t = index * delta



def integrateMethod(function, delta, X0, maximumT, order):
	# order 1 = Euler Method
	# order 2 = Runge-Kutta second-order Method
	# order 4 = Runge-Kutta fourth-order Method

	if (order == 1):
		integrateEulerMethod(function, delta, X0, maximumT)
	elif (order == 2):
		integrateRungeKutta2Method(function, delta, X0, maximumT)
	elif (order == 4):
		integrateRungeKutta4Method(function, delta, X0, maximumT)
	else:
		print ("Order not defined. Please choose 1, 2 or 4.")



function = -2*t*y**2
y0 = 1

sampleFunction = t + y

integrateMethod(sampleFunction, 0.1, 0.0, 0.4, 1)

