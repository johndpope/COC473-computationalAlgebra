from numpy.linalg import inv
import numpy as np
from math import *
from copy import deepcopy
from sympy import *
x,y,z = symbols('x y z')
init_printing(use_unicode=True)

def polynomialIntegration(numberOfpoints,interval,function):
	# If numberOfPoints is less than 5, calculates the derivative using the given interval
	# Althought, if numberOfPoints is bigger than 5, converts the interval to the Legendre-Gauss interval [-1,1] 
	delta = (interval[1]-interval[0])/(numberOfpoints-1)
	points =[]

	# Points is a list of lists
	# Each list contain the points to use when the user set N points

	for n in range(0,5):
		nPoints = []
		for n in range(0,5):
			nPoints.append(interval[0] + n*delta)
		points.append(nPoints)

	# Points extract from https://pomax.github.io/bezierinfo/legendre-gauss.html
	points.append([0.6612093864662645,-0.6612093864662645,-0.2386191860831969,0.2386191860831969,-0.9324695142031521,0.9324695142031521])
	points.append([0,0.4058451513773972,-0.4058451513773972,-0.7415311855993945,0.7415311855993945,-0.9491079123427585,0.9491079123427585])
	points.append([-0.1834346424956498,0.1834346424956498,-0.5255324099163290,0.5255324099163290,-0.7966664774136267,0.7966664774136267,-0.9602898564975363,0.9602898564975363])
	points.append([-0.8360311073266358,0.8360311073266358,-0.9681602395076261,0.9681602395076261,-0.3242534234038089,0.3242534234038089,-0.6133714327005904,0.6133714327005904])
	points.append([-0.1488743389816312,0.1488743389816312,-0.4333953941292472,0.4333953941292472,-0.6794095682990244,0.6794095682990244,-0.8650633666889845,0.8650633666889845,-0.9739065285171717,0.9739065285171717])
	
	elementIndex = numberOfpoints-1

	print("PONTOS: ",points[elementIndex])
	l = points[elementIndex][-1] - points[elementIndex][0]
	result = 0

	# Weights is a list of lists
	# Each list contain the weight to use when the user set N points
	# Weights extract from https://pomax.github.io/bezierinfo/legendre-gauss.html

	weights = [[l],
			  [l/2,l/2],
			  [l/2,(2*l)/3,l/6],
			  [l/8,(3*l)/8,(3*l)/8,l/8],
			  [(7*l)/90,(16*l)/45,(2*l)/15,(16*l)/45,(7*l/90)],
			  [0.3607615730481386,0.3607615730481386,0.4679139345726910,0.4679139345726910,0.1713244923791704,0.1713244923791704],
			  [0.4179591836734694,0.3818300505051189,0.3818300505051189,0.3818300505051189,0.2797053914892766,0.1294849661688697,0.1294849661688697],
			  [0.3626837833783620,0.3626837833783620,0.3137066458778873,0.3137066458778873,0.2223810344533745,0.2223810344533745,0.1012285362903763,0.1012285362903763],
			  [0.3302393550012598,0.1806481606948574,0.1806481606948574,0.0812743883615744,0.0812743883615744,0.3123470770400029,0.3123470770400029,0.2606106964029354,0.2606106964029354],
			  [0.2955242247147529,0.2955242247147529,0.2692667193099963,0.2692667193099963,0.2190863625159820,0.2190863625159820,0.1494513491505806,0.1494513491505806,0.0666713443086881,0.0666713443086881]]
	selectWeightArray = weights[len(points[elementIndex])-1]
	result = 0;
	for i in range(len(points[elementIndex])):
		if(numberOfpoints > 5):
			adjustTerm = (interval[-1]-interval[0])/2
			functionValue = function.subs(x,adjustTerm*points[elementIndex][i] + (interval[-1]+interval[0])/2)
			weight = selectWeightArray[i]
			result += adjustTerm*(functionValue * weight)
		else:
			functionValue = function.subs(x,points[elementIndex][i])
			weight = selectWeightArray[i]
			result += functionValue * weight

		print ("RESULTADO PARA ",i+1," PONTOS: ",result)
	return result

testFunction1 = exp(-x**2)
interval =[0,1]
print(polynomialIntegration(10,[0,1],testFunction1))
