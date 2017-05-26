import numpy as np
from math import *
from sympy import *
x,y,z,w = symbols('x y z w')
init_printing(use_unicode=True)

def polynomialIntegration(numberOfpoints,interval,function,IntegrationTerm,show = False):
	# If numberOfPoints is less than 5, calculates the derivative using the given interval
	# Although, if numberOfPoints is bigger than 5, converts the interval to the Legendre-Gauss interval [-1,1] 
	delta = (interval[1]-interval[0])/(numberOfpoints-1)
	points =[]

	# Points is a list of lists
	# Each list contain the points to use when the user set N points

	for n in range(0,5):
		nPoints = []
		for n in range(0,5):
			nPoints.append(interval[0] + n*delta)
		points.append(nPoints)

	# Points [6,10] extracted from https://pomax.github.io/bezierinfo/legendre-gauss.html

	points.append([0.6612093864662645,-0.6612093864662645,-0.2386191860831969,0.2386191860831969,-0.9324695142031521,0.9324695142031521])
	points.append([0,0.4058451513773972,-0.4058451513773972,-0.7415311855993945,0.7415311855993945,-0.9491079123427585,0.9491079123427585])
	points.append([-0.1834346424956498,0.1834346424956498,-0.5255324099163290,0.5255324099163290,-0.7966664774136267,0.7966664774136267,-0.9602898564975363,0.9602898564975363])
	points.append([-0.8360311073266358,0.8360311073266358,-0.9681602395076261,0.9681602395076261,-0.3242534234038089,0.3242534234038089,-0.6133714327005904,0.6133714327005904])
	points.append([-0.1488743389816312,0.1488743389816312,-0.4333953941292472,0.4333953941292472,-0.6794095682990244,0.6794095682990244,-0.8650633666889845,0.8650633666889845,-0.9739065285171717,0.9739065285171717])
	
	elementIndex = numberOfpoints-1
	if(show):
		print("PONTOS: ",points[elementIndex])

	l = points[elementIndex][-1] - points[elementIndex][0]
	# Weights is a list of lists
	# Each list contain the weight to use when the user set N points
	# Weights extracted from https://pomax.github.io/bezierinfo/legendre-gauss.html

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
			# Adjusts the function to work using the weights and points of Legendre-Gauss interval [-1,1]
			functionValue = function.subs(IntegrationTerm,adjustTerm*points[elementIndex][i] + (interval[-1]+interval[0])/2)
			weight = selectWeightArray[i]
			result += adjustTerm*(functionValue * weight)
		else:
			functionValue = function.subs(IntegrationTerm,points[elementIndex][i])
			weight = selectWeightArray[i]
			result += functionValue * weight
		if(show):	
			print ("RESULTADO PARA ",i+1," PONTOS: ",result)
	return result

# Question 2

testFunction1 = exp(-x**2)
testFunction2 = (1/2)*np.pi*exp((-1/2)*(x**2))
N = 10
print("RESULTADO QUESTAO 2 - FUNCAO 1 (10 PONTOS): ",polynomialIntegration(10,[0,1],testFunction1,x))
print("RESULTADO QUESTAO 2 - FUNCAO 2 (10 PONTOS): ",polynomialIntegration(10,[0,1],testFunction2,x),"\n")

# Question 3
Wn = 1
E = 0.05
Sn = 2
RAO = 1/(sqrt((1-(w/Wn)**2)**2)+(2*E*(w/Wn))**2)
testFunction3 = (RAO**2)*Sn
testFunction4 = (w**2)*(RAO**2)*Sn

m0 = polynomialIntegration(6,[0,10],testFunction3,w)
m2 = polynomialIntegration(10,[0,10],testFunction4,w)

print("RESUTADO QUESTAO 3 - m0: ",m0)
print("RESUTADO QUESTAO 3 - m2: ",m2,"\n")

# Question 4
# Change the Sn value used at question 3. 

Hs = 3
Tz = 5
Sn = ((4 * (np.pi)**3 * (Hs)**2)/(w**5 * Tz**4))*exp((-16*np.pi**3)/(w**4 * Tz**4))

m0_2 = polynomialIntegration(10,[0,10],testFunction3,w)
m2_2 = polynomialIntegration(10,[0,10],testFunction4,w)

print("RESULTADO QUESTAO 4 m0: ",m0_2)
print("RESULTADO QUESTAO 4 m2: ",m2_2,"\n")

# Question 5

testFunction5 = 2 + 2*x - x**2 + 3*x**3
aValues = [None]
ponto = 2
while(ponto <= 10):
	A = polynomialIntegration(ponto,[0,4],testFunction5,x)
	if(A == aValues[-1]):
		print("RESULTADO QUESTAO 5 - INTEGRACAO POLINOMIAL ","VALOR DE A: ", A, "OBTIDO COM ",ponto," PONTOS")
		break;
	else:
		aValues.append(A)
		ponto += 1


#TODO FIX IT
W = [0.5555555555555,0.888888889,0.5555555555555]
Z = [0.774596,0,-0.774596]
L = 4
# Uses two points of integration
x1 = 1/2*(0+4+Z[0]*L)
x2 = 1/2*(0+4+Z[1]*L)
f_x1 = testFunction5.subs(x,x1)
f_x2 = testFunction5.subs(x,x2)
A = (L/2)*(f_x1*W[0]+f_x2*W[1])
print("RESULTADO QUESTAO 5 - QUADRATURA DE GAUSS ","VALOR DE A: ", A, "OBTIDO COM ",2," PONTOS","\n")

# Question 6
print("RESULTADOS QUESTAO 6","\n")
testFunction6 = 1/(1+x**2)
interval = [0,3]

# Result using mid point rule
m = (interval[1]+interval[0])/2
A_mid = testFunction6.subs(x,m) * (interval[1]-interval[0])
print("RESULTADO DA TECNICA DO PONTO MEDIO: ", A_mid)

# Result using trapeze rule
A_t = ((testFunction6.subs(x,interval[0])+testFunction6.subs(x,interval[1]))/2.0) * (interval[1]-interval[0])
print("RESULTADO DA TECNICA DO TRAPEZIO: ", A_t)

# Result using numeric method
A = polynomialIntegration(ponto,interval,testFunction6,x)
print("RESULTADO DA RESOLUCAO NUMERICA: ", A)
