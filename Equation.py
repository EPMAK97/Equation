from matplotlib import pyplot as plt
import numpy
from numpy import log10, sqrt, power

from decimal import *
getcontext().prec = 50

# коэффициент при смешивании нефти и полимера
coeffMixture = Decimal(1)

# объемная скорость потока, м^3 / с
Q = Decimal(1.425)

# диаметр трубы, м
d = Decimal(1.067)

# КИНЕМАТИЧЕСКАЯ вязкость нефти, м^2 / с
v_n0 = Decimal(0.0000125)
v_n = v_n0

# КИНЕМАТИЧЕСКАЯ вязкость присадок, м^2 / с
v_p = Decimal(0.0017)

# какой-то епсилон, м
eps = Decimal(0.0001)

# пороговая скорость, м / с
V_por = Decimal(0.05)

# средняя скорость движения жидкости в трубе, м / с
V = Decimal(4.0) * Q / (Decimal(numpy.pi) * d * d)

# число Рейндольса, безразмерное
Re = V * d / v_n

C_LEFT_BOUND = Decimal(0.000)
C_RIGHT_BOUND = Decimal(0.015)
C_STEP = Decimal(0.0001)
CC = C_STEP

LAMBDA_STEP = Decimal(0.001)
LAMBDA_LEFT_BOUND = Decimal(LAMBDA_STEP)
LAMBDA_RIGHT_BOUND = Decimal(1.0)

def left_log(x):
    x = Decimal(x)
    tmp1 = (Decimal(2.8) * V_por) / (V * Decimal(sqrt(x)));
    tmp2 = (Decimal(1000) * CC / Decimal(5.75));
    return Decimal(power(tmp1, tmp2))

def right_log(x):
    return Decimal(2.51) / (Re * Decimal(sqrt(x))) + eps / (Decimal(3.7) * d)

def right_side(x):
    return Decimal(-2.0) * Decimal(log10(left_log(x) * right_log(x)))

def left_side(x):
    return Decimal(1) / Decimal(sqrt(x))

def brute_force_lambda():
    min_diff = 1e10
    min_lambda = -1;
    for i in numpy.arange(LAMBDA_LEFT_BOUND, LAMBDA_RIGHT_BOUND, LAMBDA_STEP):
        cur_diff = abs(left_side(i) - right_side(i))
        if (min_diff > cur_diff):
            min_diff = cur_diff
            min_lambda = i
    return min_lambda

def bin_search_lambda():
    L, R = LAMBDA_STEP, LAMBDA_RIGHT_BOUND
    coeff = right_side(L) - left_side(L)
    coeffR = right_side(R) - left_side(R)
    iterations = 0

    if coeff * coeffR < 0:
    	1
    	#print("using binsearch")
    else:
    	#print("can't use binsearch!!! using brute force")
    	return brute_force_lambda()

    while iterations < 100 and abs(L - R) > 1e-10:
        iterations += 1
        mid = (L + R) / Decimal(2.0);

        if coeff > 0:
            #print('coeff greater')
            if left_side(mid) - right_side(mid) > 0.0:
                R = mid
            else:
                L = mid
        else:
            #print('coeff less')
            if left_side(mid) - right_side(mid) < 0.0:
                R = mid
            else:
                L = mid
    return L

CC_array = []
lambda_array = []

def main():
    global CC, v_n, Re
    for CC in numpy.arange(C_LEFT_BOUND, C_RIGHT_BOUND, C_STEP):
        print(CC)
        v_n = (1 - CC / 100) * v_n + coeffMixture * CC / 100 * v_p
        Re = V * d / v_n
        CC_array.append(CC)
        lambda_array.append(bin_search_lambda())

    print_critical_point()
    plt.plot(CC_array, lambda_array, marker='o', ls='-')
    plt.xlabel('Объемная концентрация полимера')
    plt.ylabel('Коэффициент трения')
    plt.show()

def print_critical_point():
    cidx = numpy.array(lambda_array).argmax()
    print('lambda_critical = {0}'.format(lambda_array[cidx]))
    print('c_critical = {0}'.format(CC_array[cidx]))
    print('Re = {0}'.format(Re))

    for i in range(len(lambda_array)):
        if lambda_array[i] < lambda_array[0]:
            print('first lambda = {0}'.format(lambda_array[0]))
            print('first concentration with lambda less then start_lambda = {0}'.format(CC_array[i]))
            return

    print('there is no lambda less then start_lambda')

main()
