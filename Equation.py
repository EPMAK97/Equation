from matplotlib import pyplot as plt
import numpy
from numpy import log10, sqrt, power

from decimal import *
getcontext().prec = 100

E = Decimal(0.000001)

Q = Decimal(1.425)
d = Decimal(1.067)
v_n0 = Decimal(0.0000125)
v_n = v_n0
v_p = Decimal(0.0017)

eps = Decimal(0.0001)
V_por = Decimal(0.05)

Re = Decimal(4.0) * Q / (Decimal(numpy.pi) * d * v_n)

C_LEFT_BOUND = Decimal(0.005)
C_RIGHT_BOUND = Decimal(0.012)
C_STEP = Decimal(0.0001)
CC = C_STEP

LAMBDA_STEP = Decimal(0.001)
LAMBDA_LEFT_BOUND = Decimal(LAMBDA_STEP)
LAMBDA_RIGHT_BOUND = Decimal(1.0)

def left_log(x):
    x = Decimal(x)
    #return math.pow(((2.8 * V_por) / (v_n * sqrt(x))), (1000 * CC / 5.75))
    tmp1 = ((Decimal(2.8) * Decimal(power(V_por, Decimal(2.0)))) / (v_n * Decimal(sqrt(x))));
    tmp2 = (Decimal(1000) * CC / Decimal(5.75));
    return Decimal(power(tmp1, tmp2))

def right_log(x):
    return Decimal(2.51) / (Re * Decimal(sqrt(x))) + eps / Decimal(3.701)

def right_side(x):
    return Decimal(-2.0) * Decimal(log10(left_log(x) * right_log(x)))
    #return math.pow(left_log(x) * right_log(x), -2.0)

def left_side(x):
    return Decimal(1) / Decimal(sqrt(x))
    #return math.pow(10.0, 1 / sqrt(x));

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

    print("in binsearch")
    if coeff * coeffR < 0:
    	print("normal binsearch")
    else:
    	print("can't use binsearch!!! using brute force")
    	#return 1.0
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
    global CC, v_n
    for CC in numpy.arange(C_LEFT_BOUND, C_RIGHT_BOUND, C_STEP):
        print(CC)
        v_n = (1 - CC) * v_n + 1 * CC * v_p
        CC_array.append(CC)
        #bfl = brute_force_lambda()
        bsl = bin_search_lambda()
        #if abs(bfl - bsl) > 1e-2:
        lambda_array.append(bsl)
        #lambda_array.append(bfl)
        #lambda_array.append(brute_force_lambda())

    print_critical_point()
    plt.plot(CC_array, lambda_array, marker='o', ls='-')
    plt.xlabel('Bulk polymer concentration')
    plt.ylabel('Friction coeffitient')
    plt.show()


def print_critical_point():
    cidx = numpy.array(lambda_array).argmax()
    print('lambda_critical = {0}'.format(lambda_array[cidx]))
    print('c_critical = {0}'.format(CC_array[cidx]));

    for i in range(len(lambda_array)):
        if lambda_array[i] < lambda_array[0]:
            #print('first lambda = {0}'.format(lambda_array[0]))
            print('first concentration with lambda less then start_lambda = {0}'.format(CC_array[i]))
            return

    print('there is no lambda less then start_lambda')

main()


# 1/x+2*log10[ ((2.8*0.05)/(1.58975e-05* sqrt[x]))^(0.2/5.75) * (2.51/(136034.965417*sqrt[x])+0.0001/3.701)]

# 1/x + 2*math.log10(math.pow(((2.8 * V_por) / (v_n * math.sqrt(x))), 1000 * CC) * math.log10(2.51 / (Re * math.sqrt(x)) + eps / 3.701))

# -2.0 * math.log10(math.pow(((2.8 * V_por) / (v_n * math.sqrt(x))), 1000 * CC) * math.log10(2.51 / (Re * math.sqrt(x)) + eps / 3.701))

# gnu multiplay precision
