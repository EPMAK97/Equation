import math
from matplotlib import pyplot as plt
import numpy
from math import log10, sqrt

E = 0.000001

Q = 1.425
d = 1.067
v_n = 0.0000125
v_p = 0.0017

eps = 0.0001
V_por = 0.01

Re = (4 * Q) / (math.pi * d * v_n)

C_LEFT_BOUND = 0
C_RIGHT_BOUND = 1
C_STEP = 0.001
CC = C_STEP

LAMBDA_STEP = 0.0001
LAMBDA_LEFT_BOUND = LAMBDA_STEP
LAMBDA_RIGHT_BOUND = 0.25

def left_log(x):
    #return math.pow(((2.8 * V_por) / (v_n * sqrt(x))), (1000 * CC / 5.75))
    return math.pow(((2.8 * math.pow(V_por, 2.0)) / (v_n * sqrt(x))), (1000 * CC / 5.75))

def right_log(x):
    return 2.51 / (Re * sqrt(x)) + eps / 3.701

def right_side(x):
    return -2 * log10(left_log(x) * right_log(x))

def left_side(x):
    #return 1 / x
    return 1 / sqrt(x)

def brute_force_lambda():
    min_diff = 1e10
    min_lambda = -1;
    for i in numpy.arange(LAMBDA_LEFT_BOUND, LAMBDA_RIGHT_BOUND, LAMBDA_STEP):
        cur_diff = abs(left_side(i) - right_side(i))
        if (min_diff > cur_diff):
            min_diff = cur_diff
            min_lambda = i
    return min_lambda

CC_array = []
lambda_array = []

def main():
    global CC, v_n
    for CC in numpy.arange(C_LEFT_BOUND, C_RIGHT_BOUND, C_STEP):
        v_n = (1 - CC) * v_n + 1 * CC * v_p
        CC_array.append(CC)
        lambda_array.append(brute_force_lambda())

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
