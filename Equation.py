import numpy
import sys
from numpy import log10, sqrt, power
from matplotlib import pyplot as plt
from decimal import *
from PyQt5.QtWidgets import *

getcontext().prec = 30

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

def main():
    global CC, v_n, Re
    CC_array = []
    lambda_array = []
    for CC in numpy.arange(C_LEFT_BOUND, C_RIGHT_BOUND, C_STEP):
        print(CC)
        v_n = (1 - CC / Decimal(100)) * v_n + coeffMixture * CC / Decimal(100) * v_p
        Re = V * d / v_n
        CC_array.append(CC)
        lambda_array.append(bin_search_lambda())

    #print_critical_point()
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

class Example(QWidget):

    def __init__(self):
        super().__init__()

        self.initUI()


    def initUI(self):

        # Labels

        self.lbl1 = QLabel(self)
        self.lbl1.setText("Epsilon")
        self.lbl1.move(70, 44)

        self.lbl2 = QLabel(self)
        self.lbl2.setText("Q")
        self.lbl2.move(70, 64)

        self.lbl3 = QLabel(self)
        self.lbl3.setText("d")
        self.lbl3.move(70, 84)

        self.lbl4 = QLabel(self)
        self.lbl4.setText("v_p")
        self.lbl4.move(70, 104)

        self.lbl5 = QLabel(self)
        self.lbl5.setText("v_n0")
        self.lbl5.move(70, 124)

        self.lbl6 = QLabel(self)
        self.lbl6.setText("V_por")
        self.lbl6.move(70, 144)

        self.lbl7 = QLabel(self)
        self.lbl7.setText("C_LEFT_BOUND")
        self.lbl7.move(20, 164)

        self.lbl8 = QLabel(self)
        self.lbl8.setText("C_RIGHT_BOUND")
        self.lbl8.move(20, 184)

        self.lbl9 = QLabel(self)
        self.lbl9.setText("C_STEP")
        self.lbl9.move(50, 204)

        self.lbl3 = QLabel(self)
        self.lbl3.setText("LAMBDA_STEP")
        self.lbl3.move(20, 224)

        # Input dialog

        self.leEpsilon = QLineEdit(self)
        self.leEpsilon.move(130, 42)
        self.leEpsilon.setText('0.0001')

        self.leQ = QLineEdit(self)
        self.leQ.move(130, 62)
        self.leQ.setText('1.425')

        self.led = QLineEdit(self)
        self.led.move(130, 82)
        self.led.setText('1.067')

        self.lev_p = QLineEdit(self)
        self.lev_p.move(130, 102)
        self.lev_p.setText('0.0017')

        self.lev_n0 = QLineEdit(self)
        self.lev_n0.move(130, 122)
        self.lev_n0.setText('0.0000125')

        self.leV_por = QLineEdit(self)
        self.leV_por.move(130, 142)
        self.leV_por.setText('0.05')

        self.leC_LEFT_BOUND = QLineEdit(self)
        self.leC_LEFT_BOUND.move(130, 162)
        self.leC_LEFT_BOUND.setText('0.005')

        self.leC_RIGHT_BOUND = QLineEdit(self)
        self.leC_RIGHT_BOUND.move(130, 182)
        self.leC_RIGHT_BOUND.setText('0.012')

        self.leC_STEP = QLineEdit(self)
        self.leC_STEP.move(130, 202)
        self.leC_STEP.setText('0.0001')

        self.leLAMBDA_STEP = QLineEdit(self)
        self.leLAMBDA_STEP.move(130, 222)
        self.leLAMBDA_STEP.setText('0.001')

        # Buttons

        self.btn = QPushButton('graph', self)
        self.btn.move(100, 250)
        self.btn.clicked.connect(self.showDialog)

        self.setWindowTitle('Equation')
        self.resize(300, 300)
        self.center()
        self.show()

    def exception(self):
        self.msg = QMessageBox.about(self, "Warning", "Wrong data")

    def showDialog(self):

        global eps, Q, d, v_n0, v_p, V_por
        global C_LEFT_BOUND, C_RIGHT_BOUND, C_STEP
        global LAMBDA_STEP, CC, v_n, Re

        try:

            # set data

            eps = Decimal(self.leEpsilon.text())
            Q = Decimal(self.leQ.text())
            d = Decimal(self.led.text())
            v_n0 = Decimal(self.lev_n0.text())
            v_p = Decimal(self.lev_p.text())
            V_por = Decimal(self.leV_por.text())
            C_LEFT_BOUND = Decimal(self.leC_LEFT_BOUND.text())
            C_RIGHT_BOUND = Decimal(self.leC_RIGHT_BOUND.text())
            C_STEP = Decimal(self.leC_STEP.text())
            LAMBDA_STEP = Decimal(self.leLAMBDA_STEP.text())
            CC = C_STEP
            v_n = v_n0
            # средняя скорость движения жидкости в трубе, м / с
            V = Decimal(4.0) * Q / (Decimal(numpy.pi) * d * d)
            # число Рейндольса, безразмерное
            Re = V * d / v_n

            # run main function

            main()

        except InvalidOperation:
            self.exception()

    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

if __name__ == '__main__':

    app = QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())
