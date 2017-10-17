# necessary in order to compile this script to .exe
from matplotlib import use as matplotlib_use
matplotlib_use("QT5Agg")

import os
import sys
import numpy
import decimal
from numpy import log10, sqrt, power, exp, log
from decimal import Decimal

from PyQt5.QtWidgets import *
from PyQt5.QtGui import QFont

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib import pyplot as plt
from matplotlib.figure import Figure

decimal.getcontext().prec = 100

class InputField:
    def __init__(self, text, coeffSI, left_bound = None, right_bound = None):
        self.text = text
        self.coeffSI = Decimal(coeffSI)
        self.left_bound = left_bound
        self.right_bound = right_bound

inputFields = {
    'Epsilon' : InputField('Шероховатость трубопровода, мм', 1e-3, 0.001, 10),
    'Q' : InputField('Объемная скорость потока, м^3 / c', 1.0, 0.001, 1000),
    'd' : InputField('Диаметр трубопровода, мм', 1e-3, 1, 2000),
    'v_n' : InputField('Вязкость нефти, сСт', 1e-6),

    'v_p' : InputField('Вязкость присадки, сСт', 1e-6),
    'V_por' : InputField('Пороговая скорость, м / с', 1.0, 0.0001, 0.3),

    'c_left_bound' : InputField('Начальная концентрация, ppm', 1e-6, 0, 350),
    'c_right_bound' : InputField('Конечная концентрация, ppm', 1e-6, 0, 350),
    'c_points_count' : InputField('Количество точек построения', 1.0, 3, 1000),

    'concentration_coefficient' : InputField('Коэффициент Z', 1.0, 0, 2e5),
}

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()

class Equation:

    def __init__(self):
        self.Epsilon = 0
        self.Q = 0
        self.d = 0
        self.v_n = 0

        self.v_p = 0
        self.V_por = 0

        self.c_left_bound = 0
        self.c_right_bound = 0
        self.c_points_count = 0

        self.lambda_left_bound = Decimal(1e-20)
        self.lambda_right_bound = Decimal(1)
        self.lambda_bruteforce_step = Decimal(1e-3)

        self.bin_search_max_operations = 100
        self.bin_search_precision = Decimal(1e-8)

    def outOfBounds(self, value, field):
        conditionLeft = conditionRight = False
        if (field.left_bound):
            left_bound = Decimal(field.left_bound)
            conditionLeft = value < left_bound and abs(value - left_bound) > 1e-10

        if (field.right_bound):
            right_bound = Decimal(field.right_bound)
            conditionRight = value > right_bound and abs(value - right_bound) > 1e-10

        return conditionLeft or conditionRight

    def setData(self, name, textValue):
        try:
            field = inputFields[name]
            assert(field)
            value = Decimal(textValue.replace(',', '.'))
            if self.outOfBounds(value, field):
                raise ValueError()
            setattr(self, name, value * field.coeffSI)

        except AssertionError:
            raise Exception('Wrong field "{0}" in function Equation.setData() \n Сообщите программистам об ошибке'.format(name))
        except decimal.InvalidOperation:
            raise Exception('Некорректные данные в поле "{0}" \n' \
                'Допускаются цифры и точка либо запятая в качестве разделителя'.format(field.text))
        except ValueError:
            raise Exception('Значение в поле "{0}" должно лежать в пределах от {1} до {2}'.format(field.text, field.left_bound, field.right_bound))

    def left_log(self, x):
        x = Decimal(x)
        base = (Decimal(2.8) * self.V_por) / (self.V * Decimal(sqrt(x)))
        p = self.beta / Decimal(5.75)
        return Decimal(power(base, p))

    def right_log(self, x):
        return Decimal(2.51) / (self.Re * Decimal(sqrt(x))) + self.Epsilon / (Decimal(3.7) * self.d)

    def right_side(self, x):
        return Decimal(-2.0) * Decimal(log10(self.left_log(x) * self.right_log(x)))

    def left_side(self, x):
        return Decimal(1) / Decimal(sqrt(x))

    def brute_force_lambda(self):
        min_diff = 1e10
        min_lambda = -1;
        for i in numpy.arange(self.lambda_left_bound, self.lambda_right_bound, self.lambda_bruteforce_step):
            cur_diff = abs(self.left_side(i) - self.right_side(i))
            if (min_diff > cur_diff):
                min_diff = cur_diff
                min_lambda = i
        return min_lambda

    def bin_search_lambda(self):
        L, R = self.lambda_left_bound, self.lambda_right_bound
        coeffL = self.left_side(L) - self.right_side(L)
        coeffR = self.left_side(R) - self.right_side(R)

        coeffL /= abs(coeffL)
        coeffR /= abs(coeffR)
        iterations = 0

        if coeffL * coeffR > 0:
        	return self.brute_force_lambda()

        while iterations < self.bin_search_max_operations and abs(L - R) > self.bin_search_precision:
            iterations += 1
            mid = (L + R) / Decimal(2.0);

            if coeffL * (self.left_side(mid) - self.right_side(mid)) > 0.0:
            	L = mid
            else:
            	R = mid

        return L

    def showPlot(self):
        CC_array = []
        lambda_array = []
        DR_array = []

        c_step = (self.c_right_bound - self.c_left_bound) / (self.c_points_count - 1)
        self.V = Decimal(4.0) * self.Q / (Decimal(numpy.pi) * self.d * self.d)
        self.Re = self.V * self.d / self.v_n
        self.beta = 0

        lambda_0 = self.bin_search_lambda()
        it = 0
        os.system('cls')

        for c_p in numpy.arange(self.c_left_bound, self.c_right_bound + c_step, c_step):
            printProgressBar(it, self.c_points_count - 1, prefix = 'Progress:', suffix = 'Complete', length = 50)
            it += 1
            c_n = Decimal(1) - c_p

            # пересчет вязкости раствора при по формуле Вальтера
            K = Decimal(0.8)
            sum = c_n * log10(log10(K + self.v_n / inputFields['v_n'].coeffSI)) + c_p * log10(log10(K + self.v_p / inputFields['v_p'].coeffSI))
            v_r = power(Decimal(10), power(Decimal(10), sum)) - K
            v_r = v_r * inputFields['v_n'].coeffSI

            # if (c_p > 1e-8):
            #     nu = (v_r - self.v_n) / (self.v_n * c_p)
            #     tay = Decimal(8.31 * 293) / Decimal(1e6) / nu / Decimal(1000)
            #     print(tay / inputFields['v_n'].coeffSI)

            # число Рейнольдса
            self.Re = self.V * self.d / v_r
            self.beta = self.concentration_coefficient * c_p * Decimal(100) # умножаем на 100, т.к. в формуле СС в процентах

            CC_array.append(c_p / inputFields['v_n'].coeffSI)
            lambda_array.append(self.bin_search_lambda())
            DR_array.append(100 * (1 - lambda_array[-1] / lambda_0))

        fig, ax1 = plt.subplots()
        ax1.plot(CC_array, lambda_array, 'b')
        ax1.set_xlabel('Объемная концентрация полимера, ppm')
        ax1.set_ylabel('Коэффициент трения', color='b')
        ax1.tick_params('y', colors='b')

        ax2 = ax1.twinx()
        ax2.plot(CC_array, DR_array, 'r')
        ax2.set_ylabel('Гидравлическая эффективность присадки, %', color='r')
        ax2.tick_params('y', colors='r')

        ax2.grid(True)

        fig.tight_layout()
        plt.show(block=False)

class Example(QWidget):

    def __init__(self):
        super().__init__()

        self.edits = {}
        self.currentY = 0
        self.initUI()

    def addLabelEditWidget(self, name):
        label = QLabel(self)
        label.setText(inputFields[name].text)
        label.setFont(QFont("Times", 9))

        edit = QLineEdit(self)
        self.edits[name] = edit

        self.layout().addWidget(label, self.currentY, 0)
        self.layout().addWidget(edit, self.currentY, 1)
        self.currentY += 1

    def addBlockCaption(self, text):
        label = QLabel(self)
        label.setText(text)
        label.setFont(QFont("Times", 13, QFont.Bold))

        self.layout().addWidget(label, self.currentY, 0)
        self.currentY += 1

    def addFormula(self, formula, _fontsize, _x, _horizonalalignment, _verticalalignment):
        bg = self.palette().window().color()
        cl = (bg.redF(), bg.greenF(), bg.blueF())
        fig = Figure(edgecolor=cl, facecolor=cl)

        canvas = FigureCanvasQTAgg(fig)        
        self.layout().addWidget(canvas, self.currentY, 0, self.currentY, 0)
        self.currentY += 1

        fig.clear()
        fig.suptitle(
            formula,
            fontsize = _fontsize,
            x=_x,
            horizontalalignment=_horizonalalignment,
            verticalalignment=_verticalalignment)
        canvas.draw()

    def addDescription(self):
        main_formula = (
            r'$\dfrac{1}{\lambda} = '
            r'-2 \log_{10}  \left[ \left( \frac{2.8\ V_{пор}^*}{V \sqrt{\lambda}} \right) ^ {\beta\ /\ 5.75}'
            r'\left( \frac{2.51}{Re \sqrt{\lambda}} + \frac{\Delta_э}{3.7 d} \right) \right]$'
        )
        
        self.addFormula(main_formula, 14, 0.5, 'center', 'top')

        legend_formula = (
            r'$\lambda$ - коэффициент трения'
            '\n'
            r'$V_{пор}^*$ - пороговая динамическая скорость'
            '\n'
            r'$V$ - средняя скорость движения жидкости в трубопроводе'
            '\n'
            r'$Re$ - число Рейнольдса'
            '\n'
            r'$d$ - диаметр трубопровода'
            '\n'
            r'$\Delta_э$ - шероховатость трубопровода'
            '\n'
            r'$\beta = Z * C$, где C - объемная концентрация присадки'
        )

        self.addFormula(legend_formula, 11, 0, 'left', 'top')

    def initUI(self):

        grid = QGridLayout(self)
        grid.setSpacing(4)
        grid.setRowMinimumHeight(0, 60)
        grid.setRowMinimumHeight(1, 170)

        self.setLayout(grid)
        self.addDescription()

        self.addBlockCaption('Параметры трубопровода')
        self.addLabelEditWidget('Epsilon')
        self.addLabelEditWidget('Q')
        self.addLabelEditWidget('d')
        self.addLabelEditWidget('v_n')

        self.currentY += 1

        self.addBlockCaption('Параметры присадки')
        self.addLabelEditWidget('v_p')
        self.addLabelEditWidget('V_por')
        self.addLabelEditWidget('concentration_coefficient')

        self.currentY += 1

        self.addBlockCaption('Параметры графика')
        self.addLabelEditWidget('c_left_bound')
        self.addLabelEditWidget('c_right_bound')
        self.addLabelEditWidget('c_points_count')

        # # Initial values

        self.edits['Epsilon'].setText('0.1')
        self.edits['Q'].setText('1.425')
        self.edits['d'].setText('1067')
        self.edits['v_n'].setText('12.2')

        self.edits['v_p'].setText('1700')
        self.edits['V_por'].setText('0.001')
        self.edits['concentration_coefficient'].setText('3000')

        self.edits['c_left_bound'].setText('0')
        self.edits['c_right_bound'].setText('150')
        self.edits['c_points_count'].setText('30')

        # # Buttons

        btn = QPushButton('Построить график', self)
        btn.clicked.connect(self.showDialog)
        self.layout().addWidget(btn, self.currentY, 0, self.currentY, 0)

        self.setWindowTitle('Применение противотурбулентных присадок')
        self.resize(510, 400)
        self.show()
        self.center()

    def exception(self, text):
        self.msg = QMessageBox.about(self, "Ошибка", text)

    def showDialog(self):
        try:
            eq = Equation()
            for key in inputFields.keys():
                eq.setData(key, self.edits[key].text())

            eq.showPlot()

        except Exception as e:
            #self.exception(e)
            self.exception(str(e))

    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())
