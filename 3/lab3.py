import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy import stats
from sympy import *
from sys import argv, path


# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
class lab3():
    def __init__(self):
        self.x = np.ndarray((n, 2))
        self.x1 = np.ndarray(n)
        self.x2 = np.ndarray(n)
        self.u = np.ndarray(n)
        self.y = np.ndarray(n)
        with open('1/report/table_1_u_y_' + str(int(rho*100)) + '.txt', 'r') as file:
            line = file.readline().rstrip()
            for i in range(n):
                line = file.readline().rstrip().rsplit()
                self.x1[i] = float(line[1])
                self.x2[i] = float(line[2])
                self.u[i] = float(line[3])
                self.y[i] = float(line[4])

        for i in range(n):
            self.x[i] = [self.x1[i], self.x2[i]]

        self.X = np.ndarray((n, m))
        for i in range(n):
            self.X[i] = f_vector_(self.x[i])

        global theta_true
        self.theta_calc = np.ndarray(m)
        with open('2/report/table_3_th_' + str(int(rho*100)) + '.txt', 'r') as file:
            line = file.readline()
            for i in range(m):
                line = file.readline().rstrip().rsplit()
                theta_true[i] = line[0]
                self.theta_calc[i] = line[1]

        with open('2/report/table_4_dist_F_' + str(int(rho*100)) + '.txt', 'r') as file:
            line = file.readline()
            data = [float(_) for _ in file.readline().rstrip().rsplit()]
            self.dist_calc = data[0]
            dist_E = data[1]
            self.F = data[2]
            self.F_t = data[3]

    def build_conf_interval_theta(self):
        X = Matrix(self.X)
        self.matrix = ((X.T*X)**-1).tolist()
        self.t = abs(stats.t.ppf(1 - alpha/2, n - m))
        self.thera_lower = np.ndarray(m)
        self.thera_upper = np.ndarray(m)
        for i in range(m):
            param = self.theta_calc[i]
            self.thera_lower[i] = param - self.t * self.matrix[i][i]
            self.thera_upper[i] = param + self.t * self.matrix[i][i]

    def check_significance_of_parameters(self):
        self.F_t = stats.f.ppf(1 - alpha, 1, n - m)
        self.F = np.ndarray(m)
        for i in range(m):
            self.F[i] = self.theta_calc[i]**2 / \
                (self.dist_calc * self.matrix[i][i])

    def save_table_4(self):
        with open('3/report/table_4_th_' + str(int(rho*100)) + '.txt', 'w') as file:
            file.write('theta_l\ttheta\ttheta_r\tF\tF_t\tsignificance\n')
            for i in range(m):
                file.write('{:.17f}\t{:.17f}\t{:.17f}\t{:.17f}\t{:.17f}\t'.format(
                    self.thera_lower[i],
                    self.theta_calc[i],
                    self.thera_upper[i],
                    self.F[i],
                    self.F_t
                ))
                if self.F[i] < self.F_t:
                    file.write('-\n')
                else:
                    file.write('+\n')

    def check_significance_of_regression(self):
        y = Matrix(self.y)
        X = Matrix(self.X)
        th = Matrix(self.theta_calc)
        RSS = (y - X*th).T * (y - X*th)
        self.RSS = np.array(RSS).astype(np.float64)[0][0]

        self.RSS_H = 0
        yMean = np.full((n), np.mean(self.y))
        for i in range(n):
            self.RSS_H += (self.y[i] - yMean[i])**2

        q = m - 1
        self.regr_F = ((self.RSS_H - self.RSS) / q) / (self.RSS / (n-m))
        self.regr_F_t = stats.f.ppf(1 - alpha, q, n - m)
        if self.regr_F <= self.regr_F_t:
            print('Гипотеза о незначимости регрессии принимается {:f} <= {:f}'.format(self.regr_F, self.regr_F_t))
        else:
            print('Гипотеза о незначимости регрессии отвергается {:f} > {:f}'.format(self.regr_F, self.regr_F_t))

    def estimate_E_for_x1(self):
        path_to_save = '3/pics/conf_y_for_x1_' + str(int(rho*100)) + '.png'
        x1_list = np.linspace(a, b, count)
        x2 = 0  # середина D(f)
        X = Matrix(self.X)
        lower = np.ndarray(count)
        center = np.ndarray(count)
        upper = np.ndarray(count)
        for i in range(count):
            x1 = x1_list[i]
            x = [x1, x2]
            fx = f(self.theta_calc, x)
            fx_vector = Matrix(f_vector(self.theta_calc, x))
            tmp = fx_vector.T * (X.T * X)**-1 * fx_vector
            sigma = self.dist_calc * \
                (1 + np.array(tmp).astype(np.float64)[0][0])
            lower[i] = fx - self.t * sigma
            center[i] = fx
            upper[i] = fx + self.t * sigma

        plt.title('оценка M(x1), помеха ' + str(int(rho * 100)) + '%')
        plt.plot(lower)
        plt.plot(center)
        plt.plot(upper)
        plt.savefig(path_to_save)
        plt.clf()

    def estimate_E_for_x2(self):
        path_to_save = '3/pics/conf_y_for_x2_' + str(int(rho*100)) + '.png'
        x1 = 0  # середина D(f)
        x2_list = np.linspace(a, b, count)
        X = Matrix(self.X)
        lower = np.ndarray(count)
        center = np.ndarray(count)
        upper = np.ndarray(count)
        for i in range(count):
            x2 = x2_list[i]
            x = [x1, x2]
            fx = f(self.theta_calc, x)
            fx1vector = Matrix(f_vector(self.theta_calc, x))
            tmp = sqrt(fx1vector.T * (X.T * X)**-1 * fx1vector)
            sigma = sqrt(self.dist_calc) * \
                np.array(tmp).astype(np.float64)[0][0]
            lower[i] = fx - self.t * sigma
            center[i] = fx
            upper[i] = fx + self.t * sigma

        plt.title('оценка M(x2), помеха ' + str(int(rho * 100)) + '%')
        plt.plot(lower)
        plt.plot(center)
        plt.plot(upper)
        plt.savefig(path_to_save)
        plt.clf()


# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
if __name__ == "__main__":
    path.insert(1, '')
    if int(argv[2]) == 3:
        from model_3parameters import *
    else:
        from model_4parameters import *
    global rho
    rho = float(argv[1])
    print('\n\nЗапушен код 3 лабораторной работы: {:d} параметров, шум {:d}%'.format(int(argv[2]), int(rho*100)))

    l3 = lab3()
    l3.build_conf_interval_theta()
    l3.check_significance_of_parameters()
    l3.save_table_4()
    l3.check_significance_of_regression()
    l3.estimate_E_for_x1()
    l3.estimate_E_for_x2()
