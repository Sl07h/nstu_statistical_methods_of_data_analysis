import matplotlib.pyplot as plt
#import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy import stats
from sympy import *
from sys import argv, path
import math


def dotproduct(v1, v2):
    return sum((a*b) for a, b in zip(v1, v2))


def length(v):
    return math.sqrt(abs(dotproduct(v, v)))


def angle2(v1, v2):
    return math.acos(abs(dotproduct(v1, v2)) / (length(v1) * length(v2)))


def angle(v1, v2):
    return np.dot(a, b) / (np.linalg.norm(a)*np.linalg.norm(b))

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


class lab5():
    def __init__(self):
        self.x = np.ndarray((n, m))
        self.X = np.ndarray((n, m))
        self.u = np.ndarray(n)
        self.y = np.ndarray(n)

        with open('1/report/table_1_u_y_' + str(int(rho*100)) + '.txt', 'r') as file:
            line = file.readline().rstrip()
            for i in range(n):
                line = file.readline().rstrip().rsplit()
                for j in range(m):
                    self.x[i][j] = float(line[j + 1])
                self.u[i] = float(line[m + 1])
                self.y[i] = float(line[m + 2])

        for i in range(n):
            self.X[i] = f_vector_(self.x[i])

    def calc_det_of_inf_matrix(self):
        X = Matrix(self.X)
        self.inf_matr = np.array(X.T * X).astype(np.float64)
        d = np.linalg.det(self.inf_matr)
        print('det = ', d)

    # def calc_det_of_inf_matrix(self):
    #     X = Matrix(self.X)
    #     self.inf_matr_sp = X.T * X
    #     self.inf_matr_sp = self.inf_matr_sp / Trace(self.inf_matr_sp)
    #     print(self.inf_matr_sp)
    #     self.inf_matr = np.array(self.inf_matr_sp).astype(np.float64)
    #     d = np.linalg.det(self.inf_matr)
    #     print('det = ', d)

    def calc_min_eig_val(self):
        self.eig_vals = np.linalg.eigvals(self.inf_matr)
        print('lambda_min = ', np.min(self.eig_vals))

    def calc_cond_number(self):
        self.cond_number = np.max(self.eig_vals) / np.min(self.eig_vals)
        print('cond_number = ', self.cond_number)

    def calc_pair_sopr(self):
        X = self.X.transpose()
        self.R = np.ndarray((m, m))
        self.pair_sopr = -9000.0
        for i in range(m):
            for j in range(m):
                self.R[i][j] = 0.0
                for t in range(n):
                    numerator = X[i][t] * X[j][t]
                    denominator1 = 0.0
                    for t1 in range(n):
                        denominator1 += X[i][t1]**2
                    denominator2 = 0.0
                    for t2 in range(n):
                        denominator2 += X[i][t2]**2
                self.R[i][j] += numerator / np.sqrt(denominator1*denominator2)
            self.R[i][i] = 1.0

        #print(max(abs(self.R.max()), abs(self.R.min())))

        for i in range(m):
            for j in range(m):
                if i != j and (abs(self.R[i][j]) > self.pair_sopr):
                    self.pair_sopr = abs(self.R[i][j])

        # print(self.R)
        print('pair_sopr = %f' % self.pair_sopr)

    # def calc_pair_sopr(self):
    #     X = self.x.transpose()
    #     # print(X)
    #     self.R = np.ndarray((m, m))
    #     self.pair_sopr = -9000
    #     for i in range(m):
    #         for j in range(m):
    #             self.R[i][j] = angle(X[i], X[j])
    #         self.R[i][i] = 1.0

    #     for i in range(m):
    #         for j in range(m):
    #             if i != j and abs(self.R[i][j]) > self.pair_sopr:
    #                 self.pair_sopr = abs(self.R[i][j])

    #     # print(self.R)
    #     print('pair_sopr = %d' % self.pair_sopr)

    def calc_max_sopr(self):
        # R_inv = Matrix(self.R)**-1
        R_inv = np.linalg.inv(self.R)
        self.R_list = np.ndarray(m)
        # print(R_inv)
        self.max_sopr = -9000.0
        for i in range(m):
            R_i = sqrt(1.0 - (1.0 / R_inv[i][i]))
            if R_i > self.max_sopr and R_i != float('inf'):
                self.max_sopr = R_i
        #print(max(abs(self.R.max()), abs(self.R.min())))

        print('max_sopr = %f' % self.max_sopr)

    def calc_theta_ridge(self):
        X = Matrix(self.X)
        y = Matrix(self.y)
        _theta_true = Matrix(theta_true)
        points = 11
        RSS = np.ndarray(points)
        theta_norm = np.ndarray(points)
        theta_calc_ridge = np.ndarray((points, m))
        lambda_params = np.linspace(2, 5, points)
        i = 0
        for lambda_param in lambda_params:
            lambdas = Matrix(np.diag(np.full(m, lambda_param)))
            _th = ((X.T*X + lambdas) ** -1) * X.T * y
            _RSS = (y - X*_th).T * (y - X*_th)
            _theta_norm = (_th - _theta_true).norm()

            RSS[i] = np.array(_RSS).astype(np.float64)
            theta_norm[i] = np.array(_theta_norm).astype(np.float64)
            theta_calc_ridge[i] = np.array(_th.T).astype(np.float64)
            i += 1

        path_to_save = '5/pics/RSS_from_lambda.png'
        plt.title('RSS(lambda)')
        plt.plot(lambda_params, RSS)
        plt.savefig(path_to_save)
        plt.clf()

        path_to_save = '5/pics/theta_from_lambda.png'
        plt.title('theta(lambda)')
        plt.plot(lambda_params, theta_norm)
        plt.savefig(path_to_save)
        plt.clf()

    def calc_theta_main_components(self):
        X = self.X.T
        for i in range(m):
            X[i] -= np.mean(X[i])
        X_centered = Matrix(X.T)
        X_inf_matr_centered = X_centered.T * X_centered
        X_centered_np = np.array(X_inf_matr_centered).astype(np.float64)

        V = Matrix(np.diag(np.linalg.eigvals(X_centered_np)))
        y = Matrix(self.y)
        Z = X_centered * V
        b = ((Z.T*Z) ** -1) * Z.T * y
        self.theta_calc_mc = V.dot(b)
        print(self.theta_calc_mc)


# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
if __name__ == "__main__":
    path.insert(1, '')
    from model_6parameters import *

    global rho
    rho = float(argv[1])
    print('Запушен код 5 лабораторной работы: {:d} параметров, шум {:d}%'.format(
        m, int(rho*100)))
    l5 = lab5()
    l5.calc_det_of_inf_matrix()
    l5.calc_min_eig_val()
    l5.calc_cond_number()
    l5.calc_pair_sopr()
    l5.calc_max_sopr()
    l5.calc_theta_ridge()
    l5.calc_theta_main_components()
