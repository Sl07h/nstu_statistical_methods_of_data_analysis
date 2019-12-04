import matplotlib.pyplot as plt
#import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy import stats
from sympy import *
from sys import argv, path


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

    def calc_min_eig_val(self):
        self.eig_vals = np.linalg.eigvals(self.inf_matr)
        print('lambda_min = ', np.min(self.eig_vals))

    def calc_cond_number(self):
        self.cond_number = np.max(self.eig_vals) / np.min(self.eig_vals)
        print('cond_number = ', self.cond_number)

    def calc_max_sopr(self):
        X = self.X.transpose()
        self.R = np.ndarray((m, m))
        self.max_sopr = -9000
        for i in range(m):
            for j in range(m):
                self.R[i][j] = 0
                for t in range(n):
                    numerator = X[i][t] * X[j][t]
                    denominator1 = 0
                    for t1 in range(n):
                        denominator1 += X[i][t1]**2
                    denominator2 = 0
                    for t2 in range(n):
                        denominator2 += X[i][t2]**2
                self.R[i][j] += numerator / np.sqrt(denominator1*denominator2)
            self.R[i][i] = 1

        for i in range(m):
            for j in range(m):
                if i!=j and abs(self.R[i][j]) > self.max_sopr:
                    self.max_sopr = abs(self.R[i][j])

        #print(self.R)
        print('max_sopr = %d' % self.max_sopr)

    def calc_pair_sopr(self):
        # R_inv = Matrix(self.R)**-1
        R_inv = np.linalg.inv(self.R)
        self.R_list = np.ndarray(m)
        # print(R_inv)
        self.pair_sopr = -9000
        for i in range(m):
            R_i = (1.0 - 1.0 / R_inv[i][i]) ** (1/2)
            # print(R_i)
            if R_i > self.pair_sopr and R_i != float('inf'):
                self.pair_sopr = R_i
        print('pair_sopr = %d' % self.pair_sopr)

    def calc_theta_ridge(self):
        X = Matrix(self.X)
        y = Matrix(self.y)
        _theta_true = Matrix(theta_true)
        points = 26
        RSS = np.ndarray(points)
        theta_norm = np.ndarray(points)
        theta_calc_ridge = np.ndarray((points, m))
        lambda_params = np.linspace(0, 25, points)
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
        pass


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
    l5.calc_max_sopr()
    l5.calc_pair_sopr()
    l5.calc_theta_ridge()
    l5.calc_theta_main_components()
