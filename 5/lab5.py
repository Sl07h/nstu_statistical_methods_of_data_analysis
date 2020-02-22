import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import inv, det, eigvals, norm
from sys import argv, path


def cos(a, b):
    ''' косинус между углами (a,b) = |a|*|b|*cos(a,b) '''
    return np.dot(a, b) / (norm(a)*norm(b))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
class lab5():
    def __init__(self):
        ''' Выделение памяти под массивы '''
        self.x = np.ndarray((n, m))
        self.X = np.ndarray((n, m))
        self.u = np.ndarray(n)
        self.y = np.ndarray(n)
        self.R = np.ndarray((m, m))

    def read_from_file(self, filename):
        ''' Считываем данные из файла '''
        data = np.loadtxt(filename,float,skiprows=1).transpose()
        self.x = data[1:m+1].transpose()
        self.u = data[-2].transpose()
        self.y = data[-1].transpose()
        for i in range(n):
            self.X[i] = f_vector_T(self.x[i])

    def calc_det_of_inf_matrix(self):
        ''' Вычисляем определитель информационной матрицы '''
        X = self.X
        self.inf_matr = X.T @ X
        d = det(self.inf_matr)
        print('определитель информационной матрицы: ', d)

    def calc_min_eig_val(self):
        ''' Поиск минимального собственного значения '''
        self.eig_vals = eigvals(self.inf_matr)
        self.min_eig_val = np.min(self.eig_vals)
        self.max_eig_val = np.max(self.eig_vals)
        print('минимальное собственное значение: ', self.min_eig_val)

    def calc_cond_number(self):
        ''' Мера обусловленности по Нейману-Голдстейну '''
        self.cond_number = self.max_eig_val / self.min_eig_val
        print('число обусловленности: ', self.cond_number)

    def calc_max_pair_sopr(self):
        '''
        Определяем максимальную парную сопряжённость.
        3 вектора могут быть коллинеарны, а попарно нет
        Максимальный коэффициент сопряжённости не несёт на себе эффекта масштаба
        '''
        X = self.X.transpose()
        self.max_pair_sopr = -9000.0
        # строим матрицу сопряжённости
        for i in range(m):
            for j in range(m):
                self.R[i][j] = cos(X[i], X[j])
            self.R[i][i] = 1.0

        for i in range(m):
            for j in range(m):
                if i != j and (abs(self.R[i][j]) > self.max_pair_sopr):
                    self.max_pair_sopr = abs(self.R[i][j])
        print('максимальная парная сопряжённость: ', self.max_pair_sopr)

    def calc_max_sopr(self):
        ''' Максимальная сопряжённость '''
        R_inv = inv(self.R)
        self.max_sopr = -9000.0
        for i in range(m):
            R_i = np.sqrt(1.0 - (1.0 / R_inv[i][i]))
            if R_i > self.max_sopr and R_i != float('inf'):
                self.max_sopr = R_i
        print('максимальная сопряжённость: ', self.max_sopr)

    def calc_theta_ridge(self):
        '''
        Ридж-оценка
        (Мультиколлинеарные данные выглядят как хребты на карте)
        '''
        X = self.X
        y = self.y
        points = 11
        RSS = np.ndarray(points)
        theta_norm = np.ndarray(points)
        theta_calc_ridge = np.ndarray((points, m))
        lambda_params = np.linspace(1, 50, points)
        i = 0
        for lambda_param in lambda_params:
            lambdas = np.diag(np.full(m, lambda_param))
            theta = inv(X.T @ X + lambdas) @ X.T @ y
            RSS[i] = (y - X @ theta).T @ (y - X @ theta)
            theta_norm[i] = norm(theta)
            theta_calc_ridge[i] = theta.T
            i += 1

        # print(theta_calc_ridge)

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
        '''
        Метод главных компонент
        '''
        X = self.X.T
        y = self.y
        for i in range(m):
            X[i] -= np.mean(X[i])
        X_centered = X.T
        X_inf_matr_centered = X_centered.T @ X_centered

        V = np.diag(eigvals(X_inf_matr_centered))
        Z = X_centered @ V
        b = inv(Z.T @ Z) @ Z.T @ y
        self.theta_calc_mc = V.dot(b)
        print(self.theta_calc_mc)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    path.insert(1, '')
    from model_6parameters import *

    global rho
    rho = float(argv[1])
    print('5 Лабораторная работа:', m, 'параметров, шум', int(rho*100), '\n')

    l5 = lab5()
    l5.read_from_file('1/report/table_1_u_y_' + str(int(rho*100)) + '.txt')
    l5.calc_det_of_inf_matrix()
    l5.calc_min_eig_val()
    l5.calc_cond_number()
    l5.calc_max_pair_sopr()
    l5.calc_max_sopr()
    l5.calc_theta_ridge()
    l5.calc_theta_main_components()