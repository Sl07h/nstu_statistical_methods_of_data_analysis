import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy import stats
from sympy import *
from sys import argv, path


# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
class lab2():
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

    def calc_theta(self):
        X = Matrix(self.X)
        y = Matrix(self.y)
        theta_calc = ((X.T * X) ** -1) * X.T * y
        self.theta_calc = np.array(theta_calc.T).astype(np.float64)[0]

    def calc_e(self):
        X = Matrix(self.X)
        y = Matrix(self.y)
        th_calc = Matrix(self.theta_calc)
        yCalc = X * th_calc
        e_calc = y - yCalc
        self.yCalc = np.array(yCalc.T).astype(np.float64)[0]
        self.e_calc = np.array(e_calc.T).astype(np.float64)[0]

    def calc_dist(self):
        e_calc = Matrix(self.e_calc)
        dist_calc = (e_calc.T * e_calc) / (n - m)
        self.dist_calc = np.array(dist_calc.T).astype(np.float64)[0][0]

    def draw(self, title, name, elevation, azimuth):
        path_to_save = '2/pics/' + name + '_' + str(int(rho*100)) + '_' + str(elevation) + '_' + str(azimuth) + '.png'
        title += ' ' + str(int(rho * 100)) + '%'
        fig = plt.figure('1')
        ax = fig.gca(projection='3d')
        _t = np.linspace(a, b, sqrt(n))
        _x, _y = np.meshgrid(_t, _t)
        _z = f_2_parameters(theta_true, _x, _y)
        fig.subplots_adjust(bottom=-0.05, top=1, left=-0.05, right=1.05)
        ax.view_init(elevation, azimuth)
        ax.plot_surface(_x, _y, _z, zorder=2, alpha=0.2)
        ax.scatter(self.x1, self.x2, self.yCalc, c='black')
        plt.title(title, fontsize=19)
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.grid(alpha=0.4)
        plt.savefig(path_to_save)
        plt.clf()

    def check_goodness_of_fit(self):
        self.F = self.dist_calc / dist_E
        d1, d2 = n-m, f_E
        self.F_t = stats.f.ppf(1 - alpha, d1, d2 )
        print('F: ' + str(self.F) + ' F_t ' + str(self.F_t))
        if self.F <= self.F_t:
            print("Гипотеза не отвергается")
        else:
            print("Полученная модель неадекватна")

    def save_table_2(self):
        with open('2/report/table_2_u_y_' + str(int(rho*100)) + '.txt', 'w') as file:
            file.write('i\tx1\tx2\tu\tyCalc\tyMinyCalc\n')
            for i in range(n):
                x1 = self.x1[i]
                x2 = self.x2[i]
                u = self.u[i]
                y = self.y[i]
                y_minus_yCalc = float(self.y[i] - self.u[i])
                file.write('{:d}\t{:.17f}\t{:.17f}\t{:.17f}\t{:.17f}\t{:.17f}\n'.format(i, x1, x2, u, y, y_minus_yCalc))

    def save_table_3(self):
        with open('2/report/table_3_th_' + str(int(rho*100)) + '.txt', 'w') as file:
            file.write('theta_true\ttheta_calc\n')
            for i in range(m):
                file.write('{:.17f}\t{:.17f}\n'. format(theta_true[i], self.theta_calc[i]))
    
    def save_table_4(self):            
        with open('2/report/table_4_dist_F_' + str(int(rho*100)) + '.txt', 'w') as file:
            file.write('dist_calc\tdist_E\tF\tF_t\n')
            file.write('{:.17f}\t{:.17f}\t{:.17f}\t{:.17f}\n'.format(self.dist_calc, dist_E, self.F, self.F_t))
    



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    path.insert(1, '')
    if int(argv[2]) == 3:
        from model_3parameters import *
    else:
        from model_4parameters import *
    global rho
    rho = float(argv[1])
    print('\n\nЗапушен код 2 лабораторной работы: {:d} параметров, шум {:d}%'.format(int(argv[2]), int(rho*100)))

    l2 = lab2()
    l2.calc_theta()
    l2.calc_e()
    l2.calc_dist()
    l2.draw('полученная модель', 'calculated', None, None)
    l2.draw('полученная модель', 'calculated', 0, 0)
    l2.draw('полученная модель', 'calculated', 0, 90)
    l2.check_goodness_of_fit()
    l2.save_table_2()
    l2.save_table_3()
    l2.save_table_4()