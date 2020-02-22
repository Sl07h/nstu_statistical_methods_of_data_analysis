# файл параметров модели
import numpy as np


# параметры уравнения 
# th0*x1 + th1*x2 + th2*x3 + th3*x4 + th4*x5 + th5*x6
# x6 = 2*x2 + x3 + 4*x4
theta_true = [1, 1, 1, 1, 1, 1]
n, m = 30, 6              # число точек, параметров и осей
rho = 0.15                      # шум в диапазоне [0.05, 0.15] или [0.50, 0.70]
a = [-1, -1, -1, -1, -1]
b = [ 1,  1,  1,  1,  1]


def sample_x():
    x = np.ndarray((m,n))
    x[0] = np.random.uniform(low=a[0], high=b[0], size=n)
    x[1] = np.random.uniform(low=a[1], high=b[1], size=n)
    x[2] = np.random.uniform(low=a[2], high=b[2], size=n)
    x[3] = np.random.uniform(low=a[3], high=b[3], size=n)
    x[4] = np.random.uniform(low=a[4], high=b[4], size=n)
    x[5] = 2*x[0] - 3*x[3] + np.random.uniform(0, 0.01)
    return x.transpose()


def f(theta, x):
    return np.dot(theta, x)

def f_vector_T(x):
    return np.array(x)