# файл параметров модели

theta_true = [0.1, -1, 0.1, 3]  # параметры уравнения t0*x1^2 + t1*x2^2 + th2*x1 + t3
a, b, n, m = -1, 1, 100, 4      # границы, число точек и параметров
rho = 0.15                      # шум в диапазоне [0.05, 0.15] или [0.50, 0.70]
alpha = 0.05                    # уровень значимости
dist_E = 0.000491836210109302   # для theta_true = [-0.05, 2, 4]
f_E = 9000                      # вместо float('inf')
count = 11


def f_2_parameters(th, x1, x2):
    return th[0]*x1**2 + th[1]*x2**2 + th[2]*x1 + th[3]


def f(th, x_vector):
    x1 = float(x_vector[0])
    x2 = float(x_vector[1])
    return th[0]*x1**2 + th[1]*x2**2 + th[2]*x1 + th[3]


def f_vector_(x_vector):
    x1 = float(x_vector[0])
    x2 = float(x_vector[1])
    return [x1**2, x2**2, x1, 1]


def f_vector(th, x_vector):
    x1 = float(x_vector[0])
    x2 = float(x_vector[1])
    return [th[0]*x1**2, th[1]*x2**2, th[2]*x1, th[3]]
