import math
import numpy as np
import matplotlib.pyplot as plt

# --------- Неявная разностная схема ---------

a = 1  # Коэффициент теплопроводности

X = 1
T = 3

M = 10
N = 10

h = X / M  # Шаг по пространству
tau = T / N  # Шаг по времени

# Цвета для графиков
colors = ['#CD5C5C', '#FF4500', '#8A2BE2', '#DAA520', '#7FFF00',
          '#7FFFD4', '#808080', '#DC143C', '#FF1493', '#808000',
          '#FFA07A', '#BDB76B', '#4B0082', '#BC8F8F', '#228B22']


# Начально-краевые задачи:
# ---- 1 ----
# (du(x,t)/dt) = (d^2u(x,t)/dx^2) + sin(t)cos(x) + cos(x)*cos(t); x-(0,1), t-(0,3)
# u(x,0) = 0
# u(0,t) = sin(t); u(X,t) = cos(X)sin(t)

# ---- 2 ----
# (du(x,t)/dt) = (d^2u(x,t)/dx^2) + x*cos(xt) + t^2*sin(xt); x-(0,pi), t-(0, 1)
# u(x,0) = 0
# u(0,t) = 0; u(pi,t) = sin(pi*t)


# Начальное условие: u(x,0) = phi(x), x из [0,X]
def initial_condition(x):
    # phi(x)
    return 0  # --1-- и --2--


# Краевые условия 1-го рода: u(0,t) = g1(t), u(X,t) = g2(t)
def boundary_condition1(t):
    # g1(t)
    return math.sin(t)  # --1--
    # return 0  # --2--


def boundary_condition2(X, t):
    # g2(t)
    return math.cos(X) * math.sin(t)  # --1--
    # return math.sin(X*t)  # --2--


# Функция f(x,t)
def heat_source(x, t):
    # f(x, t)
    return math.sin(t) * math.cos(x) + math.cos(x) * math.cos(t)  # --1--
    # return x * math.cos(x * t) + (t ** 2) * math.sin(x * t)  # --2--


# Аналитическое решение
def analytical_solution(x, t):
    return math.cos(x) * math.sin(t)  # --1--
    # return math.sin(x * t)  # --2--


print("Шаг по x:", h, "  Шаг по t:", tau)

# Значения нулевого слоя получаем из начального условия
matrix_values = []
lst = []
for j in range(M + 1):
    lst.append(initial_condition(h * j))
matrix_values.append(lst)
# print(matrix_values)


A = (a * tau) / (h ** 2)
B = 1 + (2 * a * tau) / (h ** 2)
C = A

# Метод прогонки
for n in range(1, N + 1):
    lst = []
    alpha, betta = [], []

    # Прямая прогонка
    al = 0
    alpha.append(al)
    bt = boundary_condition1(tau * n)
    betta.append(bt)
    for j in range(1, M):
        D = tau * heat_source(h * j, tau * n) + matrix_values[n - 1][j]
        al = A / (B - C * alpha[j - 1])
        bt = (C * betta[j - 1] + D) / (B - C * alpha[j - 1])
        alpha.append(al)
        betta.append(bt)

    # Одратная прогонка
    um = boundary_condition2(X, tau * n)
    lst.append(um)
    for j in range(0, M):
        u = lst[j] * alpha[M - 1 - j] + betta[M - 1 - j]
        lst.append(u)
    lst.reverse()
    matrix_values.append(lst)

# Вывести матрицу значений
# for v in matrix_values:
#     print(v)

# Сравнение численного и точного решений для фиксированных значений t
x = [j * h for j in range(M + 1)]
for n in range(0, N + 1, 2):
    an_lst = []
    t = n * tau
    for j in range(M + 1):
        an_lst.append(analytical_solution(x[j], t))
    plt.plot(x, matrix_values[n], "--", color=colors[n % 15], marker='*')
    plt.plot(x, an_lst, color=colors[n % 15], marker='o', label=f"t={t}")
plt.suptitle('\'--*\' - численное решение;\n\'-o\' - точное решение')
plt.xlabel("x")
plt.ylabel("u(x,t)")
plt.legend()
plt.grid(True)
plt.show()

t = [n * tau for n in range(N + 1)]
xgrid, tgrid = np.meshgrid(x, t)

# Численное решение
matrix_values = np.array(matrix_values)
fig = plt.figure()
ax_3d = fig.add_subplot(projection='3d')
ax_3d.plot_surface(xgrid, tgrid, matrix_values, color='green')
ax_3d.set_xlabel('x')
ax_3d.set_ylabel('t')
ax_3d.set_zlabel('u')
plt.suptitle('Численное решение')
plt.show()

# Точное решение
zgrid = np.cos(xgrid) * np.sin(tgrid)
fig2 = plt.figure()
ax_3d2 = fig2.add_subplot(projection='3d')
ax_3d2.plot_surface(xgrid, tgrid, zgrid)
ax_3d2.set_xlabel('x')
ax_3d2.set_ylabel('t')
ax_3d2.set_zlabel('u')
plt.suptitle('Точное решение')
plt.show()

error_matrix = np.abs(matrix_values - zgrid)

# Погрешность численного решения для фиксированных значений t
for n in range(1, N + 1, 2):
    t1 = n * tau
    plt.plot(x, error_matrix[n], color=colors[n % 15], marker='+', label=f"t={t1}")
plt.suptitle('Погрешность численного решения')
plt.xlabel("x")
plt.ylabel("error")
plt.legend()
plt.grid(True)
plt.show()

# Погрешности считаются для последнего слоя (t = T)
abs_error = np.max(np.abs(matrix_values[N] - zgrid[N]))
print("Абсолютная погрешность(норма разности точного и численного решений):", abs_error)
print("Норма точного решения", np.max(np.abs(zgrid[N])))
relative_error = abs_error / (np.max(np.abs(zgrid[N])))
print("Относительная погрешность(отношение абсолютной погрешности к норме точного решения):", relative_error)
