"""

    prohl.m

    Description:
        This program is used to calculate the rotor mode, based on the transfer matrix method.
        The following code shows the calculation methods of Prohl.

    All the codes had been written by Zhizhen Dong in 2022

"""

import math
import matplotlib.pyplot as plt
import numpy as np
from pylab import mpl
from tqdm import tqdm

mpl.rcParams['font.sans-serif'] = ['SimHei']
mpl.rcParams['axes.unicode_minus'] = False

l1 = 0.12
d = 0.04
A = math.pi * d * d / 4

D = 0.5  # 盘直径
h = 0.025  # 盘厚度

a = 1
u = 0.3
rou = 7800  # 密度
E = 2.0e11  # 弹性模量
G = E / (2 * (1 + u))  # 切变模量
I = math.pi * (d ** 4) / 64  #
K1 = 2.0e7  # 弹簧刚度
v1 = 6 * E * I / (a * G * A * l1 * l1)
mi = rou * math.pi * D ** 2 / 4
Jp = mi * D ** 2 / 8
Jd = Jp / 2
Ji = Jp - Jd

L = np.array([l1, l1, l1, l1, l1, l1, l1, l1, l1, l1, l1, l1, l1, 0, 0])
M = np.array([0, mi, mi, mi, mi, mi, mi, 0, 0, 0, 0, 0, mi, mi, 0])
K = np.array([K1, 0, 0, 0, 0, 0, 0, K1, 0, 0, 0, K1, 0, 0, 0])
v = np.array([v1, v1, v1, v1, v1, v1, v1, v1, v1, v1, v1, v1, v1, 0, 0])
J = np.array([0, Ji, Ji, Ji, Ji, Ji, Ji, 0, 0, 0, 0, 0, Ji, Ji, 0])
k = 0

T = np.zeros((len(L), 4, 4))
Z = np.mat(np.zeros((4, 4)))
H = np.mat(np.zeros((4, 4)))
Tit = ['第一阶模态', '第二阶模态', '第三阶模态']

wi = []
ni = []


def prohl_T(w):
    for i in range(15):
        T[i, :, :] = np.mat([[1 + (L[i] ** 3) * (1 - v[i]) * (M[i] * w ** 2 - K[i]) / (6 * E * I),
                              L[i] + L[i] ** 2 * J[i] * w ** 2 / (2 * E * I), L[i] ** 2 / (2 * E * I),
                              L[i] ** 3 * (1 - v[i]) / (6 * E * I)],
                             [(L[i] ** 2) * (M[i] * w ** 2 - K[i]) / (2 * E * I), 1 + L[i] * J[i] * w ** 2 / (E * I),
                              L[i] / (E * I), L[i] ** 2 / (2 * E * I)],
                             [L[i] * (M[i] * w ** 2 - K[i]), J[i] * w ** 2, 1.0, L[i]],
                             [M[i] * w ** 2 - K[i], 0.0, 0.0, 1.0]])

    H = np.mat(T[0, :, :])
    for i in range(1, 15):
        H = np.mat(T[i, :, :]) * H

    return H


for w in tqdm(np.arange(0, 4000 + 0.01, 0.01)):
    Z = prohl_T(w)

    F = Z[2, 0] * Z[3, 1] - Z[2, 1] * Z[3, 0]
    if F * (-1) ** k < 0:
        wi.append(w)
        w = wi[k]
        ni.append(wi[k] * 30 / math.pi)
        k += 1

x = np.zeros(len(L))
y = np.zeros(len(L))
z = np.zeros(len(L))

if len(Tit) > k:
    tk = k
else:
    tk = len(Tit)

for i in range(tk):
    w = wi[i]
    for j in range(14):
        T[j, :, :] = np.mat([[1 + (L[j] ** 3) * (1 - v[j]) * (M[j] * w ** 2 - K[j]) / (6 * E * I),
                              L[j] + L[j] ** 2 * J[j] * w ** 2 / (2 * E * I), L[j] ** 2 / (2 * E * I),
                              L[j] ** 3 * (1 - v[j]) / (6 * E * I)],
                             [(L[j] ** 2) * (M[j] * w ** 2 - K[j]) / (2 * E * I), 1 + L[j] * J[j] * w ** 2 / (E * I),
                              L[j] / (E * I), L[j] ** 2 / (2 * E * I)],
                             [L[j] * (M[j] * w ** 2 - K[j]), J[j] * w ** 2, 1, L[j]],
                             [M[j] * w ** 2 - K[j], 0, 0, 1]])

    Z = np.mat(T[0, :, :])
    for j in range(1, 15):
        Z = np.mat(T[j, :, :]) * Z

    b = -Z[3, 0] / Z[3, 1]
    X = np.mat([[1], [b], [0], [0]])

    for n in range(1, 15):
        X = np.c_[X, np.mat(T[n - 1, :, :]) * np.mat(X[:, X.shape[1] - 1])]

    for j in range(15):
        y[j] = X[1, j]
        z[j] = X[3, j]
        x[j] = (j - 1) * l1

    # y[15]=X[1,15]
    # x[15]=1.56
    # z[15]=X[3,15]
    y = y / max(abs(y))
    z = z / max(abs(z))
    plt.subplot(len(Tit), 1, i + 1)
    plt.title(Tit[i]+' '+str(round(ni[i], 2))+'rpm')
    plt.plot(x, y, 'b-')
    plt.plot(x, z, 'r:')
    plt.xlabel('轴长 / mm')
    plt.ylabel('不平衡值 / kg')
    plt.axis([0, 1.6, -1.2, 1.2])
    plt.legend(['振型', '弯矩'])
    plt.grid()

plt.tight_layout()
plt.show()

# T11 = 1 + (L[i] ** 3) * (1 - v[i]) * (M[i] * w ** 2 - K[i]) / (6 * E * I)
# T12 = L[i] + L[i] ** 2 * J[i] * w ** 2 / (2 * E * I)
# T13 = L[i] ** 2 / (2 * E * I)
# T14 = L[i] ** 3 * (1 - v[i]) / (6 * E * I)
# T21 = (L[i] ** 2) * (M[i] * w ** 2 - K[i]) / (2 * E * I)
# T22 = 1 + L[i] * J[i] * w ** 2 / (E * I)
# T23 = L[i] / (E * I)
# T24 = L[i] ** 2 / (2 * E * I)
# T31 = L[i] * (M[i] * w ** 2 - K[i])
# T32 = J[i] * w ** 2
# T33 = 1.0
# T34 = L[i]
# T41 = M[i] * w ** 2 - K[i]
# T42 = 0.0
# T43 = 0.0
# T44 = 1.0
