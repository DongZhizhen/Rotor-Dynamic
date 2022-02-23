"""

    riccati.m

    Description:
        This program is used to calculate the rotor mode, based on the transfer matrix method.
        The following code shows the calculation methods of Riccati in three different ways.

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

u11 = np.zeros((len(L), 2, 2))
u12 = np.zeros((len(L), 2, 2))
u21 = np.zeros((len(L), 2, 2))
u22 = np.zeros((len(L), 2, 2))
S = np.zeros((len(L) + 1, 2, 2))

Z = np.zeros((len(L), 2, 2))
Tit = ['第一阶模态', '第二阶模态', '第三阶模态']  # 模态阶数选择
w_max = 500  # 最大扫频频率

wi = []
ni = []

al_name = 'rct3'  # 算法选择 可选rct1 rct2 rct3

Xl = np.zeros(len(L))
Al = np.zeros(len(L))
Ml = np.zeros(len(L))
Ql = np.zeros(len(L))
zl = np.zeros(len(L))


def riccati_T(w):
    for i in range(len(L)):
        u11[i, :, :] = np.mat([[1.0, L[i]], [0.0, 1.0]])
        u12[i, :, :] = np.mat([[L[i] * (M[i] * w ** 2 - K[i]), J[i] * w ** 2], [M[i] * w ** 2 - K[i], 0.0]])
        u21[i, :, :] = np.mat([[L[i] ** 2 / (2 * E * I), L[i] ** 3 * (1 - v[i]) / (6 * E * I)],
                               [L[i] / (E * I), L[i] ** 2 / (2 * E * I)]])
        u22[i, :, :] = np.mat([[1 + (L[i] ** 3) * (1 - v[i]) * (M[i] * w ** 2 - K[i]) / (6 * E * I),
                                L[i] + L[i] ** 2 * J[i] * w ** 2 / (2 * E * I)],
                               [(L[i] ** 2) * (M[i] * w ** 2 - K[i]) / (2 * E * I),
                                1 + L[i] * J[i] * w ** 2 / (E * I)]])

        Z[i, :, :] = np.mat(u21[i, :, :] * np.mat(S[i, :, :]) + u22[i, :, :]).I
        S[i + 1, :, :] = (u11[i, :, :] * np.mat(S[i, :, :]) + u12[i, :, :]) * np.mat(Z[i, :, :])

    return Z, S


# 原始方法，有奇点，漏根
if al_name == 'rct1':
    for w in tqdm(np.arange(0, w_max + 0.01, 0.01)):
        (Z, S) = riccati_T(w)

        F = S[len(L), 0, 0] * S[len(L), 1, 1] - S[len(L), 0, 1] * S[len(L), 1, 0]
        if F * (-1) ** k < 0:
            wi.append(w)
            w = wi[k]
            ni.append(wi[k] * 30 / math.pi)
            k += 1

# 改进方法，无奇点
elif al_name == 'rct2':
    for w in tqdm(np.arange(0, w_max + 0.01, 0.01)):
        (Z, S) = riccati_T(w)

        dlt = np.mat([[1.0, 0.0], [0.0, 1.0]])
        for i in range(15):
            dlt = dlt * np.mat(Z[i, :, :])

        F = np.linalg.det(dlt * np.mat(S[len(L), :, :]))
        if F * (-1) ** k < 0:
            wi.append(w)
            w = wi[k]
            ni.append(wi[k] * 30 / math.pi)
            k += 1

# 改进方法，无奇点，奇偶运算
elif al_name == 'rct3':
    for w in tqdm(np.arange(0, w_max + 0.01, 0.01)):
        (Z, S) = riccati_T(w)

        dlt = 1
        for i in range(15):
            dlt = dlt * np.sign(np.linalg.det(np.mat(Z[i, :, :])))

        F = np.linalg.det(dlt * np.mat(S[len(L), :, :]))
        if F * (-1) ** k < 0:
            wi.append(w)
            w = wi[k]
            ni.append(wi[k] * 30 / math.pi)
            k += 1

else:
    print('none')

if len(Tit) > k:
    tk = k
else:
    tk = len(Tit)

for i in tqdm(range(tk)):
    w = wi[i]
    (Z, S) = riccati_T(w)

    b = -S[len(L), 1, 0] / S[len(L), 1, 1]
    e = np.mat([[1], [b]])

    for n in range(len(L), 1, -1):
        e = np.c_[np.mat(Z[n - 1, :, :]) * np.mat(e[:, 0]), e]

    f = np.mat(np.zeros((2, len(L))))
    for n in range(len(L)):
        f[:, n] = np.mat(S[n, :, :]) * np.mat(e[:, n])

    zl = np.array(range(15)) * l1
    for j in range(15):
        Xl[j] = e[0, j]
        Al[j] = e[1, j]
        Ml[j] = f[0, j]
        Ql[j] = f[1, j]

    # y[15]=X[1,15]
    # x[15]=1.56
    # z[15]=X[3,15]
    Xl = Xl / max(abs(Xl))
    Ml = Ml / max(abs(Ml))
    plt.subplot(len(Tit), 1, i + 1)
    plt.title(Tit[i] + ' ' + str(round(ni[i], 2)) + 'rpm')
    plt.plot(zl, Xl, 'b-')
    plt.plot(zl, Ml, 'r:')
    plt.xlabel('轴长 / mm')
    plt.ylabel('不平衡值 / kg')
    plt.axis([0, 1.6, -1.2, 1.2])
    plt.legend(['振型', '弯矩'])
    plt.grid()

plt.tight_layout()
plt.show()
