# Algoritmo de eliminación de Gauss
import numpy as np
a = np.array([[2.0, 1.0],
              [1.0, 2.0],
            ])
print('La matriz de coeficientes a es : \n',a)

b = np.array([[1.0], [0]])
print('La matriz columna b es : \n',b)
def gaussElimin(a,b):
  n = len(b)
  # Fase de eliminacion
  for k in range(0,n-1):
    for i in range(k+1,n):
      if a[i,k] != 0.0:
        lam = a [i,k]/a[k,k]
        a[i,k+1:n] = a[i,k+1:n] - lam*a[k,k+1:n]
        b[i] = b[i] - lam*b[k]
  # Fase de sustitucion hacia atras
  for k in range(n-1,-1,-1):
    b[k] = (b[k] - np.dot(a[k,k+1:n],b[k+1:n]))/a[k,k]
  return b

print(f"Las cargas de los quarks u y d son: \n {gaussElimin(a,b)}")        

# Método de interpolación de Newton


def evalPoly(a, xData, x):  # Función que evalua polinomios de Lagrange
    n = len(xData) - 1  # Grado del polinomio
    p = a[n]
    for k in range(1, n + 1):
        p = a[n - k] + (x - xData[n - k]) * p
    return p

# Ejemplo método de Newton
def coeffts(xData, yData):
    m = len(xData)  # Número de datos
    a = yData.copy()
    for k in range(1, m):
        a[k:m] = (a[k:m] - a[k - 1]) / (xData[k:m] - xData[k - 1])
    return a
import numpy as np
import matplotlib.pyplot as plt
from math import *

xData = np.array([0.15, 2.30, 3.15, 4.85, 6.25, 7.95])
yData = np.array([4.79867, 4.49013, 4.2243, 3.47313, 2.66674, 1.51909])
coeff = coeffts(xData, yData)
x = np.arange(0, 8.5, 0.5)
plt.plot(x, evalPoly(coeff, xData, x), "r", label="Newton")
plt.plot(xData, yData, "o", label="Datos")
plt.legend()
plt.grid()
plt.show()
print("  x    yExacta        yInt       Error(%)")
print("------------------------------------------")
for i in range(len(x)):
    y = evalPoly(coeff, xData, x[i])
    yExacta = 4.8 * cos(pi * x[i] / 20)
    Error = abs(((yExacta - y) / yExacta) * 100)
    print(" %.1f  %.8f   %.8f    %.8f" % (x[i], yExacta, y, Error))