
import numpy as np

def evaluate(x, y, a=1.0, b=5.1/(4.0*np.pi**2), c=5.0/np.pi, r=6.0, s=10.0, t=1.0/(8.0*np.pi)):
  """
  https://www.sfu.ca/~ssurjano/branin.html
  x \in [-5, 10]
  y \in [0, 15]

  Input
  -----
    x: np.ndarray, shape=(:, )
    y: np.ndarray, shape=(:, )
    a = constant (optional), with default value 1
    b = constant (optional), with default value 5.1/(4*pi^2)
    c = constant (optional), with default value 5/pi
    r = constant (optional), with default value 6
    s = constant (optional), with default value 10
    t = constant (optional), with default value 1/(8*pi)
  """

  term1 = a * (y - b*x**2 + c*x - r)**2
  term2 = s*(1.0 - t) * np.cos(x)

  F = term1 + term2 + s

  return F

