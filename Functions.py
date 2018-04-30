import sympy as sy
import numpy as np

from sympy import derive_by_array

from scipy.misc import derivative
from sympy.core.multidimensional import vectorize
from scipy.linalg import norm
class Functions:

    def __init__(self):
        #self.x = symbols('x')
        self.xs = sy.symbols('xs')  # starting point of the tube
        self.xe = sy.symbols('xe')  # ending point of the tube
        #self.f = sy.Function('f')

    #x = symbols('x');

    def function(self, f, xs, xe):
        f1 = sy.lambdify([self.xe, self.xs], f, "numpy")
        #g = [sy.diff(f, x) for x in [self.xs, self.xs]]
        g = derive_by_array(f, [xs, xe])
        h = sy.hessian(f, [self.xs, self.xe])
        return f1, g, h
    #compute the energy function E (x), its gradient, and its Hessian matrix for any given x
    #@vectorize(0, 1)
    def totalEnergy(self, m, g, f1, f2):
        f = m*g * (1./2.) * sy.E * (self.xs + self.xe)
        #g = [sy.diff(f, x) for x in [self.xs, self.xs]]
        #h = sy.hessian(f, [self.xs, self.xe])
        return f
        #totalEnergy = m*g * (1/2) * math.e * (x1 + x2) + (1/2) * E * (max(0, abs(x1 - x2) - l))

    #m is an array
    def gravitationalEnergy(self, m, g, e, xs, xe):
        f = m * g * (1.0/2.0) * sy.transpose(e) * (xs + xe)
        #g = [sy.diff(f, x) for x in [self.xs, self.xs]]
            #m * g (1.0/2.0)  * transpose(e)
        #h = sy.hessian(f, [self.xs, self.xe])# should be 0;
        return f
    def elasticEnergy(self, starE, l, xs, xe):
        E = float(starE / l)
        f = (1./2.) * E *(max(0, norm(xs - xe) - l))
        #g = [sy.diff(f, x) for x in [self.xs, self.xs]]
        #h = sy.hessian(f, [self.xs, self.xe])  # should be 0;d
        return f
    def hessXL(self):


