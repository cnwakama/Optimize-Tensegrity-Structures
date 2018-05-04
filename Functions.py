import sympy as sy
import numpy as np

from sympy import derive_by_array
# it would have to have a large mu and doubling it
from scipy.misc import derivative
from sympy.core.multidimensional import vectorize
from scipy.linalg import norm
class Functions:



    def __init__(self):
        #self.x = symbols('x')
        self.xs = sy.symbols('xs')  # starting point of the tube
        self.xe = sy.symbols('xe')  # ending point of the tube
        #self.f = sy.Function('f')
        self.Eprime = 2.5 * 10 ** 6
        self.e3 = np.array([[0], [0], [1]])
        self.gravity = 9.82
        self.length = sy.symbols('l')

    #x = symbols('x');

    def function(self, f, xs, xe, l):
        func = sy.lambdify([self.xe, self.xs, self.length], f, "numpy")
        #g = [sy.diff(f, x) for x in [self.xs, self.xs]]
        g = derive_by_array(f, [self.xs, self.xe])
        grad = sy.lambdify([self.xe, self.xs, self.length], g, "numpy")
        h = sy.hessian(f, [self.xs, self.xe])
        hess = sy.lambdify([self.xe, self.xs, self.length], h, "numpy")

        f1 = func(xs,xe,l)
        g1 = grad(xs,xe,l)
        h1 = hess(xs,xe,l)

        return f1, g1, h1
    #compute the energy function E (x), its gradient, and its Hessian matrix for any given x
    #@vectorize(0, 1)
    def totalEnergy(self, m, g, f1, f2):
        f = m*g * (1./2.) * sy.E * (self.xs + self.xe)
        #g = [sy.diff(f, x) for x in [self.xs, self.xs]]
        #h = sy.hessian(f, [self.xs, self.xe])
        return f
        #totalEnergy = m*g * (1/2) * math.e * (x1 + x2) + (1/2) * E * (max(0, abs(x1 - x2) - l))

    #m is an array
    def gravitationalEnergy(self, m, g, xs, xe):
        f = m * g * (1.0/2.0) * np.transpose(self.e3) * (self.xs + self.xe)
        #g = [sy.diff(f, x) for x in [self.xs, self.xs]]
            #m * g (1.0/2.0)  * transpose(e)
        #h = sy.hessian(f, [self.xs, self.xe])# should be 0;
        return f
    def elasticEnergy(self, xs, xe):
        E = float(self.Eprime / self.length)
        f = (1./2.) * E *(max(0, norm(xs - xe) - l))
        #g = [sy.diff(f, x) for x in [self.xs, self.xs]]
        #h = sy.hessian(f, [self.xs, self.xe])  # should be 0;d
        return f


    #def model(self):







