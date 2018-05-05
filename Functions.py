import sympy as sy
import numpy as np
from sympy import Matrix

from sympy import derive_by_array

# it would have to have a large mu and doubling it
from scipy.misc import derivative
from sympy.core.multidimensional import vectorize
from scipy.linalg import norm
class Functions:



    def __init__(self, l, xss, xee, tubeLength):
        #self.x = symbols('x')
        self.xs = sy.symbols('xs')  # starting point of the tube
        self.xe = sy.symbols('xe')  # ending point of the tube
        #self.f = sy.Function('f')
        self.Eprime = 2.5 * 10 ** 6
        self.e3 = sy.symbols('e3')
        self.e33 = np.array([0,0,1])
        #self.e3 = Matrix([0, 0, 1])
        self.gravity = 9.82
        self.l = sy.symbols('l')
        self.m = sy.symbols('m')
        self.length = l
        self.xss = xss
        self.xee = xee
        self.mass = 0.28 * tubeLength


    #x = symbols('x');
    def setX(self,xs, xe):
        self.xss = xs
        self.xee = xe

    def function(self, f):
        func = sy.lambdify([self.xe, self.xs, self.l, self.m, self.e3], f, "numpy")
        #g = [sy.diff(f, x) for x in [self.xs, self.xs]]
        g = derive_by_array(f, [self.xs, self.xe])
        grad = sy.lambdify([self.xe, self.xs, self.l, self.m, self.e3], g, "numpy")
        h = derive_by_array(g, [self.xs, self.xe])
        hess = sy.lambdify([self.xe, self.xs, self.l, self.m, self.e3], h, "numpy")

        f1 = np.sum(func(self.xss,self.xee,self.length, self.mass, self.e33))
        g1 = np.sum(grad(self.xss,self.xee,self.length, self.mass, self.e33),1)
        h1 = hess(self.xss,self.xee,self.length, self.mass, self.e33)

        print g

        return (f1, g1, h1)
    #compute the energy function E (x), its gradient, and its Hessian matrix for any given x
    #@vectorize(0, 1)
    def totalEnergy(self,f1, f2):
        [f1, g1, h1] = self.function(f1)
        [f2, g2, h2] = self.function(f2)


        #f = m*g * (1./2.) * sy.E * (self.xs + self.xe)
        #g = [sy.diff(f, x) for x in [self.xs, self.xs]]
        #h = sy.hessian(f, [self.xs, self.xe])
        #return f
        #totalEnergy = m*g * (1/2) * math.e * (x1 + x2) + (1/2) * E * (max(0, abs(x1 - x2) - l))

    #m is an array
    #sy.transpose(self.e3)
    def gravitationalEnergy(self):
        f = self.mass * self.gravity * (1.0/2.0) * np.dot( self.e3, (self.xs + self.xe))
        #g = [sy.diff(f, x) for x in [self.xs, self.xs]]
            #m * g (1.0/2.0)  * transpose(e)
        #h = sy.hessian(f, [self.xs, self.xe])# should be 0;
        return f
    def elasticEnergy(self):
        E = float(self.Eprime / self.length)
        f = (1./2.) * E *(max(0, norm(self.xs - self.xe) - self.length))
        #g = [sy.diff(f, x) for x in [self.xs, self.xs]]
        #h = sy.hessian(f, [self.xs, self.xe])  # should be 0;d
        return f


    #def model(self):
#, [2, 2, 4],[2, 2, 4]
#, [4, 5, 6],[1, 20, 5]
#, 40,40
#, [2, 2, 4],[2, 2, 4], [1, 1, 1], [2, 2, 4],[2, 2, 4]
#, [2, 2, 4],[2, 2, 4], [1, 1, 1], [2, 2, 4],[2, 2, 4]
#,2,2,2,2,2
s = np.array([[1, 1, 1]], dtype=float)
e = np.array([[1, 1, 1]], dtype=float)
m = np.array(2, dtype=float)
l = np.array(1, dtype=float)
#e3 = np.array([[0, 0, 1],[0, 0, 1],[0, 0, 1], [0, 0, 1],[0, 0, 1],[0, 0, 1]], dtype=float)
f = Functions(l,s ,e ,m )

[f1, g1, h1] = f.function(f.gravitationalEnergy())
print f.gravitationalEnergy()


#print np.matmul(s,e3)
print f1
print g1
print h1
#print g1.shape
#print h1
#[m ,n] = np.size(h1)
#print m
#print n








