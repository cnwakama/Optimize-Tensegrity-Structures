import sympy as sy
import numpy as np
from sympy import Matrix

from sympy import derive_by_array

# it would have to have a large mu and doubling it
from scipy.misc import derivative
from sympy.core.multidimensional import vectorize
#from scipy.linalg import norm
from numpy.linalg import norm


class Functions:
    def __init__(self, tube, cable):
        # self.x = symbols('x')
        self.tube = tube
        self.cable = cable
        self.xs = sy.symbols('xs')  # starting point of the tube
        self.xe = sy.symbols('xe')  # ending point of the tube
        self.l = sy.symbols('l')
        self.m = sy.symbols('m')
        self.e3 = sy.symbols('e3')
        # self.f = sy.Function('f')
        self.Eprime = 2.5 * 10 ** 6
        self.e33 = np.array([0, 0, 1])
        self.gravity = 9.82
        self.length = cable[:, 6]
        self.xst = tube[:, 0:3]
        self.xet = tube[:, 3:6]
        self.xsc = cable[:, 0:3]
        self.xec = cable[:, 3:6]
        self.tubeLength = tube[:, 6]
        self.mass = 0.28 * self.tubeLength
        self.difft = self.xst - self.xet
        self.diffc = self.xsc - self.xec
        self.normt = norm(self.difft, axis=-1)
        self.normc = norm(self.diffc, axis=-1)
        self.cout1 = 0
        self.cout2 = 0

    def lagrangian(self):
        f = self.totalEnergy()
    # x = symbols('x');
    def updateT(self):
        [m,_] = self.xst.shape
        if self.cout1 == m:
            self.cout1 = 0
        else:
            self.cout1 += 1
    def updateC(self):
        [m, _] = self.xsc.shape
        if self.cout2 == m:
            self.cout2 = 0
        else:
            self.cout2 += 1

    def function(self, f):
        func = sy.lambdify([self.xe, self.xs, self.l, self.m, self.e3], f, "numpy")
        # g = [sy.diff(f, x) for x in [self.xs, self.xs]]
        g = derive_by_array(f, [self.xs, self.xe])
        grad = sy.lambdify([self.xe, self.xs, self.l, self.m, self.e3], g, "numpy")
        h = derive_by_array(g, [self.xs, self.xe])
        hess = sy.lambdify([self.xe, self.xs, self.l, self.m, self.e3], h, "numpy")

        f1 = np.sum(func(self.xss, self.xee, self.length, self.mass, self.e33))
        g1 = np.sum(grad(self.xss, self.xee, self.length, self.mass, self.e33), 1)
        h1 = hess(self.xss, self.xee, self.length, self.mass, self.e33)

        print g

        return f1, g1, h1

    def setX(self, xss, xee):
        self.xss = xss
        self.xee = xee

    def setMass(self, tubeLength):
        self.tubeLength = tubeLength

    def setLength(self, length):
        self.length = length

    # compute the energy function E (x), its gradient, and its Hessian matrix for any given x
    # @vectorize(0, 1)
    '''
        tube is a m1,3 matrix contain the xs, xe, and tube length of the tube
        cable is a m2,3 matrix contain the xs, xe, and cable length of the cable
    '''
    def totalEnergy(self, f1, tube, cable):
        Tf = 0
        Tg = 0
        Th = 0
        [f1, g1, h1] = self.function(f1)
        [f2, g2, h2] = self.elasticEnergy()
        [numTube, _] = tube.shape
        [numCable, _] = cable.shape
        self.cout1 = 0
        for x in range(numTube):
            Tf += f1
            Tg += g1
            Th += h1


            #self.setX(cable[x, 0], cable[x, 1])
            #self.setMass(cable[x, 2])

            [f1, g1, h1] = self.function(f1)

        self.cout2= 0
        for x in range(numCable):
            Tf += f2
            # fix values so grad and hession are zero
            if x == 0 or x == 1:
                #self.setX(cable[x + 1, 0], cable[x + 1, 1])
                #self.setLength(cable[x + 1, 2])
                continue
            Tg += g2
            Th += h2

            #self.setX(cable[x + 1, 0], cable[x + 1, 1])
            #self.setLength(cable[x + 1, 2])
            [f2, g2, h2] = self.elasticEnergy()




            # f = m*g * (1./2.) * sy.E * (self.xs + self.xe)
            # g = [sy.diff(f, x) for x in [self.xs, self.xs]]
            # h = sy.hessian(f, [self.xs, self.xe])
            # return f
            # totalEnergy = m*g * (1/2) * math.e * (x1 + x2) + (1/2) * E * (max(0, abs(x1 - x2) - l))

    # m is an array
    # sy.transpose(self.e3)
    def gravitationalEnergy(self):
        f = self.mass * self.gravity * (1.0 / 2.0) * np.dot(self.e3, (self.xs + self.xe))
        # g = [sy.diff(f, x) for x in [self.xs, self.xs]]
        # m * g (1.0/2.0)  * transpose(e)
        # h = sy.hessian(f, [self.xs, self.xe])# should be 0;
        return f

    def elasticEnergy(self):
        g = np.zeros((2, 1))
        h = np.zeros((2, 2))
        E = float(self.Eprime / self.length)
        # q =pow(self.xss - self.xee)
        mm = (max(0, self.norm - self.length))
        f = (1. / 2.) * E * mm
        if f == 0:
            return f,g,h

        g.fill((1. / 2.) * E * mm * self.diff / self.norm)

        m = (1/self.norm) - (np.power(self.diff, 2) / (np.power(self.norm), 3))
        h.fill((-1./2.) * E*self.length * m)
        #(-1./2.) * E*self.length *

        # g = [sy.diff(f, x) for x in [self.xs, self.xs]]
        # h = sy.hessian(f, [self.xs, self.xe])  # should be 0;d
        return f, g, h

    def constraints(self):
        num = np.matmul(np.ones((1, 3)) , np.reshape(np.transpose(self.difft[self.cout1, :]), (3, 1)))
        norm = np.reshape(self.normt[self.cout1], (1, 1))
        numb = np.power(num, 2)
        #print "V:"
        #print self.normt
        dem = np.power(norm, 3)
        e = (numb / dem)
        d = 1/norm
        #print "Checks"
        #print np.reshape(np.transpose(self.difft[self.cout1, :]), (3, 1)).shape
        #print num.shape
        #print self.difft[self.cout1, :]
        #print num
        #print norm
        #print numb
        #print dem
        #print d
        #print e
        f = np.power(norm, 2) - np.power(self.tubeLength[self.cout1], 2)
        g = np.array([num[0, 0] /self.normt[self.cout1],
                      -num[0, 0] / self.normt[self.cout1]], dtype=float)
        h = np.zeros((2, 2), float)
        #print "d-e"
        #print d-e
        de = (d - e)[0, 0]
        #print de
        h.fill(de)
        h[1, 0] = - h[1, 0]
        h[0, 1] = - h[0, 1]

        #print "Function and Gradient"
        #print f
        #print g
        #print self.normt[self.cout1]
        #print self.tubeLength[self.cout1]
        #print g.shape

        return f, g, h


        # def model(self):


# , [2, 2, 4],[2, 2, 4]
# , [4, 5, 6],[1, 20, 5]
# , 40,40
# , [2, 2, 4],[2, 2, 4], [1, 1, 1], [2, 2, 4],[2, 2, 4]
# , [2, 2, 4],[2, 2, 4], [1, 1, 1], [2, 2, 4],[2, 2, 4]
# ,2,2,2,2,2
s = np.array([[1, 1, 1]], dtype=float)
e = np.array([[1, 1, 1]], dtype=float)
m = np.array(2, dtype=float)
l = np.array(1, dtype=float)
# e3 = np.array([[0, 0, 1],[0, 0, 1],[0, 0, 1], [0, 0, 1],[0, 0, 1],[0, 0, 1]], dtype=float)
#f = Functions(l, s, e, m)

#[f1, g1, h1] = f.function(f.gravitationalEnergy())
# [f2, g2, h2] = f.function(f.elasticEnergy())
#print f.gravitationalEnergy()

# print np.matmul(s,e3)
#print f1
#print g1
#print h1
# print g1.shape
# print h1
# [m ,n] = np.size(h1)
# print m
# print n
