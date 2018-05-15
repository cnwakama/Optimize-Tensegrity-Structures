#from TrustRegion import TrustRegion
#from Functions import Functions
import numpy as np
from numpy.linalg import norm

from TrustRegion import TrustRegion

f = 0
g = 0
h = 0
c = 0
gc= 0
m = 0
mu = 0
x = np.array([[1, 2, 0], [-3, 0, 5], [-2, 0, 0], [0, 5, 3], [-1, 2, 3], [0, 2, 3], [0, -5, 3], [0, 5, 3]])
xs = np.zeros((4,3))
xe = np.zeros((4,3))
tubeLength = np.array([3.9751,3.4417, 3.4417, 2.3622])
cableLength = np.array([1.4986, 1.4859, 1.4859, 2.4257, 2.95275, 2.9972, 3.2893,
                        1.4986, 1.4859, 1.4859, 2.4257, 2.95275, 2.9972, 3.2893])
xs[0, :] = x[0]
xs[1, :] = x[2]
xs[2, :] = x[4]
xs[3, :] = x[6]
xe[0, :] = x[1]
xe[1, :] = x[3]
xe[2, :] = x[5]
xe[3, :] = x[7]
print xe.shape
print xe.shape

#tube = np.zeros((4, 3))
#print tube[:, 0].shape

#tubeLength.reshape((4,1))
print tubeLength.shape
tube = np.concatenate([xs,xe,tubeLength[:, None]], axis = 1) #the tube matrix
#tube = np.array([[xs, xe, tubeLength]])
xs = np.zeros((14,3))
xe = np.zeros((14,3))
xs[0, :] = x[0]
xs[1, :] = x[0]
xs[2, :] = x[2]
xs[3, :] = x[5]
xs[4, :] = x[5]
xs[5, :] = x[1]
xs[6, :] = x[1]

xs[7, :] = x[0]
xs[8, :] = x[0]
xs[9, :] = x[4]
xs[10, :] = x[3]
xs[11, :] = x[3]
xs[12, :] = x[2]
xs[13, :] = x[2]

xe[0, :] = x[2]
xe[1, :] = x[6]
xe[2, :] = x[7]
xe[3, :] = x[6]
xe[4, :] = x[7]
xe[5, :] = x[3]
xe[6, :] = x[2]

xe[7, :] = x[4]
xe[8, :] = x[7]
xe[9, :] = x[6]
xe[10, :] = x[7]
xe[11, :] = x[6]
xe[12, :] = x[5]
xe[13, :] = x[4]
cable = np.concatenate([xs,xe,cableLength[:, None]], axis = 1) #the cable matrix
    #np.array([xs, xe, cableLength[:, None]], dtype=object)
print tube
print cable
print tube.shape
print cable.shape

#trustRegion.trustalg(200, 50, 1/5)
xst = tube[:, 0:3]
xet = tube[:, 3:6]
difft = (xst - xet)
normt = norm(difft,  axis=-1)

print difft
print normt
print np.reshape(normt, (4, 1))
print difft.shape

trustRegion = TrustRegion(tube, cable, 20)
print "Constraint and Gradient"
print trustRegion.constraint
print trustRegion.gradConstraint
[phat, N] = trustRegion.case2()
print phat
print N

trustRegion.case3(N, phat)
