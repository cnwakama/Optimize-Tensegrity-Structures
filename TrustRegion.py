
import numpy as np

from Functions import Functions
from QR import *



class TrustRegion:
    def __init__(self, tube, cable, delta):
        self.p = 0
        self.func = 0
        self.grad = 0
        self.hess = 0
        self.constraint = np.zeros((4, 1))  # functions
        self.gradConstraint = np.zeros((4, 2))  # function
        self.function = Functions(tube, cable)
        #self.merit = m
        #self.mu = mu
        self.hessLag = 0
        self.tube = tube
        self.cable = cable
        self.delta = delta
        [self.constraint, self.gradConstraint] = self.setContraintMatrix()

    def solve_tr(self, g, B, Delta, tol, trace):
        p = 0
        at_bdry = 0
        n = len(g)
        [nrows, ncols] = np.size(B)

        if nrows != n or ncols != n:
            print('solve_tr: incompatible sizes\n')
            p = np.zeros(n, 1)
            at_bdry = 0
            return

        if Delta <= 0:
            print('solve_tr: Delta must be positive\n')
            p = np.zeros(n, 1)
            at_bdry = 0
            return

        normB = np.linalg.norm(B, 'fro')
        if normB == 0.0:
            if np.linalg.norm(g) == 0:
                p = np.zeros(n, 1)
                at_bdry = 0
            else:
                p = -g / np.linalg.norm(g, 2) * Delta
                at_bdry = 1
            return

        # Compute Geshgorin bound for lambda1
        lambda1 = 0
        for i in range(0, n):
            temp = B(i, i) - sum(abs(B[i, 1:(i - 1)]) - sum(abs(B[i, (i + 1):n])))
            if temp < lambda1:
                lambda1 = temp



        lambda_hi = np.linalg.norm(g, 2) / Delta - lambda1
        lambda_lo = 0
        lambda_old = 0
        at_bdry = 1

        p = 0

        while abs(lambda_lo - lambda_hi) > tol * normB:
            C = B + lambda_lo * np.eye(n)
            try:
                R = np.linalg.cholesky(C)
            except:
                # Cholesky factorization failed
                if trace !=0:
                    print('solve_tr: Cholesky factorization failed\n')
                lambda_old = lambda_lo
                lambda_lo = (lambda_lo + lambda_hi) / 2
                continue
            p = - np.linalg.inv(C) *g
            if trace != 0:
                print('solve_tr: ||p|| = %g, Delta = %g\n', norm(p, 2), Delta)
            if np.linalg.norm(p, 2) < Delta:
                lambda_hi = lambda_lo
                lambda_lo = (lambda_lo + lambda_old) / 2
            else:
                if np.linalg.norm(p, 2) == Delta:
                    lambda_hi = lambda_lo
                break

        if lambda_hi == 0:
            at_bdry = 0
            return

        return p, at_bdry


# f is the energy equation
# min lambda for ||grad f - grad cT * lambda||


    def minLambda(self):
        lda = problem1(self.constraint, self.grad)
        return lda

    def setContraintMatrix(self):
        [m, _] = self.tube.shape
        constraint = np.zeros((m, 1))
        gradC = np.zeros((m, 2))
        self.function.cout1 = 0
        for x in range(m):
            [f, g, h] = self.function.constraints()
            constraint[x] = f
            gradC[x, :] = g
            self.function.updateT()
        return constraint, gradC
    '''
        min ||p-hat|| for p-hat subject to c(x)+ grad c(x)p-hat
    '''
    def case2(self):
        [m, _] = self.constraint.shape
        print "Size"
        print self.gradConstraint.shape
        self.function.cout1 = 0
        p_hatM = np.zeros((m, 2))
        nM = np.zeros((m, 1))
        [p_hat, N] = problem2(self.gradConstraint, -self.constraint)
        p_hatM = p_hat
        nM = N
        '''
        for x in range(m):
            print "Size"
            #print self.gradConstraint[x,:].reshape((2,1)).shape
            #print self.gradConstraint[x, 0:2].shape
            #.reshape((2, 1))
            [p_hat, N] = problem2(self.gradConstraint[x, :].reshape((2, 1)), -self.constraint[x])
            p_hatM[x, :] = p_hat.reshape((1, 2))
            nM[x] = N
        '''
        return p_hatM, nM

    '''
    min q m(bp+Nq) subject to
    ||q|| <= (delta_2 - q_2)^1/2
    '''
    def case3(self, N, p):
        B = np.transpose(np.transpose(N) * (self.g + self.hessLag + p))
        r = np.transpose(self.g)*p + (1.0/2.0) * np.transpose(p) * self.hessLag * p
        [q, _] = self.solve_tr(B, r, self.delta)
        return q



    # min d of ||d|| subject to c(x+p)+grad c(x+p)d = 0;
    def case4(self, p):
        [m, _] = self.constraint.shape
        p_hatM = None
        dM = None
        for x in range(m):
            [d, _] = problem2(self.gradConstraint[x, :], -self.constraint[x])
            if x == 0:
                dM = d
            else:
                dM = np.concatenate([dM, d], axis=1)
        return dM
        #[c, gC, _] = self.function.function(self, self.constraint, xs, xe)



    def calrol(self, lbd):
        d = 0
        p = 0
        hL = 0
        [f, g, h] = self.function.function(self.func)
        #, xs + d + p, xe + d + p
        [f1, g1, h] = self.function.function(self.func)
        [c, gC, _] = self.function.function(self)
        [c1, gC1, _] = self.function.function(self)
        #norm(lbd, 2)
        fmu = f +  self.mu * np.linalg.norm(c)
        fu = f1 +  self.mu * np.linalg.norm(c1)

        r = (fmu - fu) / ((self.model(0, f, g, hL) + norm(lbd, 2) * norm(c)) - (self.model(p, f, g, hL) + norm(lbd, 2) *
                                                                                norm(c + gC * p)))

        return r


    def getPk(self):
        pk = 0
        r = 0
        [ph, N] = self.case2()
        q = self.case3(N,ph)
        p = 0
        self.case4()
        self.calrol()

        return pk,r


    def trustalg(self, iteration, de, n):
        pr = 0
        p = 0
        d = 0
        delta = []
        x = []
        #x.append([xs, xe])
        delta.append(de)
        for k in range(iteration):
            # Obtain pk by (approximately) solving (4.3);
            # Evaluate ppk from (4.4);
            [p, pr] = self.getPk()
            if pr < 1.0/ 4.0:
                delta.append((1. / 4.) * delta[len(delta) - 1])
            else:
                if pr > 3. / 4. and abs(np.linalg.norm(pr) - delta[len(delta) - 1]) < 0.00001:
                    delta.append(min(2 * delta[len(delta) - 1], de))
                else:
                    delta.append(delta[len(delta) - 1])
            if pr > n:
                x.append(x[len(x) - 1 + p + d])
            else:
                x.append(x[len(x) - 1])

    '''
    def model(self, p, f, g, hL):
        model = f + np.transpose(g) * p(1. / 2.) * np.transpose(p) * hL * p
        return model
    '''