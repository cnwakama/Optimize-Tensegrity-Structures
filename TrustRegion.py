from QR import *
from Functions import *
import numpy as np


class TrustRegion:
    def __init__(self, f, g, h, c, gc, m, mu):
        self.p = 0
        self.func = f
        self.grad = g
        self.hess = h
        self.contraint = c  # functions
        self.gradContraint = gc  # function
        self.function = Functions()
        self.merit = m
        self.mu = mu

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

        normB = norm(B, 'fro')
        if normB == 0.0:
            if norm(g) == 0:
                p = np.zeros(n, 1)
                at_bdry = 0
            else:
                p = -g / norm(g, 2) * Delta
                at_bdry = 1
            return

        # Compute Geshgorin bound for lambda1
        lambda1 = 0
        for i in range(0, n):
            temp = B(i, i) - sum(abs(B[i, 1:(i - 1)]) - sum(abs(B[i, (i + 1):n])))
            if temp < lambda1:
                lambda1 = temp



        lambda_hi = norm(g, 2) / Delta - lambda1
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
            if norm(p, 2) < Delta:
                lambda_hi = lambda_lo
                lambda_lo = (lambda_lo + lambda_old) / 2
            else:
                if norm(p, 2) == Delta:
                    lambda_hi = lambda_lo
                break

        if lambda_hi == 0:
            at_bdry = 0
            return

        return p, at_bdry


# f is the energy equation
# min lambda for ||grad f - grad cT * lambda||


def minLambda(self):
    lda = problem1(self.contraint, self.grad)
    return lda


'''
min ||p-hat|| for p-hat subject to c(x)+ grad c(x)p-hat
'''
def case2(self):
    [p_hat, _] = problem2(self.gradContraint, -self.contraint)
    return p_hat


# def case3(self):

# min d of ||d|| subject to c(x+p)+grad c(x+p)d = 0;
def case4(self, xs, xe, p):
    [c, gC, _] = self.function.function(self, self.contraint, xs, xe)
    [d, _] = problem2(gC, -c)
    return d


def calrol(self, xs, xe, lbd):
    d = 0
    p = 0
    hL = 0
    [f, g, h] = self.function.function(self.func, xs, xe)
    [f1, g1, h] = self.function.function(self.func, xs + d + p, xe + d + p)
    [c, gC, _] = self.function.function(self, self.contraint, xs, xe)
    [c1, gC1, _] = self.function.function(self, self.contraint, xs, xe)
    fmu = f + norm(lbd, 2) * norm(c)
    fu = f1 + norm(lbd, 2) * norm(c1)

    r = (fmu - fu) / ((self.model(0, f, g, hL) + norm(lbd, 2) * norm(c)) - (self.model(p, f, g, hL) + norm(lbd, 2) *
                                                                            norm(c + gC * p)))


def trustRegion(self, iteration, de, n, xs, xe):
    pr = 0
    p = 0
    d = 0
    # Obtain pk by (approximately) solving (4.3);
    # Evaluate ppk from (4.4);
    delta = []
    x = []
    x.append([xs, xe])
    delta.append(de)
    for k in range(iteration):
        if pr < 1. / 4.:
            delta.append((1. / 4.) * delta[len(delta) - 1])
        else:
            if pr > 3. / 4. and abs(norm(pr) - delta[len(delta) - 1]) < 0.00001:
                delta.append(min(2 * delta[len(delta) - 1], de))
            else:
                delta.append(delta[len(delta) - 1])
        if pr > n:
            x.append(x[len(x) - 1 + p + d])
        else:
            x.append(x[len(x) - 1])


def model(self, p, f, g, hL):
    model = f + np.transpose(g) * p(1. / 2.) * np.transpose(p) * hL * p
    return model
