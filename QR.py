import numpy as np
from numpy.linalg import inv

'''
% function z = problem1(B,r)
%
% Returns argmin_z ||B*z-r||
% Assumes columns of B are linearly independent
% Uses QR factorization
[m,n] = size(B);
[Q,R] = qr(B);
Q1 = Q(:,1:n);
z = R(1:n,1:n) \ (Q1'*r);
'''
def problem1(B,r):
    [m, n] = B.shape
    [q, R] = np.linalg.qr(B, mode='complete')  # maybe reduce
    q1 = q[:,0:n-1]
    z = (inv(R[0:n,0:n])) * np.transpose(q1) * r
    return z

def problem2(B, r):
    [m, n] = B.shape
    print m
    print n
    [q, R] = np.linalg.qr(B, mode='complete') #maybe reduce
    q1 = q[:, 0:n]
    q2 = q[:, n:m]
    print q
    print q1 #4x2
    print q2 #4x2
    print r #4x1
    print inv(R[0:n, 0:n]) #2x2
    t = np.transpose((inv(R[0:n, 0:n])))
    print t
    #z = np.matmul(np.matmul(q1, t).transpose(),  r)
    z = np.matmul(q1, t) * r
    #print z.shape

    return [z, q2]


    '''
        function [z,Q2] = problem2(B,r)
        % function [z,Q2] = problem2(B,r)
        %
        % Returns argmin_z ||z|| subject to B'*z = r
        % Assumes columns of B are linearly independent
        % Uses QR factorization
        [m,n] = size(B);
        [Q,R] = qr(B);
        Q1 = Q(:,1:n);
        Q2 = Q(:,n+1:m);
        z = Q1*(R(1:n,1:n)' \ r);
    '''