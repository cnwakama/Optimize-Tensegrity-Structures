import numpy as np

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
    [m, n] = np.shape(B)
    [q, r] = np.linalg.qr(B, mode='complete')  # maybe reduce
    q1 = q[:,[0,n-1]]
    z = np.inv(r[1:n,1:n]) * np.transpose(q1) *
    return z

def problem2(B, r):
    [m, n] = np.shape(B)
    [q,r] = np.linalg.qr(B, mode='complete') #maybe reduce
    q1 = q[:, 1,n]
    q2 = q[:, n+1:m]
    z = q1 * (r[1:n,1:n])
    return [z,q2]


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