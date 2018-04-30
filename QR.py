import numpy as np
class QR:

    def problem2(self, B, r):
        [m, n] = np.shape(B)
        np.linalg.qr(B, mode='complete') #maybe reduce
        return z
    % Returns
    argmin_z | | z | | subject
    to
    B'*z = r
    % Assumes
    columns
    of
    B
    are
    linearly
    independent
    % Uses
    QR
    factorization
    [m, n] = size(B);
    [Q, R] = qr(B);
    Q1 = Q(:, 1:n);
    Q2 = Q(:, n + 1:m);
    z = Q1 * (R(1:n, 1:n)' \ r);