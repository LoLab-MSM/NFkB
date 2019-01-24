import sympy as sy
L = sy.symbols('L')
M = sy.Matrix([[-1.4 - L, .439, -.4206, 0],
               [1, -39.561 - L, .4206, 0],
               [-.4206, 0, -1.4- L, .439],
               [.4206, 0, 1, -39.561 - L]])
print('Characteristic Polynomial: P(L) = ', M.det())
print('WEEE')
print('Eigenvalues: ', sy.solve(M.det()))
print('DONE')
