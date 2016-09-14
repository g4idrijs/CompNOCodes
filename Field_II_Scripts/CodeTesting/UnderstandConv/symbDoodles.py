from sympy import *
a, u, k, L = symbols("a u k L")
u =  L - k - 1 # Define k, solve for u

print(k)
print("u=",u)
print("u(0)=",u.subs({k:0}))
print("u(2L-4-a)=",u.subs({k:(2*L-4-a)}))