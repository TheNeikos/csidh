#!/usr/bin/env sage
#coding: utf8

proof.all(False)

ls = list(primes(3, 374)) + [587] # Elkies primes
p = 4 * prod(ls) - 1
base = GF(p)
A = 0
B = 1
E = EllipticCurve(base, [0, A, 0, B, 0])
R = E.random_element()
print R.xy()
print (9 * R).xy()
