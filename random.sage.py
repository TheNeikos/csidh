#coding: utf8

# This file was *autogenerated* from the file random.sage
from sage.all_cmdline import *   # import sage library

_sage_const_3 = Integer(3); _sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_4 = Integer(4); _sage_const_9 = Integer(9); _sage_const_374 = Integer(374); _sage_const_587 = Integer(587)#!/usr/bin/env sage

proof.all(False)

ls = list(primes(_sage_const_3 , _sage_const_374 )) + [_sage_const_587 ] # Elkies primes
print ls
p = _sage_const_4  * prod(ls) - _sage_const_1 
base = GF(p)
A = _sage_const_0 
B = _sage_const_1 
E = EllipticCurve(base, [_sage_const_0 , A, _sage_const_0 , B, _sage_const_0 ])
R = E.random_element()
print R.xy()
print (_sage_const_9  * R).xy()

