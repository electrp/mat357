# a is always 1, returns b, c, d
def make_lagrange_cubic(x0, x1, x2, x3, y0, y1, y2, y3) -> tuple[float, float, float]:
    # Precompute differences
    xd10 = x1 - x0
    xd20 = x2 - x0
    xd30 = x3 - x0 
    xd21 = x2 - x1
    xd31 = x3 - x1
    xd32 = x3 - x2
    
    # Precompute divisors
    l0d = -xd10 * -xd20 * -xd30
    l1d =  xd10 * -xd21 * -xd31
    l2d =  xd20 *  xd21 * -xd32
    l3d =  xd30 *  xd31 *  xd32

    # Helpers
    xmult = x0 * x1 * x2 * x3
    xaccum = x0 + x1 + x2 + x3
    x01 = x0 * x1
    x02 = x0 * x2
    x03 = x0 * x3
    x12 = x1 * x2
    x13 = x1 * x3
    x23 = x2 * x3

    # Final results 
    b = \
        (-xaccum + x0) / l0d + \
        (-xaccum + x1) / l1d + \
        (-xaccum + x2) / l2d + \
        (-xaccum + x3) / l3d 
    c = \
        (x12 + x23 + x13) / l0d + \
        (x02 + x23 + x03) / l1d + \
        (x01 + x13 + x03) / l2d + \
        (x01 + x02 + x12) / l3d
    d = \
        x12 * x3 / l0d + \
        x02 * x3 / l1d + \
        x01 * x3 / l2d + \
        x01 * x2 / l3d

    return b, c, d

def depress_cubic_const_a(b, c, d) -> tuple[float, float]:
    h1 = b / 3
    h2 = h1 * h1
    h3 = h1 * h2

    c1 = -h2 * 3 + c
    d1 = 2 * h3 - c * b / 3 + d
    return c1, d1 

# Gives c and d
def depress_cubic(a, b, c, d) -> tuple[float, float]:
    # Scale so a == 1, x = t - b/3
    b1 = b / a
    c1 = c / a
    d1 = d / a
    return depress_cubic_const_a(b1, c1, d1)
    
    # h1 = b /(3 * a)
    # h2 = h1 * h1
    # h3 = h1 * h2

    # c1 = -h2 * 3 + c / a
    # d1 = 2 * h3 - c * b / (3 * a * a) + d / a
    # return c1, d1

# https://en.wikipedia.org/wiki/Cubic_equation

""" 
Algorithm steps:

something until we have 4 points

then:
- using those 4 points, make a cubic lagrange polynomial
- depress the cubic (will)
- use the formula for roots of a depressed cubic
- pick a root
    - this is our new x point

the formula for the roots of a depressed cubic:
x**3 + p*x + q = 0  # this is a depressed cubic

let x = a - b
then x**3 = (a-b)**3 = a**3 - b**3 - c*a*b*(a-b)
x**3 = 3*a*b*x-(a**3-b**3) = 0

p = 3*a*b, -q = a**3-b**3
b = p/3a, -q = a**3-(p**3/(2+a**3))
-a**3*q = (a**3)**2 - p**3/27
(a**3)**2 + a**3*q - p**3/27

let t = a**3
t**2 +_q*t - p**3 /27 = 0
now you have a quadratic in terms of t
t = -q +/- sqrt(q**2 + (4*p**3)/27)

a**3 = -q/2 +/- sqrt(q**2 /4 + p**3 /27)
a = cube_root( )... gosh so much typing
"""

import math
def root(p_, q_):
    qot = q_/2
    big_term = math.sqrt(q_**2/4 + p_**3/27)
    return (-qot + big_term)**(1/3) + (-qot - big_term)**(1/3)


"""
Discriminant Rules:

Δ = 18abcd – 4b³d + b²c² – 4ac³ – 27a²d².
If Δ > 0, the cubic has three distinct real roots
If Δ < 0, the cubic has one real root and two non-real complex conjugate roots.
else it has a multiple root
whahhahahhahat is a multiple root?
"The multiplicity of a root is the number of occurrences of this root in the complete factorization of the polynomial, by means of the fundamental theorem of algebra." (https://en.wikipedia.org/wiki/Multiplicity_(mathematics)#Multiplicity_of_a_root_of_a_polynomial)t
If our disc is 0, because we're using depressed cubics we have one simple root and a double root
(omg furthermore) bc we're using depressed cubics if the discriminant is 0 then we havedisc(a, b, c, d) -> float:
"""

def disc(a, b, c, d) -> float:
    # Δ = 18abcd – 4b³d + b²c² – 4ac³ – 27a²d².
    return 1

# returns t1 (simple), and t2 (=t3, do)
def zero_disc(q, p) -> tuple[float, float]: # both of these roots should be real8
    return 3*q/p, -3*q/(2*p)


def depressed_disc(c,d) -> float:
    return disc(1, 0, c, d)

""" root cases
Δ > 0, the cubic has three distinct real roots
Δ < 0, the cubic has one real root and two non-real complex conjugate roots.
Δ = 0 the cubic has a multiple root (one simple and one double root)
"""

def main():
    print(depressed_disc(1,2))
    print(root(1,2))




if __name__ == "__main__":
    main()
