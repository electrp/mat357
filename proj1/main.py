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


# returns discriminant of a depressed cubic
def depressed_disc(p, q) -> float:
    """t**3 + pt + q"""
    return -(4*p**3 + 27*q**2)


def close(v_, t_, tol_=1e-12) -> float:
    return abs(v_-t_) < tol_


def sign(v_) -> int:
    if type(v_) == type((-1)**0.5):
        v_ = v_.real
    if close(v_, 0): return 0
    elif v_ > 0: return 1
    return -1


# https://en.wikipedia.org/wiki/Cubic_equation
def cardano(p_, q_) -> list[float]:
    """
    one real root (u1)**(1/3) + (u2)**(1/3)
    you have to cube root them but keep their original sign for the
    summation to get the right answer. I am not entirely sure why
    that is the case but it is passing our tests.
    """
    nqot = -q_/2
    big_term = (q_**2 /4 + p_**3 /27)**(1/2)
    u, v = nqot + big_term, nqot - big_term
    # do the cube root but get rid of the imaginary component from the sqrt if there was one
    cbrt = lambda x: sign(x) * abs(x) ** (1 / 3) if x != 0 else 0.0
    return [cbrt(u) + cbrt(v)]


import math
def trig_roots(p_, q_) -> list[float]:
    """
    2*sqrt(-p/3) * cos( 1/3*arccos(3q/2p * sqrt(-3/p)) - k2pi/3 )
    """
    l = 2*(-p_/3)**0.5

    # magic. (clipping & some algebra but huh?)
    m = 2 * (-p / 3) ** 0.5
    arg = max(-1.0, min(1.0, 3 * q / (p * m)))

    a = math.acos(arg)/3
    k = 2*math.pi/3     # k will be 0, 1 ,2 and multiplied in loop in a moment
    return [l * math.cos(a - i*k) for i in range(3)]


# returns t1 (simple), and t2 (=t3, double root)
def zero_disc_roots(p, q) -> tuple[float, float]: # both of these roots should be real8
    if p == q and p == 0: return 0
    return 3*q/p, -3*q/(2*p)


def depressed_roots(p_, q_):
    """ root cases
    Δ > 0, the cubic has three distinct real roots
    Δ < 0, the cubic has one real root and two non-real complex conjugate roots.
    Δ = 0 the cubic has a multiple root (one simple and one double root, *or* if p = q = 0, root is at 0)
    """
    if close(p_, 0) and close(q_,0):    return [0]
    elif depressed_disc(p_,q_) > 1e-10: return trig_roots(p_, q_)
    else:                               return cardano(p_, q_)




def main():
    P, Q = -5,1
    print(depressed_roots(P,Q))
    print(vibe_cubic_roots(P,Q))





if __name__ == "__main__":
    main()
