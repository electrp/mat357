from typing import Callable
import math
import random

# a is always 1, returns b, c, d
def make_lagrange_cubic(x0, x1, x2, x3, y0, y1, y2, y3) -> tuple[float, float, float, float]:
    # Precompute differences
    xd10 = x1 - x0
    xd20 = x2 - x0
    xd30 = x3 - x0 
    xd21 = x2 - x1
    xd31 = x3 - x1
    xd32 = x3 - x2

    # Precompute divisors with y
    l0_m = y0 / (-xd10 * -xd20 * -xd30)
    l1_m = y1 / (xd10 * -xd21 * -xd31)
    l2_m = y2 / (xd20 *  xd21 * -xd32)
    l3_m = y3 / (xd30 *  xd31 *  xd32)

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
    a = l0_m + l1_m + l2_m + l3_m
    b = \
        (-xaccum + x0) * l0_m + \
        (-xaccum + x1) * l1_m + \
        (-xaccum + x2) * l2_m + \
        (-xaccum + x3) * l3_m 
    c = \
        (x12 + x23 + x13) * l0_m + \
        (x02 + x23 + x03) * l1_m + \
        (x01 + x13 + x03) * l2_m + \
        (x01 + x02 + x12) * l3_m
    d = \
        -x12 * x3 * l0_m + \
        -x02 * x3 * l1_m + \
        -x01 * x3 * l2_m + \
        -x01 * x2 * l3_m

    return a, b, c, d


def depress_cubic_const_a(b, c, d) -> tuple[float, float]:
    h1 = b / 3
    h2 = h1 * h1
    h3 = h1 * h2

    c1 = -h2 * 3 + c
    d1 = 2 * h3 - c * b / 3 + d
    return c1, d1 


# Gives c, d
def depress_cubic(a, b, c, d) -> tuple[float, float]:
    # Scale so a == 1, x = t - b/3
    # b1 = b / a
    # c1 = c / a
    # d1 = d / a
    # return depress_cubic_const_a(b1, c1, d1)

    h1 = b /(3 * a)
    h2 = h1 * h1
    h3 = h1 * h2

    c1 = c - b * b / (3 * a)
    d1 = 2 * h3 - c * b / (3 * a * a) + d / a
    d1 = d + 2 * b * b * b / (27 * a * a) - (b * c) / (3 * a)
    return c1 / a, d1 / a





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


def trig_roots(p_, q_) -> list[float]:
    """
    2*sqrt(-p/3) * cos( 1/3*arccos(3q/2p * sqrt(-3/p)) - k2pi/3 )
    """
    l = 2*(-p_/3)**0.5

    # magic. (clipping & some algebra but huh?)
    m = 2 * (-p_ / 3) ** 0.5
    arg = max(-1.0, min(1.0, 3 * q_ / (p_ * m)))

    a = math.acos(arg)/3
    k = 2*math.pi/3     # k will be 0, 1 ,2 and multiplied in loop in a moment
    return [l * math.cos(a - i*k) for i in range(3)]


# returns t1 (simple), and t2 (=t3, double root)
def zero_disc_roots(p, q) -> tuple[float, float]: # both of these roots should be real8
    if p == q and p == 0: return 0
    return 3*q/p, -3*q/(2*p)


def depressed_cubic_roots(p_, q_):
    """ root cases
    Δ > 0, the cubic has three distinct real roots
    Δ < 0, the cubic has one real root and two non-real complex conjugate roots.
    Δ = 0 the cubic has a multiple root (one simple and one double root, *or* if p = q = 0, root is at 0)
    """
    if close(p_, 0) and close(q_,0):    return [0]
    elif depressed_disc(p_,q_) > 1e-10: return trig_roots(p_, q_)
    else:                               return cardano(p_, q_)

def solve_quadratic_formula(a, b, c) -> list[float]:
    b1 = -b / (2 * a)
    inner = b * b - (4 * a * c)
    if inner <= 0:
        # Return our x minima if theres no root
        return [b1]
    variant = math.sqrt(inner) / (2 * a)
    return [b1 + variant, b1 - variant]


def cubic_roots(xs, ys, tol_: float=1e-8) -> list[float]:
    # create surrogate cubic (depressed cubic from lagrange polynomial)
    cubic = make_lagrange_cubic(*(xs+ys))
    print(f"cubic: {cubic}")
    # Check to make sure it's not of a smaller degree, if we don't it can explode
    if abs(cubic[0]) < tol_:
        if abs(cubic[1]) < tol_:
            # Solve linear
            return [-cubic[3]/cubic[2]]
        # Solve quadratic
        return solve_quadratic_formula(cubic[1], cubic[2], cubic[3])
    else:
        # Cubic if it matches best
        p, q = depress_cubic(*cubic)
        print(f"depressed: {p}, {q}")
        depressed_roots = depressed_cubic_roots(p, q)
        return [v - (cubic[1] / (3 * cubic[0])) for v in depressed_roots]


def distance_from_bounds(low, high, val):
    return min(0, val - high, low - val)

# Returns a new xs and ys 
# This was a very messy way to do this, but the value selection can change convergance
# greatly. It would take far too long to find the best way to do this
def adjust_bounds(xs, ys, roots, fn) -> tuple[list[float], list[float]]:
    # Filter into negative and positive
    y_pos = [v for v in ys if v >= 0]
    y_neg = [v for v in ys if v < 0]

    # Calculate our best bounds
    y_pos_bound = None
    y_neg_bound = None

    if len(y_pos) == 0:
        y_pos_bound = max(ys)
    else:
        y_pos_bound = min(y_pos)
    if len(y_neg) == 0:
        y_neg_bound = min(ys)
    else:
        y_neg_bound = max(y_neg)

    y_pos_bound_idx = ys.index(y_pos_bound)
    y_neg_bound_idx = ys.index(y_neg_bound)
    y_pos_bound_x = xs[y_pos_bound_idx]
    y_neg_bound_x = xs[y_neg_bound_idx]
    x_min_bound = min(y_pos_bound_x, y_neg_bound_x)
    x_max_bound = max(y_pos_bound_x, y_neg_bound_x)

    # Get our remaining indices
    remaining = [0, 1, 2, 3]
    remaining.remove(y_pos_bound_idx)
    remaining.remove(y_neg_bound_idx)

    r0_y = ys[remaining[0]]
    r1_y = ys[remaining[1]]

    # Pick our closest one
    r_selected = None
    if abs(r0_y) < abs(r1_y):
        r_selected = remaining[0]
    else:
        r_selected = remaining[1]

    # Select our root (dosn't matter too much right now)
    root = roots[0]

    return \
        [y_pos_bound_x, y_neg_bound_x, xs[r_selected], root], \
        [y_pos_bound, y_neg_bound, ys[r_selected], fn(root)]


# implement this
def rootfind(f_:Callable, a_:float, b_:float, max_iter_:int=100, tol_:float=1e-12) -> float:
    """
    given sample function and range [a,b] to search for roots, find a root
    """

    """
    ...
    """


    # get four points
    xs = [a_, a_+(b_-a_)/3, a_+2*(b_-a_)/3, b_]
    ys = [f_(x) for x in xs]

    best_x = None
    best_y = None

    for i in range(max_iter_):
        # Get roots
        roots = cubic_roots(xs, ys)
        print(f"roots: {roots}")

        # Check and see if any of our roots are good
        for root in roots:
            root_y = abs(f_(root))
            if best_y == None or root_y < best_y:
                best_y = root_y
                best_x = root
        if best_y < tol_:
            print(f"steps: {i}")
            return best_x

        print(f"{xs}")
        print(f"{ys}")
        # Reset bounds and continue
        xs, ys = adjust_bounds(xs, ys, roots, f_)
        print(f"{xs}")
        print(f"{ys}")
        print("____")

    print(f"Stopped at max iterations")
    return best_x



def main():
    # import sys
    # P, Q = [float(v) for v in sys.argv[1:]]
    # print(depressed_cubic_roots(P,Q))
    def a(x):
        return 19 * math.sin(x + 2) + x/2
    l, h = -3, 1
    root = rootfind(a, l, h)
    root_v = a(root);
    print(f"{root}, {root_v}")


if __name__ == "__main__":
    main()
