from typing import Callable
import math

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


# Gives c and d
def depress_cubic(a, b, c, d) -> tuple[float, float]:
    # Prevents div by 0 cases in linear equations
    if a == 0:
        a = .000000001 # 
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



def cubic_roots(a, b, c, d, tol_: float=1e-8) -> list[float]:
    if abs(a) < tol_:
        if abs(b) < tol_:
            return [-d/c]
        return solve_quadratic_formula(b, c, d)
    else:
        p, q = depress_cubic(a, b, c, d)
        shift = -b / (3 * a)
        return [t + shift for t in depressed_cubic_roots(p, q)]

def sample_polynomial(x, a, b, c, d):
    return a*x**3 + b*x**2 + c*x + d


def rootfind(f_:Callable, a_:float, b_:float, max_iter_:int=100, tol_:float=1e-12) -> tuple[float, int, int, int, int, bool]:
    calls = 0
    def cf(x):
        nonlocal calls
        calls += 1
        return f_(x)

    cubic_steps = 0
    bisect_steps = 0
    failed = False

    # get four points
    xs = [a_, a_ + (b_ - a_) / 3, a_ + 2 * (b_ - a_) / 3, b_]
    ys = [cf(x) for x in xs]
    best = 0
    best_y = 0
    last_best = None

    for i in range(max_iter_):
        # Sort xs and ys so it goes furthest on the negative side to furthest on the positive side
        xs, ys = zip(*sorted(zip(xs, ys), key=lambda p: p[0]))
        xs, ys = list(xs), list(ys)

        a, b, c, d = make_lagrange_cubic(*(xs + ys))
        roots = cubic_roots(a, b, c, d)

        # Filter the roots to meet the following:
        #   Root is within the range of our xs
        #   Root of polynomial is within a tolerance of being a root
        #   Root is not too close to the last best root (to prevent loops)
        filtered = [r for r in roots
                     if xs[0] <= r <= xs[-1]
                     and abs(sample_polynomial(r, a, b, c, d)) < tol_
                     and (last_best is None or abs(r - last_best) > tol_)]

        if len(filtered):
            best = filtered[0]
            cubic_steps += 1
        else:
            # If no candidate roots are found we bisect the closest positive and negative points
            pos = [(x, y) for x, y in zip(xs, ys) if y > 0]
            neg = [(x, y) for x, y in zip(xs, ys) if y < 0]
            if not pos or not neg:
                failed = True
                break
            # Get closest positive and negative points and bisect
            min_pos = min(pos, key=lambda p: abs(p[0]))
            min_neg = min(neg, key=lambda p: abs(p[0]))
            best = (min_pos[0] + min_neg[0]) / 2
            bisect_steps += 1

        best_y = cf(best)

        # Check tolerance
        if abs(best_y) < tol_:
            return best, i + 1, calls, cubic_steps, bisect_steps, False

        # Replace one of the points
        xs.insert(2, best)
        ys.insert(2, best_y)

        # Remove the nearest point on both sides of the root
        pos = [(x, y) for x, y in zip(xs, ys) if y > 0]
        neg = [(x, y) for x, y in zip(xs, ys) if y < 0]
        if not pos or not neg:
            print("All points on one side, this should not happen")
            failed = True
            break
        min_pos = min(pos, key=lambda p: abs(p[0]))
        min_neg = min(neg, key=lambda p: abs(p[0]))

        # get xs and ys without the closest pos and neg
        filtered_pairs = [(x, y) for x, y in zip(xs, ys) if (x, y) != min_pos and (x, y) != min_neg]
        # Get furthest point in filtered list
        furthest = max(filtered_pairs, key=lambda p: abs(p[0]))
        # Remove furthest point
        xs.remove(furthest[0])
        ys.remove(furthest[1])

        last_best = best

    else:
        # Loop completed without break or early return — hit max iterations
        failed = True

    return best, max_iter_, calls, cubic_steps, bisect_steps, failed


if __name__ == "__main__":
    def a(x):
        return math.cos(x) - x
    l, h = 0, 2
    root, iters, calls, c_steps, b_steps = rootfind(a, l, h)
    print(f"Root: {root}, f(root)={a(root)}, iters={iters}, calls={calls}, c_steps={c_steps}, b_steps={b_steps}")