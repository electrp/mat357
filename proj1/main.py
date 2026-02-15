from typing import Callable
import math
import os
file_path = os.path.dirname(os.path.realpath(__file__))

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


def cubic_roots(xs, ys, tol_: float=1e-8) -> list[float]:
    cubic = make_lagrange_cubic(*(xs+ys))
    a, b, c, d = cubic
    if abs(a) < tol_:
        if abs(b) < tol_:
            return [-d/c]
        return solve_quadratic_formula(b, c, d)
    else:
        p, q = depress_cubic(a, b, c, d)
        shift = -b / (3 * a)
        return [t + shift for t in depressed_cubic_roots(p, q)]


def newton_rootfind(f_: Callable, a_: float, b_: float, max_iter_: int = 100, tol_: float = 1e-12, debug: bool = False) -> tuple[float, int]:
    """
    Standard Newton's method with numerical derivative, bracketed.
    Returns (root, iterations_used).
    """
    h = 1e-8
    # Start from midpoint
    x = (a_ + b_) / 2
    
    for i in range(max_iter_):
        fx = f_(x)
        if abs(fx) < tol_:
            if debug:
                print(f"  Newton converged at iteration {i}")
            return x, i + 1
        
        # Numerical derivative
        dfx = (f_(x + h) - f_(x - h)) / (2 * h)
        
        if abs(dfx) < 1e-15:
            # Derivative too small, nudge
            x = x + h * 100
            continue
        
        x_new = x - fx / dfx
        
        # Clamp to bracket
        x_new = max(a_, min(b_, x_new))
        
        # If clamped to boundary, bisect instead
        if x_new == a_ or x_new == b_:
            x_new = (a_ + b_) / 2
            # Tighten bracket
            if sign(f_(a_)) != sign(fx):
                b_ = x
            else:
                a_ = x
        
        x = x_new
    
    if debug:
        print(f"  Newton did not converge in {max_iter_} iterations, best x={x}, f(x)={f_(x):.2e}")
    return x, max_iter_


def rootfind(f_: Callable, a_: float, b_: float, max_iter_: int = 100, tol_: float = 1e-12, debug: bool = True) -> float:
    """
    Cubic-interpolation rootfinder. Also runs Newton's method as comparison.
    """
    import numpy as np

    if debug:
        import matplotlib.pyplot as plt

    # ---- Run Newton as comparison ----
    newton_root, newton_iters = newton_rootfind(f_, a_, b_, max_iter_=max_iter_, tol_=tol_, debug=debug)
    newton_val = f_(newton_root)
    if debug:
        print(f"=== NEWTON COMPARISON ===")
        print(f"  Root: {newton_root:.12f}")
        print(f"  f(root): {newton_val:.2e}")
        print(f"  Iterations: {newton_iters}")
        print(f"=========================\n")
    # ---- Our cubic method ----
    xs = [a_, a_ + (b_ - a_) / 3, a_ + 2 * (b_ - a_) / 3, b_]
    ys = [f_(x) for x in xs]

    # Validate bracket
    if sign(ys[0]) == sign(ys[-1]):
        print("Warning: no sign change detected in initial bracket")

    best = 0
    best_y = 0

    method_counts = {"cubic": 0, "newton": 0, "newton_clamped": 0, "bisection": 0}
    total_iters = 0

    for i in range(max_iter_):
        # Sort points by x to keep them ordered
        combined = sorted(zip(xs, ys), key=lambda p: p[0])
        xs = [p[0] for p in combined]
        ys = [p[1] for p in combined]

        cubic = make_lagrange_cubic(*(xs + ys))
        a_c, b_c, c_c, d_c = cubic
        roots = cubic_roots(xs, ys)

        # Filter to roots inside the current interval
        filtered = [root for root in roots if xs[0] <= root <= xs[3]]

        used_method = "cubic"
        if filtered:
            evals = [(root, f_(root)) for root in filtered]
            best, best_y = min(evals, key=lambda pair: abs(pair[1]))
        else:
            # Cubic root out of bounds — use Newton step from cubic's derivative
            best_sample_idx = min(range(4), key=lambda j: abs(ys[j]))
            x_near = xs[best_sample_idx]
            y_near = ys[best_sample_idx]

            deriv = 3 * a_c * x_near**2 + 2 * b_c * x_near + c_c

            if abs(deriv) > 1e-15:
                newton_step = x_near - y_near / deriv
                if xs[0] <= newton_step <= xs[3]:
                    best = newton_step
                    best_y = f_(best)
                    used_method = "newton"
                else:
                    best = max(xs[0], min(xs[3], newton_step))
                    if best == xs[0] or best == xs[3]:
                        for j in range(3):
                            if sign(ys[j]) != sign(ys[j + 1]):
                                best = (xs[j] + xs[j + 1]) / 2
                                break
                    best_y = f_(best)
                    used_method = "newton_clamped"
            else:
                best = (xs[0] + xs[3]) / 2
                for j in range(3):
                    if sign(ys[j]) != sign(ys[j + 1]):
                        best = (xs[j] + xs[j + 1]) / 2
                        break
                best_y = f_(best)
                used_method = "bisection"

        method_counts[used_method] += 1
        total_iters = i + 1

        if debug:
            newton_pct = 100 * (method_counts["newton"] + method_counts["newton_clamped"]) / total_iters
            print(f"--- Iter {i} ---")
            print(f"  Method: {used_method.upper()}")
            print(f"  Window: [{xs[0]:.10f}, {xs[3]:.10f}]  size={xs[3] - xs[0]:.2e}")
            print(f"  Current points:")
            for j, (xi, yi) in enumerate(zip(xs, ys)):
                print(f"    x{j}={xi:.10f}  y{j}={yi:.2e}")
            print(f"  All roots: {[f'{r:.10f}' for r in roots]}")
            print(f"  In-bounds: {[f'{r:.10f}' for r in filtered]}")
            print(f"  Chosen: x={best:.10f}  f(x)={best_y:.2e}")
            print(f"  Method stats so far: {dict(method_counts)}  Newton%: {newton_pct:.1f}%")
            print()

            fig, axes = plt.subplots(1, 2, figsize=(14, 5))
            ax1, ax2 = axes

            margin = (xs[3] - xs[0]) * 0.15
            plot_xs = np.linspace(xs[0] - margin, xs[3] + margin, 300)
            plot_ys_actual = [f_(x) for x in plot_xs]
            plot_ys_cubic = [a_c * x**3 + b_c * x**2 + c_c * x + d_c for x in plot_xs]

            ax1.plot(plot_xs, plot_ys_actual, 'b-', label='f(x)', linewidth=2)
            ax1.plot(plot_xs, plot_ys_cubic, 'r--', label='cubic fit', linewidth=1.5)
            ax1.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
            ax1.scatter(xs, ys, color='green', zorder=5, s=80, label='sample points')

            colors = {'cubic': 'purple', 'newton': 'orange', 'newton_clamped': 'darkorange', 'bisection': 'red'}
            ax1.scatter([best], [best_y], color=colors.get(used_method, 'red'), zorder=6, s=120,
                       marker='*', label=used_method)

            ax1.set_title(f'Iter {i} — {used_method.upper()}')
            ax1.legend(fontsize=8)
            ax1.set_xlabel('x')
            ax1.set_ylabel('y')

            newton_pct_str = f"{newton_pct:.1f}%"
            info_text = (
                f"Iteration: {i}\n"
                f"Window: [{xs[0]:.8f}, {xs[3]:.8f}]\n"
                f"Window size: {xs[3] - xs[0]:.2e}\n"
                f"Cubic coeffs: a={a_c:.4e}, b={b_c:.4e}, c={c_c:.4e}, d={d_c:.4e}\n"
                f"\nAll roots: {[f'{r:.8f}' for r in roots]}\n"
                f"In-bounds roots: {[f'{r:.8f}' for r in filtered]}\n"
                f"\nChosen x: {best:.10f}\n"
                f"f(x): {best_y:.2e}\n"
                f"Method: {used_method}\n"
                f"\nMethod breakdown:\n"
                f"  cubic: {method_counts['cubic']}\n"
                f"  newton: {method_counts['newton']}\n"
                f"  newton_clamped: {method_counts['newton_clamped']}\n"
                f"  bisection: {method_counts['bisection']}\n"
                f"  Newton%: {newton_pct_str}\n"
                f"\nNewton comparison:\n"
                f"  Newton root: {newton_root:.10f}\n"
                f"  Newton iters: {newton_iters}\n"
                f"\nSample points:\n"
            )
            for j, (xi, yi) in enumerate(zip(xs, ys)):
                info_text += f"  x{j}={xi:.8f}  y{j}={yi:.2e}\n"

            ax2.text(0.05, 0.95, info_text, transform=ax2.transAxes,
                    fontsize=8, verticalalignment='top', fontfamily='monospace')
            ax2.axis('off')
            ax2.set_title('Debug Info')

            plt.tight_layout()
            plt.savefig(os.path.join(file_path, f"iter_{i:02d}.png"))
            plt.close()

        if abs(best_y) < tol_:
            if debug:
                print(f"Converged at iteration {i}")
            break

        xs.append(best)
        ys.append(best_y)
        combined = sorted(zip(xs, ys), key=lambda p: p[0])
        xs = [p[0] for p in combined]
        ys = [p[1] for p in combined]

        sign_change_idx = None
        smallest_gap = float('inf')
        for j in range(len(xs) - 1):
            if sign(ys[j]) != sign(ys[j + 1]):
                gap = xs[j + 1] - xs[j]
                if gap < smallest_gap:
                    smallest_gap = gap
                    sign_change_idx = j

        if sign_change_idx is not None:
            keep = {sign_change_idx, sign_change_idx + 1}
            remaining = [j for j in range(5) if j not in keep]
            remaining.sort(key=lambda j: abs(ys[j]))
            keep.update(remaining[:2])
            keep = sorted(keep)
            xs = [xs[j] for j in keep]
            ys = [ys[j] for j in keep]
        else:
            indexed = sorted(range(5), key=lambda j: abs(ys[j]))
            keep = sorted(indexed[:4])
            xs = [xs[j] for j in keep]
            ys = [ys[j] for j in keep]

        # Redistribute straggler: if one endpoint is disproportionately far
        # from the tight cluster near the root, pull it in closer
        for j in range(len(xs) - 1):
            if sign(ys[j]) != sign(ys[j + 1]):
                sc_gap = xs[j + 1] - xs[j]
                break
        else:
            sc_gap = xs[3] - xs[0]

        redistrib_ratio = 10
        dist_left = xs[1] - xs[0]
        dist_right = xs[3] - xs[2]

        if dist_right > redistrib_ratio * sc_gap and sc_gap > 0:
            # Right endpoint is a straggler — pull it in
            xs[3] = xs[2] + 2 * sc_gap
            ys[3] = f_(xs[3])
            if debug:
                print(f"  >> Redistributed right endpoint to {xs[3]:.10f}")
        elif dist_left > redistrib_ratio * sc_gap and sc_gap > 0:
            # Left endpoint is a straggler — pull it in
            xs[0] = xs[1] - 2 * sc_gap
            ys[0] = f_(xs[0])
            if debug:
                print(f"  >> Redistributed left endpoint to {xs[0]:.10f}")
    else:
        print(f"Stopped at max iterations")
        if debug:
            print(f"Final best: x={best:.10f}, f(x)={best_y:.2e}")

    # Final summary
    if debug and total_iters > 0:
        newton_pct = 100 * (method_counts["newton"] + method_counts["newton_clamped"]) / total_iters
        print(f"\n{'='*50}")
        print(f"FINAL COMPARISON")
        print(f"{'='*50}")
        print(f"Our method:")
        print(f"  Root:       {best:.12f}")
        print(f"  f(root):    {best_y:.2e}")
        print(f"  Iterations: {total_iters}")
        print(f"  Method breakdown: {dict(method_counts)}")
        print(f"  Newton fallback%: {newton_pct:.1f}%")
        print(f"Newton's method:")
        print(f"  Root:       {newton_root:.12f}")
        print(f"  f(root):    {newton_val:.2e}")
        print(f"  Iterations: {newton_iters}")
        print(f"{'='*50}")

    print(f"Stopped at max iterations")
    if debug:
        print(f"Final best: x={best:.10f}, f(x)={best_y:.2e}")
    return best


def main():
    # import sys
    # P, Q = [float(v) for v in sys.argv[1:]]
    # print(depressed_cubic_roots(P,Q))
    def a(x):
        return 19 * math.sin(x + 2) + x/2
    l, h = -1, 3
    root = rootfind(a, l, h)
    root_v = a(root);
    print(f"{root}, {root_v}")


if __name__ == "__main__":
    main()
