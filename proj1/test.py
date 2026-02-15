import math
import random
from main import rootfind

def newton_rootfind(f, x0, max_iter=100, tol=1e-12, h=1e-8):
    x = x0
    for i in range(max_iter):
        fx = f(x)
        if abs(fx) < tol:
            return x, i
        dfx = (f(x + h) - f(x - h)) / (2 * h)
        if abs(dfx) < 1e-15:
            return x, i
        x = x - fx / dfx
    return x, max_iter

def rootfind_counted(f, a, b, max_iter=100, tol=1e-12):
    count = [0]
    def counted_f(x):
        count[0] += 1
        return f(x)
    result = rootfind(counted_f, a, b, max_iter, tol)
    return result[0], count[0]

test_functions = [
    ("sin(x)+0.2",  lambda x: math.sin(x) + 0.2, -1, 1),
    ("x^3 - x - 2", lambda x: x**3 - x - 2,       1, 2),
    ("cos(x) - x",  lambda x: math.cos(x) - x,    0, 2),
    ("e^x - 3",     lambda x: math.exp(x) - 3,    0, 2),
    ("x^5 - x - 1", lambda x: x**5 - x - 1,       1, 2),
    ("ln(x) - 1",   lambda x: math.log(x) - 1,    1, 4),
    ("tan(x) - 2x", lambda x: math.tan(x) - 2*x,  0.1, 1.2),
    ("x^2 - 2",     lambda x: x**2 - 2,           1, 2),
]

def run_comparison():
    random.seed(42)

    print(f"{'Function':<20} {'Method':<10} {'Root':<22} {'f(root)':<15} {'F-calls':<8}")
    print("-" * 80)

    newton_totals = []
    cubic_totals = []

    for name, f, lo, hi in test_functions:
        jitter_lo = lo + random.uniform(-0.1, 0.1)
        jitter_hi = hi + random.uniform(-0.1, 0.1)
        if jitter_lo >= jitter_hi:
            jitter_lo, jitter_hi = lo, hi

        x0 = (jitter_lo + jitter_hi) / 2

        c_root, c_calls = rootfind_counted(f, jitter_lo, jitter_hi)
        n_root, n_iters = newton_rootfind(f, x0)
        n_calls = n_iters * 2 + 1

        print(f"{name:<20} {'cubic':<10} {c_root:<22.15g} {f(c_root):<15.3e} {c_calls:<8}")
        print(f"{'':<20} {'newton':<10} {n_root:<22.15g} {f(n_root):<15.3e} {n_calls:<8}")
        print()

        cubic_totals.append(c_calls)
        newton_totals.append(n_calls)

    print("-" * 80)
    print(f"{'Average f-calls':<20} {'cubic':<10} {'':<22} {'':<15} {sum(cubic_totals)/len(cubic_totals):<8.1f}")
    print(f"{'':<20} {'newton':<10} {'':<22} {'':<15} {sum(newton_totals)/len(newton_totals):<8.1f}")

if __name__ == "__main__":
    run_comparison()