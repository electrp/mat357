import math
import random
from main import rootfind

def newton_rootfind(f, x0, max_iter=100, tol=1e-12, h=1e-8):
    calls = 0
    def cf(x):
        nonlocal calls
        calls += 1
        return f(x)
    x = x0
    for i in range(max_iter):
        fx = cf(x)
        if abs(fx) < tol:
            return x, i, calls
        dfx = (cf(x + h) - cf(x - h)) / (2 * h)
        if abs(dfx) < 1e-15:
            return x, i, calls
        x = x - fx / dfx
    return x, max_iter, calls

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

NUM_TRIALS = 10

def run_comparison():
    random.seed(42)

    results = {name: {"c_iters": [], "c_calls": [], "c_cubic": [], "c_bisect": [],
                       "n_iters": [], "n_calls": []}
               for name, *_ in test_functions}

    for trial in range(NUM_TRIALS):
        for name, f, lo, hi in test_functions:
            jitter_lo = lo + random.uniform(-0.1, 0.1)
            jitter_hi = hi + random.uniform(-0.1, 0.1)
            if jitter_lo >= jitter_hi:
                jitter_lo, jitter_hi = lo, hi

            x0 = (jitter_lo + jitter_hi) / 2

            try:
                c_root, c_iters, c_calls, c_cubic, c_bisect = rootfind(f, jitter_lo, jitter_hi)
            except Exception as e:
                print(f"  Cubic error on {name} trial {trial}: {e}")
                c_iters, c_calls, c_cubic, c_bisect = 100, 100, 0, 100

            n_root, n_iters, n_calls = newton_rootfind(f, x0)

            results[name]["c_iters"].append(c_iters)
            results[name]["c_calls"].append(c_calls)
            results[name]["c_cubic"].append(c_cubic)
            results[name]["c_bisect"].append(c_bisect)
            results[name]["n_iters"].append(n_iters)
            results[name]["n_calls"].append(n_calls)

    avg = lambda lst: sum(lst) / len(lst)

    # Console table
    print(f"\n{'Function':<16} {'Cubic Iters':<13} {'Newton Iters':<14} {'Cubic Calls':<13} {'Newton Calls':<14} {'Cubic':<9} {'Bisect':<9}")
    print("-" * 88)

    all_r = {"c_iters": [], "c_calls": [], "c_cubic": [], "c_bisect": [], "n_iters": [], "n_calls": []}

    for name, *_ in test_functions:
        r = results[name]
        print(f"{name:<16} {avg(r['c_iters']):<13.1f} {avg(r['n_iters']):<14.1f} {avg(r['c_calls']):<13.1f} {avg(r['n_calls']):<14.1f} {avg(r['c_cubic']):<9.1f} {avg(r['c_bisect']):<9.1f}")
        for k in all_r:
            all_r[k].extend(r[k])

    print("-" * 88)
    print(f"{'Average':<16} {avg(all_r['c_iters']):<13.1f} {avg(all_r['n_iters']):<14.1f} {avg(all_r['c_calls']):<13.1f} {avg(all_r['n_calls']):<14.1f} {avg(all_r['c_cubic']):<9.1f} {avg(all_r['c_bisect']):<9.1f}")
    print(f"\n({NUM_TRIALS} trials per function)\n")

    # LaTeX table
    print("\n% LaTeX table\n")
    print(r"\begin{table}[h]")
    print(r"\centering")
    print(r"\caption{Rootfinding comparison averaged over " + str(NUM_TRIALS) + r" trials with randomized bounds}")
    print(r"\begin{tabular}{l c c c c c c}")
    print(r"\hline")
    print(r" & \multicolumn{2}{c}{Iterations} & \multicolumn{2}{c}{Function Calls} & \multicolumn{2}{c}{Step Type} \\")
    print(r"\cline{2-3} \cline{4-5} \cline{6-7}")
    print(r"Function & Cubic & Newton & Cubic & Newton & Cubic & Bisect \\")
    print(r"\hline")
    for name, *_ in test_functions:
        r = results[name]
        print(f"${name}$ & {avg(r['c_iters']):.1f} & {avg(r['n_iters']):.1f} & {avg(r['c_calls']):.1f} & {avg(r['n_calls']):.1f} & {avg(r['c_cubic']):.1f} & {avg(r['c_bisect']):.1f} \\\\")
    print(r"\hline")
    print(f"\\textbf{{Average}} & {avg(all_r['c_iters']):.1f} & {avg(all_r['n_iters']):.1f} & {avg(all_r['c_calls']):.1f} & {avg(all_r['n_calls']):.1f} & {avg(all_r['c_cubic']):.1f} & {avg(all_r['c_bisect']):.1f} \\\\")
    print(r"\hline")
    print(r"\end{tabular}")
    print(r"\end{table}")

if __name__ == "__main__":
    run_comparison()