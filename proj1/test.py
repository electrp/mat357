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
                       "n_iters": [], "n_calls": [], "c_fails": []}
               for name, *_ in test_functions}

    for trial in range(NUM_TRIALS):
        for name, f, lo, hi in test_functions:
            jitter_lo = lo + random.uniform(-0.1, 0.1)
            jitter_hi = hi + random.uniform(-0.1, 0.1)
            if jitter_lo >= jitter_hi:
                jitter_lo, jitter_hi = lo, hi

            x0 = (jitter_lo + jitter_hi) / 2

            try:
                c_root, c_iters, c_calls, c_cubic, c_bisect, c_fail = rootfind(f, jitter_lo, jitter_hi)
            except Exception as e:
                print(f"  Cubic error on {name} trial {trial}: {e}")
                c_iters, c_calls, c_cubic, c_bisect, c_fail = 100, 100, 0, 100, True

            n_root, n_iters, n_calls = newton_rootfind(f, x0)

            results[name]["c_iters"].append(c_iters)
            results[name]["c_calls"].append(c_calls)
            results[name]["c_cubic"].append(c_cubic)
            results[name]["c_bisect"].append(c_bisect)
            results[name]["n_iters"].append(n_iters)
            results[name]["n_calls"].append(n_calls)
            results[name]["c_fails"].append(1 if c_fail else 0)

    avg = lambda lst: sum(lst) / len(lst) if lst else 0

    # Classify functions by majority step type
    cubic_dominant = []
    bisect_dominant = []
    for name, *_ in test_functions:
        r = results[name]
        if avg(r['c_cubic']) >= avg(r['c_bisect']):
            cubic_dominant.append(name)
        else:
            bisect_dominant.append(name)

    # Console table
    print(f"\n{'Function':<16} {'Cubic Iters':<13} {'Newton Iters':<14} {'Cubic Calls':<13} {'Newton Calls':<14} {'Cubic':<9} {'Bisect':<9} {'Fails':<7}")
    print("-" * 95)

    all_r = {"c_iters": [], "c_calls": [], "c_cubic": [], "c_bisect": [], "n_iters": [], "n_calls": [], "c_fails": []}
    cubic_r = {"c_iters": [], "c_calls": [], "c_cubic": [], "c_bisect": [], "n_iters": [], "n_calls": [], "c_fails": []}
    bisect_r = {"c_iters": [], "c_calls": [], "c_cubic": [], "c_bisect": [], "n_iters": [], "n_calls": [], "c_fails": []}

    for name, *_ in test_functions:
        r = results[name]
        fail_str = f"{sum(r['c_fails'])}/{NUM_TRIALS}"
        print(f"{name:<16} {avg(r['c_iters']):<13.1f} {avg(r['n_iters']):<14.1f} {avg(r['c_calls']):<13.1f} {avg(r['n_calls']):<14.1f} {avg(r['c_cubic']):<9.1f} {avg(r['c_bisect']):<9.1f} {fail_str:<7}")
        target = cubic_r if name in cubic_dominant else bisect_r
        for k in all_r:
            all_r[k].extend(r[k])
            target[k].extend(r[k])

    print("-" * 95)
    print(f"{'Overall Avg':<16} {avg(all_r['c_iters']):<13.1f} {avg(all_r['n_iters']):<14.1f} {avg(all_r['c_calls']):<13.1f} {avg(all_r['n_calls']):<14.1f} {avg(all_r['c_cubic']):<9.1f} {avg(all_r['c_bisect']):<9.1f} {sum(all_r['c_fails'])}/{len(all_r['c_fails'])}")
    print(f"{'Mainly Cubic Avg':<16} {avg(cubic_r['c_iters']):<13.1f} {avg(cubic_r['n_iters']):<14.1f} {avg(cubic_r['c_calls']):<13.1f} {avg(cubic_r['n_calls']):<14.1f} {avg(cubic_r['c_cubic']):<9.1f} {avg(cubic_r['c_bisect']):<9.1f} {sum(cubic_r['c_fails'])}/{len(cubic_r['c_fails'])}")
    print(f"{'Mainly Bisect Avg':<16} {avg(bisect_r['c_iters']):<13.1f} {avg(bisect_r['n_iters']):<14.1f} {avg(bisect_r['c_calls']):<13.1f} {avg(bisect_r['n_calls']):<14.1f} {avg(bisect_r['c_cubic']):<9.1f} {avg(bisect_r['c_bisect']):<9.1f} {sum(bisect_r['c_fails'])}/{len(bisect_r['c_fails'])}")
    print(f"\n({NUM_TRIALS} trials per function)\n")

    # LaTeX table
    print("\n% LaTeX table\n")
    print(r"\begin{table}[h]")
    print(r"\centering")
    print(r"\caption{Rootfinding comparison averaged over " + str(NUM_TRIALS) + r" trials with randomized bounds}")
    print(r"\begin{tabular}{l c c c c c c c}")
    print(r"\hline")
    print(r" & \multicolumn{2}{c}{Iterations} & \multicolumn{2}{c}{Function Calls} & \multicolumn{2}{c}{Step Type} & \\")
    print(r"\cline{2-3} \cline{4-5} \cline{6-7}")
    print(r"Function & Cubic & Newton & Cubic & Newton & Cubic & Bisect & Fails \\")
    print(r"\hline")
    for name, *_ in test_functions:
        r = results[name]
        fail_str = f"{sum(r['c_fails'])}/{NUM_TRIALS}"
        print(f"${name}$ & {avg(r['c_iters']):.1f} & {avg(r['n_iters']):.1f} & {avg(r['c_calls']):.1f} & {avg(r['n_calls']):.1f} & {avg(r['c_cubic']):.1f} & {avg(r['c_bisect']):.1f} & {fail_str} \\\\")
    print(r"\hline")
    fail_all = f"{sum(all_r['c_fails'])}/{len(all_r['c_fails'])}"
    fail_cubic = f"{sum(cubic_r['c_fails'])}/{len(cubic_r['c_fails'])}"
    fail_bisect = f"{sum(bisect_r['c_fails'])}/{len(bisect_r['c_fails'])}"
    print(f"\\textbf{{Overall Avg}} & {avg(all_r['c_iters']):.1f} & {avg(all_r['n_iters']):.1f} & {avg(all_r['c_calls']):.1f} & {avg(all_r['n_calls']):.1f} & {avg(all_r['c_cubic']):.1f} & {avg(all_r['c_bisect']):.1f} & {fail_all} \\\\")
    print(f"\\textbf{{Mainly Cubic Avg}} & {avg(cubic_r['c_iters']):.1f} & {avg(cubic_r['n_iters']):.1f} & {avg(cubic_r['c_calls']):.1f} & {avg(cubic_r['n_calls']):.1f} & {avg(cubic_r['c_cubic']):.1f} & {avg(cubic_r['c_bisect']):.1f} & {fail_cubic} \\\\")
    print(f"\\textbf{{Mainly Bisect Avg}} & {avg(bisect_r['c_iters']):.1f} & {avg(bisect_r['n_iters']):.1f} & {avg(bisect_r['c_calls']):.1f} & {avg(bisect_r['n_calls']):.1f} & {avg(bisect_r['c_cubic']):.1f} & {avg(bisect_r['c_bisect']):.1f} & {fail_bisect} \\\\")
    print(r"\hline")
    print(r"\end{tabular}")
    print(r"\end{table}")

if __name__ == "__main__":
    run_comparison()