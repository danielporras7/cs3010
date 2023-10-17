
import sys
import os.path

# Define the three root-finding functions: bisection, Newton's method, and secant method
def bisection(f, a, b, maxIter=10000, eps=1e-7):
    fa = f(a)
    fb = f(b)
    
    if fa * fb >= 0:
        print("Inadequate values for a and b.")
        return None, 0, 'fail'

    error = b - a

    for it in range(1, maxIter+1):
        error /= 2
        c = a + error
        fc = f(c)
        if abs(error) < eps or fc == 0:
            print(f"Algorithm has converged after {it} iterations!")
            return c, it, 'success'

        if fa * fc < 0:
            b = c
            fb = fc
        else:
            a = c
            fa = fc

    print("Max iterations reached without convergence...")
    return c, maxIter, 'fail'

def newton(f, derF, x, maxIter=10000, eps=1e-7, delta=1e-7):
    fx = f(x)

    for it in range(1, maxIter+1):
        fd = derF(x)

        if abs(fd) < delta:
            print("Small slope!")
            return x, it, 'fail'
        d = fx / fd
        x -= d
        fx = f(x)

        if abs(d) < eps:
            print(f"Algorithm has converged after {it} iterations!")
            return x, it, "success"

    print("Max iterations reached without convergence...")
    return x, maxIter, "fail"

def secant(f, a, b, maxIter=10000, eps=1e-7):
    fa = f(a)
    fb = f(b)
    if abs(fa) > abs(fb):
        a, b = b, a
        fa, fb = fb, fa

    for it in range(1, maxIter+1):
        if abs(fa) > abs(fb):
            a, b = b, a
            fa, fb = fb, fa

        d = (b - a) / (fb - fa)
        b = a
        fb = fa
        d = d * fa
        
        if abs(d) < eps:
            print(f"Algorithm has converged after {it} iterations!")
            return a, it, "success"

        a = a - d
        fa = f(a)

    print("Max iterations reached without convergence...")
    return a, maxIter, "fail"

def hybrid(f, derF, a, b, maxIter=10000, eps=1e-7):
    fa = f(a)
    fb = f(b)

    if fa * fb >= 0:
        print("Inadequate values for a and b.")
        return None, 0, 'fail'

    for it in range(1, maxIter+1):
        c = (a + b) / 2
        fc = f(c)

        if abs(fc) < eps:
            print(f"Algorithm has converged after {it} iterations!")
            return c, it, 'success'

        if fa * fc < 0:
            b = c
            fb = fc
        else:
            a = c
            fa = fc

        if it <= 2:
            continue

        x = newton(f, derF, c, maxIter=1, eps=eps, delta=eps)

        if x is None:
            continue

        if a < x < b:
            if abs(f(x)) < abs(fc):
                return x, it, 'success'

    print("Maximum number of iterations reached!")
    return c, maxIter, 'fail'

def read_polynomial(filename):
    with open(filename) as file:
        n = int(file.readline().strip())
        coeffs = list(map(float, file.readline().strip().split()))
    return coeffs

def write_solution(filename, root, iterations, outcome):
    with open(os.path.splitext(filename)[0] + '.sol', 'w') as file:
        file.write(f"{root} {iterations} {outcome}\n")

def evaluate_polynomial(coeffs, x):
    return sum(c * x ** i for i, c in enumerate(coeffs))

def evaluate_polynomial_derF(coeffs, x):
    result = 0
    for i, c in enumerate(coeffs):
        if i > 0:
            result+=i*c*(x**(i-1))
    return result

def main():
    args = sys.argv[1:]
    # default values
    polyfile = args[-1]
    method = 'bisection'
    maxIter = 10000
    
    # Parse arguments
    if '-newt' in args:
        method = 'newton'
        args.remove('-newt')
    elif '-sec' in args:
        method = 'secant'
        args.remove('-sec')
    elif '-hybrid' in args:
        method = 'hybrid'
        args.remove('-hybrid')
    if '-maxIter' in args:
        index = args.index('-maxIter')
        maxIter = int(args[index+1])
        args.pop(index+1)
        args.remove('-maxIter')
    if len(args) < 2 or len(args) > 4:
        print('Usage: polRoot [-newt, -sec, -hybrid] [-maxIter n] initP [initP2] polyFileName')
        return

    # Initial points
    initial_points = []
    try:
        if method in ['bisection', 'secant', 'hybrid']:
            # if the method is one of the above, we need two initial points
            # check if the we have enough arguments
            initial_points = [int(args[0]), int(args[1])]
        else:
            initial_points = [int(args[0])]
        
    except Exception:
        print('Usage: polRoot [-newt, -sec, -hybrid] [-maxIter n] initP [initP2] polyFileName')    
    
    coeffs = read_polynomial(polyfile)
    coeffs.reverse()
    root = None
    it = None
    result = None
    if method == 'newton':
        root, it, result = newton(
            lambda x: evaluate_polynomial(coeffs, x),
            lambda x: evaluate_polynomial_derF(coeffs, x),
            initial_points[0],
            maxIter=maxIter
        )
    elif method == 'secant':
        root, it, result = secant(
            lambda x: evaluate_polynomial(coeffs, x),
            initial_points[0],
            initial_points[1],
            maxIter=maxIter
        )
    elif method == 'hybrid':
        method = 'Hybrid'
        root, it, result = hybrid(
            lambda x: evaluate_polynomial(coeffs, x),
            # lambda x: sum(i * c * x ** (i - 1) for i, c in enumerate(coeffs) if i>0),
            lambda x: evaluate_polynomial_derF(coeffs, x),
            initial_points[0],
            initial_points[1],
            maxIter=maxIter
        )
    else:
        root, it, result = bisection(
            lambda x: evaluate_polynomial(coeffs, x),
            initial_points[0],
            initial_points[1],
            maxIter=maxIter
        )

    if root is not None:
        write_solution(polyfile, root, it, result)
        print(f"{method} method: root = {root:.8f}, iterations = {it}, outcome = {result}")
    else:
        write_solution(polyfile, root, it, result)
        print(f"{method} method failed to converge after {maxIter} iterations.")

if __name__ == '__main__':
    main()
