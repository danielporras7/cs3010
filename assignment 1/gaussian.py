

from typing import List

def read_input(filename):
   
    with open(filename, 'r') as f:
        n = int(f.readline().strip())
        coeff = [[float(x) for x in f.readline().strip().split()] for _ in range(n)]
        const = [float(x) for x in f.readline().strip().split()]
    return n, coeff, const

def write_output(filename, solution):
    
    with open(filename, 'w') as f:
        f.write(' '.join(str(x) for x in solution))

'''

Naive Gaussian Elimination

'''

def FwdElimination(coeff: List[List[float]], const: List[float]):
    n = len(const)
    for k in range(0, n-1):
        for i in range(k+1, n):
            mult = coeff[i][k] / coeff[k][k]
            for j in range(k, n):
                coeff[i][j] -= mult*coeff[k][j]
            const[i] -= mult*const[k]  
    return coeff, const

def BackSubst(coeff: List[List[float]], const: List[float]):
    n = len(const)
    sol = [0 for _ in range(n)] 
    sol[n-1] = const[n-1]/coeff[n-1][n-1]
    for i in reversed(range(0, n-1)):
        sum = const[i]  
        for j in range(i+1, n):
            sum -= coeff[i][j]*sol[j]
        sol[i] = sum/coeff[i][i]     
    return sol
    
def NaiveGaussian(coeff: List[List[float]], const: List[float]):
    coeff, const = FwdElimination(coeff, const)
    # solutions = BackSubst(coeff, const)
    # print("Solutions", BackSubst(coeff, const))
    return BackSubst(coeff, const)

'''

Gaussian Elimination with Scaled Partial Pivoting

'''

# Forward Elimination
def SPPFwdElimination(coeff: List[List[float]], const: List[float], ind: List[int]):
    n = len(const)
    scaling = [0]*n # vector of scaling factors
    
    for i in range(0, n):
        smax = 0
        for j in range(0, n):
            smax = max(smax, abs(coeff[i][j]))
            # print(i, j)
        scaling[i] = smax
    
    for k in range(0, n-1):
        rmax = 0
        maxind = k
        for i in range(k, n):
            r = abs(coeff[ind[i]][k] / scaling[ind[i]]) # ratio of coefficient to scaling factor
            if (r > rmax):
                rmax = r
                maxind = i
        ind[maxind], ind[k] = ind[k], ind[maxind]
        for i in range(k+1, n):
            mult = coeff[ind[i]][k] / coeff[ind[k]][k]
            for j in range(k+1, n):
              coeff[ind[i]][j] -= mult*coeff[ind[k]][j]
            const[ind[i]] -= mult*const[ind[k]]
    return coeff, const

# Back Substitution
def SPPBackSubst(coeff: List[List[float]], const: List[float], ind: List[int]):
    n = len(const)
    sol = [0 for _ in range(n)] 
    sol[n-1] = const[ind[n-1]] / coeff[ind[n-1]][n-1]
    for i in reversed(range(0, n-1)):
        sum = const[ind[i]]
        for j in range(i+1, n):
            sum -= coeff[ind[i]][j]*sol[j]
        sol[i] = sum / coeff[ind[i]][i]
    return sol

# SPP Gaussian Algorithm
def SPPGaussian(coeff: List[List[float]], const: List[float]):
    n = len(const)
    ind = []
    for i in range(0, n):
        ind.append(i)
    coeff, const = SPPFwdElimination(coeff, const, ind)
    sol = SPPBackSubst(coeff, const, ind)
    return sol

import sys
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print('Usage: gaussian [--spp] input_file')
        sys.exit(1)
    using_spp = False
    if len(sys.argv) > 2 and sys.argv[1] == '--spp':
        use_spp = True
        filename = sys.argv[2]
        print("Solving Linear System using Gaussian Elimination with Scaled Partial Pivoting\n")
    else:
        filename = sys.argv[1]
        print("Solving Linear System using Naive Gaussian Elimination\n")
    n, coefficients, constants = read_input(filename)
    solutions_x = SPPGaussian(coefficients, constants) if using_spp else NaiveGaussian(coefficients, constants)
    outfilename = filename.replace('.lin', '.sol')
    write_output(outfilename, solutions_x)
    print(f'Solution written to {outfilename}')

    '''
    The results of the Naive Gaussian Elimination method and the Gaussian 
Elimination with the Scaled Partial Pivoting method are identical. Although it is the same in this instance, 
I have tested in other systems and these two methods can yield results that differ by a small error. 
This variation may be caused by numerical errors that can happen during computations. 
Each element is divided by the pivot element in naive Gaussian elimination, which can result in accuracy loss and solution mistakes. 
By choosing the pivot element as the largest element in the column divided by the scaling factor for that row, 
Scaled Partial Pivoting attempts to reduce these inaccuracies. 
However, there are some cases showing that SPP can yield a result having a slight numerical error by a very small fraction,
 while the Naive Gaussian method gave a more accurate answer.

    '''