import sys
sys.path.append("build")
import clique_partitioning

n = 4
costs = [1, 1, 1, -1, -1, -1]
obj, bound, values = clique_partitioning.bnc(n, costs, activate_branching=False,
                                             add_some_triangles=False,
                                             separators=[("Triangle", 0, 0), ("OddWheel", 1, 0)])

print(obj, bound, values)

"""
Executing this script produces the following output

 Iter EXPND OPNND DEPTH     TIME  LP-TIME    OBJBST    OBJBND    NODBND   %I #Constr             Triangle             OddWheel 
    0     0     0     0        0        0         0         3         3  100       0    0       0       0    0       0       0 
    1     0     0     0        0        0         0       1.5       1.5   50       3    1       3       0    0       0       0 
    2     0     0     0        0        0         0         1         1  100       4    2       3       0    1       1       0 
    3     1     0     0        0        0         1         1         1  100       4    3       3       0    2       1       0 
1.0 1.0 [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]

In iteration 0, no constraints are added to the LP and the bound is 3 while all variables remain integer.
In iteration 1, three Triangle inequalities are added to the LP and the bound becomes 1.5 and only 50% of the variables remain integer.
In iteration 2, one odd wheel inequality is added to the model and the bound becomes 1.0 and all variables are integer.
In iteration 3, no additional inequality is added to the model and the algorithm terminates.
"""