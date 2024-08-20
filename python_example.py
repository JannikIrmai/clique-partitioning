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
"""