import sys
sys.path.append("build")
import clique_partitioning

n = 4
costs = [1, 1, 1, -1, -1, -1]
obj, bound, values = clique_partitioning.bnc(n, costs, activate_branching=False,
                                             add_some_triangles=False,
                                             separators=[("Triangle", 0, 0), ("OddWheel", 1, 0)])

print(obj, bound, values)
