# A cutting plane algorithm for the clique partitioning problem

## Requirements
[Gurobi](https://www.gurobi.com/documentation/current/refman/cpp_api_overview.html)

## Build

```bash
mkdir build
cd build
cmake ..
cmake --build ./
```

## Usage
```bash
./main path/to/instance.txt
```

## Example

Executing 
```bash
./main ../data/example-from-paper.txt
```
produces the following output

```
 Iter EXPND OPNND DEPTH      TIME   LP-TIME    OBJBST    OBJBND    NODBND   %I #Constr             Triangle             OddWheel 
    0     0     0     0 0.0019531 0.0019531        11        11        11  100      14    0       0       0    0       0       0 
    0     1     0     0 0.0019531 0.0019531        11        11        11  100      14    0       0       0    0       0       0 
upper bound = 11
best integer feasible objective = 11
```
