#include "gurobi_c++.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>


void test(const bool& b)
{
    if (!b)
        throw std::logic_error("Test Failed!");
}


GRBModel get_model(GRBEnv env, size_t n, std::vector<std::vector<GRBVar>>& vars, bool binary = false)
{
    GRBModel model(env);
    model.set(GRB_IntParam_OutputFlag, 0);
    model.set(GRB_IntAttr_ModelSense, -1);  // maximize

    // add one variable for each edge {i, j} (and a copy for {j,i})
    vars = std::vector<std::vector<GRBVar>>(n, std::vector<GRBVar>(n));
    
    for (size_t i = 0; i < n; ++i)
    for (size_t j = i+1; j < n; ++j)
    {
        vars[i][j] = model.addVar(0.0, 1.0, 0, binary ? GRB_BINARY : GRB_CONTINUOUS, "x_"+std::to_string(i)+"_"+std::to_string(j));
        vars[j][i] = vars[i][j];
    }

    for (size_t i = 0; i < n; ++i) 
    for (size_t j = i+1; j < n; ++j)
    for (size_t k = j+1; k < n; ++k)
    {
        model.addConstr(  vars[i][j] + vars[i][k] - vars[j][k] <= 1);
        model.addConstr(  vars[i][j] - vars[i][k] + vars[j][k] <= 1);
        model.addConstr(- vars[i][j] + vars[i][k] + vars[j][k] <= 1);
    }

    return model;
}

void test_clique_web()
{
    GRBEnv env;
    std::vector<std::vector<GRBVar>> vars;

    std::cout << " n  p  q  r facet\n";

    for (size_t n = 3; n < 20; ++n)
    for (size_t r = 0; r <= (n-3)/2; ++r)
    for (size_t q = 1; q <= (n-2*r-1) / 2; ++q)
    {
        size_t p = n - q;
        test (p - q >= 2*r + 1);
        bool facet = true;
        if (r >= 1 && q < (p-1)/2 - r)
            facet = false;
        if (p - q == 2*r + 1 && q == 1)
            facet = false;

        GRBModel model = get_model(env, n, vars, false);
        GRBLinExpr obj_func;
        for (size_t i = 0; i < p; ++i)
        for (size_t d = r+1; d <= (p-1)/2; ++d)
        {
            obj_func -= vars[i][(i+d)%p];
        }
        if (p % 2 == 0)
        {
            for (size_t i = 0; i < p / 2; ++i)
            {
                obj_func -= vars[i][(i+p/2)%p];
            }
        }
        for (size_t i = 0; i < q; ++i)
        for (size_t j = i+1; j < q; ++j)
        {
            obj_func -= vars[p+i][p+j];
        }
        for (size_t i = 0; i < p; ++i)
        for (size_t j = 0; j < q; ++j)
        {
            obj_func += vars[i][p+j];
        }
        model.setObjective(obj_func);
        
        model.optimize();
        double obj = model.get(GRB_DoubleAttr_ObjVal);

        double violation = obj - q * (r+1);
        double euclidean = std::sqrt(p * q + p * (p-1-2*r) / 2 + q * (q-1) / 2);

        test (std::abs(2 * obj - p * q) < 1e-9);


        std::cout 
            << std::setw(2) << n << " " << std::setw(2) << p << " " << std::setw(2) << q << " " << std::setw(2) << r << " "
            << std::setw(5) << (facet ? "true" : "false") << " "
            << std::setw(5) << violation << " "
            << std::setw(10) << violation/euclidean << "\n";
    }

}

int main()
{
    try {   
        test_clique_web();
    } catch (GRBException e) {
        std::cout << "Error number: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
        throw e;
    } 
    return 0;

}