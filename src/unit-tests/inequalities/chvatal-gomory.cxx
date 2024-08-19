#include <property-maps.hxx>
#include <inequalities/chvatal-gomory.hxx>
#include <inequalities/half-chorded-odd-cycle.hxx>

#include <gurobi_c++.h>

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


void test_wheel()
{
    typedef CP::EdgePropertyMap<std::vector<double>> EPM;
    
    GRBEnv env;
    std::vector<std::vector<GRBVar>> vars;

    for (size_t k = 3; k < 15; ++k)
    {
        size_t n = k+1;
        size_t m = n * (n - 1) / 2;
        std::vector<double> edge_values(m, 0.5);
        EPM edge_value_map(n, edge_values);
        for (size_t i = 0; i < k; ++i)
        {
            edge_value_map(0, i + 1) = 0.5;
            edge_value_map(i + 1, (i+1)%k + 1) = 0.0;
        }

        CP::ChvatalGomorySeparator<EPM> separator(n);
        std::vector<CP::Inequality<int>> inequalities = separator.separate_(edge_value_map);

        if (k % 2 == 1)
        {
            test(inequalities.size() >= 1); // there may be other violated inequalities beside odd wheel
            test(inequalities[0].violation() == 0.5);
            test(inequalities[0].violation() == inequalities[0].evaluate(edge_value_map) - inequalities[0].rhs());
            test(inequalities[0].edges().size() == 2 * k);
        }
        else 
            test(inequalities.size() == 0);

        // test that all found inequalities are integer feasible
        GRBModel model = get_model(env, n, vars, true);
        for (const auto& inequality : inequalities)
        {
            // set costs to find maximum violated
            GRBLinExpr obj_func;
            for (size_t i = 0; i < inequality.edges().size(); ++i)
            {
                obj_func += inequality.coefficients()[i] * vars[inequality.edges()[i][0]][inequality.edges()[i][1]];
            }
            model.setObjective(obj_func);
            model.optimize();
            test(model.get(GRB_DoubleAttr_ObjVal) == inequality.rhs());
        }   
    }
}


void test_two_chorded_cycle()
{
    typedef CP::EdgePropertyMap<std::vector<double>> EPM;
    
    GRBEnv env;
    std::vector<std::vector<GRBVar>> vars;

    for (size_t n = 5; n < 15; ++n)
    {
        size_t m = n * (n - 1) / 2;
        std::vector<double> edge_values(m, 0.5);
        EPM edge_value_map(n, edge_values);
        for (size_t i = 0; i < n; ++i)
        {
            edge_value_map(i, (i+1)%n) = 0.0;
        }

        CP::ChvatalGomorySeparator<EPM> separator(n);
        std::vector<CP::Inequality<int>> inequalities = separator.separate_(edge_value_map);

        if (n % 2 == 1)
        {
            test(inequalities.size() >= 1); // there may be other violated inequalities beside two chorded odd cycle
            test(inequalities[0].violation() == 0.5);
            test(inequalities[0].violation() == inequalities[0].evaluate(edge_value_map) - inequalities[0].rhs());
            test(inequalities[0].edges().size() == 2 * n);
        }
        else 
            test(inequalities.size() == 0);

        // test that all found inequalities are valid
        GRBModel model = get_model(env, n, vars, true);
        for (const auto& inequality : inequalities)
        {
            // set costs to find maximum violated
            GRBLinExpr obj_func;
            for (size_t i = 0; i < inequality.edges().size(); ++i)
            {
                obj_func += inequality.coefficients()[i] * vars[inequality.edges()[i][0]][inequality.edges()[i][1]];
            }
            model.setObjective(obj_func);
            model.optimize();
            test(model.get(GRB_DoubleAttr_ObjVal) == inequality.rhs());
        }   
    }
}



void test_half_chorded_cycle()
{
    typedef CP::EdgePropertyMap<std::vector<double>> EPM;
    
    GRBEnv env;
    std::vector<std::vector<GRBVar>> vars;

    for (size_t n = 5; n < 15; ++n)
    {
        size_t m = n * (n - 1) / 2;
        size_t k = n / 2;
        std::vector<double> edge_values(m);
        EPM edge_value_map(n, edge_values);
        for (size_t i = 0; i < n; ++i)
        for (size_t d = 1; d<=k; ++d)
        {
            edge_value_map(i, (i+d)%n) = (double)(k - d) / k;
        }

        CP::ChvatalGomorySeparator<EPM> separator(n);
        std::vector<CP::Inequality<int>> inequalities = separator.separate_(edge_value_map);
        
        CP::HalfChordedOddCycleSeparator<EPM> hc_separator(n);
        std::vector<CP::Inequality<int>> inequalities_hc = hc_separator.separate_(edge_value_map);

        // assert that no half chorded odd cycle inequality is violated (which would be violated by (n-3)/(n-1))
        for (const auto& inequality: inequalities)
        {
            test(inequality.violation() <= 0.5);
        }

        if (n % 2 == 1)
        {
            test(inequalities_hc.size() == 1);
            test(std::abs(inequalities_hc[0].violation() - (double)(n-3)/(n-1)) < 1e-6);
        }

        // test that all found inequalities are valid
        GRBModel model = get_model(env, n, vars, true);
        for (const auto& inequality : inequalities)
        {
            // set costs to find maximum violated
            GRBLinExpr obj_func;
            for (size_t i = 0; i < inequality.edges().size(); ++i)
            {
                obj_func += inequality.coefficients()[i] * vars[inequality.edges()[i][0]][inequality.edges()[i][1]];
            }
            model.setObjective(obj_func);
            model.optimize();
            test(model.get(GRB_DoubleAttr_ObjVal) == inequality.rhs());
        }   
    }
}




int main()
{
    test_wheel();
    test_two_chorded_cycle();
    test_half_chorded_cycle();

    return 0;
}
