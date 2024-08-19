#include <vector>
#include <property-maps.hxx>
#include <clique-partitioning-bnc.hxx>


typedef std::vector<int> VECTOR;
typedef CP::EdgePropertyMap<VECTOR> EPM;
typedef CP::BnC<EPM> BnC;


void test(const bool& b)
{
    if (!b)
        throw std::logic_error("Test Failed!");
}


void test_copy_gurobi_model()
{
    GRBEnv env;
    GRBModel model(env);
    GRBVar var0 = model.addVar(0, 1, -1, GRB_CONTINUOUS);
    GRBVar var1 = model.addVar(0, 1, -2, GRB_CONTINUOUS);
    model.update();
    GRBModel model_copy(model);
    model_copy.addConstr(var0 + var1 <= 1);

    model.optimize();
    std::cout << model.get(GRB_DoubleAttr_ObjVal) << "\n";
    model_copy.optimize();
    std::cout << model_copy.get(GRB_DoubleAttr_ObjVal) << "\n";
}



void test_3_wheel()
{
    VECTOR costs = {1, 1, 1, -1, -1, -1};
    EPM cost_map(4, costs);

    BnC bnc(cost_map);

    bnc.separator_callback().add_separator<CP::TriangleSeparator<BnC::EPM>>(0);
    bnc.verbosity = 100;
    bnc.optimize();
}


int main()
{
    try{
        test_3_wheel();
    } catch (GRBException e) {
        std::cout << "Error number: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
        throw e;
    } 

    return 0;
}