#include <inequalities/hypermetric.hxx>
#include <property-maps.hxx>

typedef CP::EdgePropertyMap<double> EPM;


void test(const bool& b)
{
    if (!b)
        throw std::logic_error("Test Failed!");
}


void test_hyper_metric_separation()
{
    size_t n = 5;
    EPM ecm(n, 1.0/3.0);

    for (size_t i = 1; i < n; ++i)
        ecm(0, i) = 2.0/3.0;

    CP::HyperMetricSeparator<EPM> sep(n);
    auto hyper_metrics = sep.separate_hyper_metric(ecm);
    
    test (hyper_metrics.size() == 6);
    for (auto h : hyper_metrics)
    {
        test (std::abs(h.violation - 1.0/3) < 1e-6);
        test (std::abs(h.compute_violation(ecm) - 1.0/3) < 1e-6);
        test (std::abs(h.violation_depth - (1.0/3) / std::sqrt(22) ) < 1e-6);
        test (std::abs(h.compute_violation_depth() - (1.0/3) / std::sqrt(22) ) < 1e-6);
        test (h.b.size() == 5);
        test (h.b[0] == -2);
        test (h.b[1] == 1);
        test (h.b[2] == 1);
        test (h.b[3] == 1);
        test (h.b[4] == 1);

    }
}


int main()
{
    test_hyper_metric_separation();

    return 0;
}