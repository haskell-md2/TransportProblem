//
// Created by ilya on 2/28/26.
//

#ifndef TRANSPORTPROBLEM_PROBLEMDATA_H
#define TRANSPORTPROBLEM_PROBLEMDATA_H



#include <vector>
#include <numeric> // for accumulate
#include "simplex/Canonical.h"
#include <stdexcept>    // std::invalid_argument

#include "TransportationData.h"

class  ProblemData {

    void set_close_start_check();
    void set_close_final_check();
    void buildNorthWest();
    void set_close();

    long long sumS;
    long long sumD;

    bool is_fictitious_supplier;
    bool is_fictitious_consumer;

public:

    int numSuppliers;
    int numConsumers;
    std::vector<int> supply;
    std::vector<int> demand;
    std::vector<std::vector<double>> costs;

    ProblemData(int num_suppliers, int num_consumers, std::vector<int> supply, std::vector<int> demand, std::vector<std::vector<double>> costs);

    bool get_is_fictitious_supplier() const;
    bool get_is_fictitious_consumer() const;

    Approximation first_approximation;
    Canonical get_canonical();
};


#endif //TRANSPORTPROBLEM_PROBLEMDATA_H