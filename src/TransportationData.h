#pragma once

#include <vector>
#include <iostream>

struct ProblemData {
    int numSuppliers;
    int numConsumers;
    std::vector<int> supply;
    std::vector<int> demand;
    std::vector<std::vector<double>> costs;
};

struct Approximation
{
    int numSuppliers;
    int numConsumers;
    std::vector<std::vector<int>> values;
    std::vector<std::vector<bool>> isEmpty;
};


// class ISolver {
// public:
//     virtual void solve(ProblemData data) = 0;
//     virtual ~ISolver() = default;
// };
