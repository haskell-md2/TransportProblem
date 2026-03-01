#pragma once

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
