#include <iostream>

#include "Potential.h"
#include "simplex/Canonical.h"
#include "simplex/SimplexSolver.h"


void print_simplex_solution(Eigen::VectorXd &simplex_solution, int num_consumers_cols) {
    for (int i = 0; i < simplex_solution.size(); i++) {
        if (i % (num_consumers_cols) == 0 ) std::cout << std::endl;
        std::cout << simplex_solution[i] << "  ";
    }
    std::cout << std::endl;
}

double calc_simplex_cost(Eigen::VectorXd &simplex_solution, std::vector<std::vector<double>> costs, int num_consumers_cols) {
    double res = 0;
    for (int i = 0; i < costs.size(); i++) {
        for (int j = 0; j < costs[0].size(); j++) {
            res += costs[i][j] * simplex_solution[i * num_consumers_cols + j];
        }
    }

    return res;
}

int main(){

    int num_consumers = 4;
    int num_suppliers = 3;

    ProblemData pd{num_suppliers, num_consumers,{27,20,43},{33,13,27,17},
    {
        {14,28,21,28},
        {10,17,15,24},
        {14,13,25,21}
    }};

    // Ожидается ответ: 1435
    Potential::getOptimalWrap(pd);
    
    std::cout << std::endl;
    std::cout << "Simplex method:" << std::endl;

    Canonical canonical = pd.get_canonical();
    Solver s(canonical);

    Eigen::VectorXd simplex_solution = s.solve();

    print_simplex_solution(simplex_solution, num_consumers);
    std::cout << calc_simplex_cost(simplex_solution, pd.costs, num_consumers) << std::endl;

    return 0;
}


// Approximation apr{3,4,
// {
//     {27,0,0,0},
//     {6,13,1,0},
//     {0,0,26,17}
// },
// {
//     {0,1,1,1},
//     {0,0,0,1},
//     {1,1,0,0}
// }};
