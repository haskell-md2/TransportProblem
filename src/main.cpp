#include <iostream>

#include "Potential.h"
#include "simplex/Canonical.h"
#include "simplex/SimplexSolver.h"


// void print_simplex_solution(Eigen::VectorXd &simplex_solution, int num_consumers_cols) {
//     for (int i = 0; i < simplex_solution.size(); i++) {
//         if (i % (num_consumers_cols) == 0 ) std::cout << std::endl;
//         std::cout << simplex_solution[i] << "  ";
//     }
//     std::cout << std::endl;
// }

// double calc_simplex_cost(Eigen::VectorXd &simplex_solution, std::vector<std::vector<double>> costs, int num_consumers_cols) {
//     double res = 0;
//     for (int i = 0; i < costs.size(); i++) {
//         for (int j = 0; j < costs[0].size(); j++) {
//             res += costs[i][j] * simplex_solution[i * num_consumers_cols + j];
//         }
//     }
//
//     return res;
// }


void delete_temp_variables(ProblemData &pb_data,  std::vector<std::vector<double>> &matrix) {
    if (pb_data.get_is_fictitious_supplier()) {
        matrix.pop_back();
        if (pb_data.get_is_fictitious_consumer()) throw std::runtime_error("delete_temp_variables");
        // apr.numSuppliers--;
    }

    if (pb_data.get_is_fictitious_consumer()) {
        for (auto &i: matrix) {
            i.pop_back();
        }
        // apr.numConsumers--;
    }
}


void solve_simplex(ProblemData& pd) {
    std::cout << std::endl;
    std::cout << "Simplex method:" << std::endl;

    Canonical canonical = pd.get_canonical();
    Solver s(canonical);

    Eigen::VectorXd simplex_solution = s.solve();

    std::vector<std::vector<double>> matrix(pd.numSuppliers, std::vector<double> (pd.numConsumers, 0));
    //
    for (int i = 0; i < pd.numSuppliers; i++) {
        for (int j = 0; j < pd.numConsumers; j++) {
            matrix[i][j] = simplex_solution[i * pd.numConsumers + j];
        }
    }

    delete_temp_variables(pd, matrix);

    double res = 0;

    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[0].size(); j++) {

            std::cout << matrix[i][j] << "  ";
            res += matrix[i][j] * pd.costs[i][j];
        }

        std::cout << std::endl;
    }

    std::cout << res << std::endl;

}


int main(int argc, char* argv[]){

    int num_consumers = 5;
    int num_suppliers = 5;

    //Ожидаем минимальную стоимость 2440
    ProblemData pd{
        5,5 , 
        {180, 220, 150,200,250},
        {170, 200, 190, 210, 220},
        {
            {7, 3, 5, 8 ,4},
            {2, 6, 4, 9, 3},
            {8, 2, 1, 6, 7},
            {4, 7, 9, 2, 5},
            {5, 8, 6, 4, 3}
        }
    };

    // if (argc > 1){
    //     if(*argv[1] == 'h'){
    //         pd = ProblemData(
    //     6,5 ,
    //     {180, 220, 150,200,250,150},
    //     {170, 200, 190, 210, 220},
    //     {
    //         {7, 3, 5, 8 ,4},
    //         {2, 6, 4, 9, 3},
    //         {8, 2, 1, 6, 7},
    //         {4, 7, 9, 2, 5},
    //         {5, 8, 6, 4, 3},
    //         {1,1000,1000,10000,1000}
    //     });
    //     num_suppliers = 6;
    //     }
    //
    // }

    Potential::getOptimalWrap(pd);
    solve_simplex(pd);

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
