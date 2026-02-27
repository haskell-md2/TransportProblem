#include <iostream>
#include "Potential.h"
#include "TransportationData.h"

using namespace std;

int main() {

    ProblemData pd{
        3, 4, 
        {27, 20, 43},
        {33, 13, 27, 17},
        {
            {14, 28, 21, 28},
            {10, 17, 15, 24},
            {14, 13, 25, 21}
        }
    };

    cout << "=== ЭТАП 1: Поиск начального приближения ===" << endl;
    
    Approximation apr = Potential::buildNorthWest(pd);
    
    Potential::printMatrix(apr.values, "Начальный план (Северо-Западный угол)");
    
    double initialCost = Potential::calculateTotalCost(pd, apr);
    cout << "Стоимость начального плана: " << initialCost << "\n\n";

    cout << "=== ЭТАП 2: Оптимизация методом потенциалов ===" << endl;
    
    // ожидается 1435
    Potential::getOptimal(pd, apr);

    return 0;
}
