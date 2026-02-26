#include <iostream>

#include "Potential.h"
#include "TransportationData.h"

int main(){

    ProblemData pd{3,4,{27,20,43},{33,13,27,17},
    {
        {14,28,21,28},
        {10,17,15,24},
        {14,13,25,21}
    }};

    Approximation apr{3,4,
    {
        {27,0,0,0},
        {6,13,1,0},
        {0,0,26,17}
    },    
    {
        {0,1,1,1},
        {0,0,0,1},
        {1,1,0,0}
    }};

    // Ожидается ответ: 1435
    Potential::getOptimal(pd,apr);

    return 0;
}