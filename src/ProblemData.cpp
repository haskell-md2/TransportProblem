//
// Created by ilya on 2/28/26.
//

#include "ProblemData.h"
#include <iostream>

ProblemData::ProblemData(int num_suppliers, int num_consumers, std::vector<int> supply, std::vector<int> demand,
                         std::vector<std::vector<double>> costs) {

    numSuppliers = num_suppliers;
    numConsumers = num_consumers;
    this->supply = supply;
    this->demand = demand;
    this->costs = costs;

    sumS = std::accumulate(supply.begin(), supply.end(), 0LL);
    sumD = std::accumulate(demand.begin(), demand.end(), 0LL);

    is_fictitious_consumer = false;
    is_fictitious_supplier = false;

    if ( sumS != sumD) set_close();

    buildNorthWest();

}

bool ProblemData::get_is_fictitious_supplier() const {
    return is_fictitious_supplier;
}

bool ProblemData::get_is_fictitious_consumer() const {
    return is_fictitious_consumer;
}

const Approximation & ProblemData::getApproximation() const {
    return first_approximation;
}

void ProblemData::set_close() {

    set_close_start_check();

    std::cout<<"sumS " <<sumS <<std::endl;
    std::cout<<"sumD " <<sumD <<std::endl;

    if (sumS > sumD) {
        // Нужно добавить фиктивного потребителя (колонку)
        long long diff = sumS - sumD;
        // Добавляем нулевую стоимость для каждого поставщика
        for (auto &row : costs) {
            row.push_back(0.0);
        }
        // Увеличиваем demand и numConsumers
        demand.push_back(static_cast<int>(diff));
        numConsumers += 1;
        is_fictitious_consumer = true;
    } else { // sumD > sumS
        long long diff = sumD - sumS;
        // Добавляем новую строку поставщика с нулевыми cost'ами длины numConsumers
        std::vector<double> newRow(numConsumers, 0.0);
        costs.push_back(std::move(newRow));
        // Увеличиваем supply и numSuppliers
        supply.push_back(static_cast<int>(diff));
        numSuppliers += 1;
        is_fictitious_supplier = true;
    }

    set_close_final_check();
}


Canonical ProblemData::get_canonical() {
    // предполагаем, что задача уже сбалансирована:
    // sum(supply) == sum(demand)

    // количество переменных
    int A_cols = numSuppliers * numConsumers;
    // убираем одну строку зависимых ограничений (обычно последнюю demand)
    int A_rows = numSuppliers + numConsumers - 1;

    // --- вектор стоимости c ---
    Eigen::VectorXd c_canonical(A_cols);
    for (int i = 0, idx = 0; i < numSuppliers; ++i) {
        for (int j = 0; j < numConsumers; ++j, ++idx) {
            c_canonical(idx) = costs[i][j];
        }
    }

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(A_rows, A_cols);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(A_rows);

    // supply-строки (0 .. m-1)
    for (int i = 0; i < numSuppliers; ++i) {
        b(i) = supply[i];
        for (int j = 0; j < numConsumers; ++j) {
            int col = i * numConsumers + j;         // index(i,j)
            A(i, col) = 1.0;
        }
    }

    // demand-строки — все кроме последней (j = 0 .. n-2)
    for (int j = 0; j < numConsumers - 1; ++j) {
        int row = numSuppliers + j;                // строка для demand j
        b(row) = demand[j];
        for (int i = 0; i < numSuppliers; ++i) {
            int col = i * numConsumers + j;
            A(row, col) = 1.0;
        }
    }

    std::vector<int> basis_indices;
    basis_indices.reserve(A_rows);
    for (int i = 0; i < numSuppliers; ++i) {
        for (int j = 0; j < numConsumers; ++j) {
            if (!first_approximation.isEmpty[i][j]) {
                int col = i * numConsumers + j;
                basis_indices.push_back(col);
            }
        }
    }

    Canonical canon(A, b, c_canonical, basis_indices);
    canon.SetOriginalVariablesCount(A_cols); // переменных = m*n
    return canon;
}

void ProblemData::buildNorthWest() {

    Approximation apr;
    apr.numSuppliers = numSuppliers;
    apr.numConsumers = numConsumers;
    apr.values.assign(apr.numSuppliers, std::vector<int>(apr.numConsumers, 0));
    apr.isEmpty.assign(apr.numSuppliers, std::vector<bool>(apr.numConsumers, true));

    std::vector<int> s = supply;
    std::vector<int> d = demand;
    int i = 0, j = 0;

    while (i < apr.numSuppliers && j < apr.numConsumers) {
        int take = std::min(s[i], d[j]);
        apr.values[i][j] = take;
        apr.isEmpty[i][j] = false;

        s[i] -= take;
        d[j] -= take;

        if (s[i] == 0 && d[j] == 0 && (i < apr.numSuppliers - 1 || j < apr.numConsumers - 1)) {
            i++;
            apr.values[i][j] = 0;
            apr.isEmpty[i][j] = false;
        } else if (s[i] == 0) {
            i++;
        } else {
            j++;
        }
    }

    first_approximation = apr;
}



void ProblemData::set_close_start_check() {
    if (numSuppliers < 0 || numConsumers < 0) {
        throw std::invalid_argument("numSuppliers/numConsumers должны быть неотрицательными");
    }
    if (static_cast<int>(supply.size()) != numSuppliers) {
        throw std::invalid_argument("Размер вектора supply не совпадает с numSuppliers");
    }
    if (static_cast<int>(demand.size()) != numConsumers) {
        throw std::invalid_argument("Размер вектора demand не совпадает с numConsumers");
    }
    // costs должен быть numSuppliers строк по numConsumers элементов в каждой
    if (static_cast<int>(costs.size()) != numSuppliers) {
        throw std::invalid_argument("Количество строк в costs не совпадает с numSuppliers");
    }
    for (const auto &row : costs) {
        if (static_cast<int>(row.size()) != numConsumers) {
            throw std::invalid_argument("Каждая строка costs должна иметь длину numConsumers");
        }
    }

    for (int v : supply) if (v < 0) throw std::invalid_argument("supply содержит отрицательное значение");
    for (int v : demand)  if (v < 0) throw std::invalid_argument("demand содержит отрицательное значение");


}

void ProblemData::set_close_final_check() {
    if (static_cast<int>(supply.size()) != numSuppliers) {
        throw std::logic_error("Внутренняя несогласованность после изменения supply");
    }
    if (static_cast<int>(demand.size()) != numConsumers) {
        throw std::logic_error("Внутренняя несогласованность после изменения demand");
    }
    if (static_cast<int>(costs.size()) != numSuppliers) {
        throw std::logic_error("Внутренняя несогласованность после изменения costs (строки)");
    }
    for (const auto &row : costs) {
        if (static_cast<int>(row.size()) != numConsumers) {
            throw std::logic_error("Внутренняя несогласованность после изменения costs (столбцы)");
        }
    }
}
