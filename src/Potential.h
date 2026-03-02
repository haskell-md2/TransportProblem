#pragma once

#include <queue>
#include <tuple>

#include "TransportationData.h"
#include "ProblemData.h"

#define INF 100000007

using namespace std;

class Potential
{
private:

    struct Cell {
        int row;
        int col;

        Cell(int r, int c) : row(r), col(c) {}

        bool operator==(const Cell& other) const {
            return row == other.row && col == other.col;
        }
    };

    struct DFSState {
        vector<Cell> path;
        vector<vector<bool>> visited;
        Cell start;
        bool found;

        DFSState(int rows, int cols, const Cell& startCell)
            : start(startCell), found(false) {
            visited.resize(rows, vector<bool>(cols, false));
            visited[startCell.row][startCell.col] = true;
            path.push_back(startCell);
        }
    };

    static bool isOptimal(ProblemData pb_data, Approximation  apr,vector<int>& u,vector<int>& v){
        size_t rows = pb_data.numSuppliers;
        size_t columns = pb_data.numConsumers;

        for(int i = 0; i < rows; i++){
            for(int j = 0; j < columns; j++){
                if(apr.isEmpty[i][j]){
                    if (v[j] + u[i] > pb_data.costs[i][j]){
                        return false; //решение, увы, неоптимальное
                    }
                }
            }
        }

        return true;
    };



    static bool calculatePotentialsBFS(
        ProblemData pb_data, Approximation apr,
        vector<int>& u,
        vector<int>& v) {

        int m = apr.numSuppliers;
        int n = apr.numConsumers;

        u.assign(m, INF);
        v.assign(n, INF);

        vector<vector<pair<int, int>>> graph(m + n);

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                if (apr.isEmpty[i][j] == false) {
                    graph[i].push_back({m + j, pb_data.costs[i][j]});
                    graph[m + j].push_back({i, pb_data.costs[i][j]});
                }
            }
        }


        queue<int> q;
        vector<bool> visited(m + n, false);

        u[0] = 0;
        q.push(0);
        visited[0] = true;

        while (!q.empty()) {
            int current = q.front();
            q.pop();

            for (auto& edge : graph[current]) {
                int next = edge.first;
                int edgeCost = edge.second;

                if (!visited[next]) {
                    visited[next] = true;

                    if (current < m) {
                        int j = next - m;
                        v[j] = edgeCost - u[current];
                        q.push(next);
                    } else {

                        int j = current - m;
                        int i = next;
                        u[i] = edgeCost - v[j];
                        q.push(next);
                    }
                }
            }
        }

        for (int i = 0; i < m; i++) {
            if (!visited[i]) return false;
        }
        for (int j = 0; j < n; j++) {
            if (!visited[m + j]) return false;
        }

        return true;
    }


    static bool findCycleDFS(Approximation& apr, int row, int col,
                            DFSState& state, bool horizontal) {

        if (state.found) return true;

        int rows = apr.numSuppliers;
        int cols = apr.numConsumers;

        if (horizontal) {

            for (int j = 0; j < cols; j++) {
                if (j == col) continue;

                bool isBasisOrStart = !apr.isEmpty[row][j] ||
                                    (row == state.start.row && j == state.start.col);

                if (isBasisOrStart) {
                    Cell next(row, j);


                    if (next.row == state.start.row && next.col == state.start.col) {
                        if (state.path.size() >= 3) {
                            state.found = true;
                            return true;
                        }
                        continue;
                    }

                    bool alreadyInPath = false;
                    for (const auto& cell : state.path) {
                        if (cell.row == next.row && cell.col == next.col) {
                            alreadyInPath = true;
                            break;
                        }
                    }

                    if (!alreadyInPath) {
                        state.visited[row][j] = true;
                        state.path.push_back(next);


                        if (findCycleDFS(apr, row, j, state, false)) {
                            return true;
                        }

                        state.path.pop_back();
                        state.visited[row][j] = false;
                    }
                }
            }
        } else {

            for (int i = 0; i < rows; i++) {
                if (i == row) continue;

                bool isBasisOrStart = !apr.isEmpty[i][col] ||
                                    (i == state.start.row && col == state.start.col);

                if (isBasisOrStart) {
                    Cell next(i, col);


                    if (next.row == state.start.row && next.col == state.start.col) {
                        if (state.path.size() >= 3) {
                            state.found = true;
                            return true;
                        }
                        continue;
                    }

                    bool alreadyInPath = false;
                    for (const auto& cell : state.path) {
                        if (cell.row == next.row && cell.col == next.col) {
                            alreadyInPath = true;
                            break;
                        }
                    }

                    if (!alreadyInPath) {
                        state.visited[i][col] = true;
                        state.path.push_back(next);

                        if (findCycleDFS(apr, i, col, state, true)) {
                            return true;
                        }

                        state.path.pop_back();
                        state.visited[i][col] = false;
                    }
                }
            }
        }
        return false;
    }


static void redistributeAlongCycle(Approximation& apr, const vector<Cell>& cycle) {
    if (cycle.size() < 4) return;


    int minValue = INF;
    for (size_t i = 1; i < cycle.size(); i += 2) {
        const Cell& cell = cycle[i];
        if (apr.values[cell.row][cell.col] < minValue) {
            minValue = apr.values[cell.row][cell.col];
        }
    }

    if (minValue == INF || minValue == 0) return;


    for (size_t i = 0; i < cycle.size(); i++) {
        const Cell& cell = cycle[i];
        if (i % 2 == 0) {

            apr.values[cell.row][cell.col] += minValue;
            if (apr.isEmpty[cell.row][cell.col]) {
                apr.isEmpty[cell.row][cell.col] = false;
            }
        } else {

            apr.values[cell.row][cell.col] -= minValue;
        }
    }

    for (size_t i = 1; i < cycle.size(); i += 2) {
        const Cell& cell = cycle[i];
        if (apr.values[cell.row][cell.col] == 0) {
            apr.isEmpty[cell.row][cell.col] = true;
            break;
        }
    }
}

    static void runRecountCycle(Approximation& apr, vector<int> &v, vector<int> &u,ProblemData &pb_data) {


        int si = -1, sj = -1;
        int maxDelta = -INF;

        for (int i = 0; i < apr.numSuppliers; i++) {
            for (int j = 0; j < apr.numConsumers; j++) {
                if (apr.isEmpty[i][j]) {

                    int delta = u[i] + v[j] - pb_data.costs[i][j];

                    if (delta > 0) {

                        if (delta > maxDelta) {
                            maxDelta = delta;
                            si = i;
                            sj = j;
                        }
                    }
                }
            }
        }



        if (si == -1 || sj == -1) return;


        apr.isEmpty[si][sj] = false;


        Cell start(si, sj);
        DFSState state(apr.numSuppliers, apr.numConsumers, start);

        if (findCycleDFS(apr, si, sj, state, true)) {
            cout << "----------------------------------------\n";
            cout << "Найдена перспективная клетка: (" << si + 1 << ", " << sj + 1 << ")\n";
            cout << "Цикл пересчета:\n";
            for (size_t k = 0; k < state.path.size(); k++) {
                const Cell& c = state.path[k];
                cout << "Пост." << c.row + 1 << " -> Потр." << c.col + 1;
                cout << (k % 2 == 0 ? " [+]" : " [-]");
                if (k < state.path.size() - 1) cout << "  ==>  ";
            }
            cout << "\n";

            int minVal = INF;
            for (size_t k = 1; k < state.path.size(); k += 2) {
                if (apr.values[state.path[k].row][state.path[k].col] < minVal) {
                    minVal = apr.values[state.path[k].row][state.path[k].col];
                }
            }
            cout << "Объем груза для переноса (тета): " << minVal << "\n";
            cout << "----------------------------------------\n";

            redistributeAlongCycle(apr, state.path);
            if (apr.values[si][sj] == 0) {
                apr.isEmpty[si][sj] = true;
            }
        } else {

            apr.isEmpty[si][sj] = true;
        }
    }


    static void delete_temp_variables(ProblemData &pb_data, Approximation &apr) {
        if (pb_data.get_is_fictitious_supplier()) {
            apr.values.pop_back();
            if (pb_data.get_is_fictitious_consumer()) throw std::runtime_error("delete_temp_variables");
            apr.numSuppliers--;
        }

        if (pb_data.get_is_fictitious_consumer()) {
            for (auto &i: apr.values) {
                i.pop_back();
            }
            apr.numConsumers--;
        }
    }

public:

    static void printMatrix(const std::vector<std::vector<int>>& matrix,
        const std::string& name = "Матрица") {

        std::cout << "\n" << name << ":" << std::endl;


        for (const auto & i : matrix) {
            for (const auto & j : i) {
                std::cout << j << "\t";
            }
            std::cout << std::endl;
        }
    }

    static double calculateTotalCost(const ProblemData& pb_data, const Approximation& apr) {
        double totalCost = 0.0;

        for (int i = 0; i < apr.numSuppliers; i++) {
            for (int j = 0; j < apr.numConsumers; j++) {
                if (!apr.isEmpty[i][j]) {
                    totalCost += apr.values[i][j] * pb_data.costs[i][j];
                }
            }
        }

        return totalCost;
    }


    static void getOptimalWrap(ProblemData &pb_data) {
        Approximation apr = pb_data.getApproximation();
        getOptimal(pb_data, apr);
    }

    static void getOptimal(ProblemData &pb_data, Approximation &apr){
        int iteration = 1;
        while(true) {
            cout << "\n=== Итерация " << iteration++ << " ===" << endl;
            vector<int> v, u;

            if(!calculatePotentialsBFS(pb_data, apr, u, v)) {
                cout << "Ошибка: Граф вырожден (расчет потенциалов не удался)!" << endl;
                return;
            }

            if(isOptimal(pb_data, apr, u, v)) {
                cout << "\nРешение оптимально." << endl;

                delete_temp_variables(pb_data, apr); // то что сбалансировали

                printMatrix(apr.values, "Оптимальный план" );
                double res = calculateTotalCost(pb_data, apr);
                cout << "Минимальная стоимость: " << res << endl;
                break;
            }

            runRecountCycle(apr, v, u, pb_data);
        }
    }

};


