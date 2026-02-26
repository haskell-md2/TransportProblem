#pragma once

#include <queue>
#include <tuple>

#include "TransportationData.h"

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
                    if (v[j] - u[i] > pb_data.costs[i][j]){
                        return false; //решение, увы, неоптимальное
                    }
                }
            }
        }

        return true; //а вот здесь - порядок
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
            // Двигаемся по горизонтали (фиксированная строка, меняем столбцы)
            for (int j = 0; j < cols; j++) {
                if (j == col) continue; // не возвращаемся в ту же клетку сразу
                
                // Клетка должна быть базисной ИЛИ это стартовая клетка
                bool isBasisOrStart = !apr.isEmpty[row][j] || 
                                    (row == state.start.row && j == state.start.col);
                
                if (isBasisOrStart) {
                    Cell next(row, j);
                    
                    
                    // Если нашли стартовую клетку - цикл замкнулся
                    if (next.row == state.start.row && next.col == state.start.col) {
                        if (state.path.size() >= 3) { // Минимум 3 шага до возврата
                            state.found = true;
                            return true;
                        }
                        continue; // Слишком короткий цикл
                    }
                    
                    // ВАЖНО: Не идем в клетки, которые уже есть в пути (кроме стартовой)
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
                        
                        // Меняем направление на вертикальное
                        if (findCycleDFS(apr, row, j, state, false)) {
                            return true;
                        }
                        
                        // Backtracking
                        state.path.pop_back();
                        state.visited[row][j] = false;
                    }
                }
            }
        } else {
            // Двигаемся по вертикали (фиксированный столбец, меняем строки)
            for (int i = 0; i < rows; i++) {
                if (i == row) continue; // не возвращаемся в ту же клетку
                
                // Клетка должна быть базисной ИЛИ это стартовая клетка
                bool isBasisOrStart = !apr.isEmpty[i][col] || 
                                    (i == state.start.row && col == state.start.col);
                
                if (isBasisOrStart) {
                    Cell next(i, col);
                    
                    
                    // Если нашли стартовую клетку - цикл замкнулся
                    if (next.row == state.start.row && next.col == state.start.col) {
                        if (state.path.size() >= 3) {
                            state.found = true;
                            return true;
                        }
                        continue;
                    }
                    
                    // ВАЖНО: Не идем в клетки, которые уже есть в пути (кроме стартовой)
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
                                           
                        // Меняем направление на горизонтальное
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
    
    // Пересчет поставок по найденному циклу
static void redistributeAlongCycle(Approximation& apr, const vector<Cell>& cycle) {
    if (cycle.size() < 4) return;
    
    // Определяем минимальное значение в НЕЧЕТНЫХ вершинах цикла
    // (вершины, в которых будем ВЫЧИТАТЬ)
    int minValue = INF;
    for (size_t i = 1; i < cycle.size(); i += 2) {
        const Cell& cell = cycle[i];
        if (apr.values[cell.row][cell.col] < minValue) {
            minValue = apr.values[cell.row][cell.col];
        }
    }
    
    if (minValue == INF || minValue == 0) return;
    
    // Перераспределяем:
    // Четные вершины (0,2,4...) - добавляем minValue
    // Нечетные вершины (1,3,5...) - вычитаем minValue
    for (size_t i = 0; i < cycle.size(); i++) {
        const Cell& cell = cycle[i];
        if (i % 2 == 0) {
            // Четные вершины - добавляем
            apr.values[cell.row][cell.col] += minValue;
            if (apr.isEmpty[cell.row][cell.col]) {
                apr.isEmpty[cell.row][cell.col] = false;
            }
        } else {
            // Нечетные вершины - вычитаем
            apr.values[cell.row][cell.col] -= minValue;
        }
    }
    
    // Проверяем, какая клетка обнулилась
    for (size_t i = 1; i < cycle.size(); i += 2) {
        const Cell& cell = cycle[i];
        if (apr.values[cell.row][cell.col] == 0) {
            apr.isEmpty[cell.row][cell.col] = true;
            break; // Только одна клетка должна стать пустой
        }
    }
}
    
    static void runRecountCycle(Approximation& apr, vector<int> &v, vector<int> &u,ProblemData &pb_data) {

        // Ищем первую пустую клетку 
        int si = -1, sj = -1;
        int maxDelta = -INF;  // Максимальная положительная разность
        
        for (int i = 0; i < apr.numSuppliers; i++) {
            for (int j = 0; j < apr.numConsumers; j++) {
                if (apr.isEmpty[i][j]) {
                    // Вычисляем оценку для пустой клетки
                    // Δij = (ui + vj) - cij
                    // В оптимальном плане должно быть ≤ 0
                    int delta = u[i] + v[j] - pb_data.costs[i][j];
                    
                    // Если оценка положительная, клетка перспективная
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
        
        // Ищем цикл DFS
        Cell start(si, sj);
        DFSState state(apr.numSuppliers, apr.numConsumers, start);
        
        if (findCycleDFS(apr, si, sj, state, true)) {
            // Нашли цикл, выполняем перераспределение
            redistributeAlongCycle(apr, state.path);
            
            // Возвращаем исходную клетку в пустое состояние
            // (она станет заполненной только если в нее добавились поставки)
            if (apr.values[si][sj] == 0) {
                apr.isEmpty[si][sj] = true;
            }
        } else {
            // Если цикл не найден, возвращаем исходное состояние
            apr.isEmpty[si][sj] = true;
        }
    }


    static void printMatrix(const std::vector<std::vector<int>>& matrix, const std::string& name = "Матрица") {
        std::cout << "\n" << name << ":" << std::endl;
        for (const auto& row : matrix) {
            for (int val : row) {
                std::cout << val << "\t";
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

public:

    static void getOptimal(ProblemData &pb_data, Approximation  &apr){

        vector<int> v;
        vector<int> u;

        calculatePotentialsBFS(pb_data,apr,u,v);

        if(isOptimal(pb_data,apr,u,v)){
            cout << "Успех" << endl;
            printMatrix(apr.values,"Оптимальный план");
            double res = calculateTotalCost(pb_data,apr);
            cout <<  res << endl;
            return;
        }

        runRecountCycle(apr,v,u,pb_data);
        getOptimal(pb_data,apr);

    }


};


