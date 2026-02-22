#include <vector>
#include <iostream>

struct ProblemData {
    int numSuppliers; // 5
    int numConsumers; // 5
    std::vector<int> supply; // Запасы поставщиков
    std::vector<int> demand; // Потребности потребителей
    std::vector<std::vector<double>> costs; // Матрица стоимостей
};

// Абстрактный класс (интерфейс) для солвера
class ISolver {
public:
    virtual void solve(ProblemData data) = 0;
    virtual ~ISolver() = default;
};
