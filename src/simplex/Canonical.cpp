#include "Canonical.h"

#include <iostream>
#include <iomanip>
#include <stdexcept>

struct Canonical::Impl
{
    Eigen::MatrixXd A;              
    Eigen::VectorXd b;              
    Eigen::VectorXd c;              
    std::vector<int> basisIndices;  
    bool minimize;                  
    int originalVariablesCount;     

    Impl(const Eigen::MatrixXd& A_,
         const Eigen::VectorXd& b_,
         const Eigen::VectorXd& c_,
         const std::vector<int>& basisIndices_,
         bool minimize_)
        : A(A_), b(b_), c(c_),
          basisIndices(basisIndices_),
          minimize(minimize_),
          originalVariablesCount(c_.size())
    {
        if (A.rows() != b.size())
        {
            throw std::invalid_argument("Размерность матрицы A и вектора b не совпадают");
        }
        if (A.cols() != c.size())
        {
            throw std::invalid_argument("Размерность матрицы A и вектора c не совпадают");
        }
        if (static_cast<int>(basisIndices.size()) != A.rows())
        {
            throw std::invalid_argument("Количество базисных индексов не совпадает с количеством строк A");
        }

        for (int idx : basisIndices)
        {
            if (idx < 0 || idx >= A.cols())
            {
                throw std::invalid_argument("Базисный индекс выходит за пределы допустимого диапазона");
            }
        }
    }
};

Canonical::Canonical(const Eigen::MatrixXd& A,
                     const Eigen::VectorXd& b,
                     const Eigen::VectorXd& c,
                     const std::vector<int>& basisIndices)
    : pimpl_(std::make_unique<Impl>(A, b, c, basisIndices, true))
{
}

Canonical::Canonical(const Canonical& other)
    : pimpl_(std::make_unique<Impl>(*other.pimpl_))
{
}

Canonical::Canonical(Canonical&& other) noexcept = default;

Canonical& Canonical::operator=(const Canonical& other)
{
    if (this != &other)
    {
        pimpl_ = std::make_unique<Impl>(*other.pimpl_);
    }
    return *this;
}

Canonical& Canonical::operator=(Canonical&& other) noexcept = default;

Canonical::~Canonical() = default;

double Canonical::Evaluate(const Eigen::VectorXd& solution) const
{
    if (solution.size() != pimpl_->c.size())
    {
        throw std::invalid_argument("Размерность решения не совпадает с количеством переменных");
    }
    
    return pimpl_->c.dot(solution);
}

void Canonical::PrintCoffMatrixWithB() {
    for (int i = 0; i < pimpl_->A.rows(); ++i) {
        for (int j = 0; j < pimpl_->A.cols(); ++j) {
            std::cout << pimpl_->A(i, j) << ' ';
        }
        std::cout << pimpl_->b[i] << '\n';
    }
}

void Canonical::PrintObjectiveCoefficients() {
    for (int i = 0; i < pimpl_->c.size(); ++i) {
        std::cout << pimpl_->c[i];
        if (i + 1 < pimpl_->c.size()) std::cout << ' ';
    }
    std::cout << '\n';
}

void Canonical::Print() const
{
    std::cout << "=== Каноническая форма задачи ЛП ===" << std::endl;
    std::cout << (pimpl_->minimize ? "Минимизировать: " : "Максимизировать: ");
    
    for (int i = 0; i < pimpl_->c.size(); ++i)
    {
        if (i > 0 && pimpl_->c[i] >= 0) std::cout << " + ";
        std::cout << pimpl_->c[i] << "*x" << (i + 1);
    }
    std::cout << std::endl << std::endl;

    std::cout << "При ограничениях (Ax = b):" << std::endl;
    for (int i = 0; i < pimpl_->A.rows(); ++i)
    {
        for (int j = 0; j < pimpl_->A.cols(); ++j)
        {
            if (j > 0 && pimpl_->A(i, j) >= 0) std::cout << " + ";
            std::cout << pimpl_->A(i, j) << "*x" << (j + 1);
        }
        std::cout << " = " << pimpl_->b[i] << std::endl;
    }

    std::cout << std::endl << "Все переменные неотрицательны: x_i >= 0" << std::endl;

    std::cout << std::endl << "Базисные переменные: ";
    for (size_t i = 0; i < pimpl_->basisIndices.size(); ++i)
    {
        if (i > 0) std::cout << ", ";
        std::cout << "x" << (pimpl_->basisIndices[i] + 1);
    }
    std::cout << std::endl;

    std::cout << "Количество исходных переменных: " << pimpl_->originalVariablesCount << std::endl;
    std::cout << "Дополнительных переменных: " << (pimpl_->c.size() - pimpl_->originalVariablesCount) << std::endl;
}

const Eigen::MatrixXd& Canonical::GetConstraintsMatrix() const
{
    return pimpl_->A;
}

const Eigen::VectorXd& Canonical::GetRightHandSide() const
{
    return pimpl_->b;
}

const Eigen::VectorXd& Canonical::GetObjectiveCoefficients() const
{
    return pimpl_->c;
}

bool Canonical::IsMaximization() const
{
    return !pimpl_->minimize;
}

const std::vector<int>& Canonical::GetBasisIndices() const
{
    return pimpl_->basisIndices;
}

int Canonical::GetOriginalVariablesCount() const
{
    return pimpl_->originalVariablesCount;
}

void Canonical::SetOriginalVariablesCount(int count)
{
    if (count <= 0 || count > pimpl_->c.size())
    {
        throw std::invalid_argument("Некорректное количество исходных переменных");
    }
    pimpl_->originalVariablesCount = count;
}

bool Canonical::IsFeasibleBasis() const
{
    Eigen::VectorXd solution = GetBasicSolution();
    
    for (int i = 0; i < solution.size(); ++i)
    {
        if (solution[i] < -1e-9)
        {
            return false;
        }
    }
    return true;
}

Eigen::VectorXd Canonical::GetBasicSolution() const
{
    Eigen::VectorXd solution = Eigen::VectorXd::Zero(pimpl_->c.size());
    
    Eigen::MatrixXd B(pimpl_->A.rows(), pimpl_->basisIndices.size());
    for (size_t i = 0; i < pimpl_->basisIndices.size(); ++i)
    {
        B.col(i) = pimpl_->A.col(pimpl_->basisIndices[i]);
    }
    
    Eigen::VectorXd basicValues = B.colPivHouseholderQr().solve(pimpl_->b);
    
    for (size_t i = 0; i < pimpl_->basisIndices.size(); ++i)
    {
        solution[pimpl_->basisIndices[i]] = basicValues[i];
    }
    
    return solution;
}

