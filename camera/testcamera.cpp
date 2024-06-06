#include <iostream>
#include <vector>
#include <Eigen/Core>
#include <chrono>

int main() {
    const int rows = 1000;
    const int cols = 1000;
     Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(rows, cols);
    // Benchmark Eigen matrix row assignment
    Eigen::MatrixXd eigenMatrix(rows, cols);
    auto startEigen = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < rows; ++i) {
        eigenMatrix.row(i) = matrix.row(i);
    }
    auto endEigen = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> durationEigen = endEigen - startEigen;
    std::cout << "Eigen matrix row assignment time: " << durationEigen.count() << " seconds" << std::endl;

    // Benchmark std::vector push_back
    std::vector<Eigen::VectorXd> vec;
    auto startVec = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < rows; ++i) {
        vec.push_back(Eigen::VectorXd::Random(cols));
    }
    auto endVec = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> durationVec = endVec - startVec;
    std::cout << "std::vector push_back time: " << durationVec.count() << " seconds" << std::endl;

    return 0;
}

