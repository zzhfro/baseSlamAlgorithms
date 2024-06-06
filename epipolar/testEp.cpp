#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include "epipolarAlgorithm.h"
using namespace Eigen;

int main() {
    
    // 打开文本文件
    std::ifstream file("/home/zzhfro/code/ImageBasedModellingEduV1.0/examples/task2/correspondences.txt");
    if (!file.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return 1;
    }

    // 创建矩阵 A 和 B
    MatrixXd A, B;
    int num_cols = 2; // 每行前两个数据存入矩阵 A
    A.resize(0, num_cols);
    B.resize(0, num_cols);
    int n_line=0;
    std::string line;
    while (std::getline(file, line)) {
     if(n_line>0)
     {
        std::istringstream iss(line);
        double value;
        VectorXd row_A(num_cols), row_B(num_cols);
        int col_index = 0;
        while (iss >> value) {
            if (col_index < num_cols) {
                row_A(col_index) = value; // 前两个数据存入矩阵 A
            } else {
                row_B(col_index - num_cols) = value; // 后两个数据存入矩阵 B
            }
            col_index++;
        }
        // 将数据添加到矩阵 A 和 B
        A.conservativeResize(A.rows() + 1, num_cols);
        B.conservativeResize(B.rows() + 1, num_cols);
        A.row(A.rows() - 1) = row_A;
        B.row(B.rows() - 1) = row_B;
    }
    n_line++;
    }

    // 关闭文件
    file.close();
    Eigen::Matrix3d f;
    int num=fundamental_matrix_ransac(A,B,f);
    std::cout<<num<<std::endl;
    std::cout<<"fundamental matrix :"<<std::endl;
    std::cout<<f<<std::endl;
    

     
}


