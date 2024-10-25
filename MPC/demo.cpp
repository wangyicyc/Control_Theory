#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions> // 包含MatrixExponential的头文件
#include "OsqpEigen/OsqpEigen.h"
#include <limits>

constexpr int n = 2;
constexpr int p = 1;
constexpr int Np = 20;
constexpr double T = 0.05;

Eigen::MatrixXd A(n, n);
Eigen::MatrixXd B(n, p);
Eigen::MatrixXd F_;
Eigen::MatrixXd G;
Eigen::MatrixXd Fai(n *Np, n);
Eigen::MatrixXd Omega(n *Np, n *Np);
Eigen::MatrixXd Gamma = Eigen::MatrixXd::Zero(n * Np, p *Np);

Eigen::MatrixXd Q(n, n);
Eigen::MatrixXd S(n, n);
Eigen::MatrixXd R(p, p);
Eigen::MatrixXd Psi(p *Np, p *Np);

Eigen::VectorXd xk(n);

Eigen::MatrixXd F;
Eigen::MatrixXd H;
void get_MatOfVariables()
{
    Omega.setZero();
    Psi.setZero();

    // 定义一个 2x2 的矩阵
    A << 0, 1,
        0, 0;
    // 定义一个 2x1 的矩阵
    B << 0, 1;
    F_ = (A * T).exp();
    //std::cout << "Here is the matrix F_:\n" << F_ << std::endl;
    // G = A.inverse() * (F_ - Eigen::MatrixXd::Identity(n,n)) * B;
    G = B * T;
    //std::cout << "Here is the matrix G:\n" << G << std::endl;
    for (int i = 0; i < Np; i++)
    {
        Fai.block<n, n>(i * n, 0) = F_.pow(i + 1);
    }
    for (int col = 0; col < Np; col++)
    {
        for (int row = col; row < Np; row++)
        {
            Gamma.block<n, p>(n * row, col) = F_.pow(row - col) * G;
        }
    }
    Q << 1, 0,
        0, 1;

    S << 1, 0,
        0, 1;
    Omega.block<n, n>(n * Np - n, n * Np - n) = S;
    for (int i = 0; i < Np - 1; i++)
    {
        Omega.block<n, n>(i * n, i * n) = Q;
    }
    R << 1;

    for (int i = 0; i < Np; i++)
    {
        Psi.block<p, p>(i * p, i * p) = R;
    }
    F = Gamma.transpose() * Omega * Fai;

    H = Gamma.transpose() * Omega * Gamma + Psi;
    // std::cout << "Here is the matrix Psi:\n" << Psi << std::endl;
    // std::cout << "Here is the matrix Gamma:\n" << Gamma << std::endl;
    // std::cout << "Here is the matrix Fai:\n" << Fai << std::endl;
    // std::cout << "Here is the matrix Omega:\n" << Omega << std::endl;
    // std::cout << "Here is the matrix F:\n" << F << std::endl;
    // std::cout << "Here is the matrix H:\n" << H << std::endl;
}

int main()
{
    xk << 2, 0;
    for (int i = 0; i < 50; i++)
    {

        get_MatOfVariables();
        Eigen::VectorXd QPSolution = -H.inverse() * F * xk;
        xk(0) = xk(1) * T;
        xk(1) = QPSolution[0] * T;
        std::cout << xk << std::endl;
        /*std::cout << "x1 = " << QPSolution[0] << std::endl
                  << "x2 = " << QPSolution[1] << std::endl
                  << "x3 = " << QPSolution[2] << std::endl;*/
    }
    return 0;
}
