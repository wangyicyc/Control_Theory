#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions> // 包含MatrixExponential的头文件
#include "OsqpEigen/OsqpEigen.h"
#include <limits>

constexpr int n = 2;
constexpr int p = 1;
constexpr int Np = 3;
constexpr double T = 0.2;
constexpr int numOfCons = (2 * n + 2 * p) * Np + 2 * n;
Eigen::Matrix<double, n, 1> xMax;
Eigen::Matrix<double, n, 1> xMin;
Eigen::Matrix<double, p, 1> uMax;
Eigen::Matrix<double, p, 1> uMin;

Eigen::VectorXd lowerBound;
Eigen::VectorXd upperBound;
Eigen::Matrix2d A(n, n);
Eigen::MatrixXd B(n, p);
Eigen::Matrix2d F_;
Eigen::MatrixXd G;
Eigen::MatrixXd Fai(n *Np, n);
Eigen::MatrixXd Omega(n *Np, n *Np);
Eigen::MatrixXd Gamma = Eigen::MatrixXd::Zero(n * Np, p *Np);

Eigen::MatrixXd Q(n, n);
Eigen::MatrixXd S(n, n);
Eigen::MatrixXd R(p, p);
Eigen::MatrixXd Psi(p *Np, p *Np);

Eigen::VectorXd x(n);

Eigen::MatrixXd F;
Eigen::MatrixXd H;

Eigen::MatrixXd Mu(2 * p + 2 * n, n);
Eigen::MatrixXd Mu_Np(2 * n, n);
Eigen::MatrixXd f(2 * p + 2 * n, p);
Eigen::MatrixXd Beta(2 * p + 2 * n, 1);
Eigen::MatrixXd Beta_N(2 * n, 1);
Eigen::MatrixXd Mu_((2 * n + 2 * p) * Np + 2 * n, n);
Eigen::MatrixXd Mu__((2 * n + 2 * p) * Np + 2 * n, n *Np);
Eigen::MatrixXd f__((2 * n + 2 * p) * Np + 2 * n, p *Np);
Eigen::MatrixXd Beta_((2 * n + 2 * p) * Np + 2 * n, 1);
Eigen::MatrixXd b;
Eigen::MatrixXd M;

void get_MatOfConstraints()
{
    // 以上是无约束得到的解，下面是有约束部分

    Eigen::MatrixXd x_low(n, 1);
    x_low << -1000, -10;
    Eigen::MatrixXd x_high(n, 1);
    x_high << 1000, 10;
    Eigen::MatrixXd U_low(p, 1);
    U_low << -10;
    Eigen::MatrixXd U_high(p, 1);
    U_high << 10;
    Mu.setZero();
    Mu.block<n, n>(2 * p, 0) = -Eigen::MatrixXd::Identity(n, n);
    Mu.block<n, n>(2 * p + n, 0) = Eigen::MatrixXd::Identity(n, n);

    Mu_Np.setZero();
    Mu_Np.block<n, n>(0, 0) = -Eigen::MatrixXd::Identity(n, n);
    Mu_Np.block<n, n>(n, 0) = Eigen::MatrixXd::Identity(n, n);

    f.setZero();
    f.block<p, p>(0, 0) = -Eigen::MatrixXd::Identity(p, p);
    f.block<p, p>(p, 0) = Eigen::MatrixXd::Identity(p, p);

    Beta.block<p, 1>(0, 0) = -U_low;
    Beta.block<p, 1>(p, 0) = U_high;
    Beta.block<n, 1>(2 * p, 0) = -x_low;
    Beta.block<n, 1>(2 * p + n, 0) = x_high;

    Beta_N.setZero();
    Beta_N.block<n, 1>(0, 0) = -x_low;
    Beta_N.block<n, 1>(n, 0) = x_high;

    Mu_.setZero();
    Mu_.block<2 * n + 2 * p, n>(0, 0) = Mu;

    Mu__.setZero();

    for (int i = 0; i < Np-1; i++)
    {
        Mu__.block<2 * n + 2 * p, n>((i + 1) * (2 * n + 2 * p), i * n) = Mu;
    }
    Mu__.block<2 * n, n>((Np) * (2 * n + 2 * p), (Np - 1) * n) = Mu_Np;

    f__.setZero();
    for (int i = 0; i < Np; i++)
    {
        f__.block<2 * n + 2 * p, p>(i * (2 * n + 2 * p), i * p) = f;
    }

    Beta_.setZero();
    for (int i = 0; i < Np; i++)
    {
        Beta_.block<2 * n + 2 * p, 1>(i * (2 * n + 2 * p), 0) = Beta;
    }
    Beta_.block<2 * n, 1>(Np * (2 * n + 2 * p), 0) = Beta_N;
    /*std::cout << "Here is the matrix Beta_:\n"
              << Beta_ << std::endl;
    std::cout << "Here is size of the matrix Beta_:\n"
              << Beta_.rows() << "x" << Beta_.cols() << std::endl;*/
    M = Mu__ * Gamma + f__;
    b = -(Mu_ + Mu__ * Fai);
    //std::cout << "Here is the matrix Mu_:\n" << Mu_ << std::endl;
    //std::cout << "Here is the matrix Mu__:\n" << Mu__ << std::endl;
    //std::cout << "Here is the matrix Mu_Np:\n" << Mu_Np << std::endl;

    std::cout << "Here is the matrix f:\n" << f << std::endl;
    std::cout << "Here is the matrix f__:\n" << f__ << std::endl;
    //std::cout << "Here is the matrix Beta:\n" << Beta << std::endl;
    //std::cout << "Here is the matrix Beta_:\n" << Beta_ << std::endl;

    //std::cout << "Here is the matrix b:\n" << b << std::endl;
    //std::cout << "Here is size of the matrix b:\n" << b.rows() << "x" << b.cols() << std::endl;
}

void get_MatOfVariables()
{
    Omega.setZero();
    Psi.setZero();
    x << 0, 0;
    // 定义一个 2x2 的矩阵
    A << 0, 1,
        -2, -3;
    // 定义一个 2x1 的矩阵
    B << 0, 1;
    F_ = (A * T).exp();
    // std::cout << "Here is the matrix F_:\n" << F_ << std::endl;
    G = A.inverse() * (F_ - Eigen::Matrix2d::Identity()) * B;
    // std::cout << "Here is the matrix G:\n" << G << std::endl;
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
    //std::cout << "Here is the matrix Psi:\n" << Psi << std::endl;
    //std::cout << "Here is the matrix Gamma:\n" << Gamma << std::endl;
    //std::cout << "Here is the matrix Fai:\n" << Fai << std::endl;
    //std::cout << "Here is the matrix Omega:\n" << Omega << std::endl;
    //std::cout << "Here is the matrix F:\n" << F << std::endl;
    //std::cout << "Here is the matrix H:\n" << H << std::endl;

}

bool initMat(OsqpEigen::Solver &solver)
{
    get_MatOfVariables();
    get_MatOfConstraints();
    // std::cout << "Here is the matrix F_:\n"<< F_ << std::endl;
    // std::cout << "Here is the matrix G:\n"<< G << std::endl;
    // std::cout << "Here is the matrix Fai:\n" << Fai << std::endl;
    // std::cout << "Here is the matrix F_ * G:\n" << F_ * G << std::endl;
    // std::cout << "Here is the matrix Gamma:\n" << Gamma << std::endl;
    /*Eigen::MatrixXd U = -H.inverse() * F * x;
    std::cout << "Here is the matrix Gamma:\n"
              << Gamma << std::endl;
    std::cout << "Here is the matrix U:\n"
              << U << std::endl;*/

    // merge inequality and equality vectors
    // std::cout << "Here is the matrix H:\n" << H << std::endl;
    Eigen::SparseMatrix<double> H_sparse = H.sparseView();
    // 检查 Hessian 矩阵是否为半正定

    // std::cout << "Here is the matrix H_sparse:\n" << H_sparse << std::endl;

    // 假设 F 是一个 Eigen::MatrixXd，x 是一个 Eigen::VectorXd
    Eigen::VectorXd gradient = (F * x).transpose(); // 将结果存储在临时变量中

    Eigen::SparseMatrix<double> M_sparse = M.sparseView();
    // std::cout << "Here is the matrix H_sparse:\n" << M_sparse.cols()<< "x" <<M_sparse.rows()<< std::endl;

    //
    Eigen::VectorXd upperBound;
    upperBound.resize(numOfCons);
    upperBound = Beta_ + b * x;
    Eigen::VectorXd lowerBound(numOfCons);
    lowerBound.setZero();
    for (int i = 0; i < numOfCons; i++)
    {
        lowerBound[i] = -1000;
    }

    // std::cout << "Here is the matrix H_sparse:\n" << b.cols() << "x" << b.rows()<< std::endl;

    solver.data()->setNumberOfVariables(Np);
    solver.data()->setNumberOfConstraints(numOfCons);
    if (!solver.data()->setHessianMatrix(H_sparse))
        return false;
    if (!solver.data()->setGradient(gradient))
    {
        return false;
    }
    if (!solver.data()->setLinearConstraintsMatrix(M_sparse))
    {
        return false;
    }
    if (!solver.data()->setLowerBound(lowerBound))
    {
        return false;
    }
    if (!solver.data()->setUpperBound(upperBound))
    {
        return false;
    }
    //std::cout << "Here is the matrix upperBound:\n" << upperBound << std::endl;
    //std::cout << "Here is size of the matrix upperBound:\n" << upperBound.rows() << "x" << upperBound.cols() << std::endl;
    return true;
}

int main()
{

    OsqpEigen::Solver solver;
    solver.settings()->setWarmStart(true);
    if (initMat(solver))
    {
        if (!solver.initSolver())
            return 1;
    }
    else
    {
        std::cout << "initilize QP solver failed" << std::endl;
        return 1;
    }

    solver.solve();
    Eigen::VectorXd QPSolution;
    QPSolution = solver.getSolution();
    std::cout << "x1 = " << QPSolution[0] << std::endl
              << "x2 = " << QPSolution[1] << std::endl
              << "x3 = " << QPSolution[2] << std::endl;
    return 0;
}
