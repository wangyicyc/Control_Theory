#include <iostream>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class KalmanFilter
{
public:
    // 构造函数
    KalmanFilter() {}

    // 初始化卡尔曼滤波器
    void init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in, MatrixXd &Q_in)
    {
        x_ = x_in;      // 状态向量
        P_ = P_in;      // 后验估计误差协方差
        P_minus = P_in; // 先验估计误差协方差
        F_ = F_in;      // 状态转移矩阵
        Q_ = Q_in;      // 过程噪声协方差
    }

    // 预测步骤
    void predict()
    {
        x_minus = F_ * x_;                            // 预测新的状态,认为x_不变
        P_minus = F_ * P_ * F_.transpose() + Q_; // 更新预测误差协方差，由上一次的P_得到这一次的P_先验估计
    }

    // 更新步骤
    void update(const VectorXd &z, const MatrixXd &H, const MatrixXd &R)
    {
        // transpose() ：矩阵转置
        MatrixXd K = P_minus * H.transpose() * (H * P_minus * H.transpose() + R).inverse(); // 计算卡尔曼增益

        // 更新状态
        // z为测量值，H为观测矩阵，R为测量噪声矩阵协方差
        x_ = x_minus + (K * (z - H * x_minus));
        int x_size = x_.size();
        MatrixXd I = MatrixXd::Identity(x_size, x_size);

        P_ = (I - K * H) * P_minus; // 根据P_的先验估计，更新估计误差协方差矩阵
    }

    // 获取当前状态
    VectorXd getState() const { return x_; }

private:
    VectorXd x_;      // 系统状态
    VectorXd x_minus; // 先验状态
    MatrixXd P_;      // 后验估计误差协方差矩阵
    MatrixXd P_minus; // 先验估计误差协方差矩阵
    MatrixXd F_;      // 状态转移矩阵，表示x_的状态随时间变化的关系
    MatrixXd Q_;      // 过程噪声矩阵协方差
};

int main()
{
    // 初始化卡尔曼滤波器参数
    VectorXd x(1);     // 初始状态
    x << 0.0;          // 假设初始位置为0
    MatrixXd P(1, 1);  // 初始估计误差协方差
    P << 1000.0;       // 较高的不确定性
    MatrixXd F(1, 1);  // 状态转移矩阵
    F << 1.0;          // 无变化
    MatrixXd H1(1, 1); // 第一个观测矩阵
    H1 << 1.0;         // 直接观测
    MatrixXd R1(1, 1); // 第一个测量噪声协方差
    R1 << 1.0;         // 测量噪声
    MatrixXd H2(1, 1); // 第二个观测矩阵
    H2 << 1.0;         // 直接观测
    MatrixXd R2(1, 1); // 第二个测量噪声协方差
    R2 << 2.5;         // 测量噪声
    MatrixXd Q(1, 1);  // 过程噪声协方差
    Q << 0.001;        // 过程噪声

    // 创建并初始化卡尔曼滤波器
    KalmanFilter kf;
    kf.init(x, P, F, Q);

    // 模拟两个传感器的测量值
    VectorXd measurement1(1), measurement2(1);
    measurement1 << 2.5; // 第一个传感器测量
    measurement2 << 3.0; // 第二个传感器测量

    // 对每个测量值执行预测和更新
    kf.predict();
    kf.update(measurement1, H1, R1);
    std::cout << "Updated state after sensor 1: " << kf.getState().transpose() << std::endl;

    kf.predict();
    kf.update(measurement2, H2, R2);
    std::cout << "Updated state after sensor 2: " << kf.getState().transpose() << std::endl;

    return 0;
}