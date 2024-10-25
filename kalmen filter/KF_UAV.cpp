#include <iostream>
#include <Eigen/Dense>

// 初始化卡尔曼滤波器参数
Eigen::Vector3d uav_pos = Eigen::Vector3d(0.0, 0.0, 0.0); // 初始状态,假设初始位置为0
double P = 1000.0;                                        // 初始估计误差协方差,较高的不确定性
double F_ = 1;                                            // 状态转移矩阵,无变化
double H1 = 1.0;                                          // 第一个观测矩阵,直接观测位置
double R1 = 1.0;                                          // 第一个测量噪声协方差
double H2 = 1.0;                                          // 第二个观测矩阵,直接观测位置
double R2 = 2.5;                                          // 第二个测量噪声协方差
double Q_ = 0.001;                                        // 过程噪声协方差

Eigen::Vector3d KalmanFilter(const Eigen::Vector3d &z1, const double &H1, const double &R1,
                             const Eigen::Vector3d &z2, const double &H2, const double &R2)
{
    // 初始化卡尔曼滤波器
    Eigen::Vector3d x_ = Eigen::Vector3d(0.0, 0.0, 0.0);      // 系统状态
    Eigen::Vector3d x_minus = Eigen::Vector3d(0.0, 0.0, 0.0); // 先验状态
    Eigen::Vector3d P_ = Eigen::Vector3d(P, P, P);            // 后验估计误差协方差矩阵
    Eigen::Vector3d P_minus;                                  // 先验估计误差协方差矩阵
    Eigen::Vector3d K;                                        // kalman gain
    // 初始化步骤结束

    // 传感器一测量值
    //  预测步骤
    for (size_t i = 0; i < 3; i++)
    {
        x_minus(i) = F_ * x_(i);           // 预测新的状态,认为x_不变
        P_minus(i) = F_ * P_(i) * F_ + Q_; // 更新预测误差协方差，由上一次的P_得到这一次的P_先验估计
    }
    // 预测步骤结束

    // 更新步骤

    // transpose() ：矩阵转置
    for (size_t i = 0; i < 3; i++)
    {
        K(i) = P_minus(i) * H1 / (H1 * P_minus(i) * H1 + R1); // 计算卡尔曼增益
        // 更新状态
        // z为测量值，H为观测矩阵，R为测量噪声矩阵协方差
        x_(i) = x_minus(i) + (K(i) * (z1(i) - H1 * x_minus(i)));
        P_(i) = (1 - K(i) * H1) * P_minus(i); // 根据P_的先验估计，更新估计误差协方差矩阵
    }
    // 更新步骤结束
    // 传感器二测量值
    // 预测步骤
    for (size_t i = 0; i < 3; i++)
    {
        x_minus(i) = F_ * x_(i);           // 预测新的状态,认为x_不变
        P_minus(i) = F_ * P_(i) * F_ + Q_; // 更新预测误差协方差，由上一次的P_得到这一次的P_先验估计
    }
    // 预测步骤结束

    // 更新步骤

    // transpose() ：矩阵转置
    for (size_t i = 0; i < 3; i++)
    {
        K(i) = P_minus(i) * H2 / (H2 * P_minus(i) * H2 + R2); // 计算卡尔曼增益
        // 更新状态
        // z为测量值，H为观测矩阵，R为测量噪声矩阵协方差
        x_(i) = x_minus(i) + (K(i) * (z2(i) - H2 * x_minus(i)));
        P_(i) = (1 - K(i) * H2) * P_minus(i); // 根据P_的先验估计，更新估计误差协方差矩阵
    }
    // 更新步骤结束
    // 输出结果
    return x_;
};

int main()
{
    // 模拟两个传感器的测量值
    for (size_t i = 0; i < 10; i++)
    {
        Eigen::Vector3d measurement1 = Eigen::Vector3d(2.5 + i, 2.5, 2.5), measurement2 = Eigen::Vector3d(3.0 + i, 3.0, 3.0);

        // 对每个测量值执行预测和更新
        Eigen::Vector3d estimate = KalmanFilter(measurement1, H1, R1, measurement2, H2, R2);
        std::cout << "Updated state after sensor: " << estimate << std::endl;
    }

    return 0;
}