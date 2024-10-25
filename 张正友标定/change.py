import numpy as np
import sympy as sp

from scipy.spatial.transform import Rotation as R

# 定义一个函数来基于像素位置、内参矩阵、四元数、相机位置及目标高度计算位置的ENU坐标
def calculate_position_numpy(
    pixel_position: np.array,
    intrinsic_matrix: np.array,
    quaternion: np.array,
    camera_position: np.array,
    target_height
):

    plane_rotation_matrix = np.linalg.inv(R.from_quat(quaternion).as_matrix())

    plane_to_cam_rotation = np.array([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, -1]])

    camera_rotation_matrix = plane_to_cam_rotation @ plane_rotation_matrix # @表示矩阵乘法

    # 计算相机坐标系到世界坐标系的转换矩阵，即求方程AX=B的解，其中A为相机坐标系到世界坐标系的转换矩阵，X为相机坐标系的位置，B为相机位置
    t_vector = np.linalg.solve(-camera_rotation_matrix.T, # .T表示矩阵转置
                               camera_position).reshape(-1, 1)
    # 水平堆叠: camera_rotation_matrix 是一个表示旋转的 3x3 矩阵，
    # 而 t_vector 是一个表示平移的 3x1 向量（即，它在三维空间中的 x、y、z 分量）。
    # 通过 hstack，这两个矩阵被水平地堆叠在一起，形成一个 3x4 的矩阵 m_matrix
    m_matrix = np.hstack((camera_rotation_matrix, t_vector))
    # kxm_matrix 实际上是相机模型矩阵（camera model matrix），它将世界坐标系中的点转换到相机图像坐标系中的点
    kxm_matrix = intrinsic_matrix @ m_matrix
    
    # 或者使用列表推导式，如果你喜欢这种风格  
    # U = np.array([[x, y, 1]]).T  # 注意这里的双括号，因为我们需要一个二维数组，然后通过转置得到列向量
    U = np.hstack((pixel_position, np.array([1]))).reshape(-1, 1)

    solution_x_y_w = np.linalg.solve(np.hstack(
        (kxm_matrix[:, :2], -U)), - kxm_matrix[:, 2:] @ np.array([target_height, 1]))

    enu_x, enu_y, _ = solution_x_y_w

    return (enu_x, enu_y)
if __name__ == '__main__':
    fx = 575.2236  # 内参数矩阵中的焦距（x方向）
    fy = 542.1198  # 内参数矩阵中的焦距（y方向）
    cx = 310  # 图像中心点x坐标
    cy = 364  # 图像中心点y坐标
    roll = np.radians(0)  # 绕x轴旋转的角度
    pitch = np.radians(0)  # 绕y轴旋转的角度
    yaw = np.radians(0)  # 绕z轴旋转的角度(角度制)
    r = R.from_euler('xyz', [roll, pitch, yaw])  # 从欧拉角生成旋转对象

    # 转换为旋转矩阵
    quaternion = r.as_quat()  # 四元数表示相机的姿态
    pixel_position = np.array([150, 150])  # 图像上的像素位置
    intrinsic_matrix = np.array([[fx, 0, cx],  # 相机的内参数矩阵
                                 [0, fy, cy],
                                 [0, 0, 1]])
    camera_position = np.array([0, 0, 200])  # 相机在世界坐标系中的位置
    target_height = 0  # 目标点的高度

    enu_x, enu_y = calculate_position_numpy(pixel_position,  # 调用函数计算像素位置对应的ENU坐标
                                      intrinsic_matrix,
                                      quaternion,
                                      camera_position,
                                      target_height)

    print(enu_x, enu_y)  # 打印计算结果