import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
from scipy.linalg import solve_discrete_are, block_diag

def get_lqr_gain(A, B, Q, R):
    """计算离散时间 LQR 增益 K 和代价矩阵 P"""
    try:
        # 求解离散代数 Riccati 方程
        P = solve_discrete_are(A, B, Q, R)
        # 计算反馈增益 K = -(R + B'PB)^-1 * B'PA
        K = -np.linalg.inv(R + B.T @ P @ B) @ B.T @ P @ A
        return K, P
    except Exception as e:
        print(f"LQR 求解失败: {e}")
        return None, None

def run_drone_example_final():
    print("\n--- 正在运行 4.3 节 无人机示例 (最终稳定版) ---")
    
    # 定义 12 维状态: [x, vx, y, vy, z, vz, roll, d_roll, pitch, d_pitch, yaw, d_yaw]
    dt = 0.1
    # 每个自由度的基础动力学
    A_i = np.array([[1.0, dt], [0.0, 1.0]])
    B_i = np.array([[0.5 * dt**2], [dt]])
    
    # 构造全维系统 (6个自由度，每个自由度都有独立控制)
    # 这样可以保证系统 100% 可控，从而 LQR 一定有解
    A = block_diag(A_i, A_i, A_i, A_i, A_i, A_i)
    B = block_diag(B_i, B_i, B_i, B_i, B_i, B_i)
    
    n, m = A.shape[0], B.shape[1] # n=12, m=6
    
    # 权重矩阵
    Q = np.eye(n) * 10
    R = np.eye(m) * 100
    O_weight = 1e6 * np.eye(n) # 偏移代价权重，很大，驱动辅助点向真实点移动
    
    # 核心：计算 K 和 P
    K, P = get_lqr_gain(A, B, Q, R)
    if K is None:
        print("致命错误：无法计算 LQR 增益。请检查 A, B 是否可控。")
        return

    # 稳态映射：(I-A)xs = B*us -> us = pinv(B) @ (I-A) @ xs
    M_theta_u = np.linalg.pinv(B) @ (np.eye(n) - A)
    
    # 初始状态与目标点
    x_curr = np.zeros(n)
    x_ref = np.zeros(n)
    x_ref[0], x_ref[2], x_ref[4] = 2.0, 3.0, 5.0 # 目标 Px=2, Py=3, Pz=5
    
    N, M_tilde = 5, 15 # 预测时域与隐式长度
    history_p = []

    for k in range(200):
        # 决策变量
        theta_a = cp.Variable(n)   # 辅助设定点
        u_v = cp.Variable((m, N + M_tilde))
        x_v = cp.Variable((n, N + M_tilde + 1))
        
        xs_a = theta_a
        us_a = M_theta_u @ theta_a
        
        # 目标函数：阶段代价 + 终端代价 + 偏移代价
        cost = cp.quad_form(xs_a - x_ref, O_weight)
        constraints = [x_v[:, 0] == x_curr]
        
        for j in range(N + M_tilde):
            constraints += [x_v[:, j+1] == A @ x_v[:, j] + B @ u_v[:, j]]
            constraints += [cp.abs(u_v[:, j]) <= 20.0] # 输入约束
            
            # 隐式终端组件的核心：在 N 步之后强制执行反馈律
            if j >= N:
                constraints += [u_v[:, j] == K @ (x_v[:, j] - xs_a) + us_a]
            
            cost += cp.quad_form(x_v[:, j] - xs_a, Q) + cp.quad_form(u_v[:, j] - us_a, R)
        
        cost += cp.quad_form(x_v[:, N + M_tilde] - xs_a, P)
        
        # 求解
        prob = cp.Problem(cp.Minimize(cost), constraints)
        prob.solve(solver=cp.ECOS)
        
        if u_v.value is None:
            print(f"Step {k}: 求解失败")
            break
            
        # 模拟一步
        x_curr = A @ x_curr + B @ u_v.value[:, 0]
        history_p.append([x_curr[0], x_curr[2], x_curr[4]]) # 记录位置
        
        if k % 10 == 0:
            print(f"Step {k}: 坐标 = {np.round(history_p[-1], 3)}")

    # 可视化
    history_p = np.array(history_p)
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(history_p[:, 0], history_p[:, 1], history_p[:, 2], 'b-o', markersize=3, label='Drone Trajectory')
    ax.scatter(x_ref[0], x_ref[2], x_ref[4], color='red', s=100, label='Target [2, 3, 5]')
    ax.set_title("12D Drone Tracking MPC with Implicit Invariant Set")
    ax.set_xlabel("X (m)"); ax.set_ylabel("Y (m)"); ax.set_zlabel("Z (m)")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    run_drone_example_final()