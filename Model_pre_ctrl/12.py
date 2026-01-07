import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
import scipy.linalg

def solve_lqr(A, B, Q, R):
    """求解离散代数Riccati方程以获取K和P"""
    P = scipy.linalg.solve_discrete_are(A, B, Q, R)
    K = -np.linalg.inv(R + B.T @ P @ B) @ B.T @ P @ A
    return K, P

# ==========================================
# 1. 系统参数定义 (根据论文 4.1 节)
# ==========================================
A = np.array([[1.1, 1], [0, 1.1]])
B = np.array([[1, 0.5], [1, 1]])
n, m = B.shape

# 约束设置
x_max = 5.0
u_max = 0.3

# 权重矩阵
Q = np.eye(n)
R = np.eye(m)
S_offset = 10 * np.eye(n)  # 偏移代价权重

# 计算 LQR 增益 K 和终端代价矩阵 P
K, P = solve_lqr(A, B, Q, R)

# 稳态映射矩阵 M_theta: [xs; us] = M_theta * theta
# 对于该系统，取 theta = xs, 则 us = pinv(B) @ (I - A) @ xs
M_theta_x = np.eye(n)
M_theta_u = np.linalg.pinv(B) @ (np.eye(n) - A)

# ==========================================
# 2. 控制器配置
# ==========================================
# 增加预测时域 N 以提高初始可行性
N = 5         
M_tilde = 4   # 论文计算得出的隐式长度
steps = 40    # 仿真步数

# ==========================================
# 3. 仿真循环
# ==========================================
x_history = np.zeros((n, steps + 1))
u_history = np.zeros((m, steps))
# 初始状态设为更合理的范围（或增大N），确保在u=0.3约束内可控
x_history[:, 0] = np.array([-2.0, 0.5]) 
x_ref = np.array([2.0, 1.0]) # 目标设定点

print("开始仿真...")

for k in range(steps):
    x_k = x_history[:, k]
    
    # 定义优化变量
    u_vars = cp.Variable((m, N + M_tilde))
    x_vars = cp.Variable((n, N + M_tilde + 1))
    theta_a = cp.Variable(n) # 辅助设定点决策变量
    
    # 辅助稳态
    xs_a = M_theta_x @ theta_a
    us_a = M_theta_u @ theta_a
    
    cost = 0
    constraints = [x_vars[:, 0] == x_k]
    
    # --- 第一阶段: 自由预测 (0 to N-1) ---
    for j in range(N):
        cost += cp.quad_form(x_vars[:, j] - xs_a, Q) + cp.quad_form(u_vars[:, j] - us_a, R)
        constraints += [x_vars[:, j+1] == A @ x_vars[:, j] + B @ u_vars[:, j]]
        constraints += [cp.abs(x_vars[:, j]) <= x_max]
        constraints += [cp.abs(u_vars[:, j]) <= u_max]
        
    # --- 第二阶段: 隐式终端段 (N to N+M_tilde-1) ---
    # 核心：锁定反馈律 u = K(x - xs_a) + us_a
    for j in range(N, N + M_tilde):
        constraints += [u_vars[:, j] == K @ (x_vars[:, j] - xs_a) + us_a]
        constraints += [x_vars[:, j+1] == A @ x_vars[:, j] + B @ u_vars[:, j]]
        constraints += [cp.abs(x_vars[:, j]) <= x_max]
        constraints += [cp.abs(u_vars[:, j]) <= u_max]
        cost += cp.quad_form(x_vars[:, j] - xs_a, Q) + cp.quad_form(u_vars[:, j] - us_a, R)
    
    # 终端代价与偏移代价
    cost += cp.quad_form(x_vars[:, N + M_tilde] - xs_a, P)
    cost += cp.quad_form(xs_a - x_ref, S_offset)
    
    # 求解 QP
    prob = cp.Problem(cp.Minimize(cost), constraints)
    prob.solve(solver=cp.ECOS)
    
    if prob.status != cp.OPTIMAL and prob.status != cp.OPTIMAL_INACCURATE:
        print(f"步数 {k}: 求解失败 ({prob.status})。尝试减小初始状态或增大预测时域。")
        break
    
    # 执行第一步控制
    u_curr = u_vars.value[:, 0]
    u_history[:, k] = u_curr
    x_history[:, k+1] = A @ x_k + B @ u_curr
    
    if k % 5 == 0:
        print(f"步数 {k}: 状态 {x_k}, 辅助设定点 {xs_a.value}")

# ==========================================
# 4. 绘图结果
# ==========================================
plt.figure(figsize=(12, 5))

# 相平面图
plt.subplot(1, 2, 1)
plt.plot(x_history[0, :k+1], x_history[1, :k+1], 'b-o', markersize=4, label='Closed-loop State')
plt.plot(x_ref[0], x_ref[1], 'r*', markersize=12, label='True Target')
plt.axhline(0, color='k', linestyle=':', alpha=0.3)
plt.axvline(0, color='k', linestyle=':', alpha=0.3)
plt.title('State Space Trajectory (Tracking MPC)')
plt.xlabel('$x_1$')
plt.ylabel('$x_2$')
plt.legend()
plt.grid(True)

# 输入约束图
plt.subplot(1, 2, 2)
plt.step(range(k), u_history[0, :k], label='$u_1$')
plt.step(range(k), u_history[1, :k], label='$u_2$')
plt.axhline(u_max, color='r', linestyle='--', alpha=0.5, label='Constraint')
plt.axhline(-u_max, color='r', linestyle='--')
plt.title('Control Inputs (Saturated at $\pm 0.3$)')
plt.xlabel('Time Step')
plt.ylabel('Control Value')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()