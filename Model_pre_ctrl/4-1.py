import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
from scipy.linalg import solve_discrete_are

def get_lqr_params(A, B, Q, R):
    P = solve_discrete_are(A, B, Q, R)
    K = -np.linalg.inv(R + B.T @ P @ B) @ B.T @ P @ A
    return K, P

# 1. 参数设置
A = np.array([[1.1, 1], [0, 1.1]])
B = np.array([[1, 0.5], [1, 1]])
n, m = B.shape
Q, R, O_weight = 10 * np.eye(n), 100 * np.eye(m), 10000 * np.eye(n)
x_limit, u_limit, lambda_val = 3.0, 2.0, 0.9

K, P = get_lqr_params(A, B, Q, R)
M_u = np.linalg.pinv(B) @ (np.eye(n) - A)

# 2. 控制器参数
N, M_tilde = 10, 5
target_refs = [np.array([1.5, -1.0]), np.array([-1.5, -1.0]), np.array([2.5, 2.5])]

# 3. 仿真准备
x_curr = np.array([-1.0, 1.0])
history_x = [] # 改为在循环中记录，更直观
history_x.append(x_curr.copy())
history_xs_a = []

print("--- 正在运行 4.1 节学术示例仿真 ---")
for i, x_ref in enumerate(target_refs):
    print(f"正在前往参考点 {i+1}: {x_ref}")
    for k in range(40): # 每个点走40步
        u_v = cp.Variable((m, N + M_tilde))
        x_v = cp.Variable((n, N + M_tilde + 1))
        theta_a = cp.Variable(n) 
        
        xs_a = theta_a
        us_a = M_u @ theta_a
        
        cost = cp.quad_form(xs_a - x_ref, O_weight)
        constraints = [x_v[:, 0] == x_curr]
        
        for j in range(N + M_tilde):
            constraints += [x_v[:, j+1] == A @ x_v[:, j] + B @ u_v[:, j]]
            constraints += [cp.abs(x_v[:, j]) <= x_limit, cp.abs(u_v[:, j]) <= u_limit]
            constraints += [cp.abs(xs_a) <= lambda_val * x_limit, cp.abs(us_a) <= lambda_val * u_limit]
            
            if j >= N:
                constraints += [u_v[:, j] == K @ (x_v[:, j] - xs_a) + us_a]
            cost += cp.quad_form(x_v[:, j] - xs_a, Q) + cp.quad_form(u_v[:, j] - us_a, R)
            
        cost += cp.quad_form(x_v[:, N + M_tilde] - xs_a, P)
        cp.Problem(cp.Minimize(cost), constraints).solve(solver=cp.ECOS)
        
        x_curr = A @ x_curr + u_v.value[:, 0] @ B.T # 更新状态
        history_x.append(x_curr.copy())
        history_xs_a.append(xs_a.value.copy())

# 4. 绘图
history_x = np.array(history_x)
history_xs_a = np.array(history_xs_a)

plt.figure(figsize=(7, 7))
plt.plot(history_x[:, 0], history_x[:, 1], 'b-o', markersize=2, label='Actual State (history_x)')
plt.plot(history_xs_a[:, 0], history_xs_a[:, 1], 'r--', label='Artificial Setpoint')
for r in target_refs: plt.plot(r[0], r[1], 'go', markersize=8)
plt.gca().add_patch(plt.Rectangle((-3,-3),6,6, fill=False, color='k', linestyle='--'))
plt.title("Academic Example Reproduction")
plt.xlabel("x1"); plt.ylabel("x2"); plt.legend(); plt.grid(True)
plt.show()
# print("仿真结束:", history_x)
print("仿真结束，最后的坐标位置为:", history_x[-1])
