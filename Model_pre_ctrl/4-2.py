import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
from scipy.linalg import block_diag, solve_discrete_are
import time 

def get_msd_system(N_carts, dt=0.02):
    m, k_s, d_c = 1.0, 2.0, 1.5
    A_ii = np.array([[0, 1], [-2*k_s/m, -d_c/m]])
    B_ii = np.array([[0], [1/m]])
    A_ij = np.array([[0, 0], [k_s/m, 0]])
    
    A_cont = np.zeros((2*N_carts, 2*N_carts))
    B_cont = np.zeros((2*N_carts, N_carts))
    
    for i in range(N_carts):
        A_cont[2*i:2*i+2, 2*i:2*i+2] = A_ii
        B_cont[2*i:2*i+2, i] = B_ii.flatten()
        if i > 0: A_cont[2*i:2*i+2, 2*(i-1):2*(i-1)+2] = A_ij
        if i < N_carts - 1: A_cont[2*i:2*i+2, 2*(i+1):2*(i+1)+2] = A_ij
            
    Ad = np.eye(2*N_carts) + A_cont * dt
    Bd = B_cont * dt
    return Ad, Bd

def run_large_scale_carts(N_carts=50):
    print(f"\n--- 正在运行 4.2 节 (N_carts: {N_carts}, 状态维度: {2*N_carts}) ---")
    
    # 离线阶段
    start_offline = time.time()
    Ad, Bd = get_msd_system(N_carts)
    n, m = Ad.shape[0], Bd.shape[1]
    
    Qi = np.diag([1.0, 1.5])
    Q = block_diag(*([Qi] * N_carts))
    R = 20 * np.eye(m)
    O_weight = 1e5 * np.eye(n)
    
    print("正在离线计算 LQR 增益与稳态映射...")
    P = solve_discrete_are(Ad, Bd, Q, R)
    K = -np.linalg.inv(R + Bd.T @ P @ Bd) @ Bd.T @ P @ Ad
    M_u = np.linalg.pinv(Bd) @ (np.eye(n) - Ad)
    offline_time = time.time() - start_offline

    # MPC 参数
    N, M_tilde = 5, 20
    x_curr = np.zeros(n)
    x_ref = np.zeros(n)
    for i in range(N_carts): x_ref[2*i] = 1.0
        
    history_r = []
    avg_pos = 0
    step_times = []
    k = 0

    print(f"开始在线仿真 (阈值: 0.99)...")
    total_online_start = time.time()
    
    while abs(avg_pos - 1.00) >= 0.01: 
        step_start = time.time()
        
        u_v = cp.Variable((m, N + M_tilde))
        x_v = cp.Variable((n, N + M_tilde + 1))
        theta_a = cp.Variable(n)
        
        xs_a = theta_a
        us_a = M_u @ theta_a
        
        cost = cp.quad_form(xs_a - x_ref, O_weight)
        constraints = [x_v[:, 0] == x_curr]
        
        for j in range(N + M_tilde):
            constraints += [x_v[:, j+1] == Ad @ x_v[:, j] + Bd @ u_v[:, j]]
            constraints += [cp.abs(u_v[:, j]) <= 5.0]
            if j >= N:
                constraints += [u_v[:, j] == K @ (x_v[:, j] - xs_a) + us_a]
            cost += cp.quad_form(x_v[:, j] - xs_a, Q) + cp.quad_form(u_v[:, j] - us_a, R)
        
        prob = cp.Problem(cp.Minimize(cost), constraints)
        prob.solve(solver=cp.OSQP, verbose=False)
        
        if u_v.value is None:
            print("求解失败")
            break
            
        x_curr = Ad @ x_curr + Bd @ u_v.value[:, 0]
        history_r.append(x_curr[::2].copy())
        avg_pos = np.mean(x_curr[::2])
        
        step_end = time.time()
        step_times.append(step_end - step_start)
        
        if k % 50 == 0:
            print(f"Step {k}, 平均位移: {avg_pos:.4f}, 当前步耗时: {step_times[-1]:.4f}s")
        k += 1

    total_online_time = time.time() - total_online_start

    # 输出统计 
    print("\n" + "="*30)
    print(f"仿真完成统计 (N_carts = {N_carts})")
    print(f"离线设计总耗时: {offline_time:.4f} s")
    print(f"在线仿真总步数: {k}")
    print(f"在线仿真总耗时: {total_online_time:.4f} s")
    print(f"在线平均单步时间: {np.mean(step_times):.4f} s")
    print("="*30)

    # 绘图
    history_r = np.array(history_r)
    plt.figure(figsize=(10, 6))
    for i in range(N_carts):
        plt.plot(history_r[:, i], alpha=0.3)
    plt.axhline(1.0, color='r', linestyle='--', label='Target')
    plt.title(f"Large-Scale System: {N_carts} Carts ({2*N_carts} States)")
    plt.xlabel("Steps"); plt.ylabel("Displacement r_i"); plt.grid(True); plt.show()

if __name__ == "__main__":
    run_large_scale_carts(N_carts=50)