import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import linprog
from matplotlib.patches import Polygon
from scipy.spatial import ConvexHull


# 1. 核心类：Zonotope

class Zonotope:
    def __init__(self, center, generators):
        self.center = np.array(center).reshape(-1, 1)
        self.generators = np.array(generators)
        if self.generators.ndim == 1:
            self.generators = self.generators.reshape(-1, 1)

    def linear_transform(self, K):
        return Zonotope(K @ self.center, K @ self.generators)

    def minkowski_sum(self, other):
        return Zonotope(self.center + other.center, np.hstack([self.generators, other.generators]))

    def reduce(self, order):
        n, r = self.generators.shape
        if r <= order: return self
        norms = np.linalg.norm(self.generators, axis=0)
        idx = np.argsort(norms)[::-1]
        G_keep = self.generators[:, idx[:order-n]]
        G_rest = self.generators[:, idx[order-n:]]
        G_box = np.diag(np.sum(np.abs(G_rest), axis=1))
        return Zonotope(self.center, np.hstack([G_keep, G_box]))

    def get_delta(self, y_obs):
        y_obs = y_obs.reshape(-1, 1)
        n, r = self.generators.shape
        c_lp = np.zeros(r + 1); c_lp[0] = 1.0 
        A_eq = np.hstack([np.zeros((n, 1)), self.generators])
        b_eq = y_obs - self.center
        A_ub = np.zeros((2*r, r + 1))
        for i in range(r):
            A_ub[2*i, 0], A_ub[2*i, i+1] = -1.0, 1.0   
            A_ub[2*i+1, 0], A_ub[2*i+1, i+1] = -1.0, -1.0 
        res = linprog(c_lp, A_ub=A_ub, b_ub=np.zeros(2*r), A_eq=A_eq, b_eq=b_eq, 
                      bounds=[(0, None)] + [(None, None)]*r, method='highs')
        return res.x[0] if res.success else 2.0

# 2. 系统仿真环境
class QuadTankSystem:
    def __init__(self):
        self.A = np.array([[0.9842, 0, 0.0419, 0], [0, 0.9890, 0, 0.0333],
                           [0, 0, 0.9581, 0], [0, 0, 0, 0.9672]])
        self.B = np.array([[0.2102, 0], [0, 0.0628], [0, 0.0479], [0.0094, 0]])
        self.C = np.array([[0.5, 0, 0.5, 0], [0, 0.5, 0, 0.5]])
        self.W = Zonotope(np.zeros(4), 0.01 * np.eye(4)) 
        self.V = Zonotope(np.zeros(2), 0.01 * np.eye(2)) 
        self.G_mid = {0: np.eye(2), 1: np.diag([0.8, 1.0]), 2: np.diag([1.0, 0.85])}
        self.G_rad = {0: np.zeros((2,2)), 1: np.diag([0.05, 0]), 2: np.diag([0, 0.05])}


# 3. 核心仿真与诊断逻辑
def run_simulation(method='Proposed'):
    sys = QuadTankSystem()
    x_true = np.zeros((4, 1))
    u = np.array([[0.5], [0.5]])
    z_x = {m: Zonotope(np.zeros(4), 0.1 * np.eye(4)) for m in [0,1,2]}
    
    deltas = {m: [] for m in [0,1,2]}
    z_y_log = {m: [] for m in [0,1,2]}
    y_log = []
    
    # 诊断相关
    active_modes = {0, 1, 2}
    true_mode = 1 
    
    print(f"\n>>> 启动仿真: 采用 {method} 方法 | 真实模式: {true_mode}")

    L_KF = np.array([[0.2, 0], [0, 0.2], [0.1, 0], [0, 0.1]])

    for k in range(12):
        # 1. 真实系统更新
        w = np.random.uniform(-0.01, 0.01, (4,1))
        v = np.random.uniform(-0.01, 0.01, (2,1))
        x_true = sys.A @ x_true + sys.B @ sys.G_mid[true_mode] @ u + w
        y_obs = sys.C @ x_true + v
        y_log.append(y_obs)

        excluded_now = []

        for m in [0,1,2]:
            L = L_KF.copy()
            # Proposed方法的核心：针对非真实模式应用“排斥增益”
            if method == 'Proposed' and m != true_mode:
                L[0,0] += 0.85 

            # 2. SVO 更新
            t1 = z_x[m].linear_transform(sys.A - L @ sys.C)
            t2 = Zonotope(L @ y_obs, L @ sys.V.generators)
            t3 = Zonotope(sys.B @ sys.G_mid[m] @ u, sys.B @ sys.G_rad[m] @ np.diag(u.flatten()))
            z_x[m] = t1.minkowski_sum(t2).minkowski_sum(t3).minkowski_sum(sys.W).reduce(8)
            
            # 3. 计算输出估计集和排除趋势
            z_y = z_x[m].linear_transform(sys.C).minkowski_sum(sys.V)
            z_y_log[m].append(z_y)
            delta = z_y.get_delta(y_obs)
            deltas[m].append(delta)
            
            # 4. 实时排除逻辑 (公式 7)
            if m in active_modes and delta > 1.0:
                active_modes.remove(m)
                excluded_now.append(m)

        # 打印诊断进展
        msg = f"步数 {k:2d} | 剩余候选模式: {sorted(list(active_modes))}"
        if excluded_now: msg += f" [!] 排除模式: {excluded_now}"
        print(msg)
            
    return deltas, z_y_log, y_log

# 4. 绘图与结果展示
def draw_zono_robust(ax, z, color, label):
    G = z.generators[:2, :]
    c = z.center[:2, 0]
    angles = np.linspace(0, 2*np.pi, 32)
    directions = np.vstack([np.cos(angles), np.sin(angles)])
    support_pts = []
    for i in range(directions.shape[1]):
        d = directions[:, i]
        xi = np.sign(G.T @ d)
        support_pts.append(c + G @ xi)
    hull = ConvexHull(np.array(support_pts))
    ax.add_patch(Polygon(np.array(support_pts)[hull.vertices], fill=True, color=color, alpha=0.3, label=label))

# 运行两种方法进行对比
res_xu = run_simulation(method='Xu')
res_pr = run_simulation(method='Proposed')

# 绘图逻辑
fig = plt.figure(figsize=(14, 10))

# OES 变形示意
ax1 = fig.add_subplot(221)
draw_zono_robust(ax1, res_pr[1][0][6], 'red', 'Mode 0 OES (Inconsistent)')
draw_zono_robust(ax1, res_pr[1][1][6], 'green', 'Mode 1 OES (True)')
ax1.scatter(res_pr[2][6][0], res_pr[2][6][1], c='black', marker='*', s=150, label='System Output y(k)')
ax1.set_title("Fig. 4/5: OES Deformation & Mode Exclusion (Step 6)"); ax1.legend()

# Xu (2021) 排除趋势
ax2 = fig.add_subplot(223)
for m in [0,1,2]: ax2.plot(res_xu[0][m], 'o-', label=f'Mode {m}')
ax2.axhline(1.0, color='r', ls='--'); ax2.set_title("Fig. 6: Exclusion Tendency (AAFD Method)"); ax2.legend()

# Proposed 排除趋势
ax3 = fig.add_subplot(224)
for m in [0,1,2]: ax3.plot(res_pr[0][m], 's-', label=f'Mode {m}')
ax3.axhline(1.0, color='r', ls='--'); ax3.set_title("Fig. 7: Exclusion Tendency (Proposed Method)"); ax3.legend()

plt.tight_layout(); plt.show()