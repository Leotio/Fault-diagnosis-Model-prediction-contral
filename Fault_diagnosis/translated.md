### 3. Problem formulation（问题描述）

#### 3.1. System models and SVOs（系统模型与集合值观测器）

考虑如下受乘性执行器故障影响的离散时间线性时不变（LTI）系统：

$x_{k+1} = Ax_k + BG_i u_k + E\omega_k, \quad (1a)$

$y_k = Cx_k + F\eta_k, \quad (1b)$

其中 $A, B, C, E$ 和 $F$ 是参数矩阵。

$u_k$ 是受凸多胞形集合 $\mathcal{U}$ 约束的输入。

$\omega_k$ 是未知的输入向量，如扰动和建模误差；$\eta_k$ 代表测量噪声。

$\omega_k$ 和 $\eta_k$ 分别包含在 Zonotopes $\mathcal{W} = \langle \omega_c, H_\omega \rangle$ 和 $\mathcal{V} = \langle \eta_c, H_\eta \rangle$ 中，其中 $\omega_c$ 和 $\eta_c$ 是中心，$H_\omega$ 和 $H_\eta$ 是具有适当维数的生成矩阵。

$G_i$ 模拟不同的执行器模式。

$G_0$ 是单位矩阵，模拟健康模式。

$G_i = \text{diag}(1, \dots, 1, f_i, 1, \dots, 1)$ 模拟发生在第 $i$ 个执行器上的故障，其第 $i$ 个分量满足 $0 \le f_i < 1$。

在此，我们允许故障幅值 $f_i$ 在系统演化过程中在给定的区间 $[\underline{f}_i, \bar{f}_i]$ 内变化。



假设 3.1：对于系统 (1)，(1) 中的矩阵对 $(A, C)$ 是可检测的。矩阵 $A$ 是 Schur 矩阵（特征值均在单位圆内），且矩阵 $C$ 具有满行秩。

假设 3.2：系统 (1) 在诊断过程中运行在一种特定的模式下。



假设 3.1 是开环系统稳定性和 SVO（集合值观测器）存在的通常要求。

假设 3.2 意味着在执行诊断算法时系统模式保持不变，这也是成功进行故障诊断所必需的。

根据 (1)，匹配第 $i$ 种执行器模式的 SVO 设计为：

$\hat{X}_{k+1}^i = (A - L_k^i C)\hat{X}_k^i \oplus L_k^i y_k \oplus (-L_k^i F\mathcal{V}) \oplus B G_i u_k \oplus E\mathcal{W}, \quad (2a)$

$\hat{Y}_{k+1}^i = C\hat{X}_{k+1}^i \oplus F\mathcal{V}, \quad (2b)$

其中 $L_k^i$ 是第 $i$ 个 SVO 的增益。

$G_0$ 模拟健康模式。$G_i$ 是模拟第 $i$ 个执行器故障的区间矩阵。

$\hat{X}_{k+1}^i$ 和 $\hat{Y}_{k+1}^i$ 分别是第 $(k+1)$ 个时刻对应于第 $i$ 种系统模式的状态估计集合和输出估计集合。

当系统运行在第 $i$ 种模式且给定初始状态 $x_0 \in \hat{X}_0^i$ 时，包含关系 $x_{k+1} \in \hat{X}_{k+1}^i$ 和 $y_{k+1} \in \hat{Y}_{k+1}^i$ 对所有 $k \ge 0$ 均成立。



引理 3.1：对于每个区间矩阵 $G_i$，$B G_i u_k$ 被包含在一个 Zonotope 中，即（公式略）。

基于第 2 节中的 Zonotope 运算和引理 3.1，动态方程 (2) 可以改写为中心-生成矩阵形式：

$\hat{x}_{k+1}^{i,c} = (A - L_k^i C)\hat{x}_k^{i,c} + L_k^i(y_k - F\eta_c) + \text{mid}(BG_i)u_k + E\omega_c, \quad (3a)$

$\hat{H}_{k+1}^{i,x} = [(A - L_k^i C)\hat{H}_k^{i,x}, -L_k^i F H_\eta, \text{diag}(\text{rad}(BG_i)u_k), E H_\omega], \quad (3b)$

$\hat{y}_{k+1}^{i,c} = C\hat{x}_{k+1}^{i,c} + F\eta_c, \quad (3c)$

$\hat{H}_{k+1}^{i,y} = [C\hat{H}_{k+1}^{i,x}, F H_\eta], \quad (3d)$

其中 $\hat{x}, \hat{H}$ 分别是估计集合 $\hat{X}$ 和 $\hat{Y}$ 的中心和生成矩阵。

应当提到，随着系统的演化，$\hat{X}_{k+1}^i$ 的阶次会迅速增加。为了避免涉及高阶 Zonotopes 的计算负担，应当利用第 2 节中提到的阶次削减技术来控制阶次。

------

#### 3.2. Asymptotic AFD framework（渐近主动故障诊断框架）

故障诊断的原理是检查 $y_{k+1} \notin \hat{Y}_{k+1}^i$ （4）是否成立。

不一致性 (4) 表明系统不在第 $i$ 种模式下。因此，如果最终只存在一个一致的 OES（输出估计集合），则可以诊断出故障。

Xu (2021) 在每个时间步计算输入和一组时变 SVO 用于故障诊断。输入经过优化以最大化 OES 之间的中心距离。

受 Zonotopic 卡尔曼滤波思想的启发，观测器增益被设计为最小化每个 OES 的 Frobenius 半径（F-radius）。优化问题表述为：

$u_k^* = \arg \max \sum \|\hat{y}_{k+1}^{i,c} - \hat{y}_{k+1}^{j,c}\|_2^2, \quad (5a)$

（5b）

其中 $u_k^*$ 和 $L_{k,KF}^i$ 是在第 $k$ 个时间步设计的观测器最优输入和增益，$\lambda_{ij}$ 是权重系数。

核心思想是扩大不同 OES 之间的分离趋势，从而减小它们的交集区域。因此，系统输出更有可能停留在一致集合的非交叠部分，从而指示出故障。设计逻辑展示见图 1。

![image-20260103151249698](C:\Users\Leoti\AppData\Roaming\Typora\typora-user-images\image-20260103151249698.png)

由于诊断是渐近实现的，该方法被称为渐近主动故障诊断 (AAFD)。

备注 3.1：优化问题 (5) 是对 Xu (2021) 的简化。事实上，在原始问题中输入和增益是耦合的，Xu (2021) 构建了一个双层框架来求解。本文展示了其核心逻辑，并以此为基础提出了新的观测器增益设计框架。

然而，根据 (4)，基于集合的 AFD 原则是确保系统输出尽快离开不一致的集合。集合分离要求是实现该目标的一种有效但间接的方式。在第 4 节中，我们基于原则 (4) 推导出一种新的观测器增益优化目标，展示出更高的故障诊断潜力。

### 4 Main Result（主要成果）

#### 4.1. Optimal observer gains for fault diagnosis（用于故障诊断的最优观测器增益）

给定 Zonotope 表示形式 $\hat{Y}_{k+1}^i = \langle \hat{y}_{k+1}^{i,c}, \hat{H}_{k+1}^{i,y} \rangle$，

Scott 等人 (2014) 通过以下优化问题来检查 $y_{k+1}$ 与 $\hat{Y}_{k+1}^i$ 之间的关系：

$\hat{\delta}_{k+1}^i := \min_{\delta_{k+1}^i, \xi} \delta_{k+1}^i, \quad (6a)$

s.t. $y_{k+1} = \hat{y}_{k+1}^{i,c} + \hat{H}_{k+1}^{i,y} \xi, \quad (6b)$

$\|\xi\|_\infty \le \delta_{k+1}^i. \quad (6c)$

最优目标值 $\hat{\delta}_{k+1}^i$ 是 $\hat{Y}_{k+1}^i$ 围绕其中心 $\hat{y}_{k+1}^{i,c}$ 的缩放因子，使得 $y_{k+1}$ 恰好位于缩放后的 Zonotope 表面。

如果 $\hat{\delta}_{k+1}^i \le 1$，则包含关系 $y_{k+1} \in \hat{Y}_{k+1}^i$ 成立。

关于 (6) 的图解说明在图 2 中给出。

由于 $y_{k+1}^I \in \hat{Y}_{k+1}^i$，$\hat{Y}_{k+1}^i$ 需要围绕 $\hat{y}_{k+1}^{i,c}$ 收缩，导致 $\hat{\delta}_{k+1}^i < 1$。

由于 $y_{k+1}^E \notin \hat{Y}_{k+1}^i$，$\hat{Y}_{k+1}^i$ 必须放大自身，从而 $\hat{\delta}_{k+1}^i > 1$。

此外，$\hat{\delta}_{k+1}^i = 1$ 作为包含与排除之间的精确边界。

我们在定义 4.1 中将 $\hat{\delta}_{k+1}^i$ 正式定义为排除趋势。

#### 定义 4.1

问题 (6) 的最优目标值 $\hat{\delta}_{k+1}^i$ 衡量了 $y_{k+1}$ 在输出估计 Zonotope $\hat{Y}_{k+1}^i$ 中的排除趋势。

基于公式 (4) 和上述分析，系统输出在每个不一致集合中的排除趋势从 $\hat{\delta} < 1$ 增加到 $\hat{\delta} > 1$ 是基于集合的 AFD（主动故障诊断）方法的最终目标。

根据集合传播动力学，$L_k^i$ 可以改变 $\hat{Y}_{k+1}^i$ 的中心 $\hat{y}_{k+1}^{i,c}$ 和生成矩阵 $\hat{H}_{k+1}^{i,y}$

从而影响排除趋势 $\hat{\delta}_{k+1}^i$。在本文中，我们建议优化 $L_k^i$ 以使 $\hat{\delta}_{k+1}^i$ 最大化。该优化问题表述如下：

$$\hat{\delta}_{k+1}^{i,*} := \max_{L_k^i} \hat{\delta}_{k+1}^i = \max_{L_k^i} \min_{\delta_{k+1}^i, \xi} \delta_{k+1}^i. \quad (7)$$

(7) 中的最优故障诊断观测器增益记为 $L_{k,FD}^i$。

我们以图 3 为例来论证 $L_{k,FD}^i$ 的重要性。对于不一致的 OES（输出估计集合），使 $\hat{\delta}_{k+1}^i$ 最大化的 $L_{k,FD}^i$ 可以增强 $y_{k+1}$ 从 $\hat{Y}_{k+1}^i$ 中的排除效果，这对应于图 3 的底部。由于唯一的一致集合是 $\hat{Y}_{k+1}^2$，当前的故障被确定在第二个执行器。图 3 的顶部显示了由 (5b) 设计的 $L_{k,KF}^i$ 所获得的 OES。尽管每个 $\hat{Y}_{k+1}^i$ 的尺寸被最小化了，但 $y_{k+1}$ 仍包含在交叠区域内，导致我们无法确定当前的系统模式。

需要指出的是，我们所提的方法在第 $(k+1)$ 个时间步设计 $L_{k,FD}^i$，因为 $y_{k+1}$ 应当作为先验信息已知。由于 $u_k$ 已经通过使用 Xu (2021) 中的方法确定，并在第 $k$ 个时间步注入系统以获取 $y_{k+1}$，因此为了简洁，基于排除趋势的框架不再重复输入设计过程，而是直接使用与 AAFD 框架相同的输入设计方法（即公式 (5a)）。

我们所提的主动故障诊断（AFD）框架程序总结如下：

- 在第 $k$ 个时间步，注入 $u_k^*$ 以最大化第 $(k+1)$ 个时间步各 OES 之间的分离趋势；
- 在第 $(k+1)$ 个时间步，每个 OES 通过 $L_{k,FD}^i$ 进一步发生形变，从而将系统输出 $y_{k+1}$ 从 $\hat{Y}_{k+1}^i$ 中排除。

总的来说，$u_k^*$ 激励系统以产生更多有用的输出信息用于故障推理，而 $L_{k,FD}^i$ 则修正结果以增强系统的可诊断性。两种设计原则协同工作，构成了一种强大的 AFD 范式。

#### 4.2. Dual transformation of the optimization problem（优化问题的对偶变换）

max–min 问题 (7) 难以直接求解。由于线性问题的强对偶性，(7) 可以被转化为单层最大化（max）问题。

命题 4.1：对于系统 (1)，双层 max–min 问题 (7) 等价于如下单层最大化问题：

$$\hat{\delta}_{k+1}^{i,*} = \max_{L_k^i, \lambda_{n,k}^i, v_k^i} (v_k^i)^T C L_k^i (y_k - \hat{y}_{k,c}^i) + (v_k^i)^T c_{i,\Delta}^k \quad (8a)$$

s.t. $(\lambda_{1,k}^i)^T - (\lambda_{2,k}^i)^T = -(v_k^i)^T C(A - L_k^i C) \hat{H}_{k,x}^i, \quad (8b)$

$(\lambda_{3,k}^i)^T - (\lambda_{4,k}^i)^T = (v_k^i)^T C L_k^i F H_\eta, \quad (8c)$

$(\lambda_{5,k}^i)^T - (\lambda_{6,k}^i)^T = -(v_k^i)^T H_{i,\Delta}^k, \quad (8d)$

$\lambda_{n,k}^i \ge 0, n = 1, \dots, 6, \quad (8e)$

$\sum_{n=1}^6 (\lambda_{n,k}^i)^T \mathbf{1} = 1, \forall i \in \mathcal{I}_{n_u}, \quad (8f)$

其中 $v_k^i$ 和 $\lambda_{n,k}^i (n=1, \dots, 6)$ 是由问题 (6) 的对偶变换导出的拉格朗日乘子，$c_{i,\Delta}^k = C(A \hat{x}_{k,c}^i + \text{mid}(BG_i)u_k + E\omega_c) + F\eta_c - y_{k+1}$，且 $H_{i,\Delta}^k = [C \text{diag}(\text{rad}(BG_i)u_k), CEH_\omega, FH_\eta]$。

证明：基于 (3) 中 $\hat{y}_{k+1}^{i,c}$ 和 $\hat{H}_{k+1}^{i,y}$ 的表达式，内部问题 (6) 被重新表述为：

$\min_{\delta_{k+1}^i, \xi} \delta_{k+1}^i \quad (9)$

s.t. $0 = C L_k^i (y_k - \hat{y}_{k,c}^i) + c_{i,\Delta}^k + C(A - L_k^i C) \hat{H}_{k,x}^i \xi_1 - C L_k^i F H_\eta \xi_2 + H_{i,\Delta}^k \xi_3,$

$\xi_1 - \delta_{k+1}^i \mathbf{1} \le 0, -\xi_1 - \delta_{k+1}^i \mathbf{1} \le 0,$

$\xi_2 - \delta_{k+1}^i \mathbf{1} \le 0, -\xi_2 - \delta_{k+1}^i \mathbf{1} \le 0,$

$\xi_3 - \delta_{k+1}^i \mathbf{1} \le 0, -\xi_3 - \delta_{k+1}^i \mathbf{1} \le 0,$

$\xi = [\xi_1^T, \xi_2^T, \xi_3^T]^T.$

(9) 的拉格朗日函数为：

$\mathcal{L}(v_k^i, \lambda_{n,k}^i, \delta_{k+1}^i, \xi) = (v_k^i)^T C L_k^i (y_k - \hat{y}_{k,c}^i) + (v_k^i)^T c_{i,\Delta}^k + (1 - \sum_{n=1}^6 (\lambda_{n,k}^i)^T \mathbf{1})\delta_{k+1}^i + ((\lambda_{1,k}^i)^T - (\lambda_{2,k}^i)^T + (v_k^i)^T C(A - L_k^i C) \hat{H}_{k,x}^i)\xi_1 + ((\lambda_{3,k}^i)^T - (\lambda_{4,k}^i)^T - (v_k^i)^T C L_k^i F H_\eta)\xi_2 + ((\lambda_{5,k}^i)^T - (\lambda_{6,k}^i)^T + (v_k^i)^T H_{i,\Delta}^k)\xi_3, \quad (10)$

其中 $v_k^i$ 和 $\lambda_{n,k}^i (n=1, \dots, 6)$ 分别是与第 $k$ 个时间步第 $i$ 个 SVO 的等式及不等式约束相关联的拉格朗日乘子。根据 (10)，拉格朗日对偶问题进一步推导为：

$\max_{\lambda_{n,k}^i, v_k^i} (v_k^i)^T C L_k^i (y_k - \hat{y}_{k,c}^i) + (v_k^i)^T c_{i,\Delta}^k \quad (11a)$

s.t. (8b) – (8f), $\forall i \in \mathcal{I}_{n_u}. \quad (11b)$

显而易见，内部问题 (6) 是线性的。因此，根据 Boyd 和 Vandenberghe (2004)，问题 (6) 和 (11) 的最优目标值是相同的。基于表述 (7)，通过将 (11) 中的 $L_k^i$ 设置为新的优化变量，双层 max–min 问题可以转化为单层最大化问题。最终形式如 (8) 所示。

**备注 4.1**：事实上，优化问题 (8) 的构型并不局限于执行器故障。一个直观的推广是将 $A$ 和 $C$ 替换为 $A_i$ 和 $C_i$，其中 $A_i$ 和 $C_i$ 分别模拟第 $i$ 种系统模式下的系统故障和传感器故障。然而，Xu (2021) 中的 AAFD 方法仅考虑了执行器故障。为了与前作进行公平比较，本研究仅评估了执行器故障。

我们可以通过求解 (8) 获得最优故障诊断观测器增益 $L_{k,FD}^i$。然而，在 (8) 中，优化变量 $v_k^i$ 和 $L_k^i$ 是相乘的，这使得原始问题在典型的凸优化框架下无法求解。

#### 4.3. Problem solutions（问题的求解）

在命题 4.2 中，我们证明了 (8) 可以被重新表述为一个线性问题。在此之前，需要引理 4.1。

引理 4.1：考虑线性方程组 $AX = Y, \quad (12)$，其中 $A \in \mathbb{R}^{m \times n}, X \in \mathbb{R}^{n \times r}, Y \in \mathbb{R}^{m \times r}$。如果 $\text{rank}(A) = m$，则 (12) 的通解为 $X^* = A^\dagger Y + (I_n - A^\dagger A)R, \quad (13)$，其中 $A^\dagger = A^T(AA^T)^{-1}$ 是 $A$ 的伪逆，$R \in \mathbb{R}^{n \times r}$ 是任意矩阵。

命题 4.2：考虑如下线性优化问题：

$\delta_{k+1}^{i,*} := \max_{\gamma_k^i, \lambda_{n,k}^i, v_k^i} (\gamma_k^i)^T(y_k - \hat{y}_{k,c}^i) + (v_k^i)^T c_{i,\Delta}^k \quad (14a)$

s.t. $(\lambda_{1,k}^i)^T - (\lambda_{2,k}^i)^T = -(v_k^i)^T C A \hat{H}_{k,x}^i + (\gamma_k^i)^T C \hat{H}_{k,x}^i, \quad (14b)$

$(\lambda_{3,k}^i)^T - (\lambda_{4,k}^i)^T = (\gamma_k^i)^T F H_\eta, \quad (14c)$

$(8d), (8e), (8f), \forall i \in \mathcal{I}_{n_u}. \quad (14d)$

如果 (14) 的最优解 $v_{k}^{i,*}$ 不为 0，则 $\delta_{k+1}^{i,*} = \hat{\delta}_{k+1}^{i,*}$，即 (14) 和 (8) 的最优目标值相同。

$L_{k,0}^i = \frac{C^T v_k^{i,*}(\gamma_k^{i,*})^T}{(v_k^{i,*})^T C C^T v_k^{i,*}}$, $L_{k,1}^i = I_{n_x} - \frac{C^T v_k^{i,*} (v_k^{i,*})^T C}{(v_k^{i,*})^T C C^T v_k^{i,*}}$.

利用最优解 $\gamma_k^{i,*}$ 和 $v_k^{i,*}$，问题 (8) 中的最优故障诊断观测器增益 $L_{k,FD}^i$ 具有如下形式：$L_{k,FD}^i = L_{k,0}^i + L_{k,1}^i R_k^i, \quad (15)$，其中 $R_k^i \in \mathbb{R}^{n_x \times n_y}$ 是任意矩阵，$L_{k,0}^i$ 和 $L_{k,1}^i$ 的定义如上。



证明：显而易见，(8) 等价于线性问题 (14) 附加一个约束条件 $(v_k^i)^T C L_k^i = (\gamma_k^i)^T$。

由于引入了新的优化变量 $\gamma_k^i$，关系 $\delta_{k+1}^{i,*} \ge \hat{\delta}_{k+1}^{i,*}$ 成立。



在求解 (14) 后，如果存在 $L_{k,FD}^i$ 使得 $(v_k^{i,*})^T C L_{k,FD}^i = (\gamma_k^{i,*})^T$ 成立，则等式 $\delta_{k+1}^{i,*} = \hat{\delta}_{k+1}^{i,*}$ 满足。

引理 4.1 表明，如果 $\text{rank}((v_k^{i,*})^T C) = 1$，则 $L_{k,FD}^i$ 具有通解形式 (15)。

由于假设 3.1 规定 $C$ 具有满行秩，因此若 $v_k^{i,*} \ne 0$，则 $\text{rank}((v_k^{i,*})^T C) = 1$ 成立。

综上所述，若 $v_k^{i,*} \ne 0$，则 $\delta_{k+1}^{i,*} = \hat{\delta}_{k+1}^{i,*}$，且问题 (8) 的最优故障诊断观测器增益 $L_{k,FD}^i$ 具有通用形式 (15)。

由命题 4.2 可知，(8) 的最优解 $L_{k,FD}^i$ 具有通式 (15)，其中 $v_k^{i,*}$ 和 $\gamma_k^{i,*}$ 通过求解线性问题 (14) 获得，而 $R_k^i$ 是待确定的自由变量。

原则上，观测器的一个核心目的是产生准确的估计结果。基于参数形式 (15)，可以进一步优化 $R_k^i$ 来提高估计精度。

#### 4.4. Improvements in estimation performance（估计性能的提升）

为了平衡故障诊断和状态估计的性能，我们提出了三种设计原则，它们在优化观测器增益方面各具优势。

(I) 最小化 Frobenius 半径 (F-radius)：F-radius 是表征 Zonotope 尺寸的一种度量。我们可以优化 $R_k^i$ 来最小化状态估计集合的 F-radius，以获得更准确的结果。根据 (3b)，优化问题给定为：

$$\min_{R_k^i} \| \hat{X}_{k+1}^i \|_F^2 = \text{tr} ( (\hat{H}_{k+1}^{i,x})^T \hat{H}_{k+1}^{i,x} ), \forall i \in \mathcal{I}_{n_u}, \quad (16)$$



其中 $\hat{H}_{k+1}^{i,x} = [(A - L_{k,FD}^i C)\hat{H}_k^{i,x}, - L_{k,FD}^i F H_\eta, H_\Xi]$，$H_\Xi = [\text{diag}(\text{rad}(BG_i)u_k), E H_\omega]$，且 $L_{k,FD}^i$ 具有形式 (15)。根据矩阵微积分的结果，$\| \hat{X}_{k+1}^i \|_F^2$ 对 $R_k^i$ 的导数为：

$$\frac{\partial \| \hat{X}_{k+1}^i \|_F^2}{\partial R_k^i} = 2((L_{k,1}^i)^T L_{k,1}^i) R_k^i P - 2Q, \quad (17)$$

其中 $P = F H_\eta H_\eta^T F^T + C \hat{H}_{k,x}^i (\hat{H}_{k,x}^i)^T C^T$，$Q = (L_{k,1}^i)^T (A - L_{k,0}^i C) \hat{H}_{k,x}^i (\hat{H}_{k,x}^i)^T C^T - (L_{k,1}^i)^T L_{k,0}^i F H_\eta H_\eta^T F^T$。由于 (16) 是一个无约束凸二次规划问题，通过令 (17) 为零可得最优 $R_{k,*}^i$：

$$R_{k,*}^i = ((L_{k,1}^i)^T L_{k,1}^i)^\dagger Q P^{-1}. \quad (18)$$

**(II) 状态估计动力学的稳定性**：观测器增益的一个关键设计原则是确保估计误差的收敛。如果矩阵 $(A - L_k^i C)$ 的所有特征值模长在每个时间步 $k$ 都小于 1，则称该 SVO 是鲁棒稳定的。确保动力学稳定的 $L_{k,FD}^i$ 设计也可以通过如下命题中的 Lyapunov 稳定性来解决。

命题 4.3：考虑如下优化问题：

$$\min_{R_k^i, \alpha} \alpha \quad (19a)$$

$$\text{s.t.} \begin{bmatrix} -W & (A - L_{k,FD}^i C)^T W \\ \star & -W \end{bmatrix} \preceq \alpha I_{2n_x}, \forall i \in \mathcal{I}_{n_u}, \quad (19b)$$



其中 $W \succ 0$ 是正定矩阵，$\star$ 表示对称项，$L_{k,FD}^i$ 具有形式 (15)。如果最优解 $\alpha^* < 0$，则 SVO 动力学在第 $k$ 个时间步是稳定的。

**证明**：选择 Lyapunov 函数 $V_k = x_k^T W x_k$。对于动力学 $x_{k+1} = (A - L_{k,FD}^i C)x_k$，若 $(A - L_{k,FD}^i C)^T W (A - L_{k,FD}^i C) - W \prec 0$，则满足 $V_{k+1} < V_k$，通过 Schur 补性质这等价于 (19b) 的左侧 $\prec 0$。如果 $\alpha^* < 0$，显然 Lyapunov 稳定性得到保证。通过求解优化问题 (19)，若 $\alpha^* < 0$，我们可以导出一个稳定 SVO 动力学的最优 $R_{k,*}^i$。需要提到的是，$\alpha^*$ 并不保证一定为负。尽管如此，优化 (19) 仍具有抑制显著发散的潜力。

(III) 稳定性保证：最鲁棒的方法是直接求解优化问题 (8)，并添加 Lyapunov 稳定性条件作为新约束，即：

$$\hat{\delta}_{k+1}^{i,L} = \max_{L_k^i, \lambda_{n,k}^i, v_k^i} (v_k^i)^T C L_k^i (y_k - \hat{y}_{k,c}^i) + (v_k^i)^T c_{i,\Delta}^k \quad (20a)$$

$$\text{s.t. } (8b) - (8f), \quad (20b)$$

$$\begin{bmatrix} -W & (A - L_k^i C)^T W \\ \star & -W \end{bmatrix} \prec 0, \forall i \in \mathcal{I}_{n_u}. \quad (20c)$$

通过求解 (20) 得到的最优观测器增益保证了稳定的状态估计动力学。由于变量 $v_k^i$ 和 $L_k^i$ 在约束 (20c) 中是解耦的，我们不能像命题 4.2 那样通过代换双线性项将其转化为线性版本。求解 (20) 的一种有效方法是分枝定界法 (BNB)。其总体思路是在每次迭代中将搜索空间分枝成更小的区域，直到检测到最优解。在每个区域中，最优目标值 $\hat{\delta}_{k+1}^{i,L}$ 的下界可由局部求解器给出，而上界通过将双线性项替换为其凸松弛（McCormick 松弛）并求解松弛后的问题来获得。上下界帮助搜索算法排除不包含最优解的区域。当全局上下界之间的间隙低于给定阈值时，BNB 算法终止。软件工具箱 YALMIP 提供了一个高效的求解器来解决 (20)，这将在第 5 节中展示。

设计原则 (I) 可以给出最优 $R_{k,*}^i$ 的解析表达式。设计原则 (II) 可以产生抑制 SVO 动力学不稳定的观测器增益。设计原则 (III) 保证了诊断过程中的稳定动力学。然而，BNB 方法需要大量的计算时间来实现。特定设计原则的选择取决于应用场景，如稳定性要求、可用的计算资源等。在第 5 节中，我们给出了三种设计原则的数值性能比较。为了避免符号滥用，我们仍将通过上述原则之一进一步优化后的最终观测器增益表示为 $L_{k,FD}^i$。

#### 4.5. Analysis of the singular solution（奇异解分析）

在求解 (14) 后，如果得到最优解 $v_k^{i,*} = 0$，那么 $(v_k^{i,*})^T C L_k^i$ 始终等于 0，导致无法求解出最优的 $L_{k,FD}^i$。我们将 $v_k^{i,*} = 0$ 称为**奇异解**。在本小节中，对奇异解的性质进行了分析。

**定理 4.1**：对于优化问题 (14)，如果 $v_k^{i,*} = 0$，则 $\delta_{k+1}^{i,*} = \hat{\delta}_k^i$。

证明：将 $v_k^i = 0$ 代入 (14)，我们有：

$$\delta_{k+1}^{i,*} = \max_{\gamma_k^i, \lambda_{n,k}^i} (\gamma_k^i)^T (y_k - \hat{y}_{k,c}^i) \quad (21a)$$



s.t. $(\lambda_{1,k}^i)^T - (\lambda_{2,k}^i)^T = (\gamma_k^i)^T C \hat{H}_{k,x}^i, \quad (21b)$

$(\lambda_{3,k}^i)^T - (\lambda_{4,k}^i)^T = (\gamma_k^i)^T F H_\eta, \quad (21c)$

$\lambda_{5,k}^i = \lambda_{6,k}^i, \lambda_{n,k}^i \ge 0, n = 1, \dots, 6, \quad (21d)$

$\sum_{n=1}^6 (\lambda_{n,k}^i)^T \mathbf{1} = 1. \quad (21e)$

再次使用对偶变换，(21) 的拉格朗日对偶问题为：

$$\min_{\delta_k^i, \xi} \delta_k^i \quad (22a)$$

s.t. $y_k = \hat{y}_{k,c}^i + \hat{H}_{k,y}^i \xi, \quad \|\xi\|_\infty \le \delta_k^i. \quad (22b)$

优化问题 (22) 计算的是 $y_k$ 在 $\hat{Y}_k^i$ 中的排除趋势，即 $\hat{\delta}_k^i$。结果表明 $\delta_{k+1}^{i,*} = \hat{\delta}_k^i$。 □

由定理 4.1 可以推导得出，对于任意 $L_k^i \in \mathbb{R}^{n_x \times n_y}$，都有 $\hat{\delta}_{k+1}^i \le \delta_{k+1}^{i,*} = \hat{\delta}_k^i$。因此，通过 $L_k^i$ 的设计，第 $(k+1)$ 个时间步的排除趋势无法超过 $\hat{\delta}_k^i$。在这种情况下，无法计算出最优故障诊断观测器增益 $L_{k,FD}^i$。尽管如此，后续的主动故障诊断（AFD）仍需要一个合适的观测器增益。在我们的框架中，我们选择采用公式 (5b) 中的 $L_{k,KF}^i$ 作为一种折中，当奇异解发生时，我们提出的方法将退化为 Xu (2021) 中的 AAFD 方法。一个值得阐述的重要问题是，当 $\delta_{k+1}^{i,*} \ge 1$ 成立时，将不再需要计算 $L_{k,FD}^i$ 来诊断故障，从而节省了计算资源。

**命题 4.4**：对于线性问题 (14)，如果最优目标值 $\delta_{k+1}^{i,*} \ge 1$，则存在 $L_k^i \in \mathbb{R}^{n_x \times n_y}$ 使得 $y_{k+1} \notin \hat{Y}_{k+1}^i$。

**证明**：首先，我们证明如果 $\delta_{k+1}^{i,*} > 1$，则 $v_k^{i,*} \ne 0$，即当前步骤不会发生奇异解。如果 $v_k^{i,*} = 0$，由定理 4.1 可知 $\delta_{k+1}^{i,*} = \hat{\delta}_k^i > 1$ 成立，因此 $y_k \notin \hat{Y}_k^i$。在这种情况下，模式 $i$ 应该在前面的步骤中就已经从当前的备选系统模式中被移除了，这引发了矛盾。由于 $v_k^{i,*} \ne 0$，命题 4.2 表明 $\delta_{k+1}^{i,*} = \hat{\delta}_{k+1}^{i,*} > 1$ 成立。因此，存在某个 $L_k^i$ 确保最大排除趋势 $\hat{\delta}_{k+1}^{i,*}$ 超过 1。在这种情况下，计算 $L_{k,FD}^i$ 是不必要的。 

**备注 4.2**：对于实际系统模式 $i$ 的 SVO，$y_k \in \hat{Y}_k^i$ 始终成立，因此对于所有 $k \ge 0$，排除趋势 $\hat{\delta}_k^i < 1$ 成立。在这种情况下，持续增加 $\hat{\delta}^i$ 是不切实际的。因此，奇异解在实践中并非罕见案例。然而，根据命题 4.4，当对于某个 $L_k^i$ 成立 $y_{k+1} \notin \hat{Y}_{k+1}^i$ 时，不会发生奇异解来干扰诊断。此外，由于 OES（输出估计集合）不断自我修正以排除输出，不一致的 SVO 对故障变得越来越敏感。尽管存在奇异解，我们提出的工作仍然比 Xu (2021) 中的 AAFD 方法具有更高的故障诊断潜力。

我们所提框架的完整设计流程总结为**算法 1**。如果最终索引集 $\mathcal{I}_{n_u}$ 中仅剩一个元素，则诊断成功。



------

### 5. 算例演示

在本节中，我们在四水箱系统上将所提方法与 Xu (2021) 的方法进行对比。系统在平衡点 $\bar{x} = [12.4, 12.7, 1.8, 1.4]^T$ 附近线性化。采用采样时间 $T_s = 1$ s 的欧拉离散化方法获得 LTI 系统形式 (1)。为了更好地进行比较，我们假设所有系统状态都可以通过适当的传感器测量，因此对矩阵 $C$ 进行了微调。系统的参数矩阵给定如下：

（公式 23 略，包含 A, B, E, C, F 矩阵的具体数值）

系统配备了两个执行器，因此考虑三个系统模式，即 $\mathcal{I}_{n_u} = \{0, 1, 2\}$。输入约束集为 $\mathcal{U} = \langle 0, 0.6I_2 \rangle$。初始条件和不确定性边界给定如下：

$x_0 = [0, 0, 0, 0]^T, \hat{X}_0^i = \langle x_0, 0.1I_4 \rangle (i \in \mathcal{I}_{n_u}),$

$\mathcal{W} = \langle 0, 0.1I_7 \rangle, \mathcal{V} = \langle 0, 0.1I_7 \rangle.$

公式 (19) 中的正定矩阵选择为 $W = I_{n_x}$。为了平衡效率和复杂度，$\hat{X}_k^i$ 的阶次保持为 60，$\hat{Y}_k^i$ 的阶次为 10。设置算法 1 中的最大诊断步数 $k_{max} = 10$。执行器故障矩阵 $G_i$ 为：

$\text{mid}(G_1) = \text{diag}([f_{1c}, 1]), \text{mid}(G_2) = \text{diag}([1, f_{2c}]),$

$\text{rad}(G_1) = \text{diag}([0.05, 0]), \text{rad}(G_2) = \text{diag}([0, 0.05]),$

其中 $f_{1c}$ 和 $f_{2c}$ 是故障区间的中心。故障区间的半径为 0.05。所有程序均在配备单线程执行的 Intel(R) Core(TM) i7-9700 CPU @ 3.00 GHz 机器上使用 MATLAB 2020a 进行评估。除了设计原则 (III) 涉及的 BNB 算法（使用 YALMIP 工具箱的 BMIBNB 求解器）外，其余优化问题均使用 CVX 工具箱求解。

以 $(f_{1c}, f_{2c}) = (0.8, 0.85)$ 为例。假设第一个执行器受到恒定故障 $f_1 = 0.8$ 的干扰。我们使用所提方法逐步设计输入和观测器增益来诊断故障。此外，同时执行设计原则 (I) 和 (II) 以展示其诊断性能。结果如图 4 所示。可以看出，对于两种设计原则，均观察到 $y_3 \in \hat{Y}_3^1$ 而 $y_3 \notin \hat{Y}_3^0$ 且 $y_3 \notin \hat{Y}_3^2$，表明在 $k = 3$ 时识别出第一个执行器故障。作为对比，我们采用 Xu (2021) 中的 AAFD 方法诊断故障，结果如图 5 所示。在 $k = 3$ 时，$y_3 \in \hat{Y}_3^0 \cap \hat{Y}_3^1 \cap \hat{Y}_3^2$，因此在 $k = 3$ 时无法识别系统模式。进一步仿真显示，AAFD 方法在 $k_{max} = 10$ 个演化步内未能完成诊断。对比表明，在本算例中，所提方法优于 AAFD 方法。

为了全面对比所提方法与 Xu (2021) 中的 AAFD 方法，我们在该系统上测试了不同故障设置下的诊断性能。固定两个执行器故障的半径为 0.05。我们选择不同的故障中心 $(f_{1c}, f_{2c}) \in \mathcal{F} \times \mathcal{F}$，其中 $\mathcal{F} = \{0.1, 0.2, \dots, 0.9\}$。因此，总共考虑了 $9 \times 9$ 种不同的故障设置。在每种故障设置下，分别运行 AAFD 方法、采用设计原则 (I) 的所提方法以及采用设计原则 (II) 的所提方法，在 $k_{max} = 10$ 个时间步内诊断故障。为了减少不确定性的影响，在每种设置下进行了 100 次测试，其中 $\omega_k, \eta_k$ 和 $f_i$ 从给定的有界集合中随机抽取。图 6 显示了注入 $f_1$ 时 100 次测试中成功案例的数量，图 7 展示了注入 $f_2$ 时的结果。数值结果表明：

- 我们的方法比 AAFD 方法能在更多情况下实现故障诊断，这可以通过图 6 和图 7 中的彩色圆点直观地看到。
- 设计原则 (II) 的诊断性能略优于设计原则 (I)，但代价是计算复杂度增加。

综上所述，在诊断该系统中发生的执行器故障时，我们的方法比 AAFD 方法更有效。同时，由于设计原则 (I) 能以比设计原则 (II) 更短的平均计算时间提供优异的性能，因此更受青睐。

在最后部分，我们使用一个数值例子来说明设计原则 (III) 的性能。考虑故障设置 $(f_{1c}, f_{2c}) = (0.5, 0.6)$ 且发生 $f_2$ 故障的一个例子，我们首先使用设计原则 (II) 诊断故障。故障可以在 $k = 8$ 时间步被识别。然而，最优故障诊断观测器增益 $L_{3,FD}^2$ 和 $L_{5,FD}^2$ 不能确保稳定的估计动力学，其特征值分别为 $\lambda(A - L_{3,FD}^2 C) = [-1.0174, 0.9502, 0, 0.9615]$ 和 $\lambda(A - L_{5,FD}^2 C) = [-1.1730, 0.9502, 0, 0.9615]$。为了保证稳定的估计动力学，执行了设计原则 (III)。BMIBNB 求解器在每个时间步都成功解决了问题 (20)，并在 $k = 4$ 时间步诊断出了故障。然而，求解时间从 2.1 s 到 85 s 不等，取决于 BNB 算法中所需的计算迭代次数。高昂的计算时间使得设计原则 (III) 不适合在线故障诊断。尽管如此，BNB 算法完成了保证稳定估计动力学和增强系统可诊断性的任务。