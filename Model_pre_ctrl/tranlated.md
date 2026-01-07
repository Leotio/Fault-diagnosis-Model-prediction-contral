**摘要**

本文提出了一种具有隐式终端组件的跟踪模型预测控制（MPC）技术。该控制器方案将辅助设定点（artificial setpoint）作为决策变量，并为依赖于该设定点的增广系统隐式地定义终端约束。在这方面，我们不再计算不变终端集，而是考虑一个扩展的预测时域，其长度可以简单地通过求解线性规划（LP）来确定界限。这种方法克服了计算不变集所需操作相关的规模限制，同时也简化了离线 MPC 的设计。所提出的控制器能够将大型系统驱动至可容许的设定点，同时保证递归可行性和收敛性。最后，通过一个学术示例、一个规模可变的质量-弹簧-阻尼系统以及一个更具现实意义的无人机案例研究对该方法进行了说明。

**1. 引言**

模型预测控制（MPC）的理论属性（如递归可行性和稳定性）通常是通过精心设计终端成本函数和终端正不变集来保证的。然而，随着系统规模的增大，这些集合的计算变得愈发困难（Gilbert & Tan, 1991; Mayne, 2013, 2014）。特别是，即使对于线性时不变系统，由于需要计算集合交集和原像集，显式确定不变集也是具有挑战性的。事实上，目前存在一些可扩展的方法来寻找这些不变集的近似值，例如：通过考虑预定义的多面体形状（Trodden, 2016）、基于线性矩阵不等式（LMI）的椭球体（Alamo, Cepeda, & Limon, 2005）、内外近似（Comelli, Olaru, & Kofman, 2024）、分形区（zonotopes）（Morato, Cunha, Santos, Normey-Rico, & Sename, 2021）或数据驱动方法（Berberich, Köhler, Müller, & Allgöwer, 2021），或者生成它们的隐式表示（Raković & Zhang, 2022, 2023; Wang & Jungers, 2020）。在这些方法中，我们对后者（Raković & Zhang, 2022, 2023）特别感兴趣，因为它通过扩展预测时域来保证最终状态属于一个不变集。

此外，一些 MPC 公式进一步使这一问题复杂化，例如跟踪 MPC（Ferramosca, Limón, Alvarado, Alamo, & Camacho, 2009; Limón, Alvarado, Alamo, & Camacho, 2008），它考虑了一个增广终端系统——扩大了其维度——以处理非固定设定点。具体而言，Limón 等人（2008）和 Ferramosca 等人（2009）提出的跟踪 MPC 公式在优化问题中加入了辅助稳态和输入作为决策变量，以放松终端约束，并通过使用修改后的成本函数确保递归可行性和对真实设定点的收敛性。虽然用于跟踪变化设定点的 MPC 控制器家族非常广泛（例如 Bemporad, Casavola, & Mosca, 1997; Garone, Di Cairano, & Kolmanovsky, 2017; Gilbert & Kolmanovsky, 2002），但我们考虑采用 Ferramosca 等人（2009）和 Limón 等人（2008）的方法来应用所提出的隐式方法，这不仅是因为其理论保证，还因为它扩大了控制器的吸引域。

本文的主要贡献在于设计了一种结合了隐式终端成分的跟踪 MPC，扩展了 Luque, Chanfreut, Limón, 和 Maestre (2024) 引入的初步工作。特别是，所提出的方法依靠用一个预定义有限长度的扩展预测时域来取代终端约束集。所提出的方法允许通过求解线性规划（LP）来获得跟踪设置下的这一长度，从而避免了集合计算，因此能够应用于任何规模的系统，代价是由于使用更长的时域而导致在线计算负担略有增加。正如将要看到的，在这种情况下使用隐式不变集并不是显而易见的，需要针对辅助设定点和真实设定点的存在，对调节问题的结果进行量身定制的改进。

所提出的方法提供了递归可行性和对真实设定点收敛的保证，这是跟踪 MPC 使用中的关键属性。此外，本文还讨论了一种针对扩展预测时域的最大允许长度受限情况的替代方法，这需要替代策略来避免终端区域的计算，同时能够选择扩展预测时域所需的延伸长度。

本文其余部分的结构如下。第 2 节介绍了所提出问题的预备知识。第 3 节介绍了具有隐式终端组件的控制器设计，以及替代方法和理论证明，其性能在第 4 节通过三个案例研究进行了说明。最后，第 5 节提供了讨论。

**符号说明**。向量 $[x^\top, u^\top]^\top$ 记为 $(x, u)$；$I_n$ 和 $0_{m \times n}$ 分别代表维度为 $n \times n$ 和 $m \times n$ 的单位矩阵和零矩阵，而 $0_n$ 和 $1_n$ 是大小为 $n \times 1$ 的零列向量和全一列向量。$\mathbb{R}, \mathbb{N}$ 分别表示实数集和自然数集。同样地，给定 $a, b \in \mathbb{N}$ 且 $a < b$，我们定义 $\mathbb{N}_{[a, b]} := \{a, a+1, \dots, b-1, b\}$，且 $\mathbb{N}_b$ 由 $\mathbb{N}_{[0, b]}$ 给出。最后，$[u_t]_{t=0}^T$ 表示对于任何给定的 $T \in \mathbb{N}$ 的向量 $[u_0^\top, u_1^\top, \dots, u_T^\top]^\top$。对于 $\mathbb{R}^n$ 中封闭非空子集 $X$ 的支撑函数 $h(X, \cdot)$，对于所有 $y \in \mathbb{R}^n$，由 $h(X, y) := \sup_x \{y^\top x : x \in X\}$ 给出。



**2. 问题设置**

在本节中，将介绍系统动力学以及所考虑的跟踪公式下的 MPC。

**2.1. 系统动力学**

考虑如下离散时间线性系统：



$$x_{k+1} = Ax_k + Bu_k, \quad (1)$$



其中 $x_k \in \mathbb{R}^n$ 和 $u_k \in \mathbb{R}^m$ 分别是系统在时刻 $k$ 的状态和输入，矩阵 $A$ 和 $B$ 具有兼容的维度，即 $A \in \mathbb{R}^{n \times n}$，$B \in \mathbb{R}^{n \times m}$，$m$ 和 $n$ 是正整数且可能不相等。同时，考虑以下约束：



$$x_k \in \mathcal{X}, \quad u_k \in \mathcal{U}, \quad \forall k \in \mathbb{N}, \quad (2)$$



其中 $\mathcal{X} \subseteq \mathbb{R}^n$ 和 $\mathcal{U} \subseteq \mathbb{R}^m$ 分别是状态和输入约束集。我们引入以下假设：

**假设 1**。对于受约束 (2) 限制的系统 (1)，以下条件成立：

- 矩阵对 $(A, B)$ 是已知且严格可稳定的。

- 约束集 $\mathcal{X}$ 和 $\mathcal{U}$ 是包含原点在其内部的凸多面体集。

- 存在一个反馈增益 $K \in \mathbb{R}^{m \times n}$ 使得矩阵 $A + BK$ 是 Schur 稳定的，且存在一个正定矩阵 $P \in \mathbb{R}^{n \times n}$ 满足：

  

  $$(A + BK)^\top P(A + BK) - P = -(Q + K^\top RK). \quad (3)$$

系统性能将通过阶段成本函数进行评估：



$$\ell(x_k, u_k, x_s, u_s) = \|x_k - x_s\|^2_Q + \|u_k - u_s\|^2_R, \quad (4)$$



其中 $Q \in \mathbb{R}^{n \times n}$ 和 $R \in \mathbb{R}^{m \times m}$ 是对称正定矩阵，$(x_s, u_s)$ 表示我们希望系统驱动至的设定点。在这方面，请注意系统的任何设定点必须满足以下方程：



$$[A - I_n \quad B] \begin{bmatrix} x_s \\ u_s \end{bmatrix} = 0_n. \quad (5)$$



因此，我们可以通过变量 $\theta \in \mathbb{R}^m$ 对对 $(x_s, u_s)$ 进行参数化，即：



$$\begin{bmatrix} x_s \\ u_s \end{bmatrix} = \underbrace{\begin{bmatrix} M_\theta^x \\ M_\theta^u \end{bmatrix}}_{M_\theta} \theta, \quad (6)$$



其中 $M_\theta$ 是 $[A - I_n \quad B]$ 零空间的合适基矩阵，由矩阵 $M_\theta^x \in \mathbb{R}^{n \times m}$ 和 $M_\theta^u \in \mathbb{R}^{m \times m}$ 组合而成（Ferramosca et al., 2009）。

**2.2. 具有显式终端组件的跟踪 MPC**

本文考虑的跟踪 MPC 公式具有以下特征（Limón et al., 2008）：

(i) 引入辅助设定点（artificial setpoint），记作 $(x_s^a, u_s^a)$，作为 MPC 问题中的优化变量。该辅助设定点由变量 $\theta^a$ 参数化，并引入了 $m$ 个新的优化变量。

(ii) 增加了一个偏移成本函数（offset cost function），用于惩罚辅助设定点与真实设定点 $(x_s^r, u_s^r)$ 之间的偏差。

(iii) 使用增广终端不变集 $\Psi_f^{tr}$，该集合针对增广系统定义：



$$\begin{bmatrix} x_{k+1} \\ \theta^a \end{bmatrix} = \underbrace{\begin{bmatrix} A + BK & BL \\ 0_{m \times n} & I_m \end{bmatrix}}_{A_{aug}} \begin{bmatrix} x_k \\ \theta^a \end{bmatrix}, \quad (7)$$



其中 $L = [-K \quad I_m]M_\theta$。在下文中，我们将区分真实设定点 $(x_s^r, u_s^r)$ 和辅助设定点 $(x_s^a, u_s^a)$。注意 (7) 考虑了系统 (1) 采用控制律：



$$u_k = u_s^a + K(x_k - x_s^a) = Kx_k + L\theta^a. \quad (8)$$

考虑到上述情况，在每个时间时刻 $k$ 求解的跟踪 MPC 问题采用如下形式：



$$V_N^*(x_k, x_s^r) = \min_{\mathbf{u}, \theta^a} V_N(x_k, x_s^r, \mathbf{u}, \theta^a) \quad (9a)$$

$$\text{s.t. } x_{0|k} = x_k, \quad (9b)$$

$$x_{j+1|k} = Ax_{j|k} + Bu_{j|k}, \quad j \in \mathbb{N}_{[0, N-1]}, \quad (9c)$$

$$x_{j+1|k} \in \mathcal{X}, \quad j \in \mathbb{N}_{[0, N-1]}, \quad (9d)$$

$$u_{j|k} \in \mathcal{U}, \quad j \in \mathbb{N}_{[0, N-1]}, \quad (9e)$$

$$\begin{bmatrix} x_s^a \\ u_s^a \end{bmatrix} = M_\theta \theta^a, \quad \begin{bmatrix} x_{N|k} \\ \theta^a \end{bmatrix} \in \Psi_f^{tr}, \quad (9f)$$



其中 $N$ 是预测时域的长度，$\mathbf{u} = [u_{j|k}]_{j=0}^{N-1}$。

在这方面，下标 $j|k$ 表示在时刻 $k$ 对时刻 $k+j$ 相应变量的预测。此外，成本函数定义为：



$$V_N(x_k, x_s^r, \mathbf{u}, \theta^a) = \sum_{j=0}^{N-1} (\|x_{j|k} - x_s^a\|^2_Q + \|u_{j|k} - u_s^a\|^2_R) + \|x_{N|k} - x_s^a\|^2_P + \|x_s^a - x_s^r\|^2_O, \quad (10)$$



其中 $O \in \mathbb{R}^{n \times n}$ 是正定矩阵，$P$ 满足 (3)。

注意，与调节 MPC 不同，预测时域内加权的是系统相对于辅助设定点的偏差，并增加了偏移成本 $\|x_s^a - x_s^r\|^2_O$ 来惩罚辅助状态参考值与真实状态参考值之间的差异。同样，终端约束 (9f) 由不变集 $\Psi_f^{tr}$ 定义，该集合是在考虑约束 (2) 的情况下针对增广系统 (7) 计算的。特别地，$\Psi_f^{tr}$ 是最大不变集的超平面近似，满足：



$$\Psi_f^{tr} \subseteq \{(x, \theta) \in \mathbb{R}^{n+m} : (x, Kx + L\theta) \in (\mathcal{X}, \mathcal{U}), \theta \in \Theta\},$$



其中



$$\Theta := \{\theta \in \mathbb{R}^m : M_\theta^x \theta \in \mathcal{X}, M_\theta^u \theta \in \mathcal{U}\}. \quad (11)$$



由于 $A_{aug}$ 具有单位特征值，集合 $\Psi_f^{tr}$ 可能不是有限确定的（Gilbert & Tan, 1991）。尽管如此，可以通过因子 $\lambda \in (0, 1)$ 缩放 $\Theta$，使得最大可容许不变集成为一个有限确定的凸多面体，记为 $\Psi_{f, \lambda}^{tr}$ (Gilbert & Tan, 1991; Limón et al., 2008)。

注意：



$$\Psi_{f, \lambda}^{tr} \subseteq \{(x, \theta) \in \mathbb{R}^{n+m} : (x, Kx + L\theta) \in (\mathcal{X}, \mathcal{U}), \theta \in \lambda\Theta\}, \quad (12)$$



且对于所有 $\begin{bmatrix} x_k \\ \theta^a \end{bmatrix} \in \Psi_{f, \lambda}^{tr}$，有 $\begin{bmatrix} Ax_k + B(Kx_k + L\theta^a) \\ \theta^a \end{bmatrix} \in \Psi_{f, \lambda}^{tr}$。

正如 Limón 等人 (2008) 和 Ferramosca 等人 (2009) 所详述的，跟踪公式 (9) 允许将系统状态驱动到任何可容许的目标设定点。然而，它需要计算增广动力学 (7) 的不变集 $\Psi_{f, \lambda}^{tr}$，其维度为 $n+m$。虽然终端成本很容易计算，但对于大型系统，构建这样的终端集变得难以处理。为了解决这个问题，本文使用隐式终端组件重新表述了 MPC 问题 (9)，即避免了显式表征 $\Psi_{f, \lambda}^{tr}$ 的需要。

**3. 具有隐式终端组件的跟踪 MPC**

在下文中，我们将 Raković 和 Zhang (2023) 中针对调节问题导出的隐式终端组件结果扩展到跟踪问题。

**3.1. 调节问题的隐式终端集**

对于调节问题，即 $(x_s^r, u_s^r) = 0_{n+m}$，终端控制律可以简单地定义为 $u_k = Kx_k$，因此终端动力学由下式给出：



$$x_{k+1} = (A + BK)x_k, \quad (13)$$



如果在 (7) 和 (8) 中将辅助设定点和真实设定点都固定在原点，也可以得到该方程。同时，给定 (2)，终端阶段的约束可以紧凑地定义为：



$$\mathcal{X}_t := \{x \in \mathbb{R}^n : x \in \mathcal{X}, Kx \in \mathcal{U}\}. \quad (14)$$



考虑到上述内容，让我们引入以下定理，它确立了最大正不变集存在的充分条件。回想一下，集合 $\Omega \subset \mathbb{R}^n$ 被定义为约束 $x \in \mathcal{X}$ 的正不变集，当且仅当 (Kerrigan, 2001; Raković & Zhang, 2022)：



$$x_k \in \Omega \Rightarrow \exists u_k \in \mathcal{U} \text{ 使得 } x_{k+1} \in \Omega, x_{k+1} \in \mathcal{X}.$$



如果一个集合是正不变的且包含 $\Omega$ 中所有的正不变集，则该集合被定义为最大正不变集 (Kerrigan, 2001)。

定理 1 (Raković 和 Zhang (2023, 定理 1 和 2))。

假设假设 1 成立。那么，系统 (13) 和约束 (14) 的最大正不变集是有限确定的，当且仅当对于某个 $M \in \mathbb{N}$，以下条件之一成立：



$$\bigcap_{j=0}^{M} (A + BK)^{-j} \mathcal{X}_t \subseteq (A + BK)^{-(M+1)} \mathcal{X}_t \quad \text{或} \quad \mathcal{X}_t \subseteq (A + BK)^{-(M+1)} \mathcal{X}_t. \quad (15)$$

根据定理 1，系统 (13) 和约束 (14) 的最大正不变集（记为 $\Psi_f$）是一个非空封闭多面体集，其内部包含原点，定义为：



$$\Psi_f = \bigcap_{j=0}^{M} (A + BK)^{-j} \mathcal{X}_t. \quad (16)$$

可以使用长度为 $M$ 的轨迹来引入另一种检查给定状态是否属于 $\Psi_f$ 的方法 (Raković & Zhang, 2023)。也就是说，$x_N \in \Psi_f$ 当且仅当序列 $[x_j]_{j=N}^{N+M}$ 满足：



$$\forall j \in \mathbb{N}_{[N, N+M]}, x_j \in \mathcal{X}_t, \quad \forall j \in \mathbb{N}_{[N, N+M-1]}, x_{j+1} = (A + BK)x_j, \quad (17)$$



其中 $M, N \in \mathbb{N}$，且 $M \ge 1$ 满足 (15)。注意如果 $M = 0$，则 $\Psi_f = \mathcal{X}$。后者作为终端组件隐式重构的基础。

设集合 $\mathcal{X}_t \subseteq \mathbb{R}^n$（见 (14)）是一个封闭多面体集，其不可约表示为：



$$\mathcal{X}_t := \{x \in \mathbb{R}^n : (C + DK)x \le \mathbf{1}_p\}, \quad (18)$$



其中矩阵对 $(C, D) \in \mathbb{R}^{p \times n} \times \mathbb{R}^{p \times m}$ 经相应定义。注意 $p$ 表示定义 $\mathcal{X}_t$ 的不等式数量。

给定 (18)，$\mathcal{X}_t$ 也可以类似地定义为：



$$\mathcal{X}_t = \{x \in \mathbb{R}^n : \forall i \in \mathbb{N}_{[1, p]}, X_i^\top x \le 1\}, \quad (19)$$



其中 $X_i^\top$ 是矩阵 $(C + DK)$ 的第 $i$ 行。根据 Raković 和 Zhang (2023)，对于所有 $i \in \mathbb{N}_{[1, p]}$，(15) 成立当且仅当以下之一成立：



$$h \left( \bigcap_{j=0}^{M} (A + BK)^{-j} \mathcal{X}_t, ((A + BK)^{M+1})^\top X_i \right) \le 1 \quad \text{或} \quad h (\mathcal{X}_t, ((A + BK)^{M+1})^\top X_i) \le 1. \quad (20)$$

(20) 中不等式的左侧可以分别为不同的 $i$ 通过以下计算简单的线性规划（LP）来求得：



$$\sup_{(x, [x_j]_{j=0}^M)} \{X_i^\top (A + BK)^{M+1} x : x_j = (A + BK)^j x, x_j \in \mathcal{X}_t, \forall j \in \mathbb{N}_M\}$$



或



$$\sup_x \{X_i^\top (A + BK)^{M+1} x : x \in \mathcal{X}_t\}.$$



因此，给定任何整数 $M$，验证 (15) 简化为求解 $p$ 个 LP。此外，可以通过搜索任何适当生成的正递增整数序列来寻找该整数 $M$。详见 Raković 和 Zhang (2022, 2023)。

**3.2. 跟踪问题的隐式终端集**

假设系统 (1) 受控制律 (8) 控制。那么，考虑到增广动力学，对于给定的 $\theta$，以下条件成立：



$$\begin{bmatrix} x_{k+1} \\ \theta \end{bmatrix} = A_{aug} \begin{bmatrix} x_k \\ \theta \end{bmatrix}, \quad (21)$$



其中该终端增广动力学的约束定义如下：



$$\mathcal{X}_{aug,t} := \{(x, \theta) \in \mathbb{R}^{n+m} : x \in \mathcal{X}, Kx + L\theta \in \mathcal{U}, \theta \in \lambda\Theta\}. \quad (22)$$



类似于 (18)，设集合 $\mathcal{X}_{aug,t} \subseteq \mathbb{R}^{n+m}$ 为一个封闭多面体，其不可约定义为：



$$\mathcal{X}_{aug,t} := \left\{(x, \theta) \in \mathbb{R}^{n+m} : \tilde{G} \begin{bmatrix} x \\ \theta \end{bmatrix} \le \mathbf{1}_{\tilde{p}}\right\}, \quad (23)$$



其中



$$\tilde{G} = \begin{bmatrix} C + DK & DL \\ 0 & \tilde{W} \end{bmatrix} \in \mathbb{R}^{\tilde{p} \times (n+m)}. \quad (24)$$



且 $\tilde{W} = (CM_\theta^x + DM_\theta^u)/\lambda$，$\tilde{p}$ 定义为 $\mathcal{X}_{aug,t}$ 中不等式的数量。给定 (23)，集合 $\mathcal{X}_{aug,t}$ 也可以定义为：



$$\mathcal{X}_{aug,t} = \left\{(x, \theta) \in \mathbb{R}^{n+m} : \forall i \in \mathbb{N}_{[1, \tilde{p}]}, X_{aug,i}^\top \begin{bmatrix} x \\ \theta \end{bmatrix} \le 1\right\}, \quad (25)$$



其中 $X_{aug,i}^\top$ 是矩阵 $\tilde{G}$ 的第 $i$ 行。此外，让我们考虑满足类似于 (15) 但适用于增广系统的条件的某个 $\tilde{M} \in \mathbb{N}$：



$$\bigcap_{j=0}^{\tilde{M}} A_{aug}^{-j} \mathcal{X}_{aug,t} \subseteq A_{aug}^{-(\tilde{M}+1)} \mathcal{X}_{aug,t} \quad \text{或} \quad \mathcal{X}_{aug,t} \subseteq A_{aug}^{-(\tilde{M}+1)} \mathcal{X}_{aug,t}, \quad (26)$$



以及以下终端集：



$$\Psi_{f, \lambda}^{tr} = \bigcap_{j=0}^{\tilde{M}} A_{aug}^{-j} \begin{bmatrix} \mathcal{X}_t \\ \lambda\Theta \end{bmatrix}. \quad (27)$$



由于对于某些 $\lambda \in (0, 1)$，跟踪问题的增广终端集 $\Psi_{f, \lambda}^{tr}$ 是有限确定的 (Gilbert & Tan, 1991; Limón et al., 2008)，因此可以确保存在一个满足方程 (26) 的有限值 $\tilde{M}$ (参见 Raković & Zhang, 2022, 推论 4)。

定理 2。

假设假设 1 成立，并考虑满足 (26) 的某些 $\tilde{M}, N \in \mathbb{N}$（且 $\tilde{M} \ge 1$）。那么，约束 $(x_N, \theta) \in \Psi_{f, \lambda}^{tr}$ 成立当且仅当存在一个序列 $[x_j]_{j=N}^{N+\tilde{M}}$ 满足：



$$\forall j \in \mathbb{N}_{[N, N+\tilde{M}]}, (x_j, \theta) \in \mathcal{X}_{aug,t}, \quad \forall j \in \mathbb{N}_{[N, N+\tilde{M}-1]}, x_{j+1} = (A + BK)x_j + BL\theta. \quad (28)$$

为了满足 (26)，对于所有 $i \in \mathbb{N}_{[1, \tilde{p}]}$，必须分别满足以下条件之一：



$$h \left( \bigcap_{j=0}^{\tilde{M}} A_{aug}^{-j} \mathcal{X}_{aug,t}, (A_{aug}^{\tilde{M}+1})^\top X_{aug,i} \right) \le 1 \quad \text{或} \quad h (\mathcal{X}_{aug,t}, (A_{aug}^{\tilde{M}+1})^\top X_{aug,i}) \le 1. \quad (29)$$



以下计算简单的 LP 允许我们检查任何 $i \in [1, \tilde{p}]$ 的不等式 (29)：



$$\sup_{(x, \theta, [x_j]_{j=0}^{\tilde{M}})} \left\{ X_{aug,i}^\top A_{aug}^{\tilde{M}+1} \begin{bmatrix} x \\ \theta \end{bmatrix} : \begin{bmatrix} x_j \\ \theta \end{bmatrix} = A_{aug}^j \begin{bmatrix} x \\ \theta \end{bmatrix} \in \mathcal{X}_{aug,t}, \forall j \in \mathbb{N}_{\tilde{M}} \right\}$$



或



$$\sup_{(x, \theta)} \left\{ X_{aug,i}^\top A_{aug}^{\tilde{M}+1} \begin{bmatrix} x \\ \theta \end{bmatrix} : (x, \theta) \in \mathcal{X}_{aug,t} \right\}.$$



因此，可以通过搜索正递增整数并为每个整数求解 $\tilde{p}$ 个 LP 来找到满足 (26) 的新参数 $\tilde{M}$。

**3.3. 具有隐式终端组件的跟踪 MPC**

本小节详述了所提出的包含跟踪和隐式终端约束的 MPC。根据定理 2，预测时域被划分为两个阶段：第一部分长度为 $N$，终端部分长度为 $\tilde{M}$。也就是说，将计算长度为 $N + \tilde{M}$ 的序列。

基于 (9)，所提出的每个步长求解的最优控制问题定义为：



$$\min_{\mathbf{u}, \theta^a} V_{N+\tilde{M}}(x_k, x_s^r, \mathbf{u}, \theta^a) \quad (30a)$$

$$\text{s.t. } x_{0|k} = x_k, \quad (30b)$$

$$x_{j+1|k} = Ax_{j|k} + Bu_{j|k}, \quad j \in \mathbb{N}_{[0, N+\tilde{M}-1]}, \quad (30c)$$

$$x_{j+1|k} \in \mathcal{X}, \quad j \in \mathbb{N}_{[0, N+\tilde{M}-1]}, \quad (30d)$$

$$u_{j|k} \in \mathcal{U}, \quad j \in \mathbb{N}_{[0, N+\tilde{M}-1]}, \quad (30e)$$

$$u_{j|k} = Kx_{j|k} + L\theta^a, \quad j \in \mathbb{N}_{[N, N+\tilde{M}-1]}, \quad (30f)$$

$$\begin{bmatrix} x_s^a \\ u_s^a \end{bmatrix} = M_\theta \theta^a, \quad (30g)$$

$$\theta^a \in \lambda\Theta. \quad (30h)$$



这里，约束 (30f) 施加了终端控制律，且仅在预测时域的第二部分（即终端阶段）考虑。此外，目标函数考虑了直到预测时间时刻 $N + \tilde{M}$ 的性能，加上终端成本和偏移成本，如 (10) 中所定义。最后，请注意扩展预测时域取代了显式终端约束 (9f)。

**备注 1**。隐式方法通过延长预测时域来避免离线计算不变集的显式表示，这意味着增加了在线计算负担。然而，由于二次规划（QP）可以在多项式时间内求解，这在大多数实际应用中可能不是限制因素。此外，一旦已知 $\tilde{M}$，也可以使用 (27) 找到闭式形式的最大不变集，从而避免在每次迭代时检查收敛条件。

备注 2。从 Raković 和 Zhang (2023) 可以推断出，跟踪问题的隐式终端成本函数可以定义如下：



$$V_{F, imp}^{tr}(x_{N|k}, x_s^a, u_s^a) = (1 - \epsilon)^{-1} \sum_{j=N}^{N+\tilde{M}-1} (\|x_{j|k} - x_s^a\|^2_Q + \|u_{j|k} - u_s^a\|^2_R),$$



其中最小化 $\epsilon \in [0, 1)$ 使得终端成本中的隐式界限变得紧致。虽然这是一个有效的选项，但按照 (10) 计算 $P$ 通常并不昂贵，且能提供更准确的剩余成本估计。因此，本文选择了后者。

**3.4. 理论属性**

在此，我们证明优化问题 (30) 的初始可行性也意味着递归可行性。此外，还证明了如果真实设定点是可容许的，系统将收敛到该设定点。

定理 3。

假设在时刻 $k$ 存在问题 (30) 的解 $(\mathbf{u}_k^*, \theta_k^*)$。那么，对于所有时刻 $t \ge k$，我们都能找到 (30) 的可行解。

（证明详情略，其核心思想是基于 $k$ 时刻的最优解构造 $k+1$ 时刻的候选解，并利用不变集的性质证明其可行性。）

**定理 4**。设 $x_0 \in \mathcal{X}_{N+\tilde{M}}$，其中 $\mathcal{X}_{N+\tilde{M}}$ 是控制器 (30) 的吸引域。那么，如果真实设定点 $x_s^r$ 是可容许的，受控系统状态 $x_k$ 当 $k$ 趋于无穷大时将收敛至 $x_s^r$。

（证明详情略，通过验证最优成本是闭环系统的 Lyapunov 函数，并利用反证法证明辅助设定点最终收敛到真实设定点。）

**3.5. 跟踪问题的替代隐式设计**

本节根据 Limón, Ferramosca, Alvarado, 和 Alamo (2018) 介绍了一种无需终端约束的稳定性预测控制器设计的替代方法。这种替代方案特别适用于那些无法或不希望将预测时域延长 $\tilde{M}$ 长度的情况。此外，对于此类系统，显式计算不变集也可能在计算上难以处理。

让我们考虑一个类似于 (10) 的终端成本函数，即 $V_F(x_{N+\bar{M}|k}, \theta^a) = \|x_{N+\bar{M}|k} - x_s^a\|^2_P$，以及 (8) 中定义的终端控制律。需要注意的是，新变量 $\bar{M}$ 与 $M$ 和 $\tilde{M}$ 不同，可以根据每个系统所需的性能或其他特定要求进行选择。同时，选择某个标量 $\alpha > 0$ 并定义集合：



$$\Psi_\alpha = \{(x, \theta) \in \mathbb{R}^{n+m} : V_F(x, \theta) \le \alpha\}, \quad (34)$$



使得 $\Psi_\alpha$ 是一个跟踪不变集。最后，定义如下 MPC 问题，它类似于 (30)，但考虑了用户定义的时域扩展，并在目标函数中加入了一个新参数 $\gamma$：



$$\min_{\mathbf{u}, \theta^a} V_{N+\bar{M}}^\gamma(x_k, x_s^r, \mathbf{u}, \theta^a) \quad (35a)$$

$$\text{s.t. } (35b)-(35h) \text{（与 (30) 中的动力学和约束相同，但长度为 } N+\bar{M})$$



特别地，$V_{N, \bar{M}}^\gamma(x_k, x_s^r, \mathbf{u}, \theta^a)$ 表示终端成本按 $\gamma$ 缩放后的成本函数，即 $\gamma V_F(x_{N+\bar{M}|k}, \theta^a)$。

对于上述控制器，当 $\gamma \ge 1$ 时，对于所有 $x \in \Upsilon_{\bar{M}, \gamma}(x_s^r)$，递归可行性和对可容许设定点的收敛性均得到保证 (Limón et al. (2018, 定理 3))。区域 $\Upsilon_{\bar{M}, \gamma}(x_s^r)$ 随着 $\bar{M}$ 或 $\gamma$ 的增加而扩大。最后，值得一提的是，这种控制器设计仅需选择 $\bar{M}, \gamma$ 和 $\alpha$ 的值，避免了对终端不变集的任何显式或隐式估计。

**4. 实例分析**

控制器的性能将通过三个具有不同维度和动力学特性的实例进行说明，仿真使用 YALMIP 建模并调用 GUROBI 求解器（Löfberg, 2004）。

**4.1. 学术示例**

第一个示例是一个低维系统，其动力学矩阵 $A$ 和 $B$ 取自 Limón 等人 (2008)。系统约束为对所有 $k \ge 0$，$\|x_k\|_\infty \le 3$ 且 $\|u_k\|_\infty \le 2$。在控制器设计中，权重矩阵设为 $Q = 10 I_2$，$R = 100 I_2$，$O = 10,000 I_2$；$N$ 设为 10；增益 $K$ 通过离散 LQR 求解获得。同样地，在 $\lambda = 0.9$ 时满足所述条件的 $\tilde{M}$ 值为 $\tilde{M} = 5$。

在 $(x_1, x_2)$ 平面上的状态轨迹如图 1 所示（$x_1$ 和 $x_2$ 为系统状态的两个分量）。我们考虑了不同的设定点（用绿点表示），其中最后一个是不可容许的。调节问题的不变集也显示在图 1 中，可以看出它包含在我们增广不变集 $\Psi_{f, \lambda}^{tr}$ 向 $(x_1, x_2)$ 平面的投影内。

此外，表 1 给出了不同 $\lambda$ 值下累计性能成本的对比。具体而言，累计性能成本通过以下性能指标计算：



$$V_{cc} = \sum_{k=0}^{T_{sim}} (\|x_k - x_{s,k}^a\|_Q^2 + \|u_k - u_{s,k}^a\|_R^2) + \sum_{k=0}^{T_{sim}} \|x_{s,k}^a - x_s^r\|_O^2,$$



其中 $(x_{s,k}^a, u_{s,k}^a)$ 表示在时刻 $k$ 计算出的辅助设定点，$T_{sim}$ 表示仿真时间步数（本例中为 210）。

显然，累计成本随着 $\lambda$ 的减小而增加。这是预料之中的，因为 $\lambda$ 缩放了集合 $\Theta$，从而限制了 $\theta^a$ 的可容许值。同样，从表 1 可以推断出一个有趣的性质：$\tilde{M}$ 随 $\lambda$ 的减小而减小，因为这导致可达设定点远离了约束边界。然而，这可能会导致系统的某些平衡点无法到达。虽然从性能角度来看，选择接近 1 的 $\lambda$ 是最理想的，但这一观察结果提供了一个新的自由度：如果我们确定了参数化系统感兴趣的真实设定点所需的最小 $\lambda$，则可以潜在地减小所需的 $\tilde{M}$。

离线计算终端时域长度 $\tilde{M}$ 所需的时间为 0.5874 秒，而显式计算最大正不变集的时间为 0.6315 秒。至于在线控制，求解问题 (30) 的平均时间为 0.0047 秒，而求解 (9)（即不使用扩展预测时域且具有显式终端集的问题）平均需要 0.0044 秒。虽然这个学术示例的简单性无法展示出显著的计算优势（正如所预期的），但它证明了所提设计在跟踪任务中的适用性，并展示了向更复杂系统（即不变集计算更具挑战性的系统）扩展的可能性。

最后，实施了第 3.5 节提出的替代设计，图 2 展示了一些关键结果。为了说明该方法的优点，我们选择了一个小于 $\tilde{M} = 5$ 的 $\bar{M}$ 值。具体而言，我们考虑 $\bar{M} = 2, \gamma = 20, \alpha = 25$。使用这些参数，得到的 $d$ 值为 5.6322。图 2 显示了针对可容许目标状态 $x_s^r = [-0.5688, -0.0523]^\top$ 的稳定区域 $\Upsilon_{\bar{M}, \gamma}(x_s^r)$。这表明如果系统的初始状态位于该区域内，则无需终端约束的跟踪 MPC (35) 将渐近稳定系统。值得注意的是，集合 $\Upsilon_{\bar{M}, \gamma}(x_s^r)$ 几乎与不变集 $\Psi_{f, \lambda}^{tr}$ 的投影一样大。

**4.2. 规模可变的质量-弹簧-阻尼系统**

本小节将所提出的 MPC 应用于 Riverso 和 Ferrari-Trecate (2012) 以及 Trodden 和 Maestre (2017) 系统的改进版本。它由通过弹簧-阻尼结构连接的多个小车组成。每个小车 $i$ 的动力学模型为：



$$\begin{bmatrix} \dot{r}_i \\ \dot{v}_i \end{bmatrix} = A_{ii} \begin{bmatrix} r_i \\ v_i \end{bmatrix} + \begin{bmatrix} 0 \\ 1/m_i \end{bmatrix} u_i + w_i, \quad (36)$$

$$w_i = \sum_{j \in \mathcal{N}_i} A_{ij} \begin{bmatrix} r_j \\ v_j \end{bmatrix}, \quad (37)$$



其中状态由小车 $i$ 偏离平衡位置的位移 $r_i$ 及其速度 $v_i$ 定义；输入 $u_i$ 代表施加在小车上的力。相邻小车之间由于连接它们的弹簧和阻尼器而存在耦合项。连续时间动力学使用零阶保持器以 0.02 秒的采样时间进行离散化。

我们通过将小车数量 $N_{carts}$ 从 3 逐步增加到 200 来模拟不同的系统规模。目标是将所有小车控制到目标设定点，同时满足状态和输入约束。权重矩阵设为 $Q_i = \text{diag}(1, 1.5), R_i = 20, \lambda = 0.9$。

如表 2 所示，显式计算系统离线不变集的时间随着规模显著增加，最终变得不可行。相比之下，使用所提隐式方法计算 $\tilde{M}$ 的时间在所有情况下都保持在可处理范围内（离线计算时间限制在 48 小时内）。在线计算时间是整个仿真过程中求解 (30) 的平均时间。随着系统规模增大，两种方法之间的在线时间差异逐渐增加，但对于所提方法而言，这从未成为限制因素。例如，对于 100 个小车的系统，显式方法的离线计算时间是所提方法的 79 倍，而在线计算时间仅比显式方法高 2.5 倍。



**4.3. 无人机**

现在采用 Beard (2008) 和 Romagnoli 等人 (2023) 描述的具有 12 个状态变量的无人机模型，来说明所提 MPC 方法对现实系统的适用性。系统状态包括三维空间中的位置 $(p_x, p_y, p_z)$、线速度 $(v_x, v_y, v_z)$，以及滚转、俯仰、偏航角及其相应的角速度。输入 $u$ 包含推力 $F$ 和扭矩 $\tau$。

成本函数定义为 $Q = 10 I_{12}, R = 100 I_4, O = 10^6 I_{12}, N = 5$。满足条件 (26) 的终端时域长度 $\tilde{M}$ 为 92。针对四旋翼无人机的线性模型进行了测试。图 3 显示了位置轨迹以及相应的位置辅助设定点和真实设定点。可以看出，无人机能够到达所有设定点。

在计算时间方面，寻找 $\tilde{M}$ 耗时 16.9558 秒，而显式计算最大不变集需要 137.1241 秒。也就是说，离线计算时间减少了 87.63%，因为寻找 $\tilde{M}$ 简化为求解简单的 LP。拟议的无人机系统被选为一个极限案例：显式计算虽然可行但计算成本昂贵。正如预期的那样，由于预测时域的延长和约束的累积，所提方法的在线计算时间略有增加。然而，在我们的仿真中，这种增加是微乎其微的：求解所提 MPC 问题 (30) 平均需要 0.0164 秒，而求解 (9) 平均需要 0.0114 秒。

**5. 结论**

本文提出了一种基于 Ferramosca 等人 (2009) 和 Limón 等人 (2008) 研究的、具有隐式终端组件的跟踪 MPC 公式。该方法避免了为定义终端约束而对系统最大正不变集进行显式表征的需要，转而使用扩展的预测时域。所提出的控制器可以进行高效设计，同时仍能受益于表征跟踪公式的辅助变量的加入。通过这种方式，它提供了一种能够以可处理的方式处理大型系统的替代方案。研究表明，对于大型系统，控制器设计的离线成本显著降低，而在线计算时间仅略有增加。最后，证明了递归可行性和收敛属性。作为未来的研究方向，将考虑将所提方法应用于非线性系统和更具一般性的约束，例如，将其自然扩展到不仅是多面体，而且是椭球体或两者交集的约束集。