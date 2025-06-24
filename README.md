# Topology Optimization via Lagrangian Differentiation　静解析版 

## 目次
1. はじめに  
2. 問題設定  
3. アルゴリズムの全体像  
4. ステップバイステップ解説  
5. 例題：片持ち梁の剛性化  
6. フローチャート（Mermaid）  
7. 事実と考察の整理  
8. 参考文献  

---

## 1. はじめに
トポロジー最適化は**材料分布を最適化し，目標性能を最大化（または最小化）**する設計手法です。  
この README では **ラグランジアンの微分**を用いた一般的フレームワークをまとめます。  

> **キーワード**: Lagrangian, adjoint method, sensitivity analysis, SIMP, MMA, H¹‐gradient

---

## 2. 問題設定（静解析）
- **設計変数** : 体積密度場 $\rho(\mathbf{x})\in[0,1]$  
- **状態変数** : 変位場 $u(\mathbf{x})$  
- **目的関数** : $F\bigl(u,\rho\bigr)$（例：特定点の変位，構造全体のコンプライアンスなど）  
- **残差（弱形式）** : $R\bigl(u,\rho;v\bigr)=0$  ※$v$ はテスト関数  
- **ラグランジアン** :  

\[
\mathcal{L}(u,\rho,\lambda)=F(u,\rho)-R(u,\rho;\lambda),
\]

ここで $\lambda$ は随伴変数（テスト関数に対応）。

---

## 3. アルゴリズムの全体像
1. **状態方程式**（構造式）を解く  
2. **随伴方程式**を解く  
3. **感度（勾配）**を評価  
4. **密度場の更新**（MMA⁄H¹勾配法など）  
5. 収束判定 → 収束していなければ 1. へ

```mermaid
flowchart TD
  %% --------  Topology Optimization (Static)  -------- %%
  A[STEP 0<br/>初期化：設計密度 ρ₀] --> B[STEP 1<br/>状態方程式を FEM で解く<br/>(構造式) u(ρ)]
  B --> C[STEP 2<br/>随伴方程式を解く<br/>(adjoint) λ]
  C --> D[STEP 3<br/>感度解析<br/>g(x) = ∂𝓛/∂ρ]

  D --> E{{STEP 4<br/>密度更新手法}}
  E -->|H¹ 勾配法| F[ρ ← ρ - α Δ⁻¹ g]
  E -->|MMA|      G[ρ ← MMA(g, 制約)]

  F --> H{収束判定<br/>||∇F|| < ε?}
  G --> H

  H -- いいえ --> B
  H -- はい  --> I[最適トポロジー<br/>ρ\* を出力]

  %% Styling %%
  classDef step fill:#f7f9fc,stroke:#202636,stroke-width:1px,rx:6,ry:6
  class A,B,C,D,E,F,G,H,I step
  class H fill:#ffe5b4
```mermaid

---

## 4. ステップバイステップ解説
### Step 1 : 状態方程式の求解  
$\displaystyle \frac{\partial \mathcal{L}}{\partial \lambda}(u,\rho,\lambda)=0  \; \Longrightarrow \; R(u,\rho;v)=0 \quad\forall v$  
→ FEM により変位 $u$ を得る。

---

### Step 2 : 随伴方程式の求解  
$\displaystyle \frac{\partial \mathcal{L}}{\partial u}(u,\rho,\lambda)=0  \; \Longrightarrow \; \text{adjoint equation for } \lambda$  
→ 同じメッシュ／剛性行列を用いて逆問題を解く。

---

### Step 3 : 感度解析  
$\displaystyle g(\mathbf{x})=\frac{\partial \mathcal{L}}{\partial \rho}(u,\rho,\lambda)$  
→ 体積密度の局所感度。  

---

### Step 4 : 密度更新  
- **H¹‐gradient**  
  \[
  \rho^{k+1} = \rho^{k} - \alpha \, \Delta^{-1} g
  \]
- **MMA**  
  ```text
  minimize   MMA sub-problem
  subject to volume constraint
