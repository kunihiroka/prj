# 初学者のための回転dq座標磁束オブザーバ解説

この文書は、誘導機の磁束オブザーバを初めて学ぶ人が、式の意味、式変形、ゲイン設計、実装時の注意点までを順番に理解できるように書いた補足資料である。既存の技術レポートは成果物と評価結果を中心にまとめている。一方、この文書では「なぜその式になるのか」を重視し、数式の展開を省略しない。

対象とするオブザーバは以下である。

| 方式 | 座標系 | 状態 | ゲイン設計 |
|---|---|---|---|
| 方式A | 回転dq座標 | 一次磁束2成分、二次磁束2成分 | 4個の実極を指定し、Sylvester方程式で一般の実数行列 $H$ を求める |
| 方式B | 回転dq座標 | 二次磁束2成分。一次磁束と一次電流は測定一次電流から代数的に再構成 | 堀ほか1986年論文5.3節に従い、$K=k_1I+k_2J$ を $\alpha,\beta$ から直接計算する |
| 方式C | 推定ロータ磁束座標 | 推定一次電流2成分、推定ロータ磁束大きさ、推定角 | SLED 2023 Appendix Aの閉形式式を使う |

方式Aは回転dq座標4状態オブザーバである。方式Bは同じ回転dq座標で動かすが、論文5.3節に忠実な二次磁束オブザーバであり、方式Aの4状態Sylvester法ではない。方式Cは状態変数と座標の取り方がさらに異なり、すべり周波数の計算式が $\hat{\psi}_{Rq}=0$ を保つためのオブザーバ構成そのものになっている。

## 1. まず「オブザーバ」とは何か

モータ制御では、すべての内部状態を直接測れるわけではない。一次電流は電流センサで測れるが、ロータ磁束は通常は直接測れない。そこで、測れる量から測れない量を推定する仕組みを使う。この推定器をオブザーバと呼ぶ。

誘導機の磁束オブザーバでは、主に以下を推定する。

| 量 | 測定できるか | オブザーバで扱う理由 |
|---|---|---|
| 一次電圧 $v_s$ | 指令値または推定値として得られる | 磁束を積分する入力になる |
| 一次電流 $i_s$ | 測定できる | 推定値との誤差を使って補正する |
| 一次磁束 $\psi_s$ | 通常は測れない | 状態として推定する |
| 二次磁束 $\psi_r$ | 通常は測れない | ベクトル制御やすべり計算に重要 |

オブザーバの基本形は非常に単純である。

```math
\text{推定値の変化}
=
\text{モータモデルによる予測}
+
\text{測定値とのずれによる補正}
```

この「測定値とのずれによる補正」の強さを決める行列がオブザーバゲイン $H$ である。

## 2. dq座標と90度回転行列

三相交流のまま式を書くと、すべての量が正弦波になる。制御やオブザーバ設計では、回転dq座標に変換して扱う。回転dq座標では、理想的な定常状態の電流や磁束を直流量として見られる。

この文書では、任意のdqベクトルを以下で表す。

```math
u=
\begin{bmatrix}
u_d\\
u_q
\end{bmatrix}
```

90度回転行列を以下で定義する。

```math
J=
\begin{bmatrix}
0 & -1\\
1 & 0
\end{bmatrix}
```

この行列をベクトルに掛けると、

```math
Ju=
\begin{bmatrix}
0 & -1\\
1 & 0
\end{bmatrix}
\begin{bmatrix}
u_d\\
u_q
\end{bmatrix}
```

行列積を1行ずつ計算すると、

```math
Ju=
\begin{bmatrix}
0\cdot u_d+(-1)\cdot u_q\\
1\cdot u_d+0\cdot u_q
\end{bmatrix}
```

したがって、

```math
Ju=
\begin{bmatrix}
-u_q\\
u_d
\end{bmatrix}
```

である。つまり $J$ はdq平面上のベクトルを反時計回りに90度回す行列である。

実装で使う磁束方程式には $-\omega J\psi$ が出る。例えば、

```math
-\omega J
\begin{bmatrix}
\psi_d\\
\psi_q
\end{bmatrix}
=
-\omega
\begin{bmatrix}
-\psi_q\\
\psi_d
\end{bmatrix}
```

なので、

```math
-\omega J\psi
=
\begin{bmatrix}
\omega\psi_q\\
-\omega\psi_d
\end{bmatrix}
```

となる。これはC実装で、d軸式に $+\omega\psi_q$、q軸式に $-\omega\psi_d$ が現れることと一致する。

## 3. 誘導機のT形定数と磁束・電流の関係

一次漏れインダクタンス、二次漏れインダクタンス、相互インダクタンスを以下で表す。

```math
L_{ls},\quad L_{lr},\quad M
```

一次自己インダクタンスと二次自己インダクタンスは、

```math
L_s=L_{ls}+M
```

```math
L_r=L_{lr}+M
```

である。

一次磁束 $\psi_s$ と二次磁束 $\psi_r$ は、一次電流 $i_s$ と二次電流 $i_r$ から以下で決まる。

```math
\psi_s=L_s i_s+M i_r
```

```math
\psi_r=M i_s+L_r i_r
```

ここで、$\psi_s,\psi_r,i_s,i_r$ はすべてdqベクトルである。つまり、例えば

```math
\psi_s=
\begin{bmatrix}
\psi_{sd}\\
\psi_{sq}
\end{bmatrix},
\qquad
i_s=
\begin{bmatrix}
i_{sd}\\
i_{sq}
\end{bmatrix}
```

である。

磁束オブザーバでは、状態として磁束を持つ。一方、測定できるのは電流である。そのため、磁束から電流を計算できる形が必要になる。

まず、磁束と電流の関係を行列で書く。

```math
\begin{bmatrix}
\psi_s\\
\psi_r
\end{bmatrix}
=
\begin{bmatrix}
L_s I_2 & M I_2\\
M I_2 & L_r I_2
\end{bmatrix}
\begin{bmatrix}
i_s\\
i_r
\end{bmatrix}
```

ここで $I_2$ は2行2列の単位行列である。

この関係を電流について解く。2つの式をもう一度書く。

```math
\psi_s=L_s i_s+M i_r
```

```math
\psi_r=M i_s+L_r i_r
```

まず、一次電流 $i_s$ を求める。1本目の式に $L_r$ を掛ける。

```math
L_r\psi_s=L_sL_r i_s+ML_r i_r
```

2本目の式に $M$ を掛ける。

```math
M\psi_r=M^2 i_s+ML_r i_r
```

この2式を引き算する。

```math
L_r\psi_s-M\psi_r
=
(L_sL_r i_s+ML_r i_r)-(M^2 i_s+ML_r i_r)
```

右辺の $ML_r i_r$ は打ち消される。

```math
L_r\psi_s-M\psi_r
=
L_sL_r i_s-M^2 i_s
```

右辺を $i_s$ でくくる。

```math
L_r\psi_s-M\psi_r
=
(L_sL_r-M^2)i_s
```

ここで、

```math
D=L_sL_r-M^2
```

と定義する。すると、

```math
L_r\psi_s-M\psi_r=D i_s
```

なので、

```math
i_s=\frac{L_r\psi_s-M\psi_r}{D}
```

同様に、二次電流は以下になる。

```math
i_r=\frac{-M\psi_s+L_s\psi_r}{D}
```

成分で書くと、一次電流は以下である。

```math
i_{sd}=\frac{L_r\psi_{sd}-M\psi_{rd}}{D}
```

```math
i_{sq}=\frac{L_r\psi_{sq}-M\psi_{rq}}{D}
```

二次電流は以下である。

```math
i_{rd}=\frac{-M\psi_{sd}+L_s\psi_{rd}}{D}
```

```math
i_{rq}=\frac{-M\psi_{sq}+L_s\psi_{rq}}{D}
```

## 4. 回転dq座標の誘導機状態方程式

磁束オブザーバの状態を以下にする。

```math
x=
\begin{bmatrix}
\psi_{sd}\\
\psi_{sq}\\
\psi_{rd}\\
\psi_{rq}
\end{bmatrix}
```

これは、一次磁束2成分と二次磁束2成分をそのまま並べた4状態である。固定座標ではなく、角速度 $\omega_k$ で回る回転dq座標上の状態である。

一次側の磁束方程式は、

```math
\dot{\psi}_s=v_s-R_s i_s-\omega_k J\psi_s
```

である。二次側は短絡されているので二次電圧は0であり、

```math
\dot{\psi}_r=-R_r i_r-(\omega_k-\omega_r)J\psi_r
```

である。ここで、

```math
\omega_r=p\omega_m
```

である。$p$ は極対数、$\omega_m$ は機械角速度である。

一次側へ3章の式を代入する。

```math
\dot{\psi}_s
=
v_s
-R_s\left(\frac{L_r\psi_s-M\psi_r}{D}\right)
-\omega_k J\psi_s
```

括弧を外す。

```math
\dot{\psi}_s
=
v_s
-\frac{R_sL_r}{D}\psi_s
+\frac{R_sM}{D}\psi_r
-\omega_k J\psi_s
```

$\psi_s$ に掛かる項と $\psi_r$ に掛かる項に分ける。

```math
\dot{\psi}_s
=
\left(-\frac{R_sL_r}{D}I_2-\omega_k J\right)\psi_s
+
\left(\frac{R_sM}{D}I_2\right)\psi_r
+
v_s
```

次に二次側へ3章の式を代入する。

```math
\dot{\psi}_r
=
-R_r\left(\frac{-M\psi_s+L_s\psi_r}{D}\right)
-(\omega_k-\omega_r)J\psi_r
```

括弧を外す。

```math
\dot{\psi}_r
=
\frac{R_rM}{D}\psi_s
-\frac{R_rL_s}{D}\psi_r
-(\omega_k-\omega_r)J\psi_r
```

$\psi_s$ に掛かる項と $\psi_r$ に掛かる項に分ける。

```math
\dot{\psi}_r
=
\left(\frac{R_rM}{D}I_2\right)\psi_s
+
\left(-\frac{R_rL_s}{D}I_2-(\omega_k-\omega_r)J\right)\psi_r
```

以上をまとめると、

```math
\dot{x}=Ax+Bv_s
```

であり、

```math
A=
\begin{bmatrix}
-\frac{R_sL_r}{D}I_2-\omega_kJ & \frac{R_sM}{D}I_2\\
\frac{R_rM}{D}I_2 & -\frac{R_rL_s}{D}I_2-(\omega_k-\omega_r)J
\end{bmatrix}
```

```math
B=
\begin{bmatrix}
I_2\\
0_{2\times2}
\end{bmatrix}
```

である。

一次電流は測定できるので、出力を

```math
y=i_s
```

とする。3章の式より、

```math
y=i_s=\frac{L_r\psi_s-M\psi_r}{D}
```

だから、

```math
y=Cx
```

```math
C=
\begin{bmatrix}
\frac{L_r}{D}I_2 & -\frac{M}{D}I_2
\end{bmatrix}
```

である。

## 5. Luenberger型フルオーダ磁束オブザーバ

モータの状態方程式は以下だった。

```math
\dot{x}=Ax+Bv_s
```

```math
y=Cx
```

これと同じ形のモデルをオブザーバ内部に持つ。ただし、オブザーバが持つ状態は真値 $x$ ではなく推定値 $\hat{x}$ である。

モデルだけで予測するなら、

```math
\dot{\hat{x}}=A\hat{x}+Bv_s
```

となる。しかし、これだけでは初期値誤差、電圧誤差、定数誤差で推定値がずれ続ける。そこで、測定電流と推定電流の差を使って補正する。

推定電流は、

```math
\hat{y}=C\hat{x}
```

である。測定電流 $y$ と推定電流 $\hat{y}$ の差は、

```math
y-\hat{y}
=
y-C\hat{x}
```

である。この差を出力誤差、またはイノベーションと呼ぶ。

補正付きのオブザーバは、

```math
\dot{\hat{x}}=A\hat{x}+Bv_s+H(y-C\hat{x})
```

である。

ここで $y-C\hat{x}$ は2成分のベクトルである。状態 $\hat{x}$ は4成分である。したがって、2成分の誤差を4成分の状態補正へ変換するために、$H$ は4行2列になる。

```math
H\in\mathbb{R}^{4\times2}
```

実装上は、

```math
H=
\begin{bmatrix}
H_{00} & H_{01}\\
H_{10} & H_{11}\\
H_{20} & H_{21}\\
H_{30} & H_{31}
\end{bmatrix}
```

である。

## 6. 誤差方程式を丁寧に導く

オブザーバ設計で最も重要なのは、推定誤差が0へ収束するかどうかである。推定誤差を以下で定義する。

```math
\tilde{x}=x-\hat{x}
```

この式を時間で微分する。

```math
\dot{\tilde{x}}
=
\frac{d}{dt}(x-\hat{x})
```

微分は引き算に分配できるので、

```math
\dot{\tilde{x}}
=
\dot{x}-\dot{\hat{x}}
```

真値モデルとオブザーバ式を代入する。

```math
\dot{\tilde{x}}
=
(Ax+Bv_s)
-
(A\hat{x}+Bv_s+H(y-C\hat{x}))
```

右辺を分解する。

```math
\dot{\tilde{x}}
=
Ax+Bv_s
-A\hat{x}
-Bv_s
-H(y-C\hat{x})
```

$+Bv_s$ と $-Bv_s$ は打ち消される。

```math
\dot{\tilde{x}}
=
Ax-A\hat{x}-H(y-C\hat{x})
```

最初の2項を $A$ でくくる。

```math
Ax-A\hat{x}=A(x-\hat{x})
```

したがって、

```math
\dot{\tilde{x}}
=
A(x-\hat{x})-H(y-C\hat{x})
```

ここで $y=Cx$ なので、

```math
y-C\hat{x}=Cx-C\hat{x}
```

右辺を $C$ でくくる。

```math
Cx-C\hat{x}=C(x-\hat{x})
```

よって、

```math
y-C\hat{x}=C(x-\hat{x})
```

$\tilde{x}=x-\hat{x}$ だから、

```math
y-C\hat{x}=C\tilde{x}
```

これを誤差の式に戻す。

```math
\dot{\tilde{x}}
=
A\tilde{x}-HC\tilde{x}
```

共通の $\tilde{x}$ でまとめる。

```math
\dot{\tilde{x}}
=
(A-HC)\tilde{x}
```

これが誤差方程式である。この式から分かる重要なことは、誤差の動きは $A-HC$ で決まるということである。

## 7. なぜ固有値を見ると収束が分かるのか

まず1次元の簡単な微分方程式を見る。

```math
\dot{e}=\lambda e
```

この解は、

```math
e(t)=e(0)e^{\lambda t}
```

である。

もし $\lambda<0$ なら、時間が進むほど $e^{\lambda t}$ は0へ近づく。したがって、

```math
\lim_{t\to\infty}e(t)=0
```

である。

もし $\lambda>0$ なら、$e^{\lambda t}$ は増加する。したがって誤差は発散する。

4状態の誤差方程式

```math
\dot{\tilde{x}}=(A-HC)\tilde{x}
```

でも考え方は同じである。行列 $A-HC$ の固有値を

```math
\lambda_1,\lambda_2,\lambda_3,\lambda_4
```

とすると、誤差は複数のモードの足し合わせとして動く。各モードは概念的に、

```math
e^{\lambda_i t}
```

の形を持つ。

したがって、すべての固有値について、

```math
\mathrm{Re}(\lambda_i)<0
```

であれば、すべての誤差モードが減衰する。このとき、

```math
\lim_{t\to\infty}\|\tilde{x}(t)\|=0
```

となり、推定誤差は0へ収束する。これを漸近安定という。

固有値の実部は収束速度を決める。例えば、

```math
\lambda=-1000
```

なら、おおまかな時定数は

```math
\tau\simeq\frac{1}{1000}=1\ \mathrm{ms}
```

である。実部がより負に大きいほど収束は速い。

固有値に虚部がある場合、

```math
\lambda=-\alpha+j\beta
```

のように書ける。このとき、誤差は

```math
e^{-\alpha t}
```

で減衰しながら、角周波数 $\beta$ で振動する。つまり、実部 $-\alpha$ は減衰、虚部 $\beta$ は振動を決める。

## 8. 極配置法の考え方

オブザーバの極とは、誤差方程式の固有値である。つまり、

```math
\lambda(A-HC)
```

がオブザーバの極である。

設計者は、オブザーバの誤差をどのくらい速く消したいかを決める。例えば、4つの極を

```math
p_1=-2200
```

```math
p_2=-2750
```

```math
p_3=-3410
```

```math
p_4=-4400
```

に置きたいとする。これは、すべての誤差モードを負の実軸上に置く設計である。

目標は、

```math
\lambda(A-HC)=\{p_1,p_2,p_3,p_4\}
```

を満たす $H$ を求めることである。

ここで注意すべきことがある。$A$ は4行4列、$C$ は2行4列、$H$ は4行2列である。したがって $HC$ は、

```math
HC
\in
\mathbb{R}^{4\times2}
\mathbb{R}^{2\times4}
=
\mathbb{R}^{4\times4}
```

である。よって、

```math
A-HC
```

は4行4列になり、4個の固有値を持つ。

ここで、誘導機の回転dq座標4状態モデルについて、$A-HC$ を実際に展開しておく。状態は、

```math
x=
\begin{bmatrix}
\psi_{sd} & \psi_{sq} & \psi_{rd} & \psi_{rq}
\end{bmatrix}^{T}
```

である。出力は一次電流

```math
y=
\begin{bmatrix}
i_{sd} & i_{sq}
\end{bmatrix}^{T}
=Cx
```

であり、

```math
C=
\begin{bmatrix}
\frac{L_r}{D} & 0 & -\frac{M}{D} & 0\\
0 & \frac{L_r}{D} & 0 & -\frac{M}{D}
\end{bmatrix}
```

である。オブザーバゲインを、

```math
H=
\begin{bmatrix}
h_{11} & h_{12}\\
h_{21} & h_{22}\\
h_{31} & h_{32}\\
h_{41} & h_{42}
\end{bmatrix}
```

と置く。まず $HC$ を計算する。$H$ は4行2列、$C$ は2行4列なので、

```math
HC=
\begin{bmatrix}
h_{11} & h_{12}\\
h_{21} & h_{22}\\
h_{31} & h_{32}\\
h_{41} & h_{42}
\end{bmatrix}
\begin{bmatrix}
\frac{L_r}{D} & 0 & -\frac{M}{D} & 0\\
0 & \frac{L_r}{D} & 0 & -\frac{M}{D}
\end{bmatrix}
```

である。1行目を例にすると、

```math
\begin{bmatrix}
h_{11} & h_{12}
\end{bmatrix}
\begin{bmatrix}
\frac{L_r}{D} & 0 & -\frac{M}{D} & 0\\
0 & \frac{L_r}{D} & 0 & -\frac{M}{D}
\end{bmatrix}
=
\begin{bmatrix}
\frac{L_r}{D}h_{11} &
\frac{L_r}{D}h_{12} &
-\frac{M}{D}h_{11} &
-\frac{M}{D}h_{12}
\end{bmatrix}
```

となる。同じ計算を4行分行うと、

```math
HC=
\begin{bmatrix}
\frac{L_r}{D}h_{11} & \frac{L_r}{D}h_{12} & -\frac{M}{D}h_{11} & -\frac{M}{D}h_{12}\\
\frac{L_r}{D}h_{21} & \frac{L_r}{D}h_{22} & -\frac{M}{D}h_{21} & -\frac{M}{D}h_{22}\\
\frac{L_r}{D}h_{31} & \frac{L_r}{D}h_{32} & -\frac{M}{D}h_{31} & -\frac{M}{D}h_{32}\\
\frac{L_r}{D}h_{41} & \frac{L_r}{D}h_{42} & -\frac{M}{D}h_{41} & -\frac{M}{D}h_{42}
\end{bmatrix}
```

である。

一方、誘導機の $A$ は、

```math
A=
\begin{bmatrix}
-\frac{R_sL_r}{D} & \omega_k & \frac{R_sM}{D} & 0\\
-\omega_k & -\frac{R_sL_r}{D} & 0 & \frac{R_sM}{D}\\
\frac{R_rM}{D} & 0 & -\frac{R_rL_s}{D} & \omega_{\mathrm{slip}}\\
0 & \frac{R_rM}{D} & -\omega_{\mathrm{slip}} & -\frac{R_rL_s}{D}
\end{bmatrix}
```

である。ここで、

```math
\omega_{\mathrm{slip}}=\omega_k-\omega_r
```

であり、

```math
\omega_r=p\omega_m
```

である。

したがって、$A-HC$ は $A$ の各要素から $HC$ の対応する要素を引けばよい。

```math
A-HC=
\begin{bmatrix}
-\frac{R_sL_r}{D}-\frac{L_r}{D}h_{11}
&
\omega_k-\frac{L_r}{D}h_{12}
&
\frac{R_sM}{D}+\frac{M}{D}h_{11}
&
\frac{M}{D}h_{12}
\\
-\omega_k-\frac{L_r}{D}h_{21}
&
-\frac{R_sL_r}{D}-\frac{L_r}{D}h_{22}
&
\frac{M}{D}h_{21}
&
\frac{R_sM}{D}+\frac{M}{D}h_{22}
\\
\frac{R_rM}{D}-\frac{L_r}{D}h_{31}
&
-\frac{L_r}{D}h_{32}
&
-\frac{R_rL_s}{D}+\frac{M}{D}h_{31}
&
\omega_{\mathrm{slip}}+\frac{M}{D}h_{32}
\\
-\frac{L_r}{D}h_{41}
&
\frac{R_rM}{D}-\frac{L_r}{D}h_{42}
&
-\omega_{\mathrm{slip}}+\frac{M}{D}h_{41}
&
-\frac{R_rL_s}{D}+\frac{M}{D}h_{42}
\end{bmatrix}
```

この4行4列行列の特性方程式

```math
\det\{sI-(A-HC)\}=0
```

の根が、誘導機磁束オブザーバの極である。つまり、ゲイン $H$ を変えると上の行列の各要素が変わり、その結果として特性方程式の根、すなわちオブザーバ極が動く。

### 8.1 $\det\{sI-(A-HC)\}$ の具体形

前節の $A-HC$ をそのまま使うと式が長くなるため、まず記号を短く置く。

```math
a=\frac{R_sL_r}{D},\qquad
b=\frac{R_sM}{D},\qquad
c=\frac{L_r}{D},\qquad
m=\frac{M}{D}
```

```math
r=\frac{R_rM}{D},\qquad
d=\frac{R_rL_s}{D}
```

また、

```math
\omega_{\mathrm{slip}}=\omega_k-\omega_r
```

とする。このとき、前節の $A-HC$ は以下の形である。

```math
A-HC=
\begin{bmatrix}
-a-ch_{11} & \omega_k-ch_{12} & b+mh_{11} & mh_{12}\\
-\omega_k-ch_{21} & -a-ch_{22} & mh_{21} & b+mh_{22}\\
r-ch_{31} & -ch_{32} & -d+mh_{31} & \omega_{\mathrm{slip}}+mh_{32}\\
-ch_{41} & r-ch_{42} & -\omega_{\mathrm{slip}}+mh_{41} & -d+mh_{42}
\end{bmatrix}
```

特性方程式は、

```math
P(s)=\det\{sI-(A-HC)\}=0
```

である。ここで $sI-(A-HC)$ を具体的に書くと、

```math
P(s)=
\det
\begin{bmatrix}
s+a+ch_{11}
&
-\omega_k+ch_{12}
&
-b-mh_{11}
&
-mh_{12}
\\
\omega_k+ch_{21}
&
s+a+ch_{22}
&
-mh_{21}
&
-b-mh_{22}
\\
-r+ch_{31}
&
ch_{32}
&
s+d-mh_{31}
&
-\omega_{\mathrm{slip}}-mh_{32}
\\
ch_{41}
&
-r+ch_{42}
&
\omega_{\mathrm{slip}}-mh_{41}
&
s+d-mh_{42}
\end{bmatrix}
```

である。この行列式を展開すると、必ず4次多項式になる。

```math
P(s)=s^4+a_3s^3+a_2s^2+a_1s+a_0
```

まず、$s^3$ 係数は比較的簡単に書ける。特性多項式 $\det(sI-F)$ の $s^3$ 係数は $-\mathrm{tr}(F)$ である。ここで、

```math
F=A-HC
```

と置くと、

```math
a_3=-\mathrm{tr}(F)
```

である。$F=A-HC$ の対角成分は、

```math
F_{11}=-a-ch_{11}
```

```math
F_{22}=-a-ch_{22}
```

```math
F_{33}=-d+mh_{31}
```

```math
F_{44}=-d+mh_{42}
```

なので、

```math
\mathrm{tr}(F)=-2a-2d-c(h_{11}+h_{22})+m(h_{31}+h_{42})
```

したがって、

```math
a_3=2a+2d+c(h_{11}+h_{22})-m(h_{31}+h_{42})
```

である。

残りの係数を完全展開すると非常に長くなる。実装や検算では、完全展開式を手で扱うより、トレースを使った以下の形で計算するほうが安全である。

```math
a_2=\frac{1}{2}\left[\{\mathrm{tr}(F)\}^2-\mathrm{tr}(F^2)\right]
```

```math
a_1=-\frac{1}{6}\left[\{\mathrm{tr}(F)\}^3-3\mathrm{tr}(F)\mathrm{tr}(F^2)+2\mathrm{tr}(F^3)\right]
```

```math
a_0=\det(-F)
```

4次なので $\det(-F)=\det(F)$ であり、したがって

```math
a_0=\det(F)
```

としてよい。

これらを使えば、特性多項式は

```math
P(s)=s^4+a_3s^3+a_2s^2+a_1s+a_0
```

として計算できる。

極配置で $H$ が正しく設計されていれば、この特性多項式は指定極

```math
p_1,\quad p_2,\quad p_3,\quad p_4
```

に対して、

```math
P(s)=(s-p_1)(s-p_2)(s-p_3)(s-p_4)
```

になる。例えば、すべての極を負の実数に置いた場合、$p_i<0$ なので、各因子は

```math
s-p_i=s+|p_i|
```

となる。したがって、指定極が左半平面にあれば、理想モデルでは推定誤差は0へ収束する。

## 9. Sylvester方程式によるHの求め方

方式Aでは、直接 $H$ を手計算で求めるのではなく、Sylvester方程式を使う。ここでは、その式がどこから来るかを丁寧に導く。方式Bは論文5.3節の $K=k_1I+k_2J$ を直接計算するため、このSylvester方程式は使わない。

誤差方程式は、

```math
\dot{\tilde{x}}=(A-HC)\tilde{x}
```

である。設計者が望む誤差方程式を、

```math
\dot{z}=Fz
```

とする。ここで $F$ は目標極を持つ4行4列行列である。

方式Aなら、

```math
F=
\begin{bmatrix}
p_1 & 0 & 0 & 0\\
0 & p_2 & 0 & 0\\
0 & 0 & p_3 & 0\\
0 & 0 & 0 & p_4
\end{bmatrix}
```

である。

実際の誤差 $\tilde{x}$ と、設計用の誤差座標 $z$ の間に、以下の変換を置く。

```math
z=T\tilde{x}
```

ここで $T$ は正則な4行4列行列である。正則とは、逆行列が存在するという意味である。

この式を時間微分する。

```math
\dot{z}=T\dot{\tilde{x}}
```

誤差方程式を代入する。

```math
\dot{z}=T(A-HC)\tilde{x}
```

括弧を展開する。

```math
\dot{z}=TA\tilde{x}-THC\tilde{x}
```

ここで、

```math
G=TH
```

と定義する。すると、

```math
\dot{z}=TA\tilde{x}-GC\tilde{x}
```

一方、設計用の目標誤差方程式は、

```math
\dot{z}=Fz
```

である。$z=T\tilde{x}$ を代入する。

```math
\dot{z}=FT\tilde{x}
```

同じ $\dot{z}$ を表しているので、

```math
TA\tilde{x}-GC\tilde{x}=FT\tilde{x}
```

この式は任意の $\tilde{x}$ について成り立ってほしい。したがって、$\tilde{x}$ を外して、

```math
TA-GC=FT
```

を得る。これを並べ替える。

```math
TA-FT=GC
```

これが本実装で使うSylvester方程式である。

```math
TA-FT=GC
```

この式を解いて $T$ が得られれば、

```math
G=TH
```

だったので、左から $T^{-1}$ を掛ける。

```math
T^{-1}G=T^{-1}TH
```

$T^{-1}T=I$ なので、

```math
T^{-1}G=H
```

したがって、

```math
H=T^{-1}G
```

である。

設計手順をまとめると以下になる。

1. 現在の動作点から $A,C$ を作る。
2. 目標極を決め、目標行列 $F$ を作る。
3. 出力誤差を4状態へ分配する設計行列 $G$ を決める。
4. Sylvester方程式 $TA-FT=GC$ を解いて $T$ を求める。
5. $H=T^{-1}G$ を計算する。
6. 得られた $H$ により、$A-HC$ の固有値が目標極になる。

重要なのは、$G$ は物理量そのものではなく、極配置計算を成立させるための設計行列であるという点である。$G$ の選び方が悪いと $T$ が特異になり、$H$ を計算できない。

## 10. 方式A: 4実極を指定する設計

方式Aでは、4つの極をすべて実数で指定する。

```math
F=
\mathrm{diag}(p_1,p_2,p_3,p_4)
```

例えば、

```math
p_1=-\omega_o
```

```math
p_2=-1.25\omega_o
```

```math
p_3=-1.55\omega_o
```

```math
p_4=-2.0\omega_o
```

のようにする。

ここで $\omega_o$ はオブザーバ帯域を表す設計値である。すべての極の実部が負なので、理想モデルでは誤差は0へ収束する。

この方式の利点は分かりやすさである。指定した4つの極がそのまま誤差方程式の極になるので、検証しやすい。

一方で、オンラインで毎周期Sylvester方程式を解くと計算負荷が重い。組込み機器では、オフラインでゲインテーブルを作る、またはより軽い方式Bや方式Cを検討するのが現実的である。

## 11. 方式B: 論文5.3節に忠実なモデル修正型二次磁束オブザーバ

方式Bは、堀ほか1986年論文5.3節の構成に合わせた方式である。方式Aのように一次磁束と二次磁束をまとめた4状態に対して $H$ を配置するのではない。動的に積分する状態は二次磁束

```math
\hat{\psi}_r=
\begin{bmatrix}
\hat{\psi}_{rd}\\
\hat{\psi}_{rq}
\end{bmatrix}
```

であり、一次磁束は測定一次電流から代数的に再構成する。

```math
\hat{\psi}_s
=
L_{\sigma s} i_s
+
\frac{M}{L_r}\hat{\psi}_r
```

ここで、

```math
L_{\sigma s}=L_s-\frac{M^2}{L_r}
```

である。したがって、方式Bでは一次電流推定値は測定一次電流に一致する。一次電流を独立に推定する4状態フルオーダ方式ではなく、測定一次電流を使って二次磁束を推定するモデル修正型オブザーバである。

論文5.3節の補正ゲインは、一般の4行2列行列 $H$ ではなく、以下の2行2列行列で表される。

```math
K=k_1I+k_2J
```

展開すると、

```math
K=
\begin{bmatrix}
k_1 & -k_2\\
k_2 & k_1
\end{bmatrix}
```

である。$k_1$ は同相成分の補正、$k_2$ は90度回転した成分の補正を表す。複素ゲインとして実装しているわけではなく、実装ではこの実数2行2列行列をそのまま使う。

論文の固定座標系では、目標極を

```math
s=-\alpha \pm j\beta
```

としたとき、論文式(33)(34)に対応するゲインは以下になる。

```math
k_1=
\frac{L_r}{M}
\left[
1-
\frac{(R_r/L_r)\alpha+\omega_r\beta}
{\alpha^2+\beta^2}
\right]
```

```math
k_2=
\frac{L_r}{M}
\left[
\frac{\omega_r\alpha-(R_r/L_r)\beta}
{\alpha^2+\beta^2}
\right]
```

ここで $\omega_r$ は電気角のロータ角速度である。$\alpha$ は誤差の減衰速度、$\beta$ は誤差が回り込む角周波数である。

本成果物ではオブザーバを回転dq座標で実行する。回転座標の角速度を $\omega_k$ とし、

```math
\Omega=\omega_r-\omega_k
```

と置く。論文の固定座標式を回転dq座標へ移すと、実装で使うゲインは以下になる。

```math
a=\frac{R_r}{L_r}
```

```math
c=\frac{M}{L_r}
```

```math
q=\beta+\omega_k
```

```math
r_1=\frac{\alpha-a}{c}
```

```math
r_2=\frac{\beta-\Omega}{c}
```

```math
k_1=
\frac{\alpha r_1+q r_2}
{\alpha^2+q^2}
```

```math
k_2=
\frac{q r_1-\alpha r_2}
{\alpha^2+q^2}
```

$\omega_k=0$ とすれば、$\Omega=\omega_r$、$q=\beta$ となるため、この式は上の論文式(33)(34)に戻る。つまり、今回の式は論文5.3節を回転dq座標に書き換えたものであり、方式AのSylvester法とは別物である。

実装では、回転dq座標の一次電圧方程式を使って二次磁束微分を計算する。$v_s$ を一次電圧、$i_s$ を測定一次電流、$J$ を90度回転行列とすると、以下の2行2列連立方程式を各周期で解く。

```math
\left(I-\frac{M}{L_r}K\right)\dot{\hat{\psi}}_r
=
\left(-aI+\Omega J\right)\hat{\psi}_r
+Ma i_s
+R_sKi_s
+L_{\sigma s}K\dot{i}_s
+L_{\sigma s}\omega_kKJi_s
+\frac{M}{L_r}\omega_kKJ\hat{\psi}_r
-Kv_s
```

右辺の各項の意味は以下である。

| 項 | 意味 |
|---|---|
| $(-aI+\Omega J)\hat{\psi}_r$ | 二次抵抗による減衰と、回転座標から見た二次磁束の回り込み |
| $Ma i_s$ | 測定一次電流からロータ磁束を作る電流モデル項 |
| $K$ を含む項 | 一次電圧方程式との差を使ってモデルを修正する補正項 |
| $L_{\sigma s}K\dot{i}_s$ | 測定一次電流の変化分による一次漏れ磁束の補正 |
| $\omega_kKJ$ を含む項 | 回転dq座標で表れる速度起電力の補正 |

この方式の特徴は、ゲイン計算も更新式も2次元の実数演算で済むことである。毎周期4状態のSylvester方程式を解かないため、方式Aより組込み実装の計算負荷は小さい。

## 12. 方式C: SLED Appendix Aの推定ロータ磁束座標オブザーバ

方式Cは、方式Aや方式Bと構成が異なる。方式Aは状態を

```math
\hat{x}=
\begin{bmatrix}
\hat{\psi}_{sd}\\
\hat{\psi}_{sq}\\
\hat{\psi}_{rd}\\
\hat{\psi}_{rq}
\end{bmatrix}
```

とする。方式Bは二次磁束2成分を動的状態とし、一次磁束は測定一次電流から再構成する。一方、方式Cでは推定ロータ磁束方向をd軸に選ぶ。したがって、

```math
\hat{\psi}_{Rq}=0
```

を常に満たす座標系を使う。

状態は、

```math
\hat{x}=
\begin{bmatrix}
\hat{i}_{sd}\\
\hat{i}_{sq}\\
\hat{\psi}_R
\end{bmatrix}
```

の3成分に見える。ただし、これはフルオーダ情報を失ったという意味ではない。推定ロータ磁束の角度は、状態ベクトルの外側で

```math
\dot{\hat{\theta}}=\omega_s
```

として積分される。したがって、推定電流2成分、推定ロータ磁束の大きさ、推定ロータ磁束角を合わせて、物理的にはフルオーダの情報を持っている。

方式Cでは、T形定数を逆Gamma形定数へ変換する。

```math
L_{\sigma}=L_s-\frac{M^2}{L_r}
```

```math
L_M=\frac{M^2}{L_r}
```

```math
R_R=R_r\left(\frac{M}{L_r}\right)^2
```

```math
R_{\sigma}=R_s+R_R
```

```math
\alpha=\frac{R_R}{L_M}
```

ただし、これは必ず逆Gamma形へ変換してから実装しなければならない、という意味ではない。T形定数の組合せとして同じ係数を直接作れる。

```math
L_s=L_{ls}+M
```

```math
L_r=L_{lr}+M
```

```math
D=L_sL_r-M^2
```

T形の物理ロータ磁束を $\phi_r$ とすると、一次磁束は

```math
\psi_s
=
\left(L_s-\frac{M^2}{L_r}\right)i_s
+
\frac{M}{L_r}\phi_r
```

である。したがって、

```math
L_{\sigma T}=\frac{D}{L_r}
```

```math
\rho=\frac{M}{L_r}
```

と置けば、

```math
\psi_s=L_{\sigma T}i_s+\rho\phi_r
```

と書ける。Appendix Aの $\hat{\psi}_R$ は、T形の物理ロータ磁束 $\hat{\phi}_R$ そのものではなく、

```math
\hat{\psi}_R=\rho\hat{\phi}_R
```

である。この関係を使えば、逆Gamma形の名前を使わずに方式Cを書ける。

T形定数だけで使う係数は以下である。

```math
R_{\Sigma T}=R_s+R_r\left(\frac{M}{L_r}\right)^2
```

```math
\alpha_T=\frac{R_r}{L_r}
```

```math
\gamma=\alpha_i-\alpha_T
```

```math
k_1=\frac{b\alpha_T}{\alpha_T^2+\hat{\omega}_m^2}
```

```math
k_2=\frac{b\hat{\omega}_m}{\alpha_T^2+\hat{\omega}_m^2}
```

T形の物理ロータ磁束 $\hat{\phi}_R$ を状態にする場合、方式Cの更新式は以下になる。

```math
L_{\sigma T}\frac{d\hat{i}_{sd}}{dt}
=
\alpha_T\rho\hat{\phi}_R
-R_{\Sigma T}\hat{i}_{sd}
+\omega_sL_{\sigma T}\hat{i}_{sq}
+u_{sd}
+L_{\sigma T}(\gamma\tilde{i}_{sd}-\hat{\omega}_m\tilde{i}_{sq})
```

```math
L_{\sigma T}\frac{d\hat{i}_{sq}}{dt}
=
-\hat{\omega}_m\rho\hat{\phi}_R
-R_{\Sigma T}\hat{i}_{sq}
-\omega_sL_{\sigma T}\hat{i}_{sd}
+u_{sq}
+L_{\sigma T}(\gamma\tilde{i}_{sq}+\hat{\omega}_m\tilde{i}_{sd})
```

```math
\frac{d\hat{\phi}_R}{dt}
=
-\alpha_T\hat{\phi}_R
+\frac{R_rM}{L_r}\hat{i}_{sd}
+(k_1\alpha_i-\gamma)\frac{D}{M}\tilde{i}_{sd}
-(\omega_s-\hat{\omega}_m)\frac{D}{M}\tilde{i}_{sq}
```

```math
\omega_s
=
\hat{\omega}_m
+
\frac{
\frac{R_rM}{L_r}i_{sq}
+k_2\alpha_i\frac{D}{M}\tilde{i}_{sd}
-\gamma\frac{D}{M}\tilde{i}_{sq}
}{
\hat{\phi}_R-\frac{D}{M}\tilde{i}_{sd}
}
```

推定誤差が0なら、

```math
\omega_s-\hat{\omega}_m
=
\frac{(R_rM/L_r)i_{sq}}{\hat{\phi}_R}
```

となる。定常状態で $\hat{\phi}_R=Mi_{sd}$ なら、

```math
\omega_s-\hat{\omega}_m
=
\frac{R_r}{L_r}\frac{i_{sq}}{i_{sd}}
```

となり、通常のT形誘導機のすべり式に一致する。

電流推定誤差を以下で定義する。

```math
\tilde{i}_{sd}=i_{sd}-\hat{i}_{sd}
```

```math
\tilde{i}_{sq}=i_{sq}-\hat{i}_{sq}
```

方式Cの更新式は以下である。

```math
L_{\sigma}\frac{d\hat{i}_{sd}}{dt}
=
\alpha\hat{\psi}_R
-R_{\sigma}\hat{i}_{sd}
+\omega_sL_{\sigma}\hat{i}_{sq}
+u_{sd}
+L_{\sigma}(\gamma\tilde{i}_{sd}-\hat{\omega}_m\tilde{i}_{sq})
```

```math
L_{\sigma}\frac{d\hat{i}_{sq}}{dt}
=
-\hat{\omega}_m\hat{\psi}_R
-R_{\sigma}\hat{i}_{sq}
-\omega_sL_{\sigma}\hat{i}_{sd}
+u_{sq}
+L_{\sigma}(\gamma\tilde{i}_{sq}+\hat{\omega}_m\tilde{i}_{sd})
```

```math
\frac{d\hat{\psi}_R}{dt}
=
-\alpha\hat{\psi}_R
+R_R\hat{i}_{sd}
+(k_1\alpha_i-\gamma)L_{\sigma}\tilde{i}_{sd}
-(\omega_s-\hat{\omega}_m)L_{\sigma}\tilde{i}_{sq}
```

座標角速度 $\omega_s$ は、

```math
\omega_s
=
\hat{\omega}_m
+
\frac{
R_R i_{sq}
+k_2\alpha_iL_{\sigma}\tilde{i}_{sd}
-\gamma L_{\sigma}\tilde{i}_{sq}
}{
\hat{\psi}_R-L_{\sigma}\tilde{i}_{sd}
}
```

である。

ここで、

```math
\gamma=\alpha_i-\alpha
```

```math
k_1=\frac{b\alpha}{\alpha^2+\hat{\omega}_m^2}
```

```math
k_2=\frac{b\hat{\omega}_m}{\alpha^2+\hat{\omega}_m^2}
```

である。

### 12.1 方式Cのゲインが出てくる理屈

方式Cのゲインは、経験的に置いた補正ゲインではない。狙いは、推定電流誤差と推定磁束誤差の誤差方程式を、設計者が意味を持って指定できる形に変形することである。

目標は次の4つである。

1. 電流推定誤差を `alpha_i` の速さで減衰させる。
2. 回転座標で発生するd/q軸の混入を打ち消す。
3. 磁束推定誤差を `b` の速さで減衰させる。
4. 推定ロータ磁束q軸成分を0に保つ。

まず、電流推定誤差を以下で定義する。

```math
\tilde{i}_s=i_s-\hat{i}_s
```

方式Cでは、この電流誤差を概念的に

```math
\frac{d\tilde{i}_s}{dt}
\simeq
-\alpha_i\tilde{i}_s
```

のように減衰させたい。ただし、推定ロータ磁束座標は回転座標であるため、d軸電流誤差とq軸電流誤差は回転速度によって混ざる。そこで、電流補正は単なる比例補正ではなく、単位行列 `I` と90度回転行列 `J` を使って

```math
H_i
=
\gamma I+\hat{\omega}_mJ
```

と置く。ここで、

```math
\gamma=\alpha_i-\alpha
```

である。成分で書けば、

```math
H_i
=
\begin{bmatrix}
\gamma & -\hat{\omega}_m \\
\hat{\omega}_m & \gamma
\end{bmatrix}
```

となる。この形にすることで、`gamma` 成分が電流誤差の減衰率を調整し、`omega_m J` 成分が回転座標で生じるd/q軸混入を補償する。

次に、磁束側のゲインを考える。方式Cでは、磁束補正の基本形を

```math
K=k_1I+k_2J
```

と置く。ここで `k1` は電流誤差と同相の補正成分、`k2` は90度回転した補正成分である。

ロータ磁束の誤差ダイナミクスには、減衰項と回転項が含まれる。主要な形は

```math
\alpha I-\hat{\omega}_mJ
```

である。方式Cでは、この回転を含む項に `K` を掛けた結果が、設計者が指定した単純な減衰 `bI` になるようにしたい。したがって、次を満たすように `K` を選ぶ。

```math
(\alpha I-\hat{\omega}_mJ)(k_1I+k_2J)=bI
```

ここで、90度回転行列は

```math
J^2=-I
```

を満たす。左辺を展開すると、

```math
(\alpha I-\hat{\omega}_mJ)(k_1I+k_2J)
=
(\alpha k_1+\hat{\omega}_m k_2)I
+
(\alpha k_2-\hat{\omega}_m k_1)J
```

である。これを `bI` にしたいので、`I` 成分と `J` 成分を比較する。

```math
\alpha k_1+\hat{\omega}_m k_2=b
```

```math
\alpha k_2-\hat{\omega}_m k_1=0
```

2本目の式から、

```math
k_2=\frac{\hat{\omega}_m}{\alpha}k_1
```

である。これを1本目に代入すると、

```math
\alpha k_1
+
\hat{\omega}_m
\frac{\hat{\omega}_m}{\alpha}k_1
=
b
```

両辺に `alpha` を掛ける。

```math
(\alpha^2+\hat{\omega}_m^2)k_1=b\alpha
```

したがって、

```math
k_1=
\frac{b\alpha}{\alpha^2+\hat{\omega}_m^2}
```

である。また、

```math
k_2=
\frac{\hat{\omega}_m}{\alpha}k_1
=
\frac{b\hat{\omega}_m}{\alpha^2+\hat{\omega}_m^2}
```

となる。これが方式Cの `k1,k2` の由来である。

ここまでで、電流誤差モードの速さ `alpha_i` と、磁束誤差モードの速さ `b` が決まる。Appendix Aの式では、これを推定ロータ磁束座標の状態方程式へ代入して整理する。その結果、逆Gamma形磁束 `psi_R` を状態にした場合、磁束状態へ入る補正は次になる。

```math
\frac{d\hat{\psi}_R}{dt}
\supset
(k_1\alpha_i-\gamma)L_{\sigma}\tilde{i}_{sd}
-
(\omega_s-\hat{\omega}_m)L_{\sigma}\tilde{i}_{sq}
```

`k1 alpha_i` は、磁束誤差を `b` で減衰させるために必要な補正から来る。`-gamma` は、電流誤差を `alpha_i` で減衰させるためにすでに電流式へ入れた補正と整合を取るために現れる。つまり、磁束状態と電流状態を別々に適当に補正しているのではなく、誤差方程式全体が狙った形になるように、両者を同時に整理した結果である。

最後に、方式Cは推定ロータ磁束d軸座標なので、

```math
\hat{\psi}_{Rq}=0
```

が座標の定義である。したがって、

```math
\frac{d\hat{\psi}_{Rq}}{dt}=0
```

を満たす必要がある。このq軸拘束式を座標角速度 `omega_s` について解いたものが、方式Cのすべり周波数式である。

したがって、方式Cでは次の3つは切り離せない。

| 要素 | 役割 |
|---|---|
| 電流補正 `H_i` | 電流推定誤差を `alpha_i` で減衰させる |
| 磁束補正 `K=k1I+k2J` | 回転を含む磁束誤差を `b` で減衰させる |
| `omega_s` 計算式 | 推定ロータ磁束q軸成分を0に保つ |

つまり方式Cのゲインは、「推定電流だけを速くするゲイン」でも「磁束だけを滑らかにするゲイン」でもない。電流誤差、磁束誤差、推定ロータ磁束座標の拘束を同時に満たすように逆算された閉形式ゲインである。

### 12.2 T形一般状態方程式としての導出

ここでは、逆Gamma形という名前を使わず、一般的なT形等価回路の定数からSLED方式の状態方程式を導く。使う定数は以下である。

```math
L_s=L_{ls}+M
```

```math
L_r=L_{lr}+M
```

```math
D=L_sL_r-M^2
```

```math
L_{\sigma}=\frac{D}{L_r}
```

```math
\rho=\frac{M}{L_r}
```

```math
\alpha=\frac{R_r}{L_r}
```

```math
c_r=\frac{R_rM}{L_r}
```

```math
R_{\Sigma}=R_s+R_r\left(\frac{M}{L_r}\right)^2
```

ここで、$\phi_r$ はT形の物理二次磁束である。一次磁束は

```math
\psi_s=L_{\sigma}i_s+\rho\phi_r
```

と書ける。

推定二次磁束d軸座標では、推定二次磁束を

```math
\hat{\phi}_r=
\begin{bmatrix}
\hat{\phi} \\
0
\end{bmatrix}
```

と置く。座標角速度を $\omega_s$、オブザーバ内部で使う電気角の回転子速度を $\hat{\omega}_m$ とする。速度センサありの評価では、$\hat{\omega}_m$ は実測した電気角速度 $\omega_m$ に等しい。この座標上で、誤差補正を入れないT形誘導機の状態方程式は次になる。

```math
L_{\sigma}\dot{i}_{sd}
=
\alpha\rho\phi
-R_{\Sigma}i_{sd}
+\omega_sL_{\sigma}i_{sq}
+u_{sd}
```

```math
L_{\sigma}\dot{i}_{sq}
=
-\hat{\omega}_m\rho\phi
-R_{\Sigma}i_{sq}
-\omega_sL_{\sigma}i_{sd}
+u_{sq}
```

```math
\dot{\phi}
=
-\alpha\phi
+c_ri_{sd}
```

q軸二次磁束は0に拘束されるため、その微分も0である。

```math
0
=
c_ri_{sq}
-
(\omega_s-\hat{\omega}_m)\phi
```

したがって、誤差補正なしの通常すべり式は

```math
\omega_s
=
\hat{\omega}_m
+
\frac{c_ri_{sq}}{\phi}
```

である。定常状態で $\phi=Mi_{sd}$ なら、

```math
\omega_s-\hat{\omega}_m
=
\frac{R_r}{L_r}\frac{i_{sq}}{i_{sd}}
```

となる。

### 12.3 SLED方式のオブザーバゲイン要素

SLED方式では、電流推定誤差を

```math
\tilde{i}_{sd}=i_{sd}-\hat{i}_{sd}
```

```math
\tilde{i}_{sq}=i_{sq}-\hat{i}_{sq}
```

と定義する。設計パラメータは、電流推定誤差の減衰率 $\alpha_i$ と、磁束推定誤差の減衰を決める $b$ である。

```math
\gamma=\alpha_i-\alpha
```

```math
k_1=\frac{b\alpha}{\alpha^2+\hat{\omega}_m^2}
```

```math
k_2=\frac{b\hat{\omega}_m}{\alpha^2+\hat{\omega}_m^2}
```

T形物理二次磁束 $\hat{\phi}$ を状態にしたSLED方式は、次の状態方程式になる。

```math
L_{\sigma}\dot{\hat{i}}_{sd}
=
\alpha\rho\hat{\phi}
-R_{\Sigma}\hat{i}_{sd}
+\omega_sL_{\sigma}\hat{i}_{sq}
+u_{sd}
+L_{\sigma}
\left(
\gamma\tilde{i}_{sd}
-\hat{\omega}_m\tilde{i}_{sq}
\right)
```

```math
L_{\sigma}\dot{\hat{i}}_{sq}
=
-\hat{\omega}_m\rho\hat{\phi}
-R_{\Sigma}\hat{i}_{sq}
-\omega_sL_{\sigma}\hat{i}_{sd}
+u_{sq}
+L_{\sigma}
\left(
\hat{\omega}_m\tilde{i}_{sd}
+\gamma\tilde{i}_{sq}
\right)
```

```math
\dot{\hat{\phi}}
=
-\alpha\hat{\phi}
+c_r\hat{i}_{sd}
+H_{\phi d}\tilde{i}_{sd}
+H_{\phi q}\tilde{i}_{sq}
```

ここで、磁束状態へ入るゲイン要素は

```math
H_{\phi d}
=
\left(k_1\alpha_i-\gamma\right)\frac{D}{M}
```

```math
H_{\phi q}
=
-(\omega_s-\hat{\omega}_m)\frac{D}{M}
```

である。

電流状態へ入るゲインを、微分方程式の右辺、つまり $\dot{\hat{i}}$ に対するゲインとして書けば、

```math
H_i
=
\begin{bmatrix}
\gamma & -\hat{\omega}_m \\
\hat{\omega}_m & \gamma
\end{bmatrix}
```

である。一方、上の状態方程式のように $L_{\sigma}\dot{\hat{i}}$ の右辺へ入れる形で書けば、

```math
L_{\sigma}H_i
=
L_{\sigma}
\begin{bmatrix}
\gamma & -\hat{\omega}_m \\
\hat{\omega}_m & \gamma
\end{bmatrix}
```

である。したがって、T形物理二次磁束を状態にした場合のオブザーバゲイン要素は、次のように整理できる。

| 入る状態 | $\tilde{i}_{sd}$ 係数 | $\tilde{i}_{sq}$ 係数 |
|---|---:|---:|
| $\dot{\hat{i}}_{sd}$ | $\gamma$ | $-\hat{\omega}_m$ |
| $\dot{\hat{i}}_{sq}$ | $\hat{\omega}_m$ | $\gamma$ |
| $\dot{\hat{\phi}}$ | $(k_1\alpha_i-\gamma)D/M$ | $-(\omega_s-\hat{\omega}_m)D/M$ |

さらに、$\hat{\phi}_q=0$ を保つため、q軸二次磁束微分を0にする。q軸側は独立した状態として積分しないが、拘束式の中には次の補正係数が入る。

| 拘束式に入る項 | $\tilde{i}_{sd}$ 係数 | $\tilde{i}_{sq}$ 係数 |
|---|---:|---:|
| $\dot{\hat{\phi}}_q=0$ | $k_2\alpha_iD/M$ | $-\gamma D/M$ |

したがって、q軸拘束式は次になる。

```math
0
=
c_ri_{sq}
+k_2\alpha_i\frac{D}{M}\tilde{i}_{sd}
-\gamma\frac{D}{M}\tilde{i}_{sq}
-(\omega_s-\hat{\omega}_m)
\left(
\hat{\phi}-\frac{D}{M}\tilde{i}_{sd}
\right)
```

これを $\omega_s$ について解くと、

```math
\omega_s
=
\hat{\omega}_m
+
\frac{
c_ri_{sq}
+k_2\alpha_i(D/M)\tilde{i}_{sd}
-\gamma(D/M)\tilde{i}_{sq}
}{
\hat{\phi}-(D/M)\tilde{i}_{sd}
}
```

となる。ここまでの式は、逆Gamma形の $L_M,R_R,R_{\sigma}$ を使わず、T形の $R_s,R_r,L_{ls},L_{lr},M$ だけで書かれている。SLED Appendix Aの式と同じ内容だが、状態をスケーリング済み磁束 $\psi_R$ ではなく、物理二次磁束 $\phi$ として書いた形である。

### 12.4 なぜこのすべり周波数式になるのか

方式Cの核心は、推定ロータ磁束座標を使うことである。この座標では、

```math
\hat{\psi}_{Rq}=0
```

でなければならない。さらに、常に0に保ちたいので、時間微分も0でなければならない。

```math
\frac{d\hat{\psi}_{Rq}}{dt}=0
```

Appendix Aの構成では、q軸ロータ磁束の微分は以下の形に整理される。

```math
\frac{d\hat{\psi}_{Rq}}{dt}
=
R_R i_{sq}
+k_2\alpha_iL_{\sigma}\tilde{i}_{sd}
-\gamma L_{\sigma}\tilde{i}_{sq}
-(\omega_s-\hat{\omega}_m)(\hat{\psi}_R-L_{\sigma}\tilde{i}_{sd})
```

これを0にしたいので、

```math
0
=
R_R i_{sq}
+k_2\alpha_iL_{\sigma}\tilde{i}_{sd}
-\gamma L_{\sigma}\tilde{i}_{sq}
-(\omega_s-\hat{\omega}_m)(\hat{\psi}_R-L_{\sigma}\tilde{i}_{sd})
```

最後の項を右辺へ移す。

```math
(\omega_s-\hat{\omega}_m)(\hat{\psi}_R-L_{\sigma}\tilde{i}_{sd})
=
R_R i_{sq}
+k_2\alpha_iL_{\sigma}\tilde{i}_{sd}
-\gamma L_{\sigma}\tilde{i}_{sq}
```

両辺を $\hat{\psi}_R-L_{\sigma}\tilde{i}_{sd}$ で割る。

```math
\omega_s-\hat{\omega}_m
=
\frac{
R_R i_{sq}
+k_2\alpha_iL_{\sigma}\tilde{i}_{sd}
-\gamma L_{\sigma}\tilde{i}_{sq}
}{
\hat{\psi}_R-L_{\sigma}\tilde{i}_{sd}
}
```

最後に $\hat{\omega}_m$ を右辺へ移す。

```math
\omega_s
=
\hat{\omega}_m
+
\frac{
R_R i_{sq}
+k_2\alpha_iL_{\sigma}\tilde{i}_{sd}
-\gamma L_{\sigma}\tilde{i}_{sq}
}{
\hat{\psi}_R-L_{\sigma}\tilde{i}_{sd}
}
```

これが方式Cの $\omega_s$ 計算式である。つまり、この式は後付けのすべり推定ではない。$\hat{\psi}_{Rq}=0$ という座標拘束を満たすための条件である。

電気角のロータ速度を

```math
\hat{\omega}_m=p\omega_m
```

と見なすと、すべり周波数は

```math
\omega_{\mathrm{slip}}=\omega_s-\hat{\omega}_m
```

である。したがって、方式Cでは

```math
\omega_{\mathrm{slip}}
=
\frac{
R_R i_{sq}
+k_2\alpha_iL_{\sigma}\tilde{i}_{sd}
-\gamma L_{\sigma}\tilde{i}_{sq}
}{
\hat{\psi}_R-L_{\sigma}\tilde{i}_{sd}
}
```

となる。

推定誤差が0の場合、

```math
\tilde{i}_{sd}=0
```

```math
\tilde{i}_{sq}=0
```

である。すると、

```math
\omega_{\mathrm{slip}}
=
\frac{R_R i_{sq}}{\hat{\psi}_R}
```

となる。これはロータ磁束基準の通常のすべり式に対応する。

## 13. 方式Cのゲイン設計パラメータ

方式Cでは、方式Aのような4行2列の一般行列 $H$ を毎周期計算しない。また、方式Bのような $K=k_1I+k_2J$ による二次磁束補正でもない。代わりに、推定ロータ磁束座標で整理された閉形式の係数を使う。

設計パラメータは主に以下である。

| 記号 | 意味 |
|---|---|
| $\alpha_i$ | 電流推定誤差の減衰率 |
| $b$ | 磁束推定誤差の減衰率 |
| $\gamma$ | 電流誤差注入の補助係数 |
| $k_1,k_2$ | 磁束誤差補正をd/q成分へ分配する係数 |

まず、

```math
\gamma=\alpha_i-\alpha
```

を計算する。ここで $\alpha=R_R/L_M$ はモータ定数で決まる値である。したがって、$\alpha_i$ を大きくすると、電流推定誤差をより速く消そうとする。

次に、

```math
k_1=\frac{b\alpha}{\alpha^2+\hat{\omega}_m^2}
```

```math
k_2=\frac{b\hat{\omega}_m}{\alpha^2+\hat{\omega}_m^2}
```

を計算する。

この2つの式の分母は同じである。

```math
\alpha^2+\hat{\omega}_m^2
```

低速では $\hat{\omega}_m$ が小さいため、分母は主に $\alpha^2$ で決まる。高速では $\hat{\omega}_m^2$ が大きくなるため、速度に応じて $k_1,k_2$ の比率が変わる。

標準実装では、

```math
b=2\zeta_{\infty}|\hat{\omega}_m|+\alpha
```

を使う。高速になるほど $|\hat{\omega}_m|$ が大きくなり、$b$ も大きくなる。これは高速域で磁束推定モードの減衰を確保するためである。

## 14. 実装時の計算順序

方式Aの1周期の処理は以下である。

1. APIから $R_s,R_r,L_{ls},L_{lr},M,T_s$ を取得する。
2. $L_s,L_r,D$ を計算する。
3. 外部から与えられた $\omega_m,\omega_{\mathrm{slip}}$ から、$\omega_r=p\omega_m$、$\omega_k=\omega_r+\omega_{\mathrm{slip}}$ を作る。
4. 現在の $\omega_k,\omega_r$ で $A,C$ を作る。
5. 目標極から $F$ を作る。
6. Sylvester方程式 $TA-FT=GC$ を解く。
7. $H=T^{-1}G$ を計算する。
8. 測定電流と推定電流の差 $y-C\hat{x}$ を計算する。
9. $\dot{\hat{x}}=A\hat{x}+Bv_s+H(y-C\hat{x})$ を計算する。
10. $\hat{x}$ を $T_s$ で積分する。

方式Bの1周期の処理は以下である。

1. APIから $R_s,R_r,L_{ls},L_{lr},M,T_s$ を取得する。
2. $L_s,L_r,L_{\sigma s}$ を計算する。
3. 外部から与えられた $\omega_m,\omega_{\mathrm{slip}}$ から、$\omega_r=p\omega_m$、$\omega_k=\omega_r+\omega_{\mathrm{slip}}$ を作る。
4. $\alpha,\beta,\omega_r,\omega_k$ から論文5.3節の回転dq座標版の $k_1,k_2$ を計算する。
5. $K=k_1I+k_2J$ を作る。
6. 測定一次電流の差分から $\dot{i}_s$ を計算する。
7. 2行2列連立方程式を解いて $\dot{\hat{\psi}}_r$ を求める。
8. $\hat{\psi}_r$ を $T_s$ で積分する。
9. $\hat{\psi}_s=L_{\sigma s}i_s+(M/L_r)\hat{\psi}_r$ で一次磁束を再構成する。
10. 磁束と電流の関係式から二次電流推定値を計算する。

方式Cの1周期の処理は以下である。

1. APIから $R_s,R_r,L_{ls},L_{lr},M,T_s$ を取得する。
2. 逆Gamma形定数 $L_{\sigma},L_M,R_R,R_{\sigma},\alpha$ を計算する。
3. 電流推定誤差 $\tilde{i}_{sd},\tilde{i}_{sq}$ を計算する。
4. $\alpha_i,b$ から $\gamma,k_1,k_2$ を計算する。
5. $\hat{\psi}_{Rq}=0$ の拘束から $\omega_s$ を計算する。
6. $\hat{i}_{sd},\hat{i}_{sq},\hat{\psi}_R$ の微分を計算する。
7. 状態を $T_s$ で積分する。
8. $\omega_s$ を積分して推定ロータ磁束角を更新する。

方式Aと方式Bでは、すべり周波数は外部から来る入力である。方式Cでは、すべり周波数はオブザーバ内部の拘束式から決まる。この違いを混同してはいけない。

## 15. パラメータ誤差がある場合の見方

ここまでの誤差方程式

```math
\dot{\tilde{x}}=(A-HC)\tilde{x}
```

は、真値モデルとオブザーバ内部モデルが一致している理想条件で成り立つ。

実際には、抵抗、インダクタンス、電圧、電流に誤差がある。すると、誤差方程式は概念的に以下の形になる。

```math
\dot{\tilde{x}}=(A-HC)\tilde{x}+d(t)
```

ここで $d(t)$ は、定数誤差、電圧誤差、電流誤差による外乱項である。

この式から分かることは2つある。

1つ目は、$A-HC$ が安定でないと、外乱がなくても誤差が発散することである。

2つ目は、$A-HC$ が安定でも、外乱 $d(t)$ があると推定誤差は完全には0にならない場合があることである。この場合、誤差は一定値または小さな振動として残る。

したがって、オブザーバ評価では以下を両方見る必要がある。

| 評価項目 | 見ているもの |
|---|---|
| 無誤差で収束するか | 設計した誤差ダイナミクスが安定か |
| 定数誤差でどれだけずれるか | 実機誤差に対する推定精度 |
| 電圧誤差でどれだけずれるか | インバータ電圧モデル誤差への感度 |
| 電流誤差でどれだけ揺れるか | センサノイズへの感度 |

## 16. 設計時の実務的な注意点

オブザーバ極を速くしすぎると、推定は速くなるが、測定ノイズや離散化誤差に敏感になる。したがって、極は速ければよいわけではない。

電流制御帯域を $\omega_{cc}$ とすると、オブザーバ帯域 $\omega_o$ は電流制御より十分速くしたくなる。しかし、実装周期、電圧誤差、センサノイズを考えると、むやみに大きくできない。

方式Aの実務上の判断は以下である。

| 設計判断 | 影響 |
|---|---|
| 極を左へ置く | 収束は速くなるが、ノイズ感度と離散化誤差が増える |
| 極を原点に近づける | 収束は遅くなるが、ノイズには鈍感になる |
| 極を重複させる | 数学的には可能でも、Sylvester方程式が特異になりやすい |
| 共役極にする | 減衰しながら回る誤差を表現できる |

方式Bの実務上の判断は以下である。

| 設計判断 | 影響 |
|---|---|
| $\alpha$ を大きくする | 二次磁束推定は速くなるが、電圧誤差と電流微分ノイズに敏感になる |
| $\beta$ を大きくする | 誤差の回り込みを速く指定できるが、高速域では $k_2$ の寄与が大きくなる |
| 電流微分を粗く扱う | $L_{\sigma s}K\dot{i}_s$ 項がノイズ源になり、磁束推定が荒れる |
| 低速・低磁束で分母保護を入れる | 数値暴走は防げるが、論文式そのものからはわずかに外れる |

方式Cの実務上の判断は以下である。

| 設計判断 | 影響 |
|---|---|
| $\alpha_i$ を大きくする | 電流推定は速くなるが、電流ノイズに敏感になる |
| $b$ を大きくする | 磁束推定は速くなるが、電圧誤差や離散化に敏感になる |
| 分母下限を大きくする | 低磁束時の数値暴走は防げるが、厳密な拘束式からは離れる |
| 推定角の積分を雑にする | dq変換の位相誤差となり、推定電流・磁束がずれる |

## 17. 文献調査から見たゲイン設計法の系譜

ここでは、同一次元磁束オブザーバのゲイン設計法を、文献調査に基づいて整理する。目的は、単に論文名を並べることではなく、組込み機器に載せるときにどの設計思想を選ぶべきかを見通せるようにすることである。

今回の調査対象は、次の問いに直接関係する文献である。

1. 誘導機のフルオーダ、または同一次元に近い磁束オブザーバであること。
2. オブザーバゲイン `H`、またはそれに相当するフィードバックゲインの設計法を扱うこと。
3. 速度でプラント行列が変わる問題を扱うこと。
4. 組込み機器で実行できる軽いオンライン計算、またはLUT化できる設計法につながること。
5. 速度センサレスや低速回生安定性の議論に関係すること。

一方、EKF、SMO、MRAS、DREMなどは誘導機磁束・速度推定として重要だが、今回の主題である「同一次元磁束オブザーバの `H` 設計」とは設計問題が異なる。そのため、これらは本命候補ではなく、周辺技術として扱う。

同一次元オブザーバの基本形は、これまで説明した

```math
\dot{\hat{x}}=A(\omega)\hat{x}+Bu+H(y-C\hat{x})
```

である。推定誤差は

```math
\dot{\tilde{x}}=(A(\omega)-HC)\tilde{x}
```

に従う。ここで重要なのは、誘導機では $A(\omega)$ が速度で変わることである。したがって、ある速度で望みどおりの極を置く $H$ を設計しても、速度が変われば

```math
A(\omega)-HC
```

の極も変わる。すべての速度で同じ減衰や同じ安定余裕を得たいなら、本来は

```math
H=H(\omega)
```

でなければならない。

しかし、一般の4状態オブザーバで毎制御周期に極配置計算を行うと、4行4列の線形代数計算が必要になる。制御周期が100 us程度の組込み機器では、この計算を毎周期行うのは重い。したがって文献上の主流は、以下のどれかである。

1. 速度依存性を単純なスカラー関数へ押し込む。
2. ゲインをオフラインで計算し、オンラインではテーブル参照にする。
3. 安定性が証明できる構造化ゲインを使う。
4. 閉形式の式を導出し、オンラインでは四則演算だけにする。

今回確認した主要文献は以下である。`本文確認` は、手元PDFまたはオンライン本文で構成・式を確認できたものを示す。`書誌確認` は、DOI・会議・巻号・タイトルまで確認したが、本文式までは未確認のものを示す。

| 系統 | 文献 | 確認状況 | この調査での意味 |
|---|---|---|---|
| 国内制御理論 | 堀, Cotter, 茅 1986 | 本文確認 | 電流/電圧モデル修正、Gopinath型、極配置を低次元化する基礎 |
| 古典FOO | Verghese/Sanders 1988 | 書誌確認 | 誘導機磁束オブザーバの古典的基礎 |
| パラメータ適応二次磁束オブザーバ | Kubota/Matsuse 1991 | 本文確認 | 誘導機固有極を `k` 倍する構造化オブザーバゲイン設計 |
| 古典速度適応FOO | Kubota/Matsuse/Nakano 1991/1993 | 書誌確認 | DSP実装を意識した速度適応FOOの出発点 |
| 古典速度適応FOO | Yang/Chin 1993 | 書誌確認 | Kubota系と同時期の速度同定・適応オブザーバ |
| FOO解析 | Hinkkanen/Luomi 2002/2004 | 書誌確認 | フルオーダ磁束オブザーバの解析・設計を正面から扱う |
| 回生安定化 | Hinkkanen/Luomi 2003/2004 | 書誌確認 | 低速回生不安定の問題設定に直結 |
| 比較研究 | Hinkkanen/Luomi 2006 | 書誌確認 | adaptive observer と inherently sensorless observer の比較 |
| 正実性 | Sangwongwanichほか 2007 | 書誌確認 | 速度推定系を正実性で設計する流れ |
| 完全安定性 | Harnefors 2007, Harnefors/Hinkkanen 2008 | 書誌確認 | 速度適応まで含む安定性解析の重要文献 |
| ゲインスケジューリング | Qu/Hinkkanen/Harnefors 2014 | 書誌確認 | オフライン設計、オンラインLUT更新の本命候補 |
| 実務的ゲイン設計 | IET 2014, ECTI-CON 2015, ICEMS 2022 | 書誌確認 | 速度収束率、簡易ゲイン、最適化でFOOゲインを設計する流れ |
| LQR/LUT | Kullick/Hackl 2018 | 要旨確認 | 重い設計をオフライン化し、オンラインはゲインスケジューリング |
| 閉形式FOO | Tiitinen/Hinkkanen/Harnefors 2023 | 本文確認 | フルオーダ速度適応オブザーバを閉形式で再設計する近年の有力解 |
| 低速回生改善 | TTE 2024, TIE 2026 | 書誌確認 | adaptive FOOの低速回生、Rsロバスト性改善の最近の流れ |
| 観測可能性 | Koteich/Duc/Maloum/Sandou 2016 | 要旨確認 | センサレスACドライブがそもそも観測困難になる条件を扱う |
| 現代適応 | Pyrkin/Bobtsov/Ortegaほか 2020 | 要旨確認 | DREM系。主題とは違うが、未知Rr・負荷を含む現代的適応推定 |

この表から分かるように、主流は一つではない。大きく分けると、次の8系統である。

1. 電流モデル/電圧モデル修正、Gopinath型、極配置を低次元化する国内・古典制御理論系。
2. DSP実装を意識した古典的な速度適応フルオーダオブザーバ系。
3. フルオーダオブザーバの解析・設計と低速回生安定化を扱うHinkkanen/Luomi系。
4. 正実性、受動性、完全安定性を使って速度適応まで保証しようとする系。
5. ゲインをオフラインで設計し、オンラインでは速度・動作点でスケジューリングする系。
6. 速度収束率、簡易式、多目的最適化など、実務的にフィードバックゲインを選ぶ系。
7. SLED 2023のように、閉形式でフルオーダオブザーバゲインを与える系。
8. DREM、SMO、EKF、観測可能性解析など、同じ問題意識を持つが `H` 設計とは別枠の周辺系。

### 17.0 ゲイン設計法の体系図

誘導機磁束オブザーバのゲイン設計法は、まず「オブザーバ誤差系の極を設計入力として直接指定するか」で大きく分けると見通しがよい。Riccati方程式、Kalman、LQR、正実性、SLED 2023型の閉形式安定化ゲインは、結果として極を動かすが、設計入力そのものは極ではない。したがって、ここでは非極配置系に分類する。

```text
誘導機磁束オブザーバのゲイン設計
|
+-- A. 極配置系
|   |
|   +-- A1. 直接極配置
|   |   |
|   |   +-- 4状態 Luenberger / Sylvester / place
|   |   |   - A - H C の極を直接指定する
|   |   |   - 任意の4極を置けるが、毎周期オンライン計算は重い
|   |   |   - 本資料の方式A
|   |   |
|   |   +-- 堀 5.3節型
|   |       - 二次磁束オブザーバの極を指定する
|   |       - k1, k2 を閉形式で計算する
|   |       - 本資料の方式B
|   |
|   +-- A2. 制約付き極配置
|       |
|       +-- Kubota / Matsuse k倍方式
|           - 誘導機固有極を k 倍する
|           - 任意極配置ではなく、誘導機固有極に対する倍率指定
|           - 計算が軽いが、k の選定が重要
|           - 本資料の方式D1
|
+-- B. 非極配置系
    |
    +-- B1. 正実性・超安定論ベース
    |   |
    |   +-- Kinpara / Koyama Riccati方式
    |   |   - Gw(s) を強正実に近づける
    |   |   - 白色雑音外乱から状態誤差までの伝達ゲインを小さくする
    |   |   - Riccati方程式から H を求める
    |   |   - 極は結果として決まるが、極を直接指定しない
    |   |
    |   +-- Yaskawa / Takase ゲインスケジューリング方式
    |       - Popov超安定論 + Kalman-Yakubovich補題に基づく
    |       - 低速センサレス安定化を重視する
    |       - K1, K2 を速度でスケジューリングする
    |       - 本資料の方式E
    |
    +-- B2. モード分離・閉形式安定化
    |   |
    |   +-- SLED 2023 Appendix A方式
    |       - 電流推定誤差を alpha_i で減衰させる
    |       - 磁束推定誤差を b で減衰させる
    |       - 速度推定系と磁束推定系を分離しやすい形にする
    |       - 閉形式でゲインを計算する
    |       - 本資料の方式C
    |
    +-- B3. 最適推定・最適制御ベース
        |
        +-- Kalman filter / EKF
        |   - 雑音共分散から推定ゲインを決める
        |   - 非線形モデルやセンサ雑音を扱いやすい
        |   - 計算負荷とチューニング負荷が大きい
        |
        +-- LQR / LQG / LUT方式
            - オフラインでRiccati/LQR設計を行う
            - 実機では速度・動作点に応じてテーブル参照する
            - 高次元モデルやLCフィルタ込みの設計で有効
```

この分類で重要なのは、極配置系かどうかは「極を評価に使うか」ではなく、「極を設計入力として直接指定するか」で決まる点である。Riccati方式やSLED方式でも、設計後に `A - H C` 相当の極を計算して安定性や応答を評価することはできる。しかし、Riccati方式の設計入力は強正実性や伝達ゲイン最小化であり、SLED方式の設計入力は $\alpha_i,b,\zeta_{\infty}$ である。したがって、これらは非極配置系として扱う。

今回の実装・調査対象に対応させると、以下になる。

| 方式 | 分類 | 設計入力 | 備考 |
|---|---|---|---|
| 方式A | 極配置系 | 4個の目標極 | Sylvester方程式で一般の $H$ を求める |
| 方式B | 極配置系 | 二次磁束オブザーバの目標極 | 堀5.3節の $k_1,k_2$ |
| 方式C | 非極配置系 | $\alpha_i,b,\zeta_{\infty}$ | SLED 2023の閉形式安定化ゲイン |
| 方式D1 | 制約付き極配置系 | 誘導機固有極の倍率 `k` | Kubota/Matsuse型 |
| 方式E | 非極配置系 | $K_1,K_2$ の速度スケジュール | 安川/Takase型 |
| Kinpara/Koyama Riccati | 非極配置系 | 正実性、伝達ゲイン最小化、Riccati重み | Riccati方程式で $H$ を求める |
| Kalman/EKF | 非極配置系 | 雑音共分散 | 推定性能重視 |
| LQR/LQG/LUT | 非極配置系 | 評価関数、重み、動作点 | オフライン設計とテーブル化向き |

### 17.1 電気学会・堀系: 極配置を低次元化する流れ

国内文献でまず重要なのは、堀・Cotter・茅の電気学会論文誌B 1986年論文である。この論文は、誘導電動機の磁束オブザーバを制御理論の立場から整理し、電流モデル修正型、電圧モデル修正型、Gopinath型最小次元オブザーバを比較している。

この系統の重要点は、誘導機のd/q軸対称性を使って、一般の行列ゲインを小さなパラメータ数へ落とすことである。代表例が

```math
K=k_1I+k_2J
```

である。展開すると

```math
K=
\begin{bmatrix}
k_1 & -k_2\\
k_2 & k_1
\end{bmatrix}
```

となる。これは、同相補正 $k_1$ と90度回転補正 $k_2$ だけで、d/q軸の対称性を保った補正を作る考え方である。

この考え方の利点は明確である。毎周期に一般の極配置問題を解かなくても、速度と設計極から $k_1,k_2$ を計算すればよい。組込み実装では非常に軽い。一方で、この形式だけで任意の4状態同一次元オブザーバの全極を自由に置けるわけではない。どの状態を動的に持つか、どの誤差を補正するか、どのモデルを修正するかを含めて、オブザーバ構成とセットで使う必要がある。

このため、堀系の結論は以下である。

| 観点 | 評価 |
|---|---|
| 計算負荷 | 小さい。$k_1,k_2$ の計算に落とせる |
| 理論の分かりやすさ | 極配置と物理モデルの関係が見えやすい |
| 注意点 | 同一次元4状態の一般ゲイン $H$ を直接得る方法とは限らない |
| 組込み候補 | 軽量実装の候補。ただし現代的な安定性評価と組み合わせたい |

### 17.2 Kubota/Matsuse系: 誘導機固有極の `k` 倍配置

Kubota/Matsuseの1991年電気学会論文「誘導電動機のパラメータ適応二次磁束オブザーバの提案とその安定性」では、一次抵抗と二次抵抗のパラメータ適応を付加するため、二次磁束だけでなく一次電流も含めて推定する完全次元オブザーバを採用している。

この論文でのオブザーバゲイン設計法は、任意の4極を自由に配置する一般極配置ではない。論文中には、広い速度範囲でパラメータ誤差の影響を低く抑えるため、オブザーバの極を回転数に応じて変化させ、誘導電動機固有の極の `k` 倍に比例させる、と明記されている。

つまり、設計パラメータは主に比例定数 `k` である。

| `k` | 意味 |
|---:|---|
| `k < 1` | オブザーバ極を誘導機固有極より遅くする |
| `k = 1` | オブザーバ極を誘導機固有極と同じにする。このとき補正ゲインは0になる |
| `k > 1` | オブザーバ極を誘導機固有極より速くする |

論文では、`k=0.5, 1.0, 1.5` を比較し、`k` により一次抵抗誤差・二次時定数誤差に対するトルク誤差が変わることを示している。重要なのは、`k` を大きくすれば常に良いわけではない点である。論文の図4では、低速域で `k` を大きくすると二次時定数誤差に対するトルク誤差が大きくなる場合が示されている。そのため、Kubota/Matsuse 1991の結論は「ゲインだけでパラメータ誤差感度を完全には消せないので、Rs/Rrの適応機構を付加する」という流れになる。

#### 17.2.1 状態方程式

論文では固定座標上で、状態を一次電流と二次磁束に取る。

```math
x
=
\begin{bmatrix}
\boldsymbol{i}_s \\
\boldsymbol{\phi}_r
\end{bmatrix}
```

```math
\dot{x}
=
\begin{bmatrix}
A_{11} & A_{12} \\
A_{21} & A_{22}
\end{bmatrix}
x
+
\begin{bmatrix}
B_1 \\
0
\end{bmatrix}
\boldsymbol{v}_s
```

各ブロックは、2軸の単位行列 `I` と90度回転行列 `J` を使って次の形になる。

```math
A_{11}
=
a_{r11}I
```

```math
A_{12}
=
a_{r12}I+a_{i12}J
=
\frac{M}{\sigma L_s L_r}
\left(
\frac{1}{\tau_r}I-\omega_rJ
\right)
```

```math
A_{21}
=
a_{r21}I
=
\frac{M}{\tau_r}I
```

```math
A_{22}
=
a_{r22}I+a_{i22}J
=
-\frac{1}{\tau_r}I+\omega_rJ
```

```math
B_1
=
\frac{1}{\sigma L_s}I
```

ここで、

```math
\sigma
=
1-\frac{M^2}{L_sL_r}
```

```math
\tau_r
=
\frac{L_r}{R_r}
```

である。

さらに、論文では

```math
A_{22}=\alpha A_{12}
```

```math
\alpha
=
-\frac{\sigma L_sL_r}{M}
```

という関係を使ってゲイン計算を簡略化している。

#### 17.2.2 オブザーバ構成

完全次元オブザーバは次の形である。

```math
\dot{\hat{x}}
=
\hat{A}\hat{x}
+
\begin{bmatrix}
\hat{B}_1 \\
0
\end{bmatrix}
\boldsymbol{v}_s
+
G(\hat{\boldsymbol{i}}_s-\boldsymbol{i}_s)
```

補正に使う出力誤差は、推定一次電流と測定一次電流の差である。ここでは論文の符号に合わせて `hat{i}_s - i_s` と書く。この符号を逆に定義する実装では、以下の `G` の符号も反転する。ゲイン行列 `G` は4行2列であり、論文ではd/q対称性を保つため、次の構造に制限している。論文中の表記は次である。

```math
G
=
\begin{bmatrix}
g_1 & g_2 & g_3 & g_4 \\
-g_2 & g_1 & -g_4 & g_3
\end{bmatrix}^{T}
```

上式は論文の表記に合わせて転置付きで書いたもので、実装上の4行2列行列として書けば、

```math
G
=
\begin{bmatrix}
g_1 & -g_2 \\
g_2 & g_1 \\
g_3 & -g_4 \\
g_4 & g_3
\end{bmatrix}
```

である。上2行が一次電流推定値への補正、下2行が二次磁束推定値への補正である。

#### 17.2.3 `k` 倍極配置によるゲイン計算

論文の狙いは、オブザーバ誤差系の極を誘導機固有極の `k` 倍に置くことである。`A_{22}=alpha A_{12}` の関係を使うと、ゲインは次の簡略式で計算できる。

```math
g_1
=
(k-1)(a_{r11}+a_{r22})
```

```math
g_2
=
(k-1)a_{i22}
```

```math
g_3
=
(k^2-1)(-\alpha a_{r11}+a_{r21})
+
\alpha(k-1)(a_{r11}+a_{r22})
```

```math
g_4
=
\alpha(k-1)a_{i22}
```

この式が、Kubota/Matsuse 1991で使われているオブザーバゲイン設計法の核心である。

なお、Kubota/Matsuse 1991の本文では、この「誘導機固有極の `k` 倍配置」によるゲイン行列の導出元として、原島・近藤・橋本・大野・井上の1989年産業応用部門大会論文「誘導機の高性能トルク制御のための離散時間オブザーバの設計」を引用している。したがって、Kubota/Matsuse 1991は、この `k` 倍極配置ゲインをパラメータ適応二次磁束オブザーバへ組み込み、`k` とパラメータ誤差感度の関係を評価した文献として読むのが正確である。

設計手順としては次のようになる。

1. モータ定数 `Rs, Rr, Ls, Lr, M` と現在の電気角速度 `omega_r` から `sigma, tau_r` を計算する。
2. `a_r11, a_r12, a_i12, a_r21, a_r22, a_i22` を計算する。
3. 設計パラメータ `k` を選ぶ。
4. 上式で `g1..g4` を計算する。
5. `G` を構成し、一次電流推定誤差 `hat{i}_s - i_s` をオブザーバに注入する。

この方法の特徴は、オンラインで一般の4次極配置問題を解かないことである。速度 `omega_r` が変わるたびに係数 `a_i12, a_i22` が変わるが、ゲイン更新はスカラー式だけで済む。したがって、組込み実装にはかなり向いている。

ただし、配置できる極は「誘導機固有極の `k` 倍」という制約付きである。任意の4極を自由に指定する方式ではない。また、論文の主眼はパラメータ適応機構の安定性証明であり、`k` 倍ゲインだけでRs/Rr誤差に対して完全ロバストになるとはしていない。むしろ、`k` の選び方だけでは低速時のパラメータ誤差感度を消しきれないため、一次抵抗・二次抵抗の適応則が必要、という構成である。

#### 17.2.4 `g1..g4` の導出

ここでは、なぜ上の `g1..g4` になるかを示す。導出の本質は、2軸量を複素数として扱い、2次の特性多項式の係数を比較することである。

まず、2軸ベクトルに対する

```math
aI+bJ
```

は、複素数

```math
a+jb
```

による掛け算と同じである。したがって、論文の状態行列ブロックを複素スカラーで

```math
a_{11}=a_{r11}
```

```math
a_{12}=a_{r12}+ja_{i12}
```

```math
a_{21}=a_{r21}
```

```math
a_{22}=a_{r22}+ja_{i22}
```

と書ける。また、ゲインの上側ブロックと下側ブロックを

```math
g=g_1+jg_2
```

```math
h=g_3+jg_4
```

とおく。

論文の補正項は `+G(hat{i}_s-i_s)` であるため、推定誤差を

```math
\tilde{x}=\hat{x}-x
```

と取ると、誤差方程式は

```math
\dot{\tilde{x}}
=
(A+GC)\tilde{x}
```

になる。複素2状態で書けば、誤差系の行列は

```math
F
=
\begin{bmatrix}
a_{11}+g & a_{12} \\
a_{21}+h & a_{22}
\end{bmatrix}
```

である。

Kubota/Matsuseの設計条件は、オブザーバ極を誘導機固有極の `k` 倍にすることである。元の誘導機行列を

```math
A_c
=
\begin{bmatrix}
a_{11} & a_{12} \\
a_{21} & a_{22}
\end{bmatrix}
```

とすると、`A_c` の固有値を `lambda_1, lambda_2`、誤差系 `F` の固有値を `k lambda_1, k lambda_2` にしたい。2次行列では、固有値の和がトレース、積が行列式なので、条件は

```math
\mathrm{tr}(F)
=
k\mathrm{tr}(A_c)
```

```math
\det(F)
=
k^2\det(A_c)
```

である。

まずトレース条件から、

```math
a_{11}+g+a_{22}
=
k(a_{11}+a_{22})
```

となる。したがって、

```math
g
=
(k-1)(a_{11}+a_{22})
```

である。実部と虚部に分けると、

```math
g_1
=
(k-1)(a_{r11}+a_{r22})
```

```math
g_2
=
(k-1)a_{i22}
```

が得られる。

次に行列式条件を使う。

```math
\det(F)
=
(a_{11}+g)a_{22}-a_{12}(a_{21}+h)
```

```math
\det(A_c)
=
a_{11}a_{22}-a_{12}a_{21}
```

なので、

```math
(a_{11}+g)a_{22}-a_{12}(a_{21}+h)
=
k^2(a_{11}a_{22}-a_{12}a_{21})
```

である。これを `h` について解くと、

```math
h
=
(k^2-1)
\left(
a_{21}
-
a_{11}\frac{a_{22}}{a_{12}}
\right)
+
g\frac{a_{22}}{a_{12}}
```

となる。ここで論文が使っている関係

```math
A_{22}=\alpha A_{12}
```

すなわち複素表現では

```math
a_{22}=\alpha a_{12}
```

を代入すると、

```math
h
=
(k^2-1)(a_{21}-\alpha a_{11})
+
\alpha g
```

となる。さらに `a11=ar11`, `a21=ar21`, `g=(k-1)(ar11+ar22)+j(k-1)ai22` を代入すれば、

```math
h
=
(k^2-1)(-\alpha a_{r11}+a_{r21})
+
\alpha(k-1)(a_{r11}+a_{r22})
+
j\alpha(k-1)a_{i22}
```

である。したがって実部と虚部から、

```math
g_3
=
(k^2-1)(-\alpha a_{r11}+a_{r21})
+
\alpha(k-1)(a_{r11}+a_{r22})
```

```math
g_4
=
\alpha(k-1)a_{i22}
```

が得られる。

以上が、Kubota/Matsuse 1991の `g1..g4` の導出である。要するに、`g1,g2` はトレース、つまり極の和を `k` 倍にするために決まり、`g3,g4` は行列式、つまり極の積を `k^2` 倍にするために決まる。

#### 17.2.5 推定二次磁束d軸座標への移植

Kubota/Matsuse型の `k` 倍極配置は、推定二次磁束d軸座標にも移植できる。ポイントは、元論文のゲイン `G` が

```math
g_1I+g_2J
```

および

```math
g_3I+g_4J
```

というd/q対称な構造を持つことである。この形は座標回転と可換なので、固定座標で計算した `g1..g4` を、推定二次磁束d軸座標上でも同じ係数として使える。

ただし、推定二次磁束d軸座標では、

```math
\hat{\phi}_{rq}=0
```

を常に満たす必要がある。したがって、q軸二次磁束状態を独立に積分するのではなく、座標角速度 $\omega_k$ を代数的に決める。

ここでは、推定誤差をKubota/Matsuse論文と同じ符号で

```math
e_d=\hat{i}_{sd}-i_{sd}
```

```math
e_q=\hat{i}_{sq}-i_{sq}
```

と定義する。推定二次磁束d軸座標での状態は

```math
\hat{x}
=
\begin{bmatrix}
\hat{i}_{sd} \\
\hat{i}_{sq} \\
\hat{\phi}
\end{bmatrix}
```

である。$\hat{\phi}$ はT形の物理二次磁束のd軸成分であり、q軸成分は0に拘束する。

このとき、Kubota/Matsuse型を移植したオブザーバは次になる。

```math
\dot{\hat{i}}_{sd}
=
a_{r11}\hat{i}_{sd}
+\omega_k\hat{i}_{sq}
+a_{r12}\hat{\phi}
+b_1u_{sd}
+g_1e_d
-g_2e_q
```

```math
\dot{\hat{i}}_{sq}
=
a_{r11}\hat{i}_{sq}
-\omega_k\hat{i}_{sd}
+a_{i12}\hat{\phi}
+b_1u_{sq}
+g_2e_d
+g_1e_q
```

```math
\dot{\hat{\phi}}
=
a_{r21}\hat{i}_{sd}
+a_{r22}\hat{\phi}
+g_3e_d
-g_4e_q
```

座標角速度 $\omega_k$ は、$\dot{\hat{\phi}}_{rq}=0$ から決まる。q軸二次磁束の微分は

```math
0
=
a_{r21}\hat{i}_{sq}
+(a_{i22}-\omega_k)\hat{\phi}
+g_4e_d
+g_3e_q
```

である。したがって、

```math
\omega_k
=
a_{i22}
+
\frac{
a_{r21}\hat{i}_{sq}
+g_4e_d
+g_3e_q
}{
\hat{\phi}
}
```

となる。ここで $a_{i22}$ は電気角の回転子速度であり、通常は

```math
a_{i22}=\hat{\omega}_m
```

である。

推定誤差が0なら、

```math
\omega_k-\hat{\omega}_m
=
\frac{a_{r21}\hat{i}_{sq}}{\hat{\phi}}
```

となる。$a_{r21}=R_rM/L_r$ なので、

```math
\omega_k-\hat{\omega}_m
=
\frac{(R_rM/L_r)\hat{i}_{sq}}{\hat{\phi}}
```

である。定常状態で $\hat{\phi}=M\hat{i}_{sd}$ なら、

```math
\omega_k-\hat{\omega}_m
=
\frac{R_r}{L_r}\frac{\hat{i}_{sq}}{\hat{i}_{sd}}
```

となり、通常の二次磁束基準すべり式に一致する。

この移植方式の重要点は次である。

| 観点 | 内容 |
|---|---|
| ゲイン計算 | 元論文の `g1..g4` をそのまま使う |
| 座標 | 推定二次磁束をd軸に取り、$\hat{\phi}_{rq}=0$ とする |
| 追加される式 | $\dot{\hat{\phi}}_{rq}=0$ から $\omega_k$ を計算する |
| 状態数 | $\hat{i}_{sd},\hat{i}_{sq},\hat{\phi}$ の3状態 + 角度積分 |
| 注意点 | `g1..g4` は任意極配置ではなく、誘導機固有極の `k` 倍配置 |

SLED方式との違いは、$\omega_k$ の式とゲインの作り方にある。SLED方式では、$\omega_s$ の分母に

```math
\hat{\phi}-(D/M)\tilde{i}_{sd}
```

が現れ、電流推定誤差を含めたq軸磁束拘束になる。一方、Kubota/Matsuse移植方式では、元の完全次元オブザーバを座標変換して $\hat{\phi}_{rq}=0$ を課すため、基本形では分母は $\hat{\phi}$ になる。

#### 17.2.6 Kubota/Matsuse/Nakano 1993との関係

Kubota, Matsuse, Nakanoの1993年IEEE論文は、DSP上で動く速度適応磁束オブザーバとしてよく参照される古典である。1991年のIAS会議版も存在する。SLED 2023の序論でも、これらの初期研究は速度適応フルオーダオブザーバの出発点として位置付けられている。

ただし、1993年IEEE論文の本文式はこの環境では直接確認できていない。したがって、このドキュメントでは、添付PDFで本文確認できたKubota/Matsuse 1991の `k` 倍極配置ゲインをKubota/Matsuse系の具体的なゲイン設計法として扱う。1993年DSP論文側のゲイン設計が同一式か、または速度適応用に変更されているかは、本文確認後に分けて記載すべきである。

この点は重要である。Kubota/Matsuse 1991については、オブザーバゲイン設計法は明確に次である。

```math
\text{observer poles}
=
k \times
\text{induction-motor natural poles}
```

一方、Kubota/Matsuse/Nakano 1993 DSP論文については、現時点では書誌と引用関係のみ確認済みであり、ゲイン設計式の本文確認は未了である。

#### 17.2.7 速度適応FOOとしての一般構造

速度適応磁束オブザーバとして見ると、この系統で重要なのは、磁束オブザーバ単体ではなく、速度適応則まで含めた構成である。すなわち、オブザーバで一次電流と二次磁束を推定し、測定電流と推定電流の差から速度推定値を更新する。

後続文献で扱われる基本形は、次のような電流誤差注入型のフルオーダオブザーバである。

```math
\dot{\hat{x}}
=
A(\hat{\omega}_r)\hat{x}
+
Bv_s
+
H(i_s-\hat{i}_s)
```

ここで `H` または `G` がオブザーバゲインである。Kubota/Matsuse 1991では、これを一般極配置ではなく `k` 倍極配置により構成している。

速度適応則まで含める場合、推定速度が真値からずれると、推定磁束ベクトルの回転速度がずれる。その結果、推定磁束と実磁束の角度がずれ、一次電流の推定誤差に「推定磁束と直交する成分」が現れる。この直交成分を使うと、速度推定誤差の符号を取り出せる。

二軸ベクトルに対して90度回転行列を

```math
J
=
\begin{bmatrix}
0 & -1 \\
1 & 0
\end{bmatrix}
```

とすると、速度適応に使う誤差信号は、代表的には次の形で表せる。

```math
\varepsilon_{\omega}
=
\tilde{i}_s^{T}J\hat{\phi}_r
```

成分で書くと、

```math
\varepsilon_{\omega}
=
\tilde{i}_{sd}\hat{\phi}_{rq}
-
\tilde{i}_{sq}\hat{\phi}_{rd}
```

である。符号は、電流誤差を `i_s - hat{i}_s` と定義するか、`hat{i}_s - i_s` と定義するかで反転する。実装時には、この符号と速度適応PIの符号を必ず合わせる必要がある。

速度推定値は、この誤差信号をPIで処理して更新する。

```math
\dot{\xi}_{\omega}
=
\varepsilon_{\omega}
```

```math
\hat{\omega}_r
=
K_{P\omega}\varepsilon_{\omega}
+
K_{I\omega}\xi_{\omega}
```

つまり、Kubota/Matsuse/Nakano系の中核は次の組み合わせである。

| 構成要素 | 内容 |
|---|---|
| 推定状態 | 一次電流2成分と二次磁束2成分 |
| 補正入力 | 測定一次電流と推定一次電流の誤差 |
| オブザーバゲイン | 電流誤差を4状態へ戻す行列 `G` または `H` |
| 速度適応信号 | 電流誤差と推定二次磁束の外積相当量 |
| 速度推定 | 速度適応信号をPI処理して `hat{omega}_r` を更新 |
| 実装負荷 | 行列ベクトル演算とスカラーPIで済むためDSPに載せやすい |

ゲイン設計法として見ると、この系統では設計対象が二つに分かれる。一つはオブザーバゲイン `G` であり、Kubota/Matsuse 1991では誘導機固有極の `k` 倍配置で求める。もう一つは速度適応PIゲイン `K_{Pomega}, K_{Iomega}` であり、これは電流誤差から得た速度誤差信号をどれだけ速く `hat{omega}_r` に反映するかを決める。

このため、Kubota/Matsuse系は「スカラー `k` で速度依存ゲインを軽く計算できる古典的方式」として重要である。一方で、任意の4極配置ではなく、`k` 倍極配置という制約付きの方式である。速度適応PIまで含めると、`G` の選び方、速度適応PIの選び方、回生運転時の安定余裕が一体になって効くため、後続研究では正実性、完全安定性、ゲインスケジューリングが議論されている。

この方式が実用上魅力的だった理由は、毎制御周期で重い最適化やRiccati方程式を解かなくてもよい点である。主な処理は、4状態オブザーバの時間更新、2成分の電流誤差計算、速度適応PIである。したがって、1990年代のDSPでも実装を狙える構造だった。

一方で、この方式をそのまま現代の設計指針として使うには注意がある。SLED 2023の序論では、初期のフルオーダ速度適応オブザーバの安定性主張には、磁束推定誤差を無視するなど不完全な仮定が含まれており、完全な安定性保証には至っていなかったと整理されている。つまり、`A(hat{omega}_r)-HC` の極を左半平面に置くだけでは、速度適応ループまで含めた安定性を保証したことにならない。実際には、

- 電流推定誤差
- 磁束推定誤差
- 速度推定誤差
- 速度適応PI
- 回生時の動作点符号

が結合した非線形な誤差系になる。低速回生領域で速度適応オブザーバが不安定になりやすいことは、後続のHinkkanen/Luomi系、Harnefors/Hinkkanen系の研究で明確に問題として扱われている。

したがって、Kubota/Matsuse/Nakano系は「実装実績のある古典的な速度適応フルオーダオブザーバ」として重要である。ただし、今回の目的である「組込み機器で軽量に実装でき、かつ回生安定性を説明・保証しやすいゲイン設計法」を選ぶなら、この方式をそのまま最終候補にするのではなく、後続の安定性解析やゲインスケジューリング方式と比較する必要がある。

この章での位置づけは次の通りである。

| 観点 | 評価 |
|---|---|
| 計算負荷 | 小さい。`g1..g4` はスカラー式で計算できる |
| 理論の分かりやすさ | `k` がオブザーバ極の速度倍率なので意味は明確 |
| 安定性 | `k` 倍ゲインだけでパラメータ誤差を完全には消せない。論文では適応則を併用 |
| ゲイン設計の見通し | 任意極配置ではなく、誘導機固有極の比例配置 |
| 組込み候補 | 軽量な基準方式として有用。低速回生や速度適応を含む保証は後続研究と比較が必要 |

### 17.3 Hinkkanen/Luomi系: フルオーダオブザーバ解析と回生安定化

HinkkanenとLuomiの系統では、2002年IECONの "Analysis and design of full-order flux observers for sensorless induction motors"、2003年IEMDCの "Stabilization of the regenerating mode..."、2004年TIEの "Analysis and Design of Full-Order Flux Observers..." と "Stabilization of Regenerating-Mode Operation..." が重要である。

この系統は、古典的な速度適応FOOをそのまま使うのではなく、フルオーダ磁束オブザーバの安定性を動作点ごとに解析し、特に回生運転で不安定化しやすい点を設計問題として扱う。今回の議論で問題になっていた「回生側だけ不安定化しやすい」という現象に最も近い文献群である。

ここで見るべきポイントは二つある。

1. 電気系のフルオーダオブザーバの極だけでなく、速度適応ループまで含む結合系を見る必要がある。
2. 低速回生では、オブザーバゲインの選び方によって安定余裕が大きく変わる。

Hinkkanen/Luomi系は、今回の目的に対して「なぜ古典FOOでは不十分か」を説明する基盤になる。ただし、組込み機器でそのまま毎周期計算できる単純な閉形式ゲインというより、安定性解析と設計条件の整理に重点がある。

| 観点 | 評価 |
|---|---|
| 計算負荷 | 文献方式による。単純な四則演算ゲインとは限らない |
| 安定性 | 低速回生を直接扱うため重要 |
| 設計の見通し | 古典FOOより強いが、実装式だけを抜き出すには注意が必要 |
| 組込み候補 | 直接実装候補というより、安定性評価の基準 |

### 17.4 Harnefors/Hinkkanen/Sangwongwanich系: 正実性、受動性、完全安定性

Hinkkanen/Luomi系の後、Harneforsは2007年に速度適応オブザーバの大域安定性に関する短報を出し、HarneforsとHinkkanenは2008年に低次元および同一次元オブザーバの完全安定性を扱っている。さらに、Sangwongwanichらは2007年に正実性に基づく速度推定設計フレームワークを示している。

この系統のキーワードは、

- 回生運転
- 低速領域
- 正実性
- 受動性
- 完全安定性
- 速度適応ループを含む安定性

である。

ここでの重要な考え方は、単に `A-HC` の極を左半平面に置くだけでは不十分な場合がある、という点である。速度推定まで含めると、磁束推定誤差、電流推定誤差、速度推定誤差が結合する。したがって、

```math
A(\omega)-H(\omega)C
```

だけでなく、速度適応則まで含んだ拡大誤差系を見る必要がある。

この系統の結論は以下である。

| 観点 | 評価 |
|---|---|
| 計算負荷 | 設計式による。単純な閉形式とは限らない |
| 安定性 | 最も重視されている。回生・低速の説明に強い |
| 設計の見通し | 正実性や受動性の理解が必要で、初学者にはやや難しい |
| 組込み候補 | 安定性条件を満たすゲインを、テーブル化または閉形式化できれば有力 |

### 17.5 Qu/Hinkkanen/Harnefors系: ゲインスケジューリング

Qu, Hinkkanen, Harneforsの2014年論文は、フルオーダオブザーバのゲインスケジューリングを扱う。これは、今回の問題意識にかなり直接対応している。

考え方はシンプルである。

1. オフラインで、速度や動作点ごとに望ましいオブザーバゲインを計算する。
2. 組込み機器にはゲインテーブルを持たせる。
3. 実行時は速度や動作点に応じてテーブル参照、必要なら補間を行う。

この方法では、オンラインで重い極配置やRiccati方程式を解く必要がない。必要な計算は、テーブル参照と補間で済む。これは組込み実装として非常に現実的である。

一方で、欠点もある。テーブルの次元が増えるとメモリを食う。速度だけでなく、磁束、トルク、弱め界磁、電圧制限状態まで含めると、テーブル設計が複雑になる。また、テーブル外挿領域の安定性をどう保証するかも課題になる。

| 観点 | 評価 |
|---|---|
| 計算負荷 | 小さい。オンラインは参照と補間 |
| メモリ | テーブル次第で増える |
| 安定性 | オフライン設計点では確認しやすい。点間補間と外挿に注意 |
| 組込み候補 | 実務的には非常に有力 |

### 17.6 実務的フィードバックゲイン設計: 速度収束率、簡易式、最適化

Crossref調査では、2014年以降に adaptive full-order observer のフィードバックゲインを実務的に設計する文献群が複数見つかった。代表例は以下である。

| 文献 | 狙い |
|---|---|
| IET Electric Power Applications 2014 | 速度収束率に基づくフィードバックゲイン設計 |
| ECTI-CON 2015 | 実用的で単純なフィードバックゲイン設計 |
| ICPE-ECCE Asia 2015 | adaptive full-order observer のロバスト性改善 |
| IPEMC-ECCE Asia 2016 | 安定性と動特性改善 |
| ICEMS 2018 | 低速回生領域向けの速度適応スキーム |
| ICEMS 2022 | 多目的最適化に基づくフィードバックゲイン設計 |
| IEEE TTE 2024 | 低速回生モードの安定性改善 |
| IEEE TIE 2026 | 固定子抵抗ロバスト性改善のためのフィードバック行列設計 |

この系統は、Hinkkanen/Harnefors系ほど理論の一般性を前面に出すというより、実機で困る具体的な問題、例えば低速回生、速度収束率、固定子抵抗誤差、動特性を改善するために、`H` や速度適応則をどう選ぶかを扱う流れである。

現時点では本文式まで確認できていない文献が多いため、このドキュメントでは最終方式として断定しない。ただし、実務上は重要な枝である。特に2024年と2026年の文献は、今回の問題意識である「回生安定性」と「定数誤差ロバスト性」に直接関係するため、実装前に本文を確認する価値が高い。

| 観点 | 評価 |
|---|---|
| 計算負荷 | 方式による。簡易式なら小さいが、最適化系はオフライン向き |
| 安定性 | 低速回生やロバスト性の改善を狙うものが多い |
| 設計の見通し | 文献ごとに異なる。本文式の確認が必要 |
| 組込み候補 | 本文式が軽ければ候補。少なくとも比較対象に入れる価値あり |

### 17.7 LQR/LUT系: 重い設計をオフラインへ逃がす流れ

KullickとHacklの2018年論文は、LCフィルタ付き誘導機を対象に、オブザーバゲインと制御ゲインをLQRでオフライン設計し、オンラインではゲインスケジューリングで更新する構成を使っている。arXivの要旨でも、オブザーバゲインと制御ゲインをLQRでオフライン計算し、オンラインでゲインスケジューリング更新する、と明記されている。

この方式は、閉形式の美しさよりも実務性を重視している。設計時にRiccati方程式などの重い計算を行い、実機ではテーブル参照にする。オブザーバ対象がLCフィルタ込みで高次元になる場合でも、同じ考え方を使える。

この方式の価値は、「厳密な閉形式を探し続けるより、オフライン設計とLUT化で十分実用になる場合がある」ことを示している点である。

| 観点 | 評価 |
|---|---|
| 計算負荷 | オンラインは小さい |
| 設計自由度 | 高い。LQR重みで調整できる |
| 理論の見通し | 重みと極の関係が直感的でない場合がある |
| 組込み候補 | テーブル設計と検証プロセスを整備できるなら有力 |

### 17.8 Tiitinen/Hinkkanen/Harnefors 2023: 閉形式フルオーダオブザーバ

SLED 2023のTiitinen, Hinkkanen, Harnefors論文は、近年の有力な到達点である。この論文は、速度適応フルオーダオブザーバを再検討し、磁束推定ダイナミクスを速度推定から分離するゲインを提案している。

この方式の特徴は、フルオーダオブザーバのゲインを、少数の設計パラメータと速度の関数として閉形式で与える点である。論文では、フルオーダオブザーバを

```math
\frac{d\hat{\psi}_s}{dt}
=
-\omega_sJ\hat{\psi}_s
-R_s\hat{i}_s
+
u_s
+
K_{\psi}\tilde{i}_s
```

```math
L_{\sigma}\frac{d\hat{i}_s}{dt}
=
(\alpha I-\hat{\omega}_mJ)\hat{\psi}_s
-L_{\sigma}(\beta I+\hat{\omega}_rJ)\hat{i}_s
+
u_s
+
K_i\tilde{i}_s
```

の形で書き、ゲインを

```math
K_{\psi}=\alpha_iL_{\sigma}K-R_sI
```

```math
K_i=L_{\sigma}[(\alpha_i-\beta)I-\hat{\omega}_rJ]
```

のように与える。ここで `alpha_i` は電流推定誤差の減衰を決める設計パラメータであり、`K` は磁束推定モードを決める行列である。

この式の `beta` はSLED論文の記号で、逆Gamma形モデルの電流方程式に出てくる係数である。論文では

```math
\beta=\frac{R_s}{L_{\sigma}}+\omega_{rb}
```

```math
\omega_{rb}
=
\left(
\frac{1}{L_M}
+
\frac{1}{L_{\sigma}}
\right)R_R
```

と定義される。ここでも重要なのは、`K_psi` と `K_i` が一般の数値極配置で得られる任意行列ではなく、モータ定数、速度推定値、少数の設計パラメータから直接計算される構造化ゲインだという点である。

この方式の強い点は、ゲインの形が閉形式であり、オンラインで重い行列方程式を解かないことである。さらに、速度適応まで含む誤差ダイナミクスの特性多項式が整理されており、設計パラメータと安定性の関係が比較的見やすい。

組込み実装では、毎周期に必要なのは主に

- 速度に依存するスカラー係数の計算
- 2行2列の `I,J` 型行列の合成
- 電流誤差注入

である。これは一般の4行4列極配置よりかなり軽い。

| 観点 | 評価 |
|---|---|
| 計算負荷 | 小さい。閉形式ゲイン |
| 安定性 | 速度適応を含む解析がある |
| 設計の見通し | `alpha_i,b` などの設計パラメータに意味がある |
| 組込み候補 | 現時点で最有力候補 |

### 17.9 周辺系: 観測可能性、DREM、SMO、EKF

今回の主題からは外れるが、文献調査上は以下も無視できない。

Koteich, Duc, Maloum, Sandouの2016年の観測可能性解析は、センサレスACドライブがどの動作条件で情報を失うかを扱う。これは `H` の作り方そのものではないが、低速・特定動作点で「どれだけ良いオブザーバゲインを作っても推定が難しい」条件を理解するために重要である。

Pyrkin, Bobtsov, Ortegaらの2020年DREM系は、未知ロータ抵抗や負荷トルクを含む磁束・速度推定を扱う。これは同一次元FOOの `H` 設計とは違う理論体系だが、パラメータ誤差を含む推定保証という観点では参考になる。

SMOやEKFも産業応用では多い。ただし、SMOは不連続補正・境界層・チャタリング、EKFは共分散設計と計算負荷が主題になるため、今回の「回転dq座標の同一次元磁束オブザーバを、軽いゲイン計算で組込みに載せる」という要求とは別枠である。

### 17.10 どの方向を選ぶべきか

今回の目的は、同一次元磁束オブザーバを組込み機器に載せることである。そのため、重要な条件は以下である。

1. オンライン計算が軽い。
2. 速度変化に対してゲインを更新できる。
3. 回生・低速領域の安定性を説明できる。
4. 設計パラメータの意味を説明できる。
5. 実機でテストする前に、シミュレーションで安定余裕を評価できる。

この条件で見ると、方向性は以下になる。

| 候補 | 採用方針 |
|---|---|
| 毎周期Sylvester極配置 | 理解用・検証用。量産組込みには重い |
| 堀5.3型 `k1,k2` | 軽量な基準方式。国内文献の流れとして重要 |
| Kubota型古典速度適応FOO | 古典実装の基準。単独採用より後続安定化設計と比較 |
| Hinkkanen/Harnefors安定化設計 | 回生・低速安定性の理論基盤 |
| Qu型ゲインスケジューリング | 実務的候補。LUT化できるなら強い |
| 速度収束率/簡易ゲイン設計 | 本文確認が必要だが、軽量実装候補として調査継続 |
| LQR/LUT | 高次元モデルや周辺回路込みには有力 |
| SLED 2023閉形式 | 現時点の本命。閉形式、軽量、安定性説明のバランスが良い |

したがって、次に実装するなら、いきなり独自式を作るより、以下の3本立てがよい。

1. SLED 2023のフルオーダ閉形式ゲインを、論文本体の式に忠実に実装する。
2. Qu/Hinkkanen/Harnefors型のゲインスケジューリングを、比較用としてLUT実装する。
3. 2014年以降の速度収束率ベース/簡易ゲイン設計/低速回生改善系から、本文式を確認できたものを追加比較する。

この3つを比較すれば、閉形式方式、LUT方式、実務的簡易ゲイン方式のどれが実機組込みに向いているかを判断できる。堀5.3型は、国内文献に基づく軽量方式として比較軸に残す価値がある。

## 18. まとめ

磁束オブザーバの理解で重要なのは、以下の流れである。

1. 磁束と電流の関係式を作る。
2. 磁束から電流を計算できるように、インダクタンス行列を逆に解く。
3. 回転dq座標の磁束方程式に電流式を代入し、状態方程式 `dx/dt = Ax + Bv_s` を作る。
4. 測定一次電流を出力として `y = Cx` を作る。
5. オブザーバを `d x_hat/dt = A x_hat + Bv_s + H(y - C x_hat)` で構成する。
6. 誤差方程式 `d x_tilde/dt = (A - HC)x_tilde` を導く。
7. `A-HC` の固有値がすべて左半平面にあれば、誤差は0へ収束する。
8. 方式Aでは、Sylvester方程式 `TA - FT = GC` を使って `H = T^{-1}G` を求める。
9. 方式Bでは、論文5.3節の `K = k1 I + k2 J` を使って二次磁束を更新し、一次磁束は測定一次電流から再構成する。
10. 方式Cでは、推定ロータ磁束座標を使い、`psi_Rq_hat = 0` の拘束から `omega_s`、つまりすべり周波数相当の量が決まる。

方式Aは、状態空間と極配置の考え方を理解するのに適している。方式Bは、論文5.3節に沿った軽量な二次磁束オブザーバとして実装しやすい。方式Cは、組込み実装の計算負荷を下げる実用的な方式として有力である。ただし、方式Cではオブザーバ構成とすべり周波数計算式が一体であるため、`omega_s` の式を単なる外部すべり推定式として扱ってはいけない。

文献調査の観点では、量産組込みへ進める場合、SLED 2023の閉形式フルオーダオブザーバと、Qu/Hinkkanen/Harnefors型のゲインスケジューリングを中心に検討するのが自然である。堀5.3型は、国内文献に基づく軽量な比較方式として残すのがよい。

## 参考文献

1. 堀 洋一, Vincent Cotter, 茅 陽一, 「誘導電動機の磁束オブザーバに関する制御理論的考察」, 電気学会論文誌B, 106巻, 11号, pp.1001-1008, 1986. [J-STAGE PDF](https://www.jstage.jst.go.jp/article/ieejpes1972/106/11/106_11_1001/_pdf)
2. 久保田 寿夫, 松瀬 貢規, 「誘導電動機のパラメータ適応二次磁束オブザーバの提案とその安定性」, 電気学会論文誌D, 111巻, 3号, pp.188-194, 1991.
3. G. C. Verghese and S. R. Sanders, "Observers for Flux Estimation in Induction Machines", IEEE Transactions on Industrial Electronics, vol. 35, no. 1, pp. 85-94, 1988.
4. G. R. Slemon, "Modelling of Induction Machines for Electric Drives", IEEE Transactions on Industry Applications, vol. 25, no. 6, pp. 1126-1131, 1989.
5. H. Kubota, K. Matsuse, and T. Nakano, "DSP-Based Speed Adaptive Flux Observer of Induction Motor", Conference Record of the 1991 IEEE Industry Applications Society Annual Meeting, 1991. [DOI](https://doi.org/10.1109/IAS.1991.178183)
6. H. Kubota, K. Matsuse, and T. Nakano, "DSP-Based Speed Adaptive Flux Observer of Induction Motor", IEEE Transactions on Industry Applications, vol. 29, no. 2, pp. 344-348, 1993. [DOI](https://doi.org/10.1109/28.216542)
7. G. Yang and T.-H. Chin, "Adaptive-Speed Identification Scheme for a Vector-Controlled Speed Sensorless Inverter-Induction Motor Drive", IEEE Transactions on Industry Applications, vol. 29, no. 4, pp. 820-825, 1993. [DOI](https://doi.org/10.1109/28.232001)
8. J. Holtz, "The Representation of AC Machine Dynamics by Complex Signal Flow Graphs", IEEE Transactions on Industrial Electronics, vol. 42, no. 3, pp. 263-271, 1995.
9. M. Hinkkanen and J. Luomi, "Analysis and Design of Full-Order Flux Observers for Sensorless Induction Motors", IEEE Transactions on Industrial Electronics, 2004. [DOI](https://doi.org/10.1109/TIE.2004.834964)
10. M. Hinkkanen and J. Luomi, "Stabilization of Regenerating-Mode Operation in Sensorless Induction Motor Drives by Full-Order Flux Observer Design", IEEE Transactions on Industrial Electronics, vol. 51, no. 6, pp. 1318-1328, 2004. [DOI](https://doi.org/10.1109/TIE.2004.837902)
11. M. Hinkkanen and J. Luomi, "Comparative Study of Adaptive and Inherently Sensorless Observers for Variable-Speed Induction-Motor Drives", IEEE Transactions on Industrial Electronics, 2006. [DOI](https://doi.org/10.1109/TIE.2005.862314)
12. L. Harnefors, "Globally Stable Speed-Adaptive Observers for Sensorless Induction Motor Drives", IEEE Transactions on Industrial Electronics, vol. 54, no. 2, pp. 1243-1245, 2007. [DOI](https://doi.org/10.1109/TIE.2007.892729)
13. S. Sangwongwanich, S. Suwankawin, S. Po-ngam, and S. Koonlaboon, "A Unified Speed Estimation Design Framework for Sensorless AC Motor Drives Based on Positive-Real Property", PCC-Nagoya, pp. 1111-1118, 2007. [DOI](https://doi.org/10.1109/PCCON.2007.373105)
14. L. Harnefors and M. Hinkkanen, "Complete Stability of Reduced-Order and Full-Order Observers for Sensorless IM Drives", IEEE Transactions on Industrial Electronics, vol. 55, no. 3, pp. 1319-1329, 2008. [DOI](https://doi.org/10.1109/TIE.2007.909077)
15. M. Hinkkanen, L. Harnefors, and J. Luomi, "Reduced-Order Flux Observers With Stator-Resistance Adaptation for Speed-Sensorless Induction Motor Drives", IEEE Transactions on Power Electronics, vol. 25, no. 5, pp. 1173-1183, 2010. [DOI](https://doi.org/10.1109/TPEL.2009.2039650)
16. Z. Qu, M. Hinkkanen, and L. Harnefors, "Gain Scheduling of a Full-Order Observer for Sensorless Induction Motor Drives", IEEE Transactions on Industry Applications, vol. 50, no. 6, pp. 3834-3845, 2014. [DOI](https://doi.org/10.1109/TIA.2014.2323482)
17. Bin Chen, Ting Wang, Wenxi Yao, Kevin Lee, Zhengyu Lu, "Speed Convergence Rate-Based Feedback Gains Design of Adaptive Full-Order Observer in Sensorless Induction Motor Drives", IET Electric Power Applications, 2014. [DOI](https://doi.org/10.1049/iet-epa.2013.0210)
18. Hongbo Wang, Wei Sun, Yong Yu, Gaolin Wang, Dianguo Xu, "Robustness Improvement for Adaptive Full Order Observer in Sensorless Induction Motor Drives", ICPE-ECCE Asia, 2015. [DOI](https://doi.org/10.1109/ICPE.2015.7167961)
19. Apirach Rattanaudompisut, Sakorn Po-Ngam, "The Practical and Simple Feedback Gains Design of an Adaptive Full-Order Observer for Speed-Sensorless Induction Motor Drives", ECTI-CON, 2015. [DOI](https://doi.org/10.1109/ECTICON.2015.7207045)
20. Zhonggang Yin, Yanqing Zhang, Xiangqian Tong, Jing Liu, Yanru Zhong, "Stability and Dynamic Performance Improvement of Adaptive Full-Order Observers in Sensorless Induction Motor Drives", IPEMC-ECCE Asia, 2016. [DOI](https://doi.org/10.1109/IPEMC.2016.7512313)
21. Mohamad Koteich, Gilles Duc, Abdelmalek Maloum, Guillaume Sandou, "Observability of Sensorless Electric Drives", arXiv:1602.04468, 2016. [arXiv](https://arxiv.org/abs/1602.04468)
22. Julian Kullick and Christoph M. Hackl, "Speed-Sensorless State Feedback Control of Induction Machines With LC Filter", arXiv:1807.11799, 2018. [arXiv](https://arxiv.org/abs/1807.11799)
23. Cheng Luo, Bo Wang, Yong Yu, Zhixin Huo, Guoqiang Zhang, Dianguo Xu, "A Speed Adaptive Scheme-Based Full-Order Observer for Sensorless Induction Motor Drives in Low-Speed Regenerating Operation Range", ICEMS, 2018. [DOI](https://doi.org/10.23919/ICEMS.2018.8549249)
24. Anton Pyrkin, Alexey Bobtsov, Alexey Vedyakov, Romeo Ortega, Anastasiia Vediakova, Madina Sinetova, "A Flux and Speed Observer for Induction Motors with Unknown Rotor Resistance and Load Torque and no Persistent Excitation Requirement", arXiv:2009.00966, 2020. [arXiv](https://arxiv.org/abs/2009.00966)
25. "Stable Adaptive Estimation for Speed-Sensorless Induction Motor Drives: A Geometric Approach", ICEM, 2020. [DOI](https://doi.org/10.1109/ICEM49940.2020.9270926)
26. Lauri Tiitinen, Marko Hinkkanen, Lennart Harnefors, "Stable and Passive Observer-Based V/Hz Control for Induction Motors", ECCE, 2022. [DOI](https://doi.org/10.1109/ECCE50734.2022.9948057)
27. Ruhan Li, Yifei Zheng, Cheng Luo, Zhijie Xu, Kai Yang, "Multi-Objective Optimization Based Feedback Gains Design of Adaptive Full-Order Observer for Induction Motor Sensorless Drive", ICEMS, 2022. [DOI](https://doi.org/10.1109/ICEMS56177.2022.9983330)
28. Lauri Tiitinen, Marko Hinkkanen, Lennart Harnefors, "Speed-Adaptive Full-Order Observer Revisited: Closed-Form Design for Induction Motor Drives", IEEE SLED, 2023. [DOI](https://doi.org/10.1109/SLED57582.2023.10261359)
29. Hongwu Chen, Jian Li, Yang Lu, Kai Yang, Linghao Wu, Zhi Liu, "Stability Improvement of Adaptive Full-Order Observer for Sensorless Induction Motor Drives in Low-Speed Regenerating Mode", IEEE Transactions on Transportation Electrification, 2024. [DOI](https://doi.org/10.1109/TTE.2023.3257056)
30. Ying Wang, Xiancheng Huang, Bo Wang, Dianguo Xu, "Feedback Matrix Design of Adaptive Full-Order Observer for Stator Resistance Robustness Improvement in Speed Sensorless Induction Motor Drives", IEEE Transactions on Industrial Electronics, 2026. [DOI](https://doi.org/10.1109/TIE.2025.3647924)
31. 本リポジトリの実装: [c/flux_observer.c](c/flux_observer.c), [c/sled23_flux_observer.c](c/sled23_flux_observer.c).
