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

### 12.1 なぜこのすべり周波数式になるのか

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

## 17. まとめ

磁束オブザーバの理解で重要なのは、以下の流れである。

1. 磁束と電流の関係式を作る。
2. 磁束から電流を計算できるように、インダクタンス行列を逆に解く。
3. 回転dq座標の磁束方程式に電流式を代入し、状態方程式 $\dot{x}=Ax+Bv_s$ を作る。
4. 測定一次電流を出力として $y=Cx$ を作る。
5. オブザーバを $\dot{\hat{x}}=A\hat{x}+Bv_s+H(y-C\hat{x})$ で構成する。
6. 誤差方程式 $\dot{\tilde{x}}=(A-HC)\tilde{x}$ を導く。
7. $A-HC$ の固有値がすべて左半平面にあれば、誤差は0へ収束する。
8. 方式Aでは、Sylvester方程式 $TA-FT=GC$ を使って $H=T^{-1}G$ を求める。
9. 方式Bでは、論文5.3節の $K=k_1I+k_2J$ を使って二次磁束を更新し、一次磁束は測定一次電流から再構成する。
10. 方式Cでは、推定ロータ磁束座標を使い、$\hat{\psi}_{Rq}=0$ の拘束から $\omega_s$、つまりすべり周波数相当の量が決まる。

方式Aは、状態空間と極配置の考え方を理解するのに適している。方式Bは、論文5.3節に沿った軽量な二次磁束オブザーバとして実装しやすい。方式Cは、組込み実装の計算負荷を下げる実用的な方式として有力である。ただし、方式Cではオブザーバ構成とすべり周波数計算式が一体であるため、$\omega_s$ の式を単なる外部すべり推定式として扱ってはいけない。

## 参考文献

1. 堀 洋一, Vincent Cotter, 茅 陽一, 「誘導電動機の磁束オブザーバに関する制御理論的考察」, 電気学会論文誌B, 106巻, 11号, pp.1001-1008, 1986. [J-STAGE PDF](https://www.jstage.jst.go.jp/article/ieejpes1972/106/11/106_11_1001/_pdf)
2. Lauri Tiitinen, Marko Hinkkanen, Lennart Harnefors, "Speed-Adaptive Full-Order Observer Revisited: Closed-Form Design for Induction Motor Drives", 2023 IEEE International Symposium on Sensorless Control for Electrical Drives (SLED), 2023.
3. 本リポジトリの実装: [c/flux_observer.c](c/flux_observer.c), [c/sled23_flux_observer.c](c/sled23_flux_observer.c).
