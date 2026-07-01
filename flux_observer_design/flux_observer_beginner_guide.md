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

### 17.2 古典速度適応FOO: Kubota/Matsuse/NakanoとYang/Chin

Kubota, Matsuse, Nakanoの1993年論文は、DSP上で動く速度適応磁束オブザーバとしてよく参照される古典である。1991年のIAS会議版も存在する。Yang/Chinの1993年論文も、同時期に速度センサレス誘導機ドライブの適応速度同定を扱っている。SLED 2023の序論でも、これらの初期研究は速度適応フルオーダオブザーバの出発点として位置付けられている。

この系統で重要なのは、磁束オブザーバ単体ではなく、速度適応則まで含めた構成である。すなわち、オブザーバで一次電流と二次磁束を推定し、測定電流と推定電流の差から速度推定値を更新する。

Kubota/Matsuse/Nakano論文そのものは固定座標系の表現で参照されることが多いが、このドキュメントでは、以降の実装方針に合わせて回転dq座標上で同じ構造を書く。回転dq座標で見ると、同一次元フルオーダオブザーバの状態は例えば次の4状態になる。

```math
\hat{x}
=
\begin{bmatrix}
\hat{i}_{sd} &
\hat{i}_{sq} &
\hat{\phi}_{rd} &
\hat{\phi}_{rq}
\end{bmatrix}^{T}
```

ここで `hat` は推定値を表す。この状態を使うと、誘導機の電気モデルは概念的に次の状態空間形で書ける。

```math
\dot{x}
=
A(\omega_r)x
+
Bv_s
```

```math
i_s
=
Cx
```

ここで、`x` は実機の電流・磁束状態、`v_s` は固定子電圧、`i_s` は測定できる固定子電流、`omega_r` はロータ電気角速度である。速度センサレスの場合、`omega_r` は未知なので、オブザーバ内部では推定速度 `hat{omega}_r` を使う。

Kubota/Matsuse/Nakano型として後続文献で扱われる基本形は、次のような電流誤差注入型のフルオーダオブザーバである。

```math
\dot{\hat{x}}
=
A(\hat{\omega}_r)\hat{x}
+
Bv_s
+
H(i_s-\hat{i}_s)
```

```math
\hat{i}_s
=
C\hat{x}
```

ここで `H` がオブザーバゲインである。補正に使う信号は、測定可能な一次電流の推定誤差である。

```math
\tilde{i}_s
=
i_s-\hat{i}_s
```

この構成の物理的な意味は単純である。推定速度が真値からずれると、推定磁束ベクトルの回転速度がずれる。その結果、推定磁束と実磁束の角度がずれ、一次電流の推定誤差に「推定磁束と直交する成分」が現れる。この直交成分を使うと、速度推定誤差の符号を取り出せる。

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
| オブザーバゲイン | 電流誤差を4状態へ戻す行列 `H` |
| 速度適応信号 | 電流誤差と推定二次磁束の外積相当量 |
| 速度推定 | 速度適応信号をPI処理して `hat{omega}_r` を更新 |
| 実装負荷 | 行列ベクトル演算とスカラーPIで済むためDSPに載せやすい |

ゲイン設計法として見ると、この系統では設計対象が二つに分かれる。一つはオブザーバゲイン `H` であり、これは電流推定誤差を使って4状態の推定誤差をどれだけ速く減衰させるかを決める。もう一つは速度適応PIゲイン `K_{Pomega}, K_{Iomega}` であり、これは電流誤差から得た速度誤差信号をどれだけ速く `hat{omega}_r` に反映するかを決める。

このため、Kubota/Matsuse/Nakano系は「フルオーダオブザーバと速度適応則を組み合わせる基本構成」としては非常に重要である。一方で、今回欲しい「全動作点で所望の極を明示的に配置し、そのゲインを組込み機器で軽く計算する式」という意味では、Hori 5.3方式やSLED 2023方式ほど直接的ではない。特に、`H` の選び方、速度適応PIの選び方、回生運転時の安定余裕が一体になって効くため、後続研究では正実性、完全安定性、ゲインスケジューリングが議論されている。

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
| 計算負荷 | 小さい。4状態更新、電流誤差注入、速度適応PIが中心 |
| 理論の分かりやすさ | 構成は分かりやすいが、速度適応を含む安定性証明は別問題 |
| 安定性 | 古典的な基本構造は有用。ただし低速回生安定性は後続研究で補強が必要 |
| ゲイン設計の見通し | `H` と閉ループ極・安定余裕の関係は、SLED 2023方式ほど直接的ではない |
| 組込み候補 | 比較対象・基準実装として有用。単独採用より後続安定化設計と組み合わせたい |

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

ここでの重要な考え方は、単に $A-HC$ の極を左半平面に置くだけでは不十分な場合がある、という点である。速度推定まで含めると、磁束推定誤差、電流推定誤差、速度推定誤差が結合する。したがって、

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

のように与える。ここで $\alpha_i$ は電流推定誤差の減衰を決める設計パラメータであり、$K$ は磁束推定モードを決める行列である。

この式の $\beta$ はSLED論文の記号で、逆Gamma形モデルの電流方程式に出てくる係数である。論文では

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

と定義される。ここでも重要なのは、$K_{\psi}$ と $K_i$ が一般の数値極配置で得られる任意行列ではなく、モータ定数、速度推定値、少数の設計パラメータから直接計算される構造化ゲインだという点である。

この方式の強い点は、ゲインの形が閉形式であり、オンラインで重い行列方程式を解かないことである。さらに、速度適応まで含む誤差ダイナミクスの特性多項式が整理されており、設計パラメータと安定性の関係が比較的見やすい。

組込み実装では、毎周期に必要なのは主に

- 速度に依存するスカラー係数の計算
- 2行2列の $I,J$ 型行列の合成
- 電流誤差注入

である。これは一般の4行4列極配置よりかなり軽い。

| 観点 | 評価 |
|---|---|
| 計算負荷 | 小さい。閉形式ゲイン |
| 安定性 | 速度適応を含む解析がある |
| 設計の見通し | $\alpha_i,b$ などの設計パラメータに意味がある |
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
| 堀5.3型 $k_1,k_2$ | 軽量な基準方式。国内文献の流れとして重要 |
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
3. 回転dq座標の磁束方程式に電流式を代入し、状態方程式 $\dot{x}=Ax+Bv_s$ を作る。
4. 測定一次電流を出力として $y=Cx$ を作る。
5. オブザーバを $\dot{\hat{x}}=A\hat{x}+Bv_s+H(y-C\hat{x})$ で構成する。
6. 誤差方程式 $\dot{\tilde{x}}=(A-HC)\tilde{x}$ を導く。
7. $A-HC$ の固有値がすべて左半平面にあれば、誤差は0へ収束する。
8. 方式Aでは、Sylvester方程式 $TA-FT=GC$ を使って $H=T^{-1}G$ を求める。
9. 方式Bでは、論文5.3節の $K=k_1I+k_2J$ を使って二次磁束を更新し、一次磁束は測定一次電流から再構成する。
10. 方式Cでは、推定ロータ磁束座標を使い、$\hat{\psi}_{Rq}=0$ の拘束から $\omega_s$、つまりすべり周波数相当の量が決まる。

方式Aは、状態空間と極配置の考え方を理解するのに適している。方式Bは、論文5.3節に沿った軽量な二次磁束オブザーバとして実装しやすい。方式Cは、組込み実装の計算負荷を下げる実用的な方式として有力である。ただし、方式Cではオブザーバ構成とすべり周波数計算式が一体であるため、$\omega_s$ の式を単なる外部すべり推定式として扱ってはいけない。

文献調査の観点では、量産組込みへ進める場合、SLED 2023の閉形式フルオーダオブザーバと、Qu/Hinkkanen/Harnefors型のゲインスケジューリングを中心に検討するのが自然である。堀5.3型は、国内文献に基づく軽量な比較方式として残すのがよい。

## 参考文献

1. 堀 洋一, Vincent Cotter, 茅 陽一, 「誘導電動機の磁束オブザーバに関する制御理論的考察」, 電気学会論文誌B, 106巻, 11号, pp.1001-1008, 1986. [J-STAGE PDF](https://www.jstage.jst.go.jp/article/ieejpes1972/106/11/106_11_1001/_pdf)
2. G. C. Verghese and S. R. Sanders, "Observers for Flux Estimation in Induction Machines", IEEE Transactions on Industrial Electronics, vol. 35, no. 1, pp. 85-94, 1988.
3. G. R. Slemon, "Modelling of Induction Machines for Electric Drives", IEEE Transactions on Industry Applications, vol. 25, no. 6, pp. 1126-1131, 1989.
4. H. Kubota, K. Matsuse, and T. Nakano, "DSP-Based Speed Adaptive Flux Observer of Induction Motor", Conference Record of the 1991 IEEE Industry Applications Society Annual Meeting, 1991. [DOI](https://doi.org/10.1109/IAS.1991.178183)
5. H. Kubota, K. Matsuse, and T. Nakano, "DSP-Based Speed Adaptive Flux Observer of Induction Motor", IEEE Transactions on Industry Applications, vol. 29, no. 2, pp. 344-348, 1993. [DOI](https://doi.org/10.1109/28.216542)
6. G. Yang and T.-H. Chin, "Adaptive-Speed Identification Scheme for a Vector-Controlled Speed Sensorless Inverter-Induction Motor Drive", IEEE Transactions on Industry Applications, vol. 29, no. 4, pp. 820-825, 1993. [DOI](https://doi.org/10.1109/28.232001)
7. J. Holtz, "The Representation of AC Machine Dynamics by Complex Signal Flow Graphs", IEEE Transactions on Industrial Electronics, vol. 42, no. 3, pp. 263-271, 1995.
8. M. Hinkkanen and J. Luomi, "Analysis and Design of Full-Order Flux Observers for Sensorless Induction Motors", IEEE Transactions on Industrial Electronics, 2004. [DOI](https://doi.org/10.1109/TIE.2004.834964)
9. M. Hinkkanen and J. Luomi, "Stabilization of Regenerating-Mode Operation in Sensorless Induction Motor Drives by Full-Order Flux Observer Design", IEEE Transactions on Industrial Electronics, vol. 51, no. 6, pp. 1318-1328, 2004. [DOI](https://doi.org/10.1109/TIE.2004.837902)
10. M. Hinkkanen and J. Luomi, "Comparative Study of Adaptive and Inherently Sensorless Observers for Variable-Speed Induction-Motor Drives", IEEE Transactions on Industrial Electronics, 2006. [DOI](https://doi.org/10.1109/TIE.2005.862314)
11. L. Harnefors, "Globally Stable Speed-Adaptive Observers for Sensorless Induction Motor Drives", IEEE Transactions on Industrial Electronics, vol. 54, no. 2, pp. 1243-1245, 2007. [DOI](https://doi.org/10.1109/TIE.2007.892729)
12. S. Sangwongwanich, S. Suwankawin, S. Po-ngam, and S. Koonlaboon, "A Unified Speed Estimation Design Framework for Sensorless AC Motor Drives Based on Positive-Real Property", PCC-Nagoya, pp. 1111-1118, 2007. [DOI](https://doi.org/10.1109/PCCON.2007.373105)
13. L. Harnefors and M. Hinkkanen, "Complete Stability of Reduced-Order and Full-Order Observers for Sensorless IM Drives", IEEE Transactions on Industrial Electronics, vol. 55, no. 3, pp. 1319-1329, 2008. [DOI](https://doi.org/10.1109/TIE.2007.909077)
14. M. Hinkkanen, L. Harnefors, and J. Luomi, "Reduced-Order Flux Observers With Stator-Resistance Adaptation for Speed-Sensorless Induction Motor Drives", IEEE Transactions on Power Electronics, vol. 25, no. 5, pp. 1173-1183, 2010. [DOI](https://doi.org/10.1109/TPEL.2009.2039650)
15. Z. Qu, M. Hinkkanen, and L. Harnefors, "Gain Scheduling of a Full-Order Observer for Sensorless Induction Motor Drives", IEEE Transactions on Industry Applications, vol. 50, no. 6, pp. 3834-3845, 2014. [DOI](https://doi.org/10.1109/TIA.2014.2323482)
16. Bin Chen, Ting Wang, Wenxi Yao, Kevin Lee, Zhengyu Lu, "Speed Convergence Rate-Based Feedback Gains Design of Adaptive Full-Order Observer in Sensorless Induction Motor Drives", IET Electric Power Applications, 2014. [DOI](https://doi.org/10.1049/iet-epa.2013.0210)
17. Hongbo Wang, Wei Sun, Yong Yu, Gaolin Wang, Dianguo Xu, "Robustness Improvement for Adaptive Full Order Observer in Sensorless Induction Motor Drives", ICPE-ECCE Asia, 2015. [DOI](https://doi.org/10.1109/ICPE.2015.7167961)
18. Apirach Rattanaudompisut, Sakorn Po-Ngam, "The Practical and Simple Feedback Gains Design of an Adaptive Full-Order Observer for Speed-Sensorless Induction Motor Drives", ECTI-CON, 2015. [DOI](https://doi.org/10.1109/ECTICON.2015.7207045)
19. Zhonggang Yin, Yanqing Zhang, Xiangqian Tong, Jing Liu, Yanru Zhong, "Stability and Dynamic Performance Improvement of Adaptive Full-Order Observers in Sensorless Induction Motor Drives", IPEMC-ECCE Asia, 2016. [DOI](https://doi.org/10.1109/IPEMC.2016.7512313)
20. Mohamad Koteich, Gilles Duc, Abdelmalek Maloum, Guillaume Sandou, "Observability of Sensorless Electric Drives", arXiv:1602.04468, 2016. [arXiv](https://arxiv.org/abs/1602.04468)
21. Julian Kullick and Christoph M. Hackl, "Speed-Sensorless State Feedback Control of Induction Machines With LC Filter", arXiv:1807.11799, 2018. [arXiv](https://arxiv.org/abs/1807.11799)
22. Cheng Luo, Bo Wang, Yong Yu, Zhixin Huo, Guoqiang Zhang, Dianguo Xu, "A Speed Adaptive Scheme-Based Full-Order Observer for Sensorless Induction Motor Drives in Low-Speed Regenerating Operation Range", ICEMS, 2018. [DOI](https://doi.org/10.23919/ICEMS.2018.8549249)
23. Anton Pyrkin, Alexey Bobtsov, Alexey Vedyakov, Romeo Ortega, Anastasiia Vediakova, Madina Sinetova, "A Flux and Speed Observer for Induction Motors with Unknown Rotor Resistance and Load Torque and no Persistent Excitation Requirement", arXiv:2009.00966, 2020. [arXiv](https://arxiv.org/abs/2009.00966)
24. "Stable Adaptive Estimation for Speed-Sensorless Induction Motor Drives: A Geometric Approach", ICEM, 2020. [DOI](https://doi.org/10.1109/ICEM49940.2020.9270926)
25. Lauri Tiitinen, Marko Hinkkanen, Lennart Harnefors, "Stable and Passive Observer-Based V/Hz Control for Induction Motors", ECCE, 2022. [DOI](https://doi.org/10.1109/ECCE50734.2022.9948057)
26. Ruhan Li, Yifei Zheng, Cheng Luo, Zhijie Xu, Kai Yang, "Multi-Objective Optimization Based Feedback Gains Design of Adaptive Full-Order Observer for Induction Motor Sensorless Drive", ICEMS, 2022. [DOI](https://doi.org/10.1109/ICEMS56177.2022.9983330)
27. Lauri Tiitinen, Marko Hinkkanen, Lennart Harnefors, "Speed-Adaptive Full-Order Observer Revisited: Closed-Form Design for Induction Motor Drives", IEEE SLED, 2023. [DOI](https://doi.org/10.1109/SLED57582.2023.10261359)
28. Hongwu Chen, Jian Li, Yang Lu, Kai Yang, Linghao Wu, Zhi Liu, "Stability Improvement of Adaptive Full-Order Observer for Sensorless Induction Motor Drives in Low-Speed Regenerating Mode", IEEE Transactions on Transportation Electrification, 2024. [DOI](https://doi.org/10.1109/TTE.2023.3257056)
29. Ying Wang, Xiancheng Huang, Bo Wang, Dianguo Xu, "Feedback Matrix Design of Adaptive Full-Order Observer for Stator Resistance Robustness Improvement in Speed Sensorless Induction Motor Drives", IEEE Transactions on Industrial Electronics, 2026. [DOI](https://doi.org/10.1109/TIE.2025.3647924)
30. 本リポジトリの実装: [c/flux_observer.c](c/flux_observer.c), [c/sled23_flux_observer.c](c/sled23_flux_observer.c).
