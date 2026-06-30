# 回転dq座標フルオーダ磁束オブザーバ設計・評価

本資料は、指定された誘導機定数を用いて、回転dq座標上のフルオーダ磁束オブザーバを設計し、一次磁束・二次磁束、および一次電流・二次電流の推定値が真値へ収束することを確認した結果をまとめたものである。

今回の方針は、実数4状態モデルをそのまま扱う正攻法の極配置である。d軸/q軸の回転対称性を使った2極相当の省略設計は使わず、オブザーバゲインは一般の実数行列

$$
H \in \mathbb{R}^{4\times 2}
$$

として設計する。したがって、4状態に対応する4個のオブザーバ極をすべて明示的に配置する。

使用した誘導機定数は以下である。

| 記号 | 値 | 説明 |
|---|---:|---|
| $R_s$ | $0.00762\ \Omega$ | 一次抵抗 |
| $R_r$ | $0.008041\ \Omega$ | 二次抵抗 |
| $L_{ls}$ | $0.0000419\ \mathrm{H}$ | 一次漏れインダクタンス |
| $L_{lr}$ | $0.0000419\ \mathrm{H}$ | 二次漏れインダクタンス |
| $M=L_m$ | $0.0001608\ \mathrm{H}$ | 相互インダクタンス |
| $p$ | $4$ | 極対数 |

合成インダクタンスは以下で定義する。

$$
L_s=L_{ls}+M,\qquad L_r=L_{lr}+M,\qquad D=L_sL_r-M^2
$$

実装と評価は [scripts/run_flux_observer_evaluation.py](scripts/run_flux_observer_evaluation.py) にまとめている。C言語実装は [c/flux_observer.c](c/flux_observer.c) と [c/flux_observer.h](c/flux_observer.h) に置いた。PDF版は [flux_observer_design_report.pdf](flux_observer_design_report.pdf) に出力する。

## 1. 概要

目的は、一次電圧、一次電流、ロータ機械角速度、および滑り角周波数から、回転dq座標上の一次磁束と二次磁束を推定することである。

状態変数は一次磁束と二次磁束であり、実数4状態として以下のように定義する。

$$
x=
\begin{bmatrix}
\psi_{sd} & \psi_{sq} & \psi_{rd} & \psi_{rq}
\end{bmatrix}^{T}
$$

補正に使う測定出力は一次電流である。

$$
y=
\begin{bmatrix}
i_{sd} & i_{sq}
\end{bmatrix}^{T}
$$

二次電流は直接測定しない。一次磁束・二次磁束の推定値から、オブザーバの出力として計算する。

評価条件は以下とした。

| 条件 | 速度 | トルク | $i_{sd}$ | $i_{sq}$ | $\omega_{\mathrm{slip}}$ |
|---|---:|---:|---:|---:|---:|
| 力行高速 | $5000\ \mathrm{r/min}$ | $+220\ \mathrm{Nm}$ | $673.610\ \mathrm{A}$ | $+426.722\ \mathrm{A}$ | $+25.130\ \mathrm{rad/s}$ |
| 回生高速 | $5000\ \mathrm{r/min}$ | $-220\ \mathrm{Nm}$ | $673.610\ \mathrm{A}$ | $-426.722\ \mathrm{A}$ | $-25.130\ \mathrm{rad/s}$ |
| 回生低速 | $1000\ \mathrm{r/min}$ | $-220\ \mathrm{Nm}$ | $673.610\ \mathrm{A}$ | $-426.722\ \mathrm{A}$ | $-25.130\ \mathrm{rad/s}$ |

回転dq座標の角速度は指定どおり以下とした。

$$
\omega_k=p\omega_m+\omega_{\mathrm{slip}}
$$

## 2. 磁束オブザーバの構成

### 2.1 電流と磁束の関係

一次磁束、二次磁束、一次電流、二次電流をそれぞれ $\psi_s,\psi_r,i_s,i_r$ とする。d軸・q軸それぞれで、磁束と電流は以下の関係を満たす。

$$
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
$$

したがって、磁束から電流を計算する式は以下である。

$$
i_s=\frac{L_r\psi_s-M\psi_r}{D}
$$

$$
i_r=\frac{-M\psi_s+L_s\psi_r}{D}
$$

この式により、オブザーバは推定一次電流と推定二次電流を出力する。

### 2.2 回転dq座標の状態方程式

90度回転行列を以下で定義する。

$$
J=
\begin{bmatrix}
0 & -1\\
1 & 0
\end{bmatrix}
$$

一次電圧を $v_s=[v_{sd}\ v_{sq}]^T$、二次電圧を短絡として $v_r=0$ とする。回転dq座標の磁束方程式は以下である。

$$
\dot{\psi}_s=v_s-R_s i_s-\omega_k J\psi_s
$$

$$
\dot{\psi}_r=-R_r i_r-(\omega_k-\omega_r)J\psi_r
$$

ここで、ロータ電気角速度は以下である。

$$
\omega_r=p\omega_m
$$

電流式を代入すると、実数4状態の状態方程式は以下になる。

$$
\dot{x}=Ax+Bv_s
$$

$$
A=
\begin{bmatrix}
-\frac{R_sL_r}{D}I_2-\omega_kJ & \frac{R_sM}{D}I_2\\
\frac{R_rM}{D}I_2 & -\frac{R_rL_s}{D}I_2-(\omega_k-\omega_r)J
\end{bmatrix}
$$

$$
B=
\begin{bmatrix}
I_2\\
0_{2\times2}
\end{bmatrix}
$$

一次電流出力は以下である。

$$
y=Cx
$$

$$
C=
\begin{bmatrix}
\frac{L_r}{D}I_2 & -\frac{M}{D}I_2
\end{bmatrix}
$$

二次電流出力は以下で計算する。

$$
i_r=C_r x
$$

$$
C_r=
\begin{bmatrix}
-\frac{M}{D}I_2 & \frac{L_s}{D}I_2
\end{bmatrix}
$$

### 2.3 オブザーバ式

本資料で用いるLuenberger型オブザーバは以下である。

$$
\dot{\hat{x}}=A_o\hat{x}+B_o v_{s,m}+H\left(y_m-C_o\hat{x}\right)
$$

ここで、$A_o,B_o,C_o$ はオブザーバ内部で使用するモータ定数から作った行列であり、$v_{s,m},y_m$ は測定電圧と測定一次電流である。定数誤差評価では、真値モデルには基準定数を使い、オブザーバ内部行列には誤差付き定数を使った。

ゲイン $H$ は一般の実数 $4\times2$ 行列とする。

$$
H=
\begin{bmatrix}
H_{00} & H_{01}\\
H_{10} & H_{11}\\
H_{20} & H_{21}\\
H_{30} & H_{31}
\end{bmatrix}
$$

この設計では、d軸/q軸が独立、または完全対称であるという仮定でゲインを減らさない。速度干渉項を含んだ実4状態モデルに対し、4個の極をすべて配置する。

## 3. オブザーバゲイン設計法

### 3.1 極配置法の基本

連続時間の線形状態方程式を以下とする。

$$
\dot{x}=Ax+Bu
$$

$$
y=Cx
$$

オブザーバを以下で構成する。

$$
\dot{\hat{x}}=A\hat{x}+Bu+H(y-C\hat{x})
$$

推定誤差を以下で定義する。

$$
\tilde{x}=x-\hat{x}
$$

この時間微分を取ると、

$$
\dot{\tilde{x}}=\dot{x}-\dot{\hat{x}}
$$

である。真値モデルとオブザーバ式を代入すると、同じ入力 $Bu$ は差し引きで消え、さらに $y=Cx$ より以下を得る。

$$
\dot{\tilde{x}}=(A-HC)\tilde{x}
$$

つまり、オブザーバの収束速度と振動性は $A-HC$ の固有値で決まる。この固有値を設計者が指定した位置へ動かすことを、オブザーバの極配置と呼ぶ。

状態数が4なので、今回のオブザーバ極も4個である。4状態すべての誤差を設計対象にするため、目標極は

$$
p_1,p_2,p_3,p_4
$$

の4個を指定する。

### 3.2 今回の目標極

今回の基準設計では、観測帯域を

$$
\omega_o=2200\ \mathrm{rad/s}
$$

とし、目標極を以下に置いた。

$$
p_1=-\omega_o
$$

$$
p_2=-1.25\omega_o
$$

$$
p_3=-1.55\omega_o
$$

$$
p_4=-2.0\omega_o
$$

数値では以下である。

| 極 | 値 |
|---|---:|
| $p_1$ | $-2200\ \mathrm{rad/s}$ |
| $p_2$ | $-2750\ \mathrm{rad/s}$ |
| $p_3$ | $-3410\ \mathrm{rad/s}$ |
| $p_4$ | $-4400\ \mathrm{rad/s}$ |

$\omega_o=2200\ \mathrm{rad/s}$ は、電流制御帯域 $\omega_{cc}=1000\ \mathrm{rad/s}$ より速く、かつ測定ノイズや離散化誤差を過度に増幅しない範囲として選んだ。1.25、1.55、2.0の倍率は、重根を避け、4個の極を適度に分離して数値条件を悪化させにくくするための設計例である。

制御周期を $T_s=100\ \mu\mathrm{s}$ とすると、最速極 $4400\ \mathrm{rad/s}$ はサンプリング角周波数 $2\pi/T_s=62832\ \mathrm{rad/s}$ の約7.0%であり、離散化に対して余裕がある。

C APIでは、標準設定として `FluxObserver_SetPolePlacement(&observer, 2200.0f, 2.0f)` を用いる。個別に4個の極を指定したい場合は `FluxObserver_SetObserverPoles()` を使う。

### 3.3 Sylvester方程式による4極配置

今回の実装は、実4状態の行列 $A_o,C_o$ をそのまま使う。目標誤差ダイナミクスを以下で定義する。

$$
F=\mathrm{diag}(p_1,p_2,p_3,p_4)
$$

出力誤差を4状態へ配分する行列を以下とする。

$$
G=
\begin{bmatrix}
1 & 0\\
0 & 1\\
1 & 0\\
0 & 1
\end{bmatrix}
$$

この $G$ は設計上の計算行列であり、物理量を直接意味するものではない。各状態が少なくとも一方の電流誤差信号から補正を受けるように選んでいる。

次のSylvester方程式を解く。

$$
T A_o - F T = G C_o
$$

$T$ が正則であれば、オブザーバゲインを以下で計算する。

$$
H=T^{-1}G
$$

このとき、

$$
T(A_o-HC_o)=F T
$$

が成り立つ。したがって、

$$
A_o-HC_o=T^{-1}FT
$$

となり、$A_o-HC_o$ は $F$ と相似である。相似な行列は同じ固有値を持つため、オブザーバ誤差行列の4個の極は指定した $p_1,p_2,p_3,p_4$ に一致する。

実装では各制御周期で以下の順に計算する。

1. APIからモータ定数と制御周期を取得する。
2. 入力された $\omega_m$ と $\omega_{\mathrm{slip}}$ から $\omega_r=p\omega_m$ と $\omega_k=\omega_r+\omega_{\mathrm{slip}}$ を計算する。
3. 現在の速度条件で $A_o,B_o,C_o$ を作る。
4. 目標極から $F$ を作る。
5. $T A_o - F T = G C_o$ を実数16元連立一次方程式として解く。
6. $H=T^{-1}G$ を解く。
7. 得られた $H$ でオブザーバを1ステップ積分する。

Python実装では `observer_H_gain_by_pole_placement()` がこの計算を行う。C実装では `fo_observer_H()` が同じ計算を行う。複素数表現による省略やd/q対称性によるゲイン削減は使っていない。

![observer pole map](figures/observer_pole_map.png)

## 4. オブザーバの安定性の証明

まず、定数誤差、電圧誤差、電流誤差がない場合を考える。このとき真値モデルとオブザーバ内部モデルは一致する。

$$
\dot{x}=Ax+Bv_s
$$

$$
\dot{\hat{x}}=A\hat{x}+Bv_s+H(y-C\hat{x})
$$

推定誤差 $\tilde{x}=x-\hat{x}$ の時間微分は以下である。

$$
\dot{\tilde{x}}=\dot{x}-\dot{\hat{x}}
$$

上の2式を代入すると、同じ入力 $Bv_s$ は消える。

$$
\dot{\tilde{x}}=A(x-\hat{x})-H(y-C\hat{x})
$$

無誤差時は $y=Cx$ であるため、

$$
y-C\hat{x}=C(x-\hat{x})=C\tilde{x}
$$

である。したがって、

$$
\dot{\tilde{x}}=(A-HC)\tilde{x}
$$

となる。

3.3節の設計により、$A-HC$ は $F$ と相似であり、固有値は以下である。

$$
\lambda(A-HC)=\{p_1,p_2,p_3,p_4\}
$$

今回の設計では4個すべての極が負の実数である。

$$
\mathrm{Re}(p_i)<0\qquad (i=1,2,3,4)
$$

したがって、固定速度の凍結動作点では推定誤差は指数的に0へ収束する。すなわち、適当な正定数 $K,\alpha$ が存在し、以下を満たす。

$$
\|\tilde{x}(t)\|\le K e^{-\alpha t}\|\tilde{x}(0)\|
$$

定数誤差、電圧誤差、電流誤差がある場合は、誤差方程式に外乱項 $d(t)$ が加わる。

$$
\dot{\tilde{x}}=(A_o-HC_o)\tilde{x}+d(t)
$$

$A_o-HC_o$ がHurwitz、すなわち全固有値の実部が負であり、$d(t)$ が有界であれば、推定誤差も有界に保たれる。この場合、誤差は完全には0にならず、定数誤差や測定誤差に応じた定常誤差へ収束する。

なお、本評価は一定速度・一定滑りの凍結動作点で実施した。速度急変を含む一般の線形時変系としての大域安定性を主張するには、ゲインスケジューリング速度や共通Lyapunov関数の確認が別途必要である。

## 5. 誤差無し/有の評価結果

### 5.1 評価方法

評価スクリプトは以下で実行できる。

```powershell
python .\flux_observer_design\scripts\run_flux_observer_evaluation.py
```

生成物は以下である。

| ファイル | 内容 |
|---|---|
| [data/evaluation_summary.csv](data/evaluation_summary.csv) | 全評価ケースの数値結果 |
| [figures/nominal_waveform_5000rpm_motoring.png](figures/nominal_waveform_5000rpm_motoring.png) | 5000 r/min力行の無誤差波形 |
| [figures/nominal_waveform_5000rpm_regen.png](figures/nominal_waveform_5000rpm_regen.png) | 5000 r/min回生の無誤差波形 |
| [figures/nominal_waveform_1000rpm_regen.png](figures/nominal_waveform_1000rpm_regen.png) | 1000 r/min回生の無誤差波形 |
| [figures/nominal_convergence.png](figures/nominal_convergence.png) | 3動作点の収束誤差 |
| [figures/parameter_error_sweep.png](figures/parameter_error_sweep.png) | 定数誤差感度 |
| [figures/sensor_error_summary.png](figures/sensor_error_summary.png) | 電圧・電流誤差感度 |

シミュレーション条件は以下である。本評価はオブザーバ単体の評価であり、電流制御器、速度制御器、トルク制御器、滑り演算器を含む閉ループ制御系の評価ではない。

| 項目 | 値 |
|---|---:|
| シミュレーション時間 | $0.12\ \mathrm{s}$ |
| 積分刻み | $10\ \mu\mathrm{s}$ |
| 初期推定誤差 | 一次・二次磁束に大きな初期誤差を付与 |
| 誤差評価窓 | 最終20%区間のRMS |

定数誤差は、$R_s,R_r,M,L_{ls},L_{lr}$ を1つずつ以下の刻みで変化させた。

$$
\pm 5\%,\qquad \pm 10\%,\qquad \pm 20\%
$$

電圧誤差と電流誤差は以下とした。

| 誤差種別 | 条件 |
|---|---|
| 電圧ゲイン誤差 | $\pm 5\%$ |
| 電圧オフセット | $+1\ \mathrm{V}$ on d軸, $-1\ \mathrm{V}$ on q軸 |
| 電流ゲイン誤差 | $\pm 1\%$ |
| 電流オフセット | $+1\ \mathrm{A}$ on d軸, $-1\ \mathrm{A}$ on q軸 |
| 電流ノイズ | $1\ \mathrm{A_{rms}}$ on d/q軸 |

### 5.2 無誤差時の収束

無誤差では、3動作点すべてで一次磁束・二次磁束、および一次電流・二次電流の推定値が真値へ収束した。

| 条件 | 二次磁束RMS誤差 | 一次電流RMS誤差 | 二次磁束1%収束時間 | 一次電流1A収束時間 |
|---|---:|---:|---:|---:|
| 5000 r/min, +220 Nm | $6.52\times10^{-13}\%$ | $1.17\times10^{-12}\ \mathrm{A}$ | $2.75\ \mathrm{ms}$ | $3.64\ \mathrm{ms}$ |
| 5000 r/min, -220 Nm | $5.34\times10^{-13}\%$ | $8.51\times10^{-13}\ \mathrm{A}$ | $2.82\ \mathrm{ms}$ | $3.71\ \mathrm{ms}$ |
| 1000 r/min, -220 Nm | $1.74\times10^{-12}\%$ | $2.45\times10^{-13}\ \mathrm{A}$ | $3.50\ \mathrm{ms}$ | $3.65\ \mathrm{ms}$ |

以下に、各動作点での一次磁束、二次磁束、一次電流、二次電流の真値とオブザーバ推定値を示す。破線が推定値であり、初期推定誤差を付与した後、各成分が真値へ収束している。

![5000 r/min motoring nominal waveform](figures/nominal_waveform_5000rpm_motoring.png)

![5000 r/min regeneration nominal waveform](figures/nominal_waveform_5000rpm_regen.png)

![1000 r/min regeneration nominal waveform](figures/nominal_waveform_1000rpm_regen.png)

収束誤差を対数軸で表示した結果を以下に示す。

![nominal convergence](figures/nominal_convergence.png)

### 5.3 定数誤差の評価

全ての定数誤差ケースで発散は発生しなかった。5000 r/min, -220 Nm条件で、二次磁束推定誤差が最大となった各定数のケースは以下である。

| 誤差対象 | 最悪ケース | 二次磁束RMS誤差 | 一次電流RMS誤差 |
|---|---:|---:|---:|
| $R_s$ | $-20\%$ | $0.655\%$ | $0.094\ \mathrm{A}$ |
| $R_r$ | $-20\%$ | $0.405\%$ | $1.410\ \mathrm{A}$ |
| $M$ | $-20\%$ | $5.742\%$ | $2.651\ \mathrm{A}$ |
| $L_{ls}$ | $-20\%$ | $7.502\%$ | $0.941\ \mathrm{A}$ |
| $L_{lr}$ | $-20\%$ | $2.583\%$ | $0.181\ \mathrm{A}$ |

今回の範囲では、二次磁束推定は $L_{ls}$ と $M$ の誤差に比較的敏感であり、$R_s$ と $R_r$ の単独誤差に対しては比較的鈍感であった。ただし、これはオブザーバ単体を凍結動作点で評価した結果である。$R_r$ 誤差が滑り演算やベクトル制御の軸角に与える影響は、この評価には含めていない。

![parameter error sweep](figures/parameter_error_sweep.png)

### 5.4 電圧誤差・電流誤差の評価

5000 r/min, -220 Nm条件の結果を以下に示す。全ケースで発散は発生しなかった。

| ケース | 二次磁束RMS誤差 | 一次電流RMS誤差 | 安定性 |
|---|---:|---:|---|
| 無誤差 | $0.000\%$ | $0.000\ \mathrm{A}$ | 安定 |
| 電圧ゲイン $+5\%$ | $7.823\%$ | $0.764\ \mathrm{A}$ | 安定 |
| 電圧ゲイン $-5\%$ | $7.823\%$ | $0.764\ \mathrm{A}$ | 安定 |
| 電圧オフセット | $0.763\%$ | $0.101\ \mathrm{A}$ | 安定 |
| 電流ゲイン $+1\%$ | $0.647\%$ | $7.929\ \mathrm{A}$ | 安定 |
| 電流ゲイン $-1\%$ | $0.647\%$ | $7.929\ \mathrm{A}$ | 安定 |
| 電流オフセット | $0.115\%$ | $1.403\ \mathrm{A}$ | 安定 |
| 電流ノイズ $1\ \mathrm{A_{rms}}$ | $0.020\%$ | $0.285\ \mathrm{A}$ | 安定 |

電圧ゲイン誤差は磁束スケールに直接影響するため、今回のセンサ誤差条件では最も二次磁束誤差が大きかった。一方、電流ゲイン誤差は一次電流出力の真値比較ではゲイン差がそのまま残るため、一次電流RMS誤差が大きく見えるが、磁束推定誤差は1%未満であった。

![sensor error summary](figures/sensor_error_summary.png)

### 5.5 結論

本設計では、回転dq座標の一次・二次磁束を状態とし、一次電流誤差で補正するフルオーダ磁束オブザーバを構成した。ゲイン $H$ は一般の実数 $4\times2$ 行列として扱い、各動作点の凍結モデルで4個のオブザーバ誤差極を以下へ配置した。

$$
-2200,\quad -2750,\quad -3410,\quad -4400\ \mathrm{rad/s}
$$

無誤差では、一次・二次磁束および一次・二次電流の推定値は真値へ収束した。誤差あり評価では、指定範囲の定数誤差、電圧誤差、電流誤差に対して発散は見られず、安定性は維持された。精度面では、電圧ゲイン誤差、一次漏れインダクタンス誤差、相互インダクタンス誤差が二次磁束推定精度に対して支配的であった。

## 付録A. C言語実装

C言語実装は以下に追加した。

| ファイル | 内容 |
|---|---|
| [c/flux_observer.h](c/flux_observer.h) | 公開API、入出力構造体、オブザーバ状態構造体 |
| [c/flux_observer.c](c/flux_observer.c) | 回転dq座標フルオーダ磁束オブザーバ本体 |

実装は `float` を使用し、動的メモリ確保は行わない。モータ定数および制御周期は、プラットフォーム側APIから取得する前提である。オブザーバは各 `FluxObserver_Step()` 呼び出しでAPIを呼び、最新の $R_s,R_r,L_{ls},L_{lr},M,p,T_s$ を取得する。

プラットフォーム側は、次のコールバックを用意する。

```c
static int GetMotorConfig(void *user, FluxObserverMotorConfig *config)
{
    (void)user;
    config->rs_ohm = 0.00762f;
    config->rr_ohm = 0.008041f;
    config->lls_h = 0.0000419f;
    config->llr_h = 0.0000419f;
    config->lm_h = 0.0001608f;
    config->pole_pairs = 4u;
    config->control_period_s = 100.0e-6f;
    return 0;
}
```

初期化は以下で行う。

```c
FluxObserver observer;
FluxObserverApi api;

api.user = NULL;
api.get_motor_config = GetMotorConfig;
FluxObserver_Init(&observer, api);
FluxObserver_SetPolePlacement(&observer, 2200.0f, 2.0f);
FluxObserver_ResetFlux(&observer, 0.0f, 0.0f, 0.0f, 0.0f);
```

4個の極を直接指定する場合は以下を使う。

```c
FluxObserver_SetObserverPoles(
    &observer,
    -2200.0f,
    -2750.0f,
    -3410.0f,
    -4400.0f);
```

制御周期ごとの呼び出しは以下である。

```c
FluxObserverInput input;
FluxObserverOutput output;
FluxObserverStatus status;

input.vsd_v = vd;
input.vsq_v = vq;
input.isd_a = id_meas;
input.isq_a = iq_meas;
input.omega_m_rad_s = omega_m;
input.omega_slip_rad_s = omega_slip;

status = FluxObserver_Step(&observer, &input, &output);
if (status != FLUX_OBSERVER_OK) {
    /* handle API, parameter, or singularity error */
}
```

`FluxObserverOutput` には以下を出力する。

| 出力 | 内容 |
|---|---|
| `psi_sd_wb`, `psi_sq_wb` | 推定一次磁束 |
| `psi_rd_wb`, `psi_rq_wb` | 推定二次磁束 |
| `isd_hat_a`, `isq_hat_a` | 推定一次電流 |
| `ird_hat_a`, `irq_hat_a` | 推定二次電流 |
| `omega_r_rad_s`, `omega_k_rad_s` | API定数と入力から計算した電気角速度 |
| `H[4][2]` | 実数4状態のオブザーバゲイン行列 |

C実装のコンパイル確認は以下で行った。

```powershell
gcc -std=c99 -Wall -Wextra -pedantic -fsyntax-only .\flux_observer_design\c\flux_observer.c
```
