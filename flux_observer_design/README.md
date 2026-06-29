# 回転dq座標フルオーダ磁束オブザーバ設計・評価

本資料は、指定された誘導機定数を用いて、回転dq座標上のフルオーダ磁束オブザーバを設計し、推定一次磁束・二次磁束、および出力一次電流・二次電流が真値へ収束することを確認した結果をまとめたものである。オブザーバゲインは極配置法でオンライン計算する。

使用した誘導機定数は以下である。

| 記号 | 値 | 説明 |
|---|---:|---|
| \(R_s\) | \(0.00762\ \Omega\) | 一次抵抗 |
| \(R_r\) | \(0.008041\ \Omega\) | 二次抵抗 |
| \(L_{ls}\) | \(0.0000419\ \mathrm{H}\) | 一次漏れインダクタンス |
| \(L_{lr}\) | \(0.0000419\ \mathrm{H}\) | 二次漏れインダクタンス |
| \(M=L_m\) | \(0.0001608\ \mathrm{H}\) | 相互インダクタンス |
| \(p\) | \(4\) | 極対数 |

合成インダクタンスは

$$
L_s = L_{ls}+M,\qquad L_r=L_{lr}+M,\qquad D=L_sL_r-M^2
$$

である。

## 1. 概要

目的は、一次電圧、一次電流、回転速度、およびdq座標の回転角速度から、一次磁束と二次磁束を推定することである。オブザーバの状態変数は

$$
x=
\begin{bmatrix}
\psi_s\\
\psi_r
\end{bmatrix}
$$

とする。ここで \(\psi_s=\psi_{sd}+j\psi_{sq}\)、\(\psi_r=\psi_{rd}+j\psi_{rq}\) は回転dq座標上の複素表現である。これは実数状態

$$
\begin{bmatrix}
\psi_{sd} & \psi_{sq} & \psi_{rd} & \psi_{rq}
\end{bmatrix}^T
$$

と等価である。

オブザーバの測定出力は一次電流であり、

$$
i_s=i_{sd}+ji_{sq}
$$

を補正入力に使う。二次電流 \(i_r=i_{rd}+ji_{rq}\) は直接測定しないが、推定磁束からオブザーバ出力として計算する。

評価条件は以下とした。

| 条件 | 速度 | トルク | \(i_{sd}\) | \(i_{sq}\) | \(\omega_{\mathrm{slip}}\) |
|---|---:|---:|---:|---:|---:|
| 力行高速 | \(5000\ \mathrm{r/min}\) | \(+220\ \mathrm{Nm}\) | \(673.610\ \mathrm{A}\) | \(+426.722\ \mathrm{A}\) | \(+25.130\ \mathrm{rad/s}\) |
| 回生高速 | \(5000\ \mathrm{r/min}\) | \(-220\ \mathrm{Nm}\) | \(673.610\ \mathrm{A}\) | \(-426.722\ \mathrm{A}\) | \(-25.130\ \mathrm{rad/s}\) |
| 回生低速 | \(1000\ \mathrm{r/min}\) | \(-220\ \mathrm{Nm}\) | \(673.610\ \mathrm{A}\) | \(-426.722\ \mathrm{A}\) | \(-25.130\ \mathrm{rad/s}\) |

dq座標の回転角速度は指定どおり

$$
\omega_k=p\omega_m+\omega_{\mathrm{slip}}
$$

とした。

実装と評価は [scripts/run_flux_observer_evaluation.py](scripts/run_flux_observer_evaluation.py) にまとめている。PDF版は [flux_observer_design_report.pdf](flux_observer_design_report.pdf) に出力した。

## 2. 磁束オブザーバの構成

### 2.1 電流と磁束の関係

d軸・q軸それぞれで、磁束と電流は

$$
\begin{bmatrix}
\psi_s\\
\psi_r
\end{bmatrix}
=
\begin{bmatrix}
L_s & M\\
M & L_r
\end{bmatrix}
\begin{bmatrix}
i_s\\
i_r
\end{bmatrix}
$$

を満たす。したがって、磁束から電流を計算する式は

$$
i_s = \frac{L_r\psi_s-M\psi_r}{D}
$$

$$
i_r = \frac{-M\psi_s+L_s\psi_r}{D}
$$

である。オブザーバでは、この式を用いて推定一次電流と推定二次電流を出力する。

### 2.2 回転dq座標の状態方程式

一次電圧を \(v_s=v_{sd}+jv_{sq}\)、二次電圧を短絡として \(v_r=0\) とする。回転dq座標の状態方程式は

$$
\dot{\psi}_s
=
v_s
-R_s i_s
-j\omega_k\psi_s
$$

$$
\dot{\psi}_r
=
-R_r i_r
-j(\omega_k-\omega_r)\psi_r
$$

である。ここで

$$
\omega_r=p\omega_m
$$

は電気角のロータ速度である。

電流式を代入すると、

$$
\dot{x}=Ax+Bv_s
$$

となる。ただし、

$$
A=
\begin{bmatrix}
-R_sL_r/D-j\omega_k & R_sM/D\\
R_rM/D & -R_rL_s/D-j(\omega_k-\omega_r)
\end{bmatrix}
$$

$$
B=
\begin{bmatrix}
1\\
0
\end{bmatrix}
$$

である。一次電流出力は

$$
i_s=C_sx,\qquad
C_s=
\begin{bmatrix}
L_r/D & -M/D
\end{bmatrix}
$$

である。

### 2.3 オブザーバ式

オブザーバは一次電流誤差で補正する。

$$
\dot{\hat{x}}
=
A_o\hat{x}
B_o v_{s,m}
l\left(i_{s,m}-C_{s,o}\hat{x}\right)
$$

ここで添字 \(o\) はオブザーバ内部で使用する定数を表し、添字 \(m\) は測定値を表す。定数誤差評価では \(A_o,C_{s,o}\) に誤差付き定数を使用し、真値モデルには基準定数を使用した。

複素ゲインを

$$
l=
\begin{bmatrix}
l_s\\
l_r
\end{bmatrix}
$$

とすると、実数4状態のゲイン行列は

$$
L_{\mathrm{real}}=
\begin{bmatrix}
\mathrm{Re}(l_s) & -\mathrm{Im}(l_s)\\
\mathrm{Im}(l_s) & \mathrm{Re}(l_s)\\
\mathrm{Re}(l_r) & -\mathrm{Im}(l_r)\\
\mathrm{Im}(l_r) & \mathrm{Re}(l_r)
\end{bmatrix}
$$

である。実装では複素2状態で計算しているが、これは上記の4実状態オブザーバと等価である。

## 3. オブザーバゲイン設計法

オブザーバゲインは、各制御周期で現在の \(\omega_r,\omega_k\) を用いて \(A_o\) を更新し、極配置法によりオンライン計算する。

目標極を

$$
p_1=-\omega_o,\qquad p_2=-1.55\omega_o
$$

とした。本評価では

$$
\omega_o=2200\ \mathrm{rad/s}
$$

を用いたため、目標極は

$$
p_1=-2200,\qquad p_2=-3410
$$

である。

閉じたオブザーバ誤差行列は

$$
A_L=A_o-lC_{s,o}
$$

である。2状態複素モデルでは、所望特性多項式を

$$
\chi_d(s)=(s-p_1)(s-p_2)
$$

とし、

$$
\mathrm{tr}(A_L)=p_1+p_2
$$

$$
\det(A_L)=p_1p_2
$$

を満たすように \(l\) を決定する。

ランク1更新の性質より、

$$
\mathrm{tr}(A_o-lC_{s,o})
=
\mathrm{tr}(A_o)-C_{s,o}l
$$

$$
\det(A_o-lC_{s,o})
=
\det(A_o)-C_{s,o}\operatorname{adj}(A_o)l
$$

である。したがって \(l\) は次の2元連立一次方程式で得られる。

$$
\begin{bmatrix}
C_{s,o}\\
C_{s,o}\operatorname{adj}(A_o)
\end{bmatrix}
l
=
\begin{bmatrix}
\mathrm{tr}(A_o)-(p_1+p_2)\\
\det(A_o)-p_1p_2
\end{bmatrix}
$$

この式を各サンプルで解けば、速度に依存するオンライン極配置オブザーバになる。Python実装では `observer_gain_by_pole_placement()` がこの計算を行う。

![observer pole map](figures/observer_pole_map.png)

## 4. オブザーバの安定性の証明

定数誤差、電圧誤差、電流誤差がない場合を考える。このとき真値モデルとオブザーバ内部モデルは一致し、

$$
\dot{x}=Ax+Bv_s
$$

$$
\dot{\hat{x}}=A\hat{x}+Bv_s+l(i_s-C_s\hat{x})
$$

である。推定誤差を

$$
\tilde{x}=x-\hat{x}
$$

と定義すると、

$$
\dot{\tilde{x}}
=
(A-lC_s)\tilde{x}
$$

となる。

3章の設計により、複素2状態の誤差行列 \(A-lC_s\) の固有値は \(p_1,p_2\) に配置される。今回は

$$
\mathrm{Re}(p_1)<0,\qquad \mathrm{Re}(p_2)<0
$$

であるため、任意の初期誤差に対して

$$
\|\tilde{x}(t)\|
\le
K e^{-\alpha t}\|\tilde{x}(0)\|
$$

を満たす \(K>0,\alpha>0\) が存在する。したがって、固定速度条件ではオブザーバ誤差は指数安定であり、一次磁束・二次磁束推定値は真値へ収束する。

実数4状態表現では、固有値は複素モデルの極とその共役で構成される。今回の目標極は実数負極であるため、4実状態モデルでも全固有値は左半平面にある。

定数誤差や測定誤差がある場合、誤差方程式は

$$
\dot{\tilde{x}}
=
(A_o-lC_{s,o})\tilde{x}
+d(t)
$$

となる。ここで \(d(t)\) は定数ずれ、電圧誤差、電流誤差、ノイズによる有界外乱である。\(A_o-lC_{s,o}\) がHurwitzであれば、推定誤差は有界入力有界状態となる。したがって、誤差ありでは真値への完全収束ではなく、誤差要因に応じた定常推定誤差へ収束する。

なお、速度が急変する一般の線形時変系に対する大域安定性は、共通Lyapunov関数またはゲインスケジューリング速度の制限を別途確認する必要がある。本評価は指定条件に合わせ、速度一定の凍結動作点で実施した。

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

シミュレーション条件は以下である。

| 項目 | 値 |
|---|---:|
| シミュレーション時間 | \(0.12\ \mathrm{s}\) |
| 積分刻み | \(10\ \mu\mathrm{s}\) |
| 初期推定誤差 | 一次・二次磁束に大きな初期誤差を付与 |
| 誤差評価窓 | 最終20%区間のRMS |

定数誤差は、\(R_s,R_r,M,L_{ls},L_{lr}\) を1つずつ

$$
\pm 5\%,\quad \pm 10\%,\quad \pm 20\%
$$

変化させた。電圧誤差と電流誤差は以下とした。

| 誤差種別 | 条件 |
|---|---|
| 電圧ゲイン誤差 | \(\pm 5\%\) |
| 電圧オフセット | \(+1\ \mathrm{V}\) on d軸, \(-1\ \mathrm{V}\) on q軸 |
| 電流ゲイン誤差 | \(\pm 1\%\) |
| 電流オフセット | \(+1\ \mathrm{A}\) on d軸, \(-1\ \mathrm{A}\) on q軸 |
| 電流ノイズ | \(1\ \mathrm{A_{rms}}\) on d/q軸 |

### 5.2 無誤差時の収束

無誤差では、3動作点すべてで一次磁束・二次磁束、および一次電流・二次電流の推定値が真値へ収束した。

| 条件 | 二次磁束RMS誤差 | 一次電流RMS誤差 | 二次磁束1%収束時間 | 一次電流1A収束時間 |
|---|---:|---:|---:|---:|
| 5000 r/min, +220 Nm | \(2.06\times10^{-13}\%\) | \(2.01\times10^{-13}\ \mathrm{A}\) | \(2.76\ \mathrm{ms}\) | \(3.65\ \mathrm{ms}\) |
| 5000 r/min, -220 Nm | \(5.95\times10^{-13}\%\) | \(8.04\times10^{-13}\ \mathrm{A}\) | \(2.81\ \mathrm{ms}\) | \(3.70\ \mathrm{ms}\) |
| 1000 r/min, -220 Nm | \(2.53\times10^{-12}\%\) | \(4.58\times10^{-13}\ \mathrm{A}\) | \(3.53\ \mathrm{ms}\) | \(3.69\ \mathrm{ms}\) |

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
| \(R_s\) | \(+20\%\) | \(0.655\%\) | \(0.101\ \mathrm{A}\) |
| \(R_r\) | \(-20\%\) | \(0.441\%\) | \(1.585\ \mathrm{A}\) |
| \(M\) | \(-20\%\) | \(5.576\%\) | \(3.899\ \mathrm{A}\) |
| \(L_{ls}\) | \(-20\%\) | \(7.456\%\) | \(1.296\ \mathrm{A}\) |
| \(L_{lr}\) | \(-20\%\) | \(2.589\%\) | \(0.234\ \mathrm{A}\) |

今回の範囲では、二次磁束推定は \(L_{ls}\) と \(M\) の誤差に比較的敏感であり、\(R_s\) と \(R_r\) の単独誤差に対しては比較的鈍感であった。ただし \(R_r\) 誤差は二次電流および滑り推定を使う制御系では別途影響が大きくなる可能性がある。

![parameter error sweep](figures/parameter_error_sweep.png)

### 5.4 電圧誤差・電流誤差の評価

5000 r/min, -220 Nm条件の結果を以下に示す。全ケースで発散は発生しなかった。

| ケース | 二次磁束RMS誤差 | 一次電流RMS誤差 | 安定性 |
|---|---:|---:|---|
| 無誤差 | \(0.000\%\) | \(0.000\ \mathrm{A}\) | 安定 |
| 電圧ゲイン \(+5\%\) | \(7.773\%\) | \(1.201\ \mathrm{A}\) | 安定 |
| 電圧ゲイン \(-5\%\) | \(7.773\%\) | \(1.201\ \mathrm{A}\) | 安定 |
| 電圧オフセット | \(0.763\%\) | \(0.118\ \mathrm{A}\) | 安定 |
| 電流ゲイン \(+1\%\) | \(0.636\%\) | \(7.918\ \mathrm{A}\) | 安定 |
| 電流ゲイン \(-1\%\) | \(0.636\%\) | \(7.918\ \mathrm{A}\) | 安定 |
| 電流オフセット | \(0.113\%\) | \(1.404\ \mathrm{A}\) | 安定 |
| 電流ノイズ \(1\ \mathrm{A_{rms}}\) | \(0.016\%\) | \(0.265\ \mathrm{A}\) | 安定 |

電圧ゲイン誤差は磁束スケールに直接影響するため、今回のセンサ誤差条件では最も二次磁束誤差が大きかった。一方、電流ゲイン誤差は、一次電流出力の真値比較ではゲイン差がそのまま残るため一次電流RMS誤差が大きく見えるが、磁束推定誤差は1%未満であった。

![sensor error summary](figures/sensor_error_summary.png)

### 5.5 結論

本設計では、回転dq座標の一次・二次磁束を状態とし、一次電流誤差で補正するフルオーダ磁束オブザーバを構成した。極配置法により、各動作点の凍結モデルでオブザーバ誤差極を \(-2200\) および \(-3410\ \mathrm{rad/s}\) に配置した。

無誤差では、一次・二次磁束および一次・二次電流の推定値は真値へ収束した。誤差あり評価では、指定範囲の定数誤差、電圧誤差、電流誤差に対して発散は見られず、安定性は維持された。精度面では、電圧ゲイン誤差、一次漏れインダクタンス誤差、相互インダクタンス誤差が二次磁束推定精度に対して支配的であった。


