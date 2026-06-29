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

オブザーバは一次電流誤差で補正する。設計上の状態は実数4状態

$$
x_R=
\begin{bmatrix}
\psi_{sd} & \psi_{sq} & \psi_{rd} & \psi_{rq}
\end{bmatrix}^T
$$

とし、測定出力は

$$
y_R=
\begin{bmatrix}
i_{sd} & i_{sq}
\end{bmatrix}^T
$$

とする。2.2節の空間ベクトル状態方程式を実部・虚部に分解すると、実数状態方程式

$$
\dot{x}_R=A_Rx_R+B_Rv_R
$$

$$
y_R=C_Rx_R
$$

が得られる。ここで

$$
v_R=
\begin{bmatrix}
v_{sd} & v_{sq}
\end{bmatrix}^T
$$

である。

本資料で用いるオブザーバは

$$
\dot{\hat{x}}_R
=
A_{R,o}\hat{x}_R
+B_{R,o}v_{R,m}
+H\left(y_{R,m}-C_{R,o}\hat{x}_R\right)
$$

である。ここで、オブザーバゲインは

$$
H\in\mathbb{R}^{4\times2}
$$

である。添字 \(o\) はオブザーバ内部で使用する定数を表し、添字 \(m\) は測定値を表す。定数誤差評価では \(A_{R,o},B_{R,o},C_{R,o}\) に誤差付き定数を使用し、真値モデルには基準定数を使用した。

実装では、計算を簡潔にするため空間ベクトル表記の補助量

$$
h=
\begin{bmatrix}
h_s\\
h_r
\end{bmatrix}
$$

を用いる。これは正式な設計ゲインを複素数で置き換えるものではなく、次の実数ゲイン行列 \(H\) と一対一に対応する省略表現である。

$$
H=
\begin{bmatrix}
\mathrm{Re}(h_s) & -\mathrm{Im}(h_s)\\
\mathrm{Im}(h_s) & \mathrm{Re}(h_s)\\
\mathrm{Re}(h_r) & -\mathrm{Im}(h_r)\\
\mathrm{Im}(h_r) & \mathrm{Re}(h_r)
\end{bmatrix}
$$

誘導機のdq軸は、今回の線形モデルではd軸/q軸で同じ抵抗・インダクタンスを持つ。そのため、上式のような回転対称な \(H\) は自然な構造である。d軸/q軸で非対称な飽和やセンサ特性を入れる場合は、上記構造に制限せず、一般の実数 \(4\times2\) ゲイン \(H\) を直接設計する。

## 3. オブザーバゲイン設計法

オブザーバゲイン \(H\) は、各制御周期で現在の \(\omega_r,\omega_k\) を用いて \(A_{R,o}\) を更新し、極配置法によりオンライン計算する。

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
A_H=A_{R,o}-HC_{R,o}
$$

である。実数4状態表現では、今回の実数目標極に対して

$$
\lambda(A_H)=\{p_1,p_1,p_2,p_2\}
$$

となるように \(H\) を決定する。

Python実装では、同じ極配置を空間ベクトル表記で計算している。これは実数4状態設計と等価であり、計算後に上式で \(H\) へ戻せる。空間ベクトル表記では、閉じた誤差行列は

$$
A_h=A_o-hC_{s,o}
$$

であり、所望特性多項式を

$$
\chi_d(s)=(s-p_1)(s-p_2)
$$

とする。したがって、

$$
\mathrm{tr}(A_h)=p_1+p_2
$$

$$
\det(A_h)=p_1p_2
$$

を満たすように補助量 \(h\) を決定する。

ランク1更新の性質より、

$$
\mathrm{tr}(A_o-hC_{s,o})
=
\mathrm{tr}(A_o)-C_{s,o}h
$$

$$
\det(A_o-hC_{s,o})
=
\det(A_o)-C_{s,o}\operatorname{adj}(A_o)h
$$

である。したがって \(h\) は次の2元連立一次方程式で得られる。

$$
\begin{bmatrix}
C_{s,o}\\
C_{s,o}\operatorname{adj}(A_o)
\end{bmatrix}
h
=
\begin{bmatrix}
\mathrm{tr}(A_o)-(p_1+p_2)\\
\det(A_o)-p_1p_2
\end{bmatrix}
$$

この式を各サンプルで解き、得られた \(h\) から実数ゲイン \(H\) を構成すれば、速度に依存するオンライン極配置オブザーバになる。Python実装では `observer_h_gain_by_pole_placement()` がこの計算を行う。

![observer pole map](figures/observer_pole_map.png)

## 4. オブザーバの安定性の証明

定数誤差、電圧誤差、電流誤差がない場合を考える。このとき真値モデルとオブザーバ内部モデルは一致し、

$$
\dot{x}_R=A_Rx_R+B_Rv_R
$$

$$
\dot{\hat{x}}_R=A_R\hat{x}_R+B_Rv_R+H(y_R-C_R\hat{x}_R)
$$

である。推定誤差を

$$
\tilde{x}_R=x_R-\hat{x}_R
$$

と定義すると、

$$
\dot{\tilde{x}}_R
=
(A_R-HC_R)\tilde{x}_R
$$

となる。

3章の設計により、実数4状態の誤差行列 \(A_R-HC_R\) の固有値は \(p_1,p_1,p_2,p_2\) に配置される。今回は

$$
\mathrm{Re}(p_1)<0,\qquad \mathrm{Re}(p_2)<0
$$

であるため、任意の初期誤差に対して

$$
\|\tilde{x}_R(t)\|
\le
K e^{-\alpha t}\|\tilde{x}_R(0)\|
$$

を満たす \(K>0,\alpha>0\) が存在する。したがって、固定速度条件ではオブザーバ誤差は指数安定であり、一次磁束・二次磁束推定値は真値へ収束する。

空間ベクトル表記で計算した \(h\) は、2.3節の変換により実数ゲイン \(H\) と等価である。したがって、実装上の空間ベクトル計算を用いても、安定性は上記の実数4状態オブザーバ \(A_R-HC_R\) の固有値で判断できる。

定数誤差や測定誤差がある場合、誤差方程式は

$$
\dot{\tilde{x}}_R
=
(A_{R,o}-HC_{R,o})\tilde{x}_R
+d(t)
$$

となる。ここで \(d(t)\) は定数ずれ、電圧誤差、電流誤差、ノイズによる有界外乱である。\(A_{R,o}-HC_{R,o}\) がHurwitzであれば、推定誤差は有界入力有界状態となる。したがって、誤差ありでは真値への完全収束ではなく、誤差要因に応じた定常推定誤差へ収束する。

なお、速度が急変する一般の線形時変系に対する大域安定性は、共通Lyapunov関数またはゲインスケジューリング速度の制限を別途確認する必要がある。本評価は指定条件に合わせ、速度一定の凍結動作点で実施した。

## 5. C言語実装

C言語実装は以下に追加した。

| ファイル | 内容 |
|---|---|
| [c/flux_observer.h](c/flux_observer.h) | 公開API、入出力構造体、オブザーバ状態構造体 |
| [c/flux_observer.c](c/flux_observer.c) | 回転dq座標フルオーダ磁束オブザーバ本体 |

実装は `float` を使用し、動的メモリ確保は行わない。モータ定数および制御周期は、プラットフォーム側APIから取得する前提である。オブザーバは各 `FluxObserver_Step()` 呼び出しでAPIを呼び、最新の

$$
R_s,
R_r,
L_{ls},
L_{lr},
M,
p,
T_s
$$

を取得する。

### 5.1 APIの構成

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
FluxObserver_SetPolePlacement(&observer, 2200.0f, 1.55f);
FluxObserver_ResetFlux(&observer, 0.0f, 0.0f, 0.0f, 0.0f);
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

ここで、`omega_m_rad_s` は機械角速度、`omega_slip_rad_s` は滑り角周波数である。C実装内部ではAPIから取得した極対数 `pole_pairs` を使い、

$$
\omega_r=p\omega_m
$$

$$
\omega_k=\omega_r+\omega_{\mathrm{slip}}
$$

を計算する。

### 5.2 C実装の出力

`FluxObserverOutput` には以下を出力する。

| 出力 | 内容 |
|---|---|
| `psi_sd_wb`, `psi_sq_wb` | 推定一次磁束 |
| `psi_rd_wb`, `psi_rq_wb` | 推定二次磁束 |
| `isd_hat_a`, `isq_hat_a` | 推定一次電流 |
| `ird_hat_a`, `irq_hat_a` | 推定二次電流 |
| `omega_r_rad_s`, `omega_k_rad_s` | API定数と入力から計算した電気角速度 |
| `h[4][2]` | 実数4状態のオブザーバゲイン行列 \(H\) |

C実装は空間ベクトルの補助量 \(h\) を内部計算に使うが、公開出力は実数4状態の \(H\) である。`h[4][2]` の並びは

$$
H=
\begin{bmatrix}
H_{00} & H_{01}\\
H_{10} & H_{11}\\
H_{20} & H_{21}\\
H_{30} & H_{31}
\end{bmatrix}
$$

であり、行は \(\psi_{sd},\psi_{sq},\psi_{rd},\psi_{rq}\)、列は \(i_{sd}\) 誤差、\(i_{sq}\) 誤差に対応する。

C実装のコンパイル確認は以下で行った。

```powershell
gcc -std=c99 -Wall -Wextra -pedantic -c .\flux_observer_design\c\flux_observer.c
```

## 6. 誤差無し/有の評価結果

### 6.1 評価方法

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

### 6.2 無誤差時の収束

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

### 6.3 定数誤差の評価

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

### 6.4 電圧誤差・電流誤差の評価

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

### 6.5 結論

本設計では、回転dq座標の一次・二次磁束を状態とし、一次電流誤差で補正するフルオーダ磁束オブザーバを構成した。極配置法により、各動作点の凍結モデルでオブザーバ誤差極を \(-2200\) および \(-3410\ \mathrm{rad/s}\) に配置した。

無誤差では、一次・二次磁束および一次・二次電流の推定値は真値へ収束した。誤差あり評価では、指定範囲の定数誤差、電圧誤差、電流誤差に対して発散は見られず、安定性は維持された。精度面では、電圧ゲイン誤差、一次漏れインダクタンス誤差、相互インダクタンス誤差が二次磁束推定精度に対して支配的であった。
