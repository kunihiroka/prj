# 5000 r/min 回生時に電流制御が不安定化するメカニズム解析

## 結論

5000 r/min、220 Nm級の回生で電流制御が不安定化する主因は、電流PI単体の高周波安定性ではなく、以下の低周波結合ループである。

```text
軸ずれ Δδ
→ dq電流誤差 Δed, Δeq
→ PI電圧 Δvd, Δvq
→ 実 Δiq, Δφ
→ 実すべり変動 Δωslip,real
→ Δδdot
→ 角度積分で Δδ
```

この軸ずれループは数Hz帯にゲイン交差を持つ。力行と回生でゲインは大きく変わらないが、回生では位相が $-180^\circ$ 方向へ寄るため、安定余裕が小さくなる。今回の縮約ループでは、

| 条件 | 0 dB交差周波数 | 位相 |
|---|---:|---:|
| 力行 +220 Nm | 3.52 Hz | +40.1 deg |
| 回生 -220 Nm | 3.42 Hz | -143.7 deg |

である。

さらに、全閉ループをMIMO return ratioで見ると、回生だけ

$$
\det(I+L_i(j\omega))
$$

が4.5 Hz付近で原点へ接近する。

| 条件 | $\min|\det(I+L_i)|$ | 周波数 |
|---|---:|---:|
| 力行 +220 Nm | 1.000 | 10000 Hz |
| 回生 -220 Nm | 0.332 | 4.51 Hz |

したがって、「回生だけ不安定化する」現象は、軸ずれループの位相悪化として直感的に説明でき、MIMO特性方程式の周波数応答で定量的にも再現できる。

実務上の対策は、滑り計算を `id_ref` 起点の磁束だけに依存させないことである。`id_feedback` または実磁束を滑り計算に入れると、実磁束変動が滑り指令へ戻り、危険な数Hz帯の低周波モードが抑制される。

## 1. 対象条件と現象

対象条件は以下である。

| 項目 | 値 |
|---|---:|
| $R_s$ | 0.00762 ohm |
| $R_r$ | 0.008041 ohm |
| $L_{ls}$ | 0.0000419 H |
| $L_{lr}$ | 0.0000419 H |
| $L_m$ | 0.0001583 H |
| 極対数 | 4 |
| 速度 | 5000 r/min |
| 電流振幅 | $i_d=|i_q|=673.6$ A |
| 電流PI | $C_i(s)=K_p(1+K_i/s)$ |
| $K_p$ | $\sigma L_s\omega_{cc}$ |
| $K_i$ | $R_s/(\sigma L_s)$ |
| $\omega_{cc}$ | 1000 rad/s |

標準設定では、磁束推定と滑り計算は指令ベースで作られる。

$$
\hat{\phi}
=
\mathrm{LPF}(L_m i_d^*)
$$

$$
\omega_{\mathrm{slip,cmd}}
=
\frac{R_r L_m}{L_r}
\frac{i_q^*}{\hat{\phi}}
$$

つまり、実際の $i_d$、$i_q$、実ロータ磁束が変動しても、その変動は滑り指令に直接戻らない。

現象としては、同じ速度・同じ電流振幅でも、力行より回生の方が低周波振動を起こしやすい。連続系モデルでPIを

$$
C_i(s)=K_p\left(1+\frac{K_i}{s}\right)
$$

にした場合、回生側で数Hz帯の発散または低減衰振動が出る。

![regen waveform](./continuous_foc_regen_pi_kp_1_ki_over_s_2s.png)

対照として、力行側は同じPI構成でも安定しやすい。

![motoring waveform](./continuous_foc_motoring_pi_kp_1_ki_over_s_2s.png)

## 2. 軸ずれループによる原因説明

制御角と実ロータ磁束角の差を

$$
\delta
=
\theta_{\mathrm{ctrl}}
-
\theta_{\mathrm{flux}}
$$

とする。角速度差から、

$$
\dot{\delta}
=
\omega_{\mathrm{slip,cmd}}
-
\omega_{\mathrm{slip,real}}
$$

である。実ロータ磁束基準の滑り周波数は、

$$
\omega_{\mathrm{slip,real}}
=
K\frac{i_{q,\mathrm{actual}}}{\phi_{\mathrm{actual}}},
\qquad
K=\frac{R_r L_m}{L_r}
$$

である。

### 2.1 軸ずれから電流誤差への変換

微小な軸ずれがあると、測定dq電流は小角近似で

$$
\begin{aligned}
\Delta i_{d,\mathrm{meas}}
&\simeq
I_{q,\mathrm{actual}}\Delta\delta \\
\Delta i_{q,\mathrm{meas}}
&\simeq
-I_{d,\mathrm{actual}}\Delta\delta
\end{aligned}
$$

となる。電流誤差は

$$
e_d=i_d^*-i_{d,\mathrm{meas}},
\qquad
e_q=i_q^*-i_{q,\mathrm{meas}}
$$

なので、

$$
\begin{aligned}
\Delta e_d
&\simeq
-I_{q,\mathrm{actual}}\Delta\delta \\
\Delta e_q
&\simeq
+I_{d,\mathrm{actual}}\Delta\delta
\end{aligned}
$$

である。したがって

$$
M_{\delta e}
=
\begin{bmatrix}
-I_{q,\mathrm{actual}} \\
I_{d,\mathrm{actual}}
\end{bmatrix}
$$

と書ける。ここが力行と回生で符号が変わる第一の場所である。

| 条件 | $I_q$ | $\Delta\delta\rightarrow\Delta e_d$ |
|---|---:|---|
| 力行 | $I_q>0$ | 逆相 |
| 回生 | $I_q<0$ | 同相 |

この符号反転により、同じ軸ずれでも、PIが作る $v_d$ の向きが力行と回生で反転する。

![axis mixing vector diagram](./axis_mixing_motoring_regen_vector_diagram.png)

### 2.2 実すべり変動への変換

実滑りの小信号は、動作点まわりの一次近似で

$$
\Delta\omega_{\mathrm{slip,real}}
=
\frac{K}{\Phi}\Delta i_{q,\mathrm{actual}}
-
\frac{K I_{q,\mathrm{actual}}}{\Phi^2}
\Delta\Phi_{\mathrm{actual}}
$$

となる。これを

$$
S_{\omega}
=
\begin{bmatrix}
0 &
\frac{K}{\Phi} &
-\frac{K I_{q,\mathrm{actual}}}{\Phi^2}
\end{bmatrix}
$$

と書けば、

$$
\Delta\omega_{\mathrm{slip,real}}
=
S_{\omega}
\begin{bmatrix}
\Delta i_d \\
\Delta i_q \\
\Delta\Phi
\end{bmatrix}
$$

である。ここでも $I_q$ の符号が入る。

### 2.3 縮約軸ずれループ

電流PIを

$$
C_i(j\omega)
=
K_p\left(1+\frac{K_i}{j\omega}\right)
$$

モータの小信号応答を

$$
\begin{bmatrix}
\Delta i_d \\
\Delta i_q \\
\Delta\Phi
\end{bmatrix}
=
G_m(j\omega)
\begin{bmatrix}
\Delta v_d \\
\Delta v_q
\end{bmatrix}
$$

とすると、軸ずれから軸ずれへ戻る縮約ループは

$$
L_{\delta}(j\omega)
=
\frac{
-
S_{\omega}
G_m(j\omega)
C_i(j\omega)
M_{\delta e}
}{j\omega}
$$

である。最後の $1/(j\omega)$ は $\Delta\dot{\delta}$ から $\Delta\delta$ への角度積分で、常に $-90^\circ$ を加える。

力行と回生の $L_{\delta}$ は以下である。

![L delta motoring regen compare](./l_delta_motoring_regen_compare.png)

重要なのは、ゲイン差ではなく位相差である。0 dB交差周波数は力行と回生でほぼ同じだが、回生だけ $-180^\circ$ 側に寄る。

| 運転 | 0 dB交差周波数 | $\angle L_\delta$ | $-1+j0$ への距離 |
|---|---:|---:|---:|
| 力行 +220 Nm | 3.519 Hz | +40.1 deg | 0.982 |
| 回生 -220 Nm | 3.422 Hz | -143.7 deg | 0.623 |

### 2.4 位相の足し上げ

3.42 Hz付近で、単位軸ずれ $\Delta\delta=1$ を入れた場合の位相を順に追う。

力行では以下になる。

| 段 | 位相 |
|---|---:|
| $\Delta e_d$ | +180 deg |
| $\Delta e_q$ | 0 deg |
| $C_i$ | -78 deg |
| $\Delta v_d$ | +102 deg |
| $\Delta v_q$ | -78 deg |
| $G_m$ 後の $\Delta i_q$ | -76 deg |
| $G_m$ 後の $\Delta\Phi$ | -90 deg |
| $\Delta\omega_{\mathrm{slip,real}}$ | -51 deg |
| $\Delta\dot{\delta}=-\Delta\omega_{\mathrm{slip,real}}$ | +129 deg |
| $\Delta\delta=\Delta\dot{\delta}/j\omega$ | +39 deg |

回生では以下になる。

| 段 | 位相 |
|---|---:|
| $\Delta e_d$ | 0 deg |
| $\Delta e_q$ | 0 deg |
| $C_i$ | -78 deg |
| $\Delta v_d$ | -78 deg |
| $\Delta v_q$ | -78 deg |
| $G_m$ 後の $\Delta i_q$ | +102 deg |
| $G_m$ 後の $\Delta\Phi$ | -90 deg |
| $\Delta\omega_{\mathrm{slip,real}}$ | +126 deg |
| $\Delta\dot{\delta}=-\Delta\omega_{\mathrm{slip,real}}$ | -54 deg |
| $\Delta\delta=\Delta\dot{\delta}/j\omega$ | -144 deg |

したがって、回生だけ位相余裕が小さくなる理由は、以下の順序で説明できる。

1. 回生では $\Delta\delta\rightarrow\Delta e_d$ が同相になる。
2. そのため $\Delta v_d$ と $\Delta v_q$ がほぼ同相で出る。
3. 高速回転座標の $G_m$ により、$\Delta i_q$ が力行とほぼ反対側の位相に出る。
4. $\Delta i_q$ と $\Delta\Phi$ から作る $\Delta\omega_{\mathrm{slip,real}}$ の位相が回生側で $+126^\circ$ になる。
5. $\Delta\dot{\delta}=-\Delta\omega_{\mathrm{slip,real}}$ と角度積分 $1/(j\omega)$ により、最終的な $L_\delta$ が $-144^\circ$ となる。

つまり、回生だけが悪い理由は「回生のゲインが大きい」ことではなく、$I_q$ 符号反転を起点にして、軸ずれループの戻り位相が $-1+j0$ 方向へ回ることである。

### 2.5 $\Delta i_q$ と $\Delta\Phi$ の寄与

$\Delta i_q$ への寄与を

$$
\Delta i_q
=
G_{qd}\Delta v_d
+
G_{qq}\Delta v_q
$$

に分けると以下になる。

![iq contribution bode](./iq_contribution_bode.png)

3.42 Hz付近では、

| 運転 | 経路 | 振幅 | 位相 |
|---|---|---:|---:|
| 力行 | $G_{qd}\Delta v_d$ | 700 A/rad | -67.8 deg |
| 力行 | $G_{qq}\Delta v_q$ | 311 A/rad | -94.9 deg |
| 力行 | 合成 | 987 A/rad | -76.1 deg |
| 回生 | $G_{qd}\Delta v_d$ | 758 A/rad | +111.8 deg |
| 回生 | $G_{qq}\Delta v_q$ | 311 A/rad | +78.3 deg |
| 回生 | 合成 | 1030 A/rad | +102.3 deg |

回生では $v_d$ 経由と $v_q$ 経由が相殺せず、ほぼ同方向に加算される。

$\Delta\Phi$ への寄与は以下である。

![phi contribution bode](./phi_contribution_bode.png)

3.42 Hz付近では、

| 運転 | 経路 | 振幅 | 位相 |
|---|---|---:|---:|
| 力行 | $G_{\phi d}\Delta v_d$ | 0.0273 Wb/rad | -100.9 deg |
| 力行 | $G_{\phi q}\Delta v_q$ | 0.0787 Wb/rad | -86.4 deg |
| 力行 | 合成 | 0.105 Wb/rad | -90.1 deg |
| 回生 | $G_{\phi d}\Delta v_d$ | 0.0320 Wb/rad | -95.5 deg |
| 回生 | $G_{\phi q}\Delta v_q$ | 0.0818 Wb/rad | -87.2 deg |
| 回生 | 合成 | 0.114 Wb/rad | -89.5 deg |

$\Delta\Phi$ 応答自体は力行と回生で大きく変わらない。力行/回生差を強く作っているのは、$\Delta i_q$ の位相と、滑り式に含まれる $I_q$ 符号である。

### 2.6 高速ほど悪化する理由

速度を振ると、回生側の交差位相は高回転ほど $-180^\circ$ 側へ寄る。

![L delta speed sweep](./l_delta_speed_sweep_bode.png)

| 運転 | 速度 | 0 dB交差周波数 | 位相 |
|---|---:|---:|---:|
| 力行 | 1000 r/min | 18.6 Hz | +98.6 deg |
| 力行 | 3000 r/min | 4.91 Hz | +53.4 deg |
| 力行 | 5000 r/min | 3.52 Hz | +40.1 deg |
| 力行 | 6000 r/min | 3.15 Hz | +36.3 deg |
| 回生 | 1000 r/min | 9.47 Hz | -134.5 deg |
| 回生 | 3000 r/min | 4.62 Hz | -135.4 deg |
| 回生 | 5000 r/min | 3.42 Hz | -143.8 deg |
| 回生 | 6000 r/min | 3.09 Hz | -146.6 deg |

この縮約ループだけでは不安定を最終判定しない。ナイキスト上では5000 r/min条件でも $-1+j0$ を包囲しないためである。ただし、「高速回生ほど危険な方向へ移動する」ことは明確に説明できる。

## 3. MIMO特性方程式による不安定周波数応答の再現

$L_\delta$ は、現象を理解するために軸ずれループを1本へ縮約したモデルである。実際の電流制御は $d/q$ の2入力2出力系なので、安定性はMIMO return ratioで見る必要がある。

PI入力でループを切り、PI出力電圧から検出dq電流までを

$$
\Delta i_{\mathrm{meas}}(s)
=
G_i(s)\Delta v_{\mathrm{PI}}(s)
$$

と定義する。$G_i(s)$ には、モータ磁束状態、磁束推定、滑り計算、非干渉項、制御角、dq座標変換を含める。

![MIMO current-loop return ratio block diagram](./mimo_return_ratio_block_diagram.png)

電流PIを

$$
C_i(s)=K_p\left(1+\frac{K_i}{s}\right)
$$

とすると、MIMO return ratioは

$$
L_i(s)=C_i(s)G_i(s)
$$

である。閉ループの自由振動条件は

$$
(I+L_i(s))\Delta i=0
$$

であり、非ゼロ解が存在する条件は

$$
\det(I+L_i(s))=0
$$

である。したがって、MIMOの特性方程式は

$$
\boxed{\det(I+L_i(s))=0}
$$

である。

`det(I+L_i)` を直接ナイキスト表示する場合、危険点は $-1$ ではなく原点である。

### 3.1 力行/回生比較

5000 r/min、`id_ref` 起点滑り、磁束フィードバック非干渉あり条件で、力行と回生を比較すると以下になる。

![motoring vs regen MIMO current-loop return ratio](./full_current_loop_motoring_regen_compare.png)

| 運転状態 | $\min|\det(I+L_i)|$ | 周波数 | 判定 |
|---|---:|---:|---|
| 力行 +220 Nm | 1.000 | 10000 Hz | 低周波の危険な谷なし |
| 回生 -220 Nm | 0.332 | 4.51 Hz | 4.5 Hz付近で原点へ接近 |

これは、縮約 $L_\delta$ で示した「回生だけ位相余裕が小さい」という説明と整合する。MIMO特性方程式でも、回生だけ数Hz帯に危険な谷が現れる。

### 3.2 滑り計算の磁束入力による差

滑り計算に使う磁束を `id_ref` 起点から、`id_feedback` 起点、さらに実ロータ磁束へ近づけると、低周波の危険な谷は小さくなる。

![full current-loop return ratio](./full_current_loop_return_ratio.png)

| 滑り計算の磁束 | $\min|\det(I+L_i)|$ | 周波数 |
|---|---:|---:|
| $i_{d,\mathrm{ref}}$ 起点 | 0.332 | 4.52 Hz |
| $i_{d,\mathrm{feedback}}$ 起点 | 0.879 | 5.59 Hz |
| 実ロータ磁束 | 1.000 | 10000 Hz |

`id_ref` 起点では、実磁束が変動しても滑り指令へ戻らない。一方、`id_feedback` または実磁束を使うと、実磁束変動が滑り指令側へ戻るため、

$$
\omega_{\mathrm{slip,cmd}}
\simeq
\omega_{\mathrm{slip,real}}
$$

に近づく。これは位相進み補償ではなく、滑り誤差を作る低周波ゲインを構造的に下げる効果である。

### 3.3 速度依存

`id_ref` 起点滑りのまま速度を振ると、高速ほど低周波の谷が深くなる。

![full current-loop speed sweep idref](./full_current_loop_speed_sweep_idref.png)

| 速度 | $\min|\det(I+L_i)|$ | 周波数 |
|---:|---:|---:|
| 1000 r/min | 1.000 | 10000 Hz |
| 3000 r/min | 1.000 | 10000 Hz |
| 5000 r/min | 0.332 | 4.51 Hz |
| 6000 r/min | 0.128 | 3.93 Hz |

これは「速度が高いほど不安定化しやすい」という観測と一致する。

### 3.4 電流PIゲイン依存

`id_ref` 起点滑り、5000 r/min回生で電流PIの比例ゲインを振ると以下になる。

![full current-loop kp sweep idref](./full_current_loop_kp_sweep_idref.png)

| $K_p$倍率 | $\min|\det(I+L_i)|$ | 周波数 |
|---:|---:|---:|
| 0% | 1.000 | 0.100 Hz |
| 50% | 0.188 | 2.61 Hz |
| 100% | 0.332 | 4.51 Hz |
| 300% | 1.002 | 10000 Hz |

実験で「電流制御ゲインを0、または3000 rad/s相当にすると安定化する」という傾向は、この図で説明できる。中途半端なゲインでは数Hz帯のMIMO特性方程式が原点へ近づき、危険な低周波モードを作る。一方、ゲイン0ではPIが軸ずれループを励起せず、非常に高いゲインでは低周波の電流偏差を強く抑え込む。

### 3.5 固有値による確認

全閉ループ状態を線形化して固有値を見ると、非干渉なしの連続系では回生だけ右半平面極が出る。

![full closed-loop eigenvalues](./full_closed_loop_eigenvalues.png)

代表値は以下である。

| 条件 | 支配極 | 周波数 | 判定 |
|---|---:|---:|---|
| 力行 / no decoupling / id_ref slip | $-16.35 \pm j15.53$ | 2.47 Hz | 安定 |
| 回生 / no decoupling / id_ref slip | $+0.863 \pm j20.89$ | 3.33 Hz | 不安定 |
| 回生 / no decoupling / id_feedback slip | $-3.70 \pm j19.65$ | 3.13 Hz | 安定化 |
| 回生 / no decoupling / actual flux slip | $-2.29 \pm j15.54$ | 2.47 Hz | 安定化 |

非干渉を含めた簡略モデルでは右半平面極までは出ないが、回生側の減衰は力行より小さく、高速ほど安定余裕が小さい。

![full closed-loop speed torque sweep](./full_closed_loop_speed_torque_sweep.png)

## 4. 対策

### 4.1 滑り計算を指令磁束だけに依存させない

最も本質的な対策は、滑り計算の磁束を `id_ref` 起点だけにしないことである。

推奨順は以下。

| 対策 | 効果 |
|---|---|
| 実ロータ磁束を滑り計算に使う | 最も直接的。実すべりと指令すべりが一致しやすい |
| $i_d$ フィードバックから磁束推定する | 実機で実装しやすい。低周波滑り誤差ゲインを下げる |
| $i_q$ フィードバックを滑り分子に混ぜる | $i_q$ ずれによる滑り誤差を減らす |

`id_feedback` は位相を進ませる補償ではない。磁束推定器は

$$
\frac{\hat{\phi}}{i_{d,\mathrm{fb}}}
=
\frac{L_m}{1+sT_r}
$$

であり、単体では一次遅れである。それでも安定化する理由は、実 $i_d$ 変動を滑り指令側へ戻すことで、

$$
\omega_{\mathrm{slip,cmd}}
-
\omega_{\mathrm{slip,real}}
$$

の低周波ゲインを下げるためである。

### 4.2 電流PI設計

`id_ref` 起点滑りを使い続ける場合、問題となる数Hz帯で電流偏差を十分に抑える必要がある。低周波電流プラントを

$$
P_{\mathrm{current}}(s)
\simeq
\frac{1}{\sigma L_s s+R_s}
$$

と近似すると、問題周波数 $\omega_{\mathrm{bad}}$ で

$$
\left|
C_i(j\omega_{\mathrm{bad}})
P_{\mathrm{current}}(j\omega_{\mathrm{bad}})
\right|
\ge
L_{\mathrm{req}}
$$

を満たすようにする。

低周波でI項支配なら、

$$
C_i(j\omega)
\simeq
\frac{K_pK_i}{j\omega}
$$

なので、

$$
K_i
\ge
\frac{
L_{\mathrm{req}}
\omega_{\mathrm{bad}}
R_s
}{K_p}
$$

となる。ただし、PIゲインだけで押さえ込む方法は、モータ定数誤差、電圧制限、検出遅れ、ノイズに弱い。設計上は、滑り計算側にフィードバックを入れて低周波正帰還経路を構造的に弱める方が妥当である。

### 4.3 実機/Simulinkでの確認方法

不安定動作点では、閉ループのまま長時間のsinestream/chirpを入れると動作点が崩れる。そのため、測定時は電流PIループを開いて、PI出力電圧相当へ小信号を注入する。

```text
vd_PI = vdPI0 + d_inj
vq_PI = vqPI0 + q_inj

入力: d_inj, q_inj
出力: id_meas, iq_meas
```

これにより

$$
G_i(j\omega)
=
\begin{bmatrix}
i_d/v_d & i_d/v_q \\
i_q/v_d & i_q/v_q
\end{bmatrix}
$$

を測定できる。その後、MATLAB上で

$$
L_i(j\omega)
=
C_i(j\omega)G_i(j\omega)
$$

を作り、

$$
\det(I+L_i(j\omega))
$$

を評価すればよい。

## 5. まとめ

本現象の説明は以下の流れで整理できる。

1. 回生では $I_q<0$ となり、軸ずれからd軸電流誤差への符号が力行と反転する。
2. その結果、PIが作る $v_d/v_q$ の位相関係が変わる。
3. 高速回転座標のモータ応答 $G_m$ を通ると、回生側では $\Delta i_q$ と $\Delta\Phi$ から作られる実すべり変動が、軸ずれを戻す方向から外れる。
4. 角度積分 $1/(j\omega)$ がさらに $-90^\circ$ を加え、回生側の $L_\delta$ は $-180^\circ$ 方向へ近づく。
5. MIMO特性方程式 $\det(I+L_i)=0$ で見ると、回生だけ4.5 Hz付近に低周波の谷が出る。
6. 滑り計算に実磁束または $i_d$ フィードバックを入れると、この低周波の谷が浅くなり、安定方向へ移る。

したがって、問題の本質は「電流PIが悪い」ではなく、**指令ベース滑りと実磁束の切り離しにより、電流PIが軸ずれ/滑りの低周波ループを励起すること**である。

## 付録A. 滑り周波数小信号式の導出

実滑りを

$$
\omega_{\mathrm{slip,real}}
=
K\frac{i_q}{\Phi}
$$

とする。動作点を

$$
i_q=I_q,
\qquad
\Phi=\Phi_0
$$

とし、微小変動を

$$
i_q=I_q+\Delta i_q,
\qquad
\Phi=\Phi_0+\Delta\Phi
$$

と置く。

多変数テイラー展開の一次項までを取ると、

$$
\Delta\omega_{\mathrm{slip,real}}
\simeq
\left.
\frac{\partial}{\partial i_q}
K\frac{i_q}{\Phi}
\right|_0
\Delta i_q
+
\left.
\frac{\partial}{\partial \Phi}
K\frac{i_q}{\Phi}
\right|_0
\Delta\Phi
$$

である。偏微分は

$$
\left.
\frac{\partial}{\partial i_q}
K\frac{i_q}{\Phi}
\right|_0
=
\frac{K}{\Phi_0}
$$

$$
\left.
\frac{\partial}{\partial \Phi}
K\frac{i_q}{\Phi}
\right|_0
=
-
\frac{K I_q}{\Phi_0^2}
$$

なので、

$$
\Delta\omega_{\mathrm{slip,real}}
\simeq
\frac{K}{\Phi_0}\Delta i_q
-
\frac{K I_q}{\Phi_0^2}\Delta\Phi
$$

を得る。

ここで重要なのは、必要なのは

$$
K\Delta\left(\frac{i_q}{\Phi}\right)
$$

であって、

$$
K\frac{\Delta i_q}{\Delta\Phi}
$$

ではないことである。後者は変化量どうしの比であり、滑り周波数の小信号変化ではない。

## 付録B. 数値線形化と周波数応答

非線形モデルを

$$
\dot{x}=f(x,u)
$$

とする。動作点 $(x_0,u_0)$ まわりで

$$
x=x_0+\Delta x,
\qquad
u=u_0+\Delta u
$$

と置き、一次テイラー展開すると

$$
\Delta\dot{x}
=
A\Delta x+B\Delta u
$$

となる。ここで

$$
A=
\left.
\frac{\partial f}{\partial x}
\right|_{x_0,u_0},
\qquad
B=
\left.
\frac{\partial f}{\partial u}
\right|_{x_0,u_0}
$$

である。出力を

$$
\Delta y=C\Delta x+D\Delta u
$$

とすれば、周波数応答は

$$
G(j\omega)
=
C(j\omega I-A)^{-1}B+D
$$

で求められる。

今回のMIMO return ratioでは、入力をPI出力電圧、出力を検出dq電流として

$$
G_i(j\omega)
=
\frac{
\begin{bmatrix}
\Delta i_d \\
\Delta i_q
\end{bmatrix}
}{
\begin{bmatrix}
\Delta v_{d,\mathrm{PI}} \\
\Delta v_{q,\mathrm{PI}}
\end{bmatrix}
}
$$

を求め、その後

$$
L_i(j\omega)
=
C_i(j\omega)G_i(j\omega)
$$

として評価した。
