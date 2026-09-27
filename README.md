# matlab-simscape-tracer

液体の流れと一緒にトレーサー濃度（パッシブスカラー）を運ぶ、Simscape用のカスタムドメインとコンポーネントです。
同じ式をMATLABだけで解く試作と、その検証スクリプトも含みます。

## 状態

| 対象 | 状態 |
|---|---|
| `+mylib/*.ssc`（Simscape版） | **未検証**。Simscapeのない環境で作成したため、`ssc_build` もシミュレーションも未実施 |
| `tracer_*.m`（MATLAB版） | `tracer_verify.m` で理論値との一致を確認済み |

## Simscapeライブラリ（`+mylib`）

| ファイル | 内容 |
|---|---|
| `tracer_liquid.ssc` | ドメイン。Across：圧力 p、濃度 c（質量分率）／Through：質量流量 mdot、トレーサー質量流量 mdot_c。物性はドメインのパラメータ |
| `tracer_fluid_properties.ssc` | 物性（ρ₀、基準圧力、β、動粘度 ν、拡散係数 D_c）を設定するブロック |
| `tracer_orifice.ssc` | Local Restriction (IL) と同じ式。Re_cr から Δp_cr を計算、面積比の補正、圧力回復 |
| `tracer_pipe.ssc` | 配管の摩擦（層流 Hagen–Poiseuille／遷移／乱流 Haaland）。体積は持たない |
| `tracer_volume.ssc` | 一定体積。密度は圧力の関数、トレーサーは完全混合 |
| `tracer_reservoir.ssc` | 圧力と濃度を固定する境界 |
| `tracer_injector.ssc` | 物理信号で注入流量を与える注入源 |

トレーサーは風上差分（tanh で平滑化）と拡散項で運びます。

### ビルド（Simscapeが必要）

```matlab
cd <このリポジトリ>
ssc_build mylib   % mylib_lib.slx が生成される
mylib_lib
```

抵抗系（オリフィス・配管）と体積を交互に並べて使います。

## MATLAB版の試作と検証

| ファイル | 内容 |
|---|---|
| `tracer_proto.m` | 最初の試作（固定物性、オリフィスのみ） |
| `tracer_study.m` | 入口条件の比較、分割数 N の影響、RTD（槽列モデルとの比較） |
| `tracer_net.m` | 拡張版の試作。圧力依存の密度、オリフィス、配管摩擦、拡散。ode15s で DAE として解く |
| `tracer_verify.m` | 拡張版の検証（オリフィス、摩擦係数、質量保存、RTD、流量ゼロでの拡散） |

```matlab
tracer_study    % study1〜3_*.png を出力
tracer_verify   % verify1〜4_*.png を出力
```
