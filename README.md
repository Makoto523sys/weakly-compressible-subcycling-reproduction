# weakly-compressible-subcycling-reproduction

ローカルCodexへの引継ぎは [docs/HANDOFF_LOCAL_JA.md](docs/HANDOFF_LOCAL_JA.md) を最初に読む。実行場所については、この最新のローカル方針を優先する。

Juliaによる独立再現です。**4.2節の気泡上昇速度・振動抑制を検証済み。4.4節の多孔質は未実装です。**
精度と制約を含む[最終レビュー](docs/BUBBLE_REVIEW_JA.md)を参照してください。
2026-09-07の実測・未達事項・実行ログは [docs/LOCAL_VALIDATION_JA.md](docs/LOCAL_VALIDATION_JA.md) を参照してください。

対象論文：Shu Yamashita, Shintaro Matsushita, Tetsuya Suekane,
[Weakly Compressible Subcycling for Accelerating Simulations of Surface-Tension-Dominated Incompressible Two-Phase Flows](https://arxiv.org/abs/2608.23110), arXiv:2608.23110v1.

目的は4.2節の二次元気泡上昇で手法を検証した後、4.4節の多孔質媒体内二相流へ進むことです。
著者の公式実装ではなく、本文・付録をもとにした独立実装です。

## 現在の状態

| 項目 | 状態 |
|---|---|
| Juliaソルバー | Julia 1.12.7で読込・実行確認 |
| 標準非圧縮／EPP単独／subcycling | 3モードを記述済み |
| 保存性・圧力補正等のテスト | 既存21件＋前処理・substep丸め・再開・メモリ保護の検証、合計49件が合格 |
| 図5の参照曲線 | SVGから抽出済み、出典・校正式・ハッシュを保存 |
| 気泡上昇の論文比較 | 3方式256×512・0.07秒、512格子・時間感度のレビュー完了。NRMSE約0.18–0.21%、細格子への速度感度約0.82% |
| 多孔質媒体 | 図10・本文の条件整理済み。ソルバー未実装 |
| GitHub | 気泡の検証結果を指定先へ公開する段階 |

ローカルではJulia 1.12.7、Core Ultra 7 265、約31 GiB RAMを確認しました。
静止気泡には圧力差約1.3%の誤差と寄生流が残ります。原因診断と、上昇問題の格子・時間・保存性検証を含めて再現範囲を判断しました。

## 実行方法

検証済み環境はJulia 1.12.7で、Manifest.tomlもこの版で固定しています。
互換性宣言は1.10以上ですが、他のJulia版は未検証です。実行依存は標準ライブラリのみです。
以下の単体テストと短時間計算はローカルで実行済みです。既存の出力名は再使用できません。

```bash
julia --project=. test/runtests.jl
julia --project=. scripts/run_bubble.jl standard 32 1 0.001 results/smoke-standard
julia --project=. scripts/run_bubble.jl subcycling 32 30 0.001 results/smoke-proposed
```

テスト・小規模計算を確認した後、論文と同じ256×512格子・0.07秒で比較します。
weak単独は約500秒、subcyclingは約697秒で0.07秒まで完走しました（同時実行負荷あり）。
標準法はmultigrid前処理で約2150秒でした。前処理と同時実行負荷が異なるため、性能比の比較には使いません。
標準法には任意の第6引数 `multigrid` で検証済みのPCG前処理を選べます。

```bash
julia --project=. scripts/run_bubble.jl standard 256 1 0.07 results/bubble-standard multigrid
julia --project=. scripts/run_bubble.jl weak 256 1 0.07 results/bubble-weak
julia --project=. scripts/run_bubble.jl subcycling 256 30 0.07 results/bubble-proposed
julia --project=. scripts/compare_bubble.jl results/bubble-standard reference/bubble_standard.csv reference/bubble_standard.toml results/standard-comparison.toml
julia --project=. scripts/compare_bubble.jl results/bubble-proposed reference/bubble_standard.csv reference/bubble_standard.toml results/proposed-comparison.toml
```

引数は `MODE NX FACTOR T_END OUTPUT`、`NY=2NX`、SI単位です。
任意の第6引数は `jacobi`（既定）または `multigrid` です。物理条件と真の圧力残差基準は共通です。
既存の出力ディレクトリは上書きしません。

出力は `history.csv`、ParaViewで読める `field_*.vtk`、実行条件・Juliaバージョン・
ソースハッシュ・終了状態を記録した `run.toml` です。計算終了だけでは検証合格になりません。
自動比較の合格は認定ではありません。この版の総合判断は最終レビューに記録しています。

実測の[3方式重ね合わせ](evidence/local-20260907/bubble-overlay.png)、
[判定結果・限界](docs/LOCAL_VALIDATION_JA.md)、
[標準法の暫定判定](evidence/local-20260907/standard-acceptance.toml)、
[subcyclingの暫定判定](evidence/local-20260907/proposed-acceptance.toml)を保存しています。

## 実装と論文との差

- ACDI相輸送、相流束と整合する運動量輸送、WENO3、局所化CSFを記述しています。
- 小ステップで相・運動量の流束発散と力を時間積分し、主ステップを初期状態から再構成して圧力補正します。
- 論文のhypre FlexGMRES/PFMGに対し、この版は零平均空間でのJacobi前処理PCGです。性能の直接比較には使えません。
- 初期拡散界面、WENOの微小定数、速度最大値の評価、重力と壁面圧力の扱いには実装上の選択があります。
  詳細は [docs/IMPLEMENTATION.md](docs/IMPLEMENTATION.md) を参照してください。
- 相分率を輸送後に強制クリップして保存誤差を隠す処理はありません。

## 再現判定と次の作業

[docs/VALIDATION.md](docs/VALIDATION.md) に判定基準を示します。
4.4節は気泡上昇の検証後に実装します。円柱を単なる格子マスクに置き換えたり、
接触角を省略した計算を論文再現と扱うことはしません。

`.github/workflows/validate.yml` は、Juliaを利用できるGitHub Actionsで
テストと短時間計算を実施するための未実行ワークフローです。
最新のユーザー指示により、気泡上昇の検証・問題解決後に公開し、その後に多孔質流れへ進みます。
論文規模の計算は手動起動のみで、各ジョブに時間上限を設けます。

アップロード先：[Makoto523sys/weakly-compressible-subcycling-reproduction](https://github.com/Makoto523sys/weakly-compressible-subcycling-reproduction)
