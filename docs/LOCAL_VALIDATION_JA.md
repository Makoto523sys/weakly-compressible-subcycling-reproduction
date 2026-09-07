# ローカル実行記録（2026-09-07、進行中）

対象は [arXiv:2608.23110v1](https://arxiv.org/html/2608.23110v1) の4.2節。
**2026-09-08追記：4.2節の限定した再現を確認。[最終レビュー](BUBBLE_REVIEW_JA.md)を優先する。**
以下は検証の経過記録であり、「実行中」「保留」は当該時点の記述。多孔質は未実装。
引継ぎ時の「Juliaなし・全計算未実行」は現在には当てはまらない。

## 環境と保全

- Julia 1.12.7、Linux/WSL2、Core Ultra 7 265、20論理CPU、RAM約31 GiB。
- 権限付きnvidia-smiでRTX A400（4 GiB）を確認。今回の計算ではGPUを使用していない。
- 適用されるAGENTS.mdは見つからなかった。同梱ファイルは元のSHA256一覧と全件一致した。
- 元ソルバーは `evidence/local-20260907/source/` に保存。
- リモートmainは `fe301db516317568a7984acedb08591626088666`、追跡ファイルはMIT LICENSEのみ。
  `.remote-review`へ読取fetchした。公開前には再fetchが必要。
- 実行は `JULIA_DEPOT_PATH="$PWD/.julia_depot:/home/uh488080/.julia" OPENBLAS_NUM_THREADS=1`。
  シェルは `set -o pipefail` と `tee logs/<名前>.log` を使用。

## 実行した検証

1. 初回 `using WeaklyCompressibleSubcycling` 成功：`logs/01-load.log`。
2. 元のテスト21件成功：`logs/02-unit-tests.log`。
3. 静止二相気泡64/128/256格子、重力のみ0、他の物性・半径・領域は気泡問題と同一。
   p=u=0から0.001秒。これは4.2節の上昇計算とは別の検証問題。
4. 32×64、0.001秒のstandard/weak/subcycling完走：`logs/04-smoke-*.log`。
   weak/subcyclingの相分率上限は1.00294。完走だけで合格としなかった。
5. 256×512、0.07秒のstandard・weak・subcyclingがすべて完走。
6. multigrid前処理の対称性・正定値性、圧力解、3ステップの元実装との一致を追加。
   全27件成功：`logs/12-tests-final-mg.log`。
7. 丸め誤差による余分なsubstepを旧コードで再現し修正。チェックポイントの完全状態保存・再開も追加。
   合計48件成功：`logs/24-tests-checkpoint.log`。
8. 粗格子の密行列が巨大になる格子をmultigridが拒否するメモリ保護を追加。
   現在の合計49件成功：`logs/27-tests-coarse-memory-guard.log`。

## 静止気泡の実測

理論圧力差は2Dのsigma/R=57.6 Pa。内側r<R/2と外側r>2Rの平均圧力差を測定。
速度欄は `hypot(max(abs(uf)),max(abs(vf)))` という上界で、真の最大ベクトル速度とは異なる。

| NX | 最終圧力差 Pa | 相対誤差 | 0.001秒までの最大速度上界 m/s |
|---|---:|---:|---:|
| 64 | 59.92214 | 4.0315% | 0.00607001 |
| 128 | 58.68193 | 1.8783% | 0.00331559 |
| 256 | 58.39319 | 1.3771% | 0.00546389 |
| 512 | 58.33124 | 1.2695% | 0.00802508 |

気体面積の相対変化は約2e-16以下。256格子で時間刻み半減後の最大速度は0.00546147 m/s、
圧力rtolを1e-9から1e-11に変更しても0.00546389 m/s。
格子細分化で寄生流が単調減少するとは確認できていない。0.001秒は定常到達の保証ではない。
256格子を0.01秒へ延長すると最終圧力差58.38307 Pa（1.3595%誤差）、
最終速度上界0.00532224 m/s、全区間の最大0.00923862 m/s、面積変化1.30e-13となった。
`static-fine`の旧実行は終了時にソースハッシュを採る記録不備があったため検証根拠から除外し、
起動時ハッシュとソース保存を修正して`static-fine-fixed`で再実行完了（上表512行）。

`scripts/audit_surface_balance.jl` は初回投影の力の不釣り合いを切り分ける診断。
256格子で、論文CSFの投影後加速度上界8.0322 m/s²、解析的な一定曲率でも7.2122 m/s²、
離散ポテンシャル勾配では1.69e-6 m/s²となった。
後者は**診断専用で、本計算には採用していない**。論文式32–33の中点評価は厳密な離散勾配ではない。
実際、H(phi)=3phi²−2phi³に対し
`6*mean(phi)*(1-mean(phi))*Δphi = ΔH + (Δphi)^3/2`。
epsilon/dx=1を維持する格子細分化では、この誤差の単純な消失を仮定できない。

## 全区間を完了した気泡計算

すべて論文の領域・配置・境界条件・物性・零初期速度圧力を維持。初期拡散界面等の解釈は
`IMPLEMENTATION.md` に残す。参照CSVは図5の抽出曲線で、計算結果ではない。

| 計算 | 状態 | 速度比較 |
|---|---|---|
| weak、256×512 | 0.07秒完了、8628ステップ、約500秒 | 対応するweak参照へNRMSE 0.1925%、最大誤差0.001710 m/s |
| subcycling factor30、256×512 | 0.07秒完了、288主ステップ/8628 substeps、約697秒 | 標準参照へNRMSE 0.1965%、最大誤差0.0009975 m/s |
| standard、256×512 | multigrid版0.07秒完了、8628ステップ、約2150秒。Jacobi版は一致確認後に途中停止 | 標準参照へNRMSE 0.1812%、最大誤差0.0010004 m/s |
| subcycling、128×256 | 0.07秒完了 | 256格子との差RMS 0.004139、最大0.006244 m/s |
| subcycling factor15、256×512 | 修正版0.07秒完了、576主ステップ/8628 substeps | factor30との差RMS 0.0001340、最大0.0002512 m/s |
| subcycling、512×1024 | 旧版は途中停止、丸め修正版で再実行中 | 格子依存性は未確定 |

standardとsubcycling factor30は既存の速度・面積保存・相分率範囲・発散の暫定チェックを全て通過した。
判定ファイルは `evidence/local-20260907/standard-acceptance.toml` と `proposed-acceptance.toml`、状態は
`velocity_checks_passed_requires_review`。これを再現認定とはしない。
weakの相分率最大は1.005413、最大発散188.14 s^-1。弱圧縮性を許すモードであり、
投影後非圧縮条件の合格とは区別する。相分率超過の事実は隠さない。
standardの最大発散は4.80e-7 s^-1、相対気体面積変化1.51e-13。
subcyclingはそれぞれ1.21e-5 s^-1、4.44e-16。
0.03–0.07秒の標準参照曲線からの速度偏差のpeak-to-peakはweakで0.07207 m/s、
subcyclingで0.001012 m/s。これは周波数フィルタで分離した音響振幅ではない。

## 最適化と計算費用

プロファイルで圧力解法が支配的だった。追加した前処理は保存的な2×2集約、対称Jacobi平滑化、
零平均粗格子解を用いたV-cycle。元の細格子圧力方程式、真の残差確認、rtolは変更しない。
256×512の初期右辺でJacobi 1638反復/0.996秒、multigrid 89反復/0.241秒
（JITを除く一回の測定、他の計算と同時実行）。これは全計算の速度比でもH100との比較でもない。
共通の初期区間で元実装とmultigrid版の上昇速度差は約7e-12 m/s以下。

旧subcyclingは例えばfactor15に対して16個目のほぼ零時間幅のsubstepを作った。
残時間が安定限界を64機械イプシロン以内で上回る場合、残時間全体を最後のsubstepとする修正を追加。
旧データを上書きせず、`*-fixed`で時間・格子感度を再実行している。
factor15の修正前後は全区間の速度差がRMS 1.32e-9、最大2.86e-9 m/sだった。
128→256格子の速度差は参照最大速度のRMS約2.3%で、512格子の完了前には格子不確かさを認定できない。

別ディレクトリ`evidence/local-20260907/kernel-candidate`で境界チェック除去と8-thread化も試験した。
256格子3ステップで全状態がbitwise一致し、相・運動量流束の単体測定は約1.8倍だった。
これは実験候補であり、`src/`の正式ソルバーや現在の計算には採用していない。

初期静止気泡の標準法124ステップ/45.35秒からの0.07秒見積りは約53分だった。
実際の上昇計算では反復数と同時実行負荷が変化するため、この見積りは保証値ではない。
実測時間の比較には異なる圧力前処理と同時実行負荷が混在するため、性能再現の根拠にはしない。
多孔質の0.18秒には約677061毛管substepsが必要。
weakの実測500秒をセル数とsubstep数だけで外挿すると約78時間であり、
cut-cell/固体境界処理や圧力解法を含まない粗い見積りである。現時点で多孔質の実測値はない。

## 再実行と確認

```bash
julia --project=. test/runtests.jl
julia --project=. scripts/validate_static.jl results/static-new
julia --project=. scripts/validate_static.jl results/static-half-new 256 0.001 0.5 1e-9
julia --project=. scripts/audit_surface_balance.jl
julia --project=. scripts/benchmark_pressure.jl
julia --project=. scripts/run_bubble.jl standard 256 1 0.07 results/standard-new multigrid
julia --project=. scripts/run_bubble.jl subcycling 256 15 0.07 results/time-new multigrid
python3 scripts/review_local.py
python3 scripts/review_sensitivity.py
python3 scripts/collect_evidence.py
```

Python図生成にはNumPyとMatplotlibを使う。ソルバー本体はJulia標準ライブラリのみ。
`review_local.py` は途中データも描くが、比較区間・完了状態・入力ハッシュをJSONに記録する。
初期結果は元ソースを各結果の`src/`に保存している。比較処理はそのハッシュを実行メタデータと照合する。
新規ランナーは自動でソースを保存し、`snapshots.csv`にVTKと時刻の対応を記録する。
`collect_evidence.py` は完了した履歴・設定・実行ソースを集めてハッシュを照合する。
旧ソースによる再実行には汎用の `scripts/replay_run.jl` を追加した：

```bash
julia --project=evidence/local-20260907/runs/bubble-weak-initial scripts/replay_run.jl evidence/local-20260907/runs/bubble-weak-initial/run.toml results/weak-replay-new
```

旧weakソースで32格子の短時間計算を再実行し、元のhistory.csvとバイト単位で一致した。
修正版は各スナップショット時に`checkpoint.jls`を保存する。
第7引数にそのパスを指定すると、Julia版・ソースハッシュ・全設定・方式・factorが同じ場合だけ再開する。
再開先は新しい出力ディレクトリが必要。再開先の履歴は途中時刻から始まるため、単独で全区間比較に使えない。
この機能追加前の実行にはチェックポイントはない。VTKだけから完全状態を再構成した扱いにはしない。
作業ツリーのソースが更新された後は、実行に保存されたプロジェクトを使える：

```bash
julia --project=results/bubble-proposed-512-fixed scripts/run_bubble.jl subcycling 512 30 0.07 results/bubble-proposed-512-resumed multigrid results/bubble-proposed-512-fixed/checkpoint.jls
```

同じ実行がまだ稼働している間に再開して重複計算しない。再開後は旧・新履歴の時間範囲と
チェックポイント時刻の整合を確認して結合する必要がある。

## 未達成

512格子の全時間区間の評価、静止気泡の残留寄生流の総合評価は未完了。
このため4.2節の総合合格は保留。4.4節のcut-cell、immersed boundary、液体側150°接触角、
侵入経路・圧力比較は未実装・未実行。最新ユーザー指示により、気泡検証・問題解決後に公開し、
続けて多孔質へ進む。公開について追加確認は不要。

修正版512格子計算は継続しており、保存済みチェックポイントがある。
`scripts/watch_refinement.py` が完了後に比較JSON・図・根拠の収集を更新する。
この監視処理は再現認定、多孔質計算開始、GitHub公開を行わない。最新の計算完了状態は
`results/bubble-proposed-512-fixed/run.toml`、監視ログは `logs/38-refinement-watcher.log` を参照する。

標準法の時間刻み半減（256×512、factor0.5、0.07秒）の追加計算を開始：
`results/bubble-standard-half`、`logs/40-standard-half.log`。完了前のため時間誤差は未判定。


## 追加の切り分け（公開前レビュー）

初期CSFの投影だけを検査する `scripts/audit_static_width.jl` を実行した。
界面厚さを物理値epsilon=7.8125e-5 mで固定すると、NX=128/256/512で
投影後加速度の上界が4.46138/1.13263/0.284697 m/s²となり、ほぼ二次で減少する。
一方、epsilon=dxの条件では4.46138/8.03220/15.94620 m/s²。
この診断は物理厚さを固定した場合の離散化の切り分けであり、論文条件の置換ではない。
圧力残差の欄は絶対L2値で、相対停止条件の分母と単位を省いて良否を解釈しない。

時間依存性の追加計算としてsubcycling factor7.5も開始した（logs/42-factor7p5.log）。
速度ノルムの実装選択を検査するため、独立した診断ソースで接線成分を面へ補間し、
各面の速度大きさの最大をGammaに使う256格子計算を開始した（logs/43-face-norm.log）。
この選択を正式ソルバーに採用したわけではない。結果は完了後に比較する。
