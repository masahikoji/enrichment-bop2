# GitHub 更新手順（Mac・2026-10-25版）

今回のZIPは、前回版の更新版です。**前回版と今回版の両方を重ねる必要はありません。** 同梱の解析コードの主な計算アルゴリズムは変更せず、原著BOP2の検証結果・S5感度解析のCSV・文書を追加しました。

## 更新

ZIPを展開すると `enrichment-bop2-github-update-2026-10-25` フォルダができます。ターミナルで、既存リポジトリへその中身をコピーしてください。

```bash
cd ~/Downloads
git clone https://github.com/masahikoji/enrichment-bop2.git
cp -R enrichment-bop2-github-update-2026-10-25/. enrichment-bop2/
cd enrichment-bop2
python3 scripts/verify_reference_results.py
git status
```

既にclone済みの場合は既存の作業ツリーで `git pull` してからコピーします。差分を確認後、

```bash
git add .
git commit -m "Add revision analyses and verify original BOP2 numerical calibration"
git push origin main
```

必要に応じて、`git diff --stat` や `git diff README.md` で差分を確認してください。GitHubへの直接書き込みはこのパッケージ作成段階では行っていません。

## Mac上でのR検証

```bash
Rscript R/validation/01_reproduce_bop2_boundaries.R
Rscript run_revision.R r1-benchmarks
Rscript run_revision.R r1-reverse
Rscript run_revision.R r2-claims
Rscript run_revision.R r2-adjusted
```

原著BOP2の検証結果は `results/revision/bop2_validation/` に出力されます（Git追跡外）。S5の基準結果はgamma上限8です。**gamma上限12は探索的感度解析であり、Mainに記載した主探索仕様とは区別します。** 大規模な複雑なエンドポイントのシミュレーションは時間がかかるため、必要に応じて別途実行してください。

## 公開上の注意

生の患者データおよび原著論文のPDFは含めていません。投稿原稿・Supplement・査読回答書との整合性を最終確認し、ライセンスを適切に設定してください。詳細は `PUBLIC_RELEASE_CHECKLIST.md` を参照してください。
