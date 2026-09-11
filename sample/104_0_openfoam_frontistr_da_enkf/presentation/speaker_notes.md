# 発表者ノート

所属・発表者名・学会名は未指定。表紙を発表前に設定する。

## 01. 少数の温度・変位観測から
温度分布と発熱量を推定できるか

本資料はページ数を制限しない詳細版。学会名・発表者・所属は未指定のため記入欄としている。実測試験ではなく数値的な双子実験であることを最初に明示する。

出典：対象：sample/104_0_openfoam_frontistr_da_enkf

## 02. 第0部　予報と観測を統合する理論



出典：Särkkä & Svensson, Bayesian Filtering and Smoothing (2023)

## 03. まず、何を推定する問題なのか

ここでは加法ノイズの状態空間モデルを導入する。一般の非加法モデルにもベイズフィルタは拡張できる。

出典：Särkkä & Svensson, Bayesian Filtering and Smoothing (2023)

## 04. 記号と単位：発熱量Qと誤差共分散を混同しない

文献で過程ノイズ共分散をQと表す場合が多いため、ここでは混同防止のためΩを使用する。

出典：Särkkä & Svensson, Bayesian Filtering and Smoothing (2023)

## 05. 分散と共分散が表すもの

以下のPは推定誤差についての共分散。EnKFではメンバーのばらつきでこれを近似する。

出典：Särkkä & Svensson, Bayesian Filtering and Smoothing (2023)

## 06. 逐次推定：予報してから、新しい観測で更新する



出典：Särkkä & Svensson, Bayesian Filtering and Smoothing (2023)

## 07. ベイズ更新：どの状態が観測を説明しやすいか



出典：Särkkä & Svensson, Bayesian Filtering and Smoothing (2023)

## 08. 予報分布：前の事後分布をモデルで伝える



出典：Särkkä & Svensson, Bayesian Filtering and Smoothing (2023)

## 09. 第0部A　最適内挿法（OI）



出典：Bouttier & Courtier, Data Assimilation Concepts and Methods, §4–7 (ECMWF)

## 10. OIの前提：背景と観測は誤差を持つ

ガウス分布を仮定すれば後の解は事後平均かつMAP。ガウス性なしでも線形不偏推定の最小分散として同じゲインを導ける。

出典：Bouttier & Courtier, Data Assimilation Concepts and Methods, §4–7 (ECMWF)

## 11. OI導出①：背景と観測からのずれを評価する



出典：Bouttier & Courtier, Data Assimilation Concepts and Methods, §4–7 (ECMWF)

## 12. OI導出②：Jを微分して0と置く



出典：Bouttier & Courtier, Data Assimilation Concepts and Methods, §4–7 (ECMWF)

## 13. OI導出③：背景との差δxで書き直す



出典：Bouttier & Courtier, Data Assimilation Concepts and Methods, §4–7 (ECMWF)

## 14. OI導出④：観測空間の逆行列へ変形する



出典：Bouttier & Courtier, Data Assimilation Concepts and Methods, §4–7 (ECMWF)

## 15. OIの更新式とゲイン



出典：Bouttier & Courtier, Data Assimilation Concepts and Methods, §4–7 (ECMWF)

## 16. OI導出⑤：更新後の誤差共分散



出典：Bouttier & Courtier, Data Assimilation Concepts and Methods, §4–7 (ECMWF)

## 17. 「最適」の意味：線形不偏推定の誤差を最小化する

ガウス性がない場合、これが全ての非線形推定を含むMMSE解とは限らない。

出典：Bouttier & Courtier, Data Assimilation Concepts and Methods, §4–7 (ECMWF)

## 18. OIの数値例：温度予報20 ℃、観測22 ℃



出典：Bouttier & Courtier, Data Assimilation Concepts and Methods, §4–7 (ECMWF)

## 19. OIは何を省略しているか



出典：Bouttier & Courtier, Data Assimilation Concepts and Methods, §4–7 (ECMWF)

## 20. 第0部B　カルマンフィルタ（KF）



出典：Kalman (1960); Särkkä & Svensson (2023)

## 21. KFの線形状態空間モデル



出典：Kalman (1960); Särkkä & Svensson (2023)

## 22. KF導出①：平均の予報



出典：Kalman (1960); Särkkä & Svensson (2023)

## 23. KF導出②：共分散の予報

誤差を推定値−真値で定義する。真値には過程ノイズηが加わるため予報誤差には−ηが入るが、共分散の寄与は＋Ωとなる。

出典：Kalman (1960); Särkkä & Svensson (2023)

## 24. KF導出③：予報を背景としてOIの式を適用



出典：Kalman (1960); Särkkä & Svensson (2023)

## 25. KFの数値例①：1回目の予報と更新



出典：Kalman (1960); Särkkä & Svensson (2023)

## 26. KFの数値例②：2回目は前回の結果から再開



出典：Kalman (1960); Särkkä & Svensson (2023)

## 27. 非線形・高次元になると何が難しいか

メモリは共分散1枚の概算であり、行列積やソルバ全体の使用量ではない。

出典：Kalman (1960); Särkkä & Svensson (2023)

## 28. 第0部C　アンサンブルカルマンフィルタ（EnKF）



出典：Evensen (2003); 104/dacore/enkf.py

## 29. EnKF導出①：各メンバーを予報する



出典：Evensen (2003); 104/dacore/enkf.py

## 30. EnKF導出②：KFのPを標本共分散で近似する



出典：Evensen (2003); 104/dacore/enkf.py

## 31. EnKF導出③：予報観測の偏差を作る



出典：Evensen (2003); 104/dacore/enkf.py

## 32. EnKF導出④：線形観測ならKFの式に一致する構造



出典：Evensen (2003); 104/dacore/enkf.py

## 33. EnKF導出⑤：確率的EnKFのメンバー更新



出典：Evensen (2003); 104/dacore/enkf.py

## 34. なぜ観測摂動を入れるのか：共分散の展開

決定論的平方根EnKFは、摂動を使わず別の変換で目標共分散を実現する。単にεを削除する方法とは違う。

出典：Evensen (2003); 104/dacore/enkf.py

## 35. 理論の列ベクトルを、NumPyの行配列へ直す



出典：Evensen (2003); 104/dacore/enkf.py

## 36. 低ランク・インフレーション・局所化



出典：Evensen (2003); 104/dacore/enkf.py

## 37. 未知のQも状態に入れる：拡大状態



出典：Evensen (2003); 104/dacore/enkf.py

## 38. 第0部D　粒子フィルタ（PF）



出典：Särkkä & Svensson (2023); 104/dacore/pf.py（実装との相違は本文参照）

## 39. PF導出①：確率分布を重み付き粒子で近似



出典：Särkkä & Svensson (2023); 104/dacore/pf.py（実装との相違は本文参照）

## 40. PF導出②：提案分布から粒子を生成する



出典：Särkkä & Svensson (2023); 104/dacore/pf.py（実装との相違は本文参照）

## 41. PF導出③：ブートストラップPFでは遷移密度が消える



出典：Särkkä & Svensson (2023); 104/dacore/pf.py（実装との相違は本文参照）

## 42. PF導出④：ガウス観測誤差から尤度を作る



出典：Särkkä & Svensson (2023); 104/dacore/pf.py（実装との相違は本文参照）

## 43. PFの数値例①：5つの温度候補に尤度を付ける



出典：Särkkä & Svensson (2023); 104/dacore/pf.py（実装との相違は本文参照）

## 44. PFの数値例②：重み付き平均とESS



出典：Särkkä & Svensson (2023); 104/dacore/pf.py（実装との相違は本文参照）

## 45. PF導出⑤：系統リサンプリングを具体的に追う

この1回の再標本化後の平均は21.4℃。再標本化前の重み付き平均と毎回ぴったり同じになるわけではない。

出典：Särkkä & Svensson (2023); 104/dacore/pf.py（実装との相違は本文参照）

## 46. PFの弱点：重みの退化と粒子の重複



出典：Särkkä & Svensson (2023); 104/dacore/pf.py（実装との相違は本文参照）

## 47. 104のPF参考コードと、標準理論の差

本作業は資料の改訂であり、PFソルバコードや過去の解析結果を変更しない。pf_updateは一様事前重みからの一回更新としては理解できる。

出典：Särkkä & Svensson (2023); 104/dacore/pf.py（実装との相違は本文参照）

## 48. 4手法の違いを、同じ軸で比較する



出典：Särkkä & Svensson, Bayesian Filtering and Smoothing (2023)

## 49. 「最適」「非線形」「非ガウス」の意味を整理する



出典：Särkkä & Svensson, Bayesian Filtering and Smoothing (2023)

## 50. 104の説明を読む順番：式→行列→ファイル



出典：Särkkä & Svensson, Bayesian Filtering and Smoothing (2023)

## 51. 理論の参考資料と本資料の位置付け



出典：Särkkä & Svensson, Bayesian Filtering and Smoothing (2023)

## 52. 第1部　OpenFOAM ＋ FrontISTR でデータ同化

章区切り

出典：104ケースの設定・実装を基に作成

## 53. 出発点：内部の温度と、実際の発熱量がわからない

まず対象と困りごとを説明する。センサのない内部温度と、設定値だけでは確定しない発熱量を知りたい。この動機に対してデータ同化を使う。実験計測は今回実施しておらず、観測を想定した数値実験である。

出典：README.md; docs/02_beginner_guide.md

## 54. 解析対象：片側から加熱する中空円筒

形状寸法は102_1のrun_thermal_expansion.pyの定数に基づく。実ソルバの固体セル数は104の保存結果で確認した。

出典：102_1/python/run_thermal_expansion.py; 104/results/openfoam_fem_enkf_summary.yaml

## 55. 104で確かめたいこと：外した予測を修正できるか

何をもって成功とするかを先に示す。温度の観測点だけでなく全セルの誤差を評価する。Qは直接観測せず推定する。本検証は短時間・限定条件での成立確認。

出典：daof/of_fem_twin.py; results/openfoam_fem_enkf_summary.yaml

## 56. なぜ変位も使うのか：温度とは別の手がかりを得たい

研究の狙いと今回確認できた範囲を区別する。変位追加で精度が何倍良くなったとは言えない。今回試したのは温度と変位を合わせて用いる仕組みで、比較実験は次段階である。

出典：README.md; docs/00_temp_displacement_da.md; fem/fem_obs.py

## 57. OpenFOAMで何を計算したか：固体の温度場

chtMultiRegionFoamは固体伝導＋流体NS＋輻射の連成。出力は固体温度場。

出典：sample/102_0（CHT）

## 58. 流体の支配方程式（OpenFOAM）

流体はNS＋浮力＋エネルギー。浮力項ρgが自然対流を生む。

出典：sample/102_0（CHT・流体）

## 59. OpenFOAMは流体も解く：自然対流の流速分布

流体もNSで解いており、自然対流が温度場に効く。

出典：sample/102_0 の流体U(0-600s)

## 60. FrontISTRで何を計算したか：熱膨張の変位

FrontISTRは温度から変形（変位）を計算するソルバ。

出典：sample/102_1

## 61. 連成前進モデルの全体像（どのソルバが何を解くか）

OpenFOAMが熱流体＋輻射、FrontISTRが変形。温度→変位の一方向連成。

出典：104の連成(daof・fem)

## 62. 構造解析の物性と支持条件

拘束指定はそれぞれ自由度1〜3、2〜3、3を拘束する。底面全固定と要約しない。今回は物性と支持条件の同時推定は行っていない。

出典：102_1/config/material_properties_steel.yaml; python/fistr_case.py: write_cnt()

## 63. 観測ベクトル：温度2成分と変位2成分

温度プローブ座標はhot=(0.032,0,0.05025)m、cold=(-0.032,0,0.05025)m。観測ノイズは仮定した値で実機評価値ではない。

出典：openfoam/da_openfoam_config.yaml; fem/fem_obs.py; 102_1/run_one_time()

## 64. どの点を観測に、どの点を確認に使ったか

観測=温度2点＋変位2点。未観測のmid/top/coreで効果を検証する。

出典：104のnode_locations

## 65. 状態ベクトル：温度場とQを一緒に更新

拡大状態によるパラメータ推定。Qを観測しているわけではない。温度・変位への影響を通じて間接推定する。

出典：daof/of_fem_twin.py: n_aug, iQ, Z; of_case.set_solid_state()

## 66. 観測演算子：温度場から予報観測を作る

既存docsでは非線形観測演算子と呼んでいるが、固定補間・一定物性・線形静解析なら温度差から変位への写像は線形またはアフィンに扱える。ソルバを通すこと自体は非線形性の証明ではない。今回は行列を明示しない評価方式として説明する。

出典：daof/of_fem_twin.py: member_obs(); fem/fem_obs.py: displacement_obs()

## 67. 平均・偏差の理論式を確認してから、数値を代入する

理論式を省略せず、次の数値例に対応させる。

出典：dacore/enkf.py: enkf_update(); 理論章のEnKF導出

## 68. 平均の具体例①：5メンバーを行に並べる

説明用に選んだ仮想データ。104の保存結果ではない。温度2列とQに縮小し、インフレーションは1として計算する。数値はensemble_numeric_example.pyで再現できる。実コードのインフレーション1.05を適用する場合は偏差と共分散およびゲインが変わる。

出典：説明用仮想例：presentation/ensemble_numeric_example.py ／ 実装：dacore/enkf.py

## 69. 平均の具体例②：列ごとの合計を5で割る

説明用に選んだ仮想データ。104の保存結果ではない。温度2列とQに縮小し、インフレーションは1として計算する。数値はensemble_numeric_example.pyで再現できる。実コードのインフレーション1.05を適用する場合は偏差と共分散およびゲインが変わる。

出典：説明用仮想例：presentation/ensemble_numeric_example.py ／ 実装：dacore/enkf.py

## 70. 偏差の具体例：各行から同じ平均ベクトルを引く

説明用に選んだ仮想データ。104の保存結果ではない。温度2列とQに縮小し、インフレーションは1として計算する。数値はensemble_numeric_example.pyで再現できる。実コードのインフレーション1.05を適用する場合は偏差と共分散およびゲインが変わる。

出典：説明用仮想例：presentation/ensemble_numeric_example.py ／ 実装：dacore/enkf.py

## 71. 予報観測も同じ計算：Yf → ybar → dY

説明用に選んだ仮想データ。104の保存結果ではない。温度2列とQに縮小し、インフレーションは1として計算する。数値はensemble_numeric_example.pyで再現できる。実コードのインフレーション1.05を適用する場合は偏差と共分散およびゲインが変わる。

出典：説明用仮想例：presentation/ensemble_numeric_example.py ／ 実装：dacore/enkf.py

## 72. 観測偏差dY：温度と変位を列ごとに中心化する

説明用に選んだ仮想データ。104の保存結果ではない。温度2列とQに縮小し、インフレーションは1として計算する。数値はensemble_numeric_example.pyで再現できる。実コードのインフレーション1.05を適用する場合は偏差と共分散およびゲインが変わる。

出典：説明用仮想例：presentation/ensemble_numeric_example.py ／ 実装：dacore/enkf.py

## 73. 偏差を1.05倍して、過度な収縮を緩和する

inflationという名前だけで共分散を1.05倍すると説明しない。実装の倍率と一般の論文における記号定義を区別する。

出典：openfoam/da_openfoam_config.yaml: inflation; dacore/enkf.py

## 74. 共分散の理論式：成分の和を行列積へ書き直す

理論式を省略せず、次の数値例に対応させる。

出典：dacore/enkf.py: enkf_update(); 理論章のEnKF導出

## 75. 共分散①：転置すると、状態ごとの偏差が1行になる

説明用に選んだ仮想データ。104の保存結果ではない。温度2列とQに縮小し、インフレーションは1として計算する。数値はensemble_numeric_example.pyで再現できる。実コードのインフレーション1.05を適用する場合は偏差と共分散およびゲインが変わる。

出典：説明用仮想例：presentation/ensemble_numeric_example.py ／ 実装：dacore/enkf.py

## 76. 共分散②：C_zyの全要素を計算する

説明用に選んだ仮想データ。104の保存結果ではない。温度2列とQに縮小し、インフレーションは1として計算する。数値はensemble_numeric_example.pyで再現できる。実コードのインフレーション1.05を適用する場合は偏差と共分散およびゲインが変わる。

出典：説明用仮想例：presentation/ensemble_numeric_example.py ／ 実装：dacore/enkf.py

## 77. 共分散③：C_yyは、dY同士の内積から作る

説明用に選んだ仮想データ。104の保存結果ではない。温度2列とQに縮小し、インフレーションは1として計算する。数値はensemble_numeric_example.pyで再現できる。実コードのインフレーション1.05を適用する場合は偏差と共分散およびゲインが変わる。

出典：説明用仮想例：presentation/ensemble_numeric_example.py ／ 実装：dacore/enkf.py

## 78. 共分散④：対角と非対角を1要素ずつ手計算する

説明用に選んだ仮想データ。104の保存結果ではない。温度2列とQに縮小し、インフレーションは1として計算する。数値はensemble_numeric_example.pyで再現できる。実コードのインフレーション1.05を適用する場合は偏差と共分散およびゲインが変わる。

出典：説明用仮想例：presentation/ensemble_numeric_example.py ／ 実装：dacore/enkf.py

## 79. 観測誤差共分散Rは、想定した測定誤差

Kとmmの混在は整合した共分散と単位で扱う。変位の単位をµmへ変更するときはYf、y、Rをすべて対応して変換する必要がある。

出典：daof/of_fem_twin.py: R = np.diag(...)

## 80. ゲインと更新の理論式：なぜsolveを使うのか

理論式を省略せず、次の数値例に対応させる。

出典：dacore/enkf.py: enkf_update(); 理論章のEnKF導出

## 81. ゲイン①：観測誤差Rを足してSを作る

説明用に選んだ仮想データ。104の保存結果ではない。温度2列とQに縮小し、インフレーションは1として計算する。数値はensemble_numeric_example.pyで再現できる。実コードのインフレーション1.05を適用する場合は偏差と共分散およびゲインが変わる。

出典：説明用仮想例：presentation/ensemble_numeric_example.py ／ 実装：dacore/enkf.py

## 82. ゲイン②：S @ K.T = C_zy.T を解く

説明用に選んだ仮想データ。104の保存結果ではない。温度2列とQに縮小し、インフレーションは1として計算する。数値はensemble_numeric_example.pyで再現できる。実コードのインフレーション1.05を適用する場合は偏差と共分散およびゲインが変わる。

出典：説明用仮想例：presentation/ensemble_numeric_example.py ／ 実装：dacore/enkf.py

## 83. 更新①：member_00の「観測−予報」を計算する

説明用に選んだ仮想データ。104の保存結果ではない。温度2列とQに縮小し、インフレーションは1として計算する。数値はensemble_numeric_example.pyで再現できる。実コードのインフレーション1.05を適用する場合は偏差と共分散およびゲインが変わる。

出典：説明用仮想例：presentation/ensemble_numeric_example.py ／ 実装：dacore/enkf.py

## 84. 更新②：全メンバーにゲインを掛けて補正する

説明用に選んだ仮想データ。104の保存結果ではない。温度2列とQに縮小し、インフレーションは1として計算する。数値はensemble_numeric_example.pyで再現できる。実コードのインフレーション1.05を適用する場合は偏差と共分散およびゲインが変わる。

出典：説明用仮想例：presentation/ensemble_numeric_example.py ／ 実装：dacore/enkf.py

## 85. 実データで手計算：変位1点がメンバー全員を真値へ戻す

実ログ5メンバーの手計算で本番の同化を再現。ゲイン0.827K/µm。FrontISTR=変位予報係、補正=EnKF。

出典：run_enkf_600.log 実測値／docs/13_displacement_correction_by_hand.md

## 86. 同化後の書き戻しが、次の予報に反映される

流体場は各メンバーが持つものを継続。全状態の整合性を厳密に投影する処理はない。クリッピングは更新分布を変え得るため、将来は適用頻度の記録も有用。

出典：daof/of_fem_twin.py; daof/of_case.py: set_solid_state()

## 87. Q更新の理論式：4成分の残差の重み付き和

理論式を省略せず、次の数値例に対応させる。

出典：dacore/enkf.py: enkf_update(); 理論章のEnKF導出

## 88. Qの更新①：C_zyとKの「Qの行」を取り出す

説明用に選んだ仮想データ。104の保存結果ではない。温度2列とQに縮小し、インフレーションは1として計算する。数値はensemble_numeric_example.pyで再現できる。実コードのインフレーション1.05を適用する場合は偏差と共分散およびゲインが変わる。

出典：説明用仮想例：presentation/ensemble_numeric_example.py ／ 実装：dacore/enkf.py

## 89. Qの更新②：4成分の寄与を足してQを修正する

説明用に選んだ仮想データ。104の保存結果ではない。温度2列とQに縮小し、インフレーションは1として計算する。数値はensemble_numeric_example.pyで再現できる。実コードのインフレーション1.05を適用する場合は偏差と共分散およびゲインが変わる。

出典：説明用仮想例：presentation/ensemble_numeric_example.py ／ 実装：dacore/enkf.py

## 90. フォルダ構成（実ソルバ版）

実ソルバ版=daof+fem+openfoam。EnKF解析はdacore/enkf.pyを共有。

出典：104のフォルダ

## 91. 計算の流れ：温度→変位→観測→補正

温度2点は温度場から直接、変位2点はFrontISTR経由でEnKFへ。

出典：104の連成構成

## 92. 5ケースを順番に前進し、そろってから同化する

メンバー並列の実装ではない。スレッド数とアンサンブル数は別の概念。真値は先に2区間計算して合成観測を用意する。

出典：daof/of_fem_twin.py: for i, m in enumerate(members); daof/of_case.py: run_window()

## 93. 何をしたか②：5通りの予測を、観測で2回補正した

EnKFの数式の前に操作を説明する。5ケースの平均と共分散が補正方向を決める。同化後の各ケースを次の予報へ使う。詳細な計算式は後半の補足資料に残す。

出典：daof/of_fem_twin.py; dacore/enkf.py

## 94. 計算回数とコストの構造

これは同化ドライバの呼び出し回数。準備や可視化用の追加解析は除外。壁時計時間の総計は今回集計していないため、リアルタイム性や何倍高速といった主張を避ける。

出典：daof/of_fem_twin.py; openfoam/run_fem_enkf/member_*/log.run_*

## 95. EnKFへ渡す配列のサイズ

例えばZf[:,100]は5ケースにおけるセル100の温度。1ケースの全セルを平均する操作と明確に区別する。

出典：daof/of_fem_twin.py: Z, Yf; dacore/enkf.py: enkf_update()

## 96. 温度場RMSEは、2回の同化で大きく低下

CSVの温度場RMSEのみで図を再描画した。既存図にある観測セルRMSEの重複系列は掲載しない。曲線の間の時刻で同じ誤差だったことは意味しない。

出典：results/openfoam_fem_enkf_history.csv（発表用に再描画）

## 97. 発熱量Qは、真値15 Wの近傍へ

帯はQ.std()によるメンバー平均±表示用標準偏差で、平均の信頼区間ではない。単一条件での結果であり、任意の熱入力履歴の同定が保証されたわけではない。

出典：results/openfoam_fem_enkf_history.csv; summary.yaml（発表用に再描画）

## 98. 3次元温度分布：同化後と真値を比較

左は初期アンサンブル平均を表示したもので、自由予報を20秒まで走らせた結果ではない。左と中央・右で色範囲が違う。中央と右は共通色範囲。

出典：run/plot_da_temperature_field3d.py; fields/ensmean_t0.npy, ensmean_t20.npy, truth_final.npy

## 99. 予報変位：次の区間ではメンバーの散らばりが縮小

値は保存ログの丸め値。破線は真値ではなくノイズを加えた合成観測。左右は別スケールである。既存変位図の0秒はゼロ代入で実計算ではないため、この図では使用しない。

出典：openfoam/run_fem_enkf.log（予報値の丸め値から再描画）

## 100. 補正した平均温度場から、変位分布も再計算

左は初期平均温度からの変位。中央と右は20秒の同化後平均温度・真値温度からの解析である。変形は見やすく拡大され、左と中央・右の色範囲は異なる。変位分布の定量RMSEは未集計。

出典：run/plot_da_displacement_field.py

## 101. 結果の数値を一表で確認する

有効サンプルサイズESSの列は互換用のNaNであり、このEnKF結果の評価指標として使用しない。統計の定義を揃えた表を掲載する。

出典：results/openfoam_fem_enkf_history.csv

## 102. 第1部まとめ：正確だが計算が重い

実ソルバ版は正確だが重い。全時間・多数試行には軽いROMが要る。

出典：results/openfoam_fem_enkf_summary

## 103. 第2部　ROM にした

章区切り

出典：104ケースの設定・実装を基に作成

## 104. ROM章の目的：軽いモデルで、何を比較するのか

ROM章では計算A〜Dを区別して説明する。温度のみの既定計算と、別ドライバの変位追加比較では観測ノイズの生成順序が異なる。同じ乱数シードでも同じノイズ系列にはならない。

出典：dacore/calibrate.py; run/run_rom_fem.py; compare_rom_disp.py; make_da_gifs.py

## 105. ROM章の図の地図：何を入力し、何と比較するか

ROM章では計算A〜Dを区別して説明する。温度のみの既定計算と、別ドライバの変位追加比較では観測ノイズの生成順序が異なる。同じ乱数シードでも同じノイズ系列にはならない。

出典：計算A〜Dの生成プログラムに対応

## 106. 何の方程式を、どうROMにしたか

同じ熱収支をセル2万→5点に粗くした軽い版。C・K・hは校正で決定。

出典：104のdacore

## 107. ROMの発想（図解）：重いPDEを5点のODEに縮約し校正する

重いPDEを『代表5点＋熱の通り道』のODEに縮約し、本物に校正する。

出典：104のdacore

## 108. ROM：温度を5つの代表ノードで近似する

5ノードが全3次元温度場を任意条件で高精度に再現すると主張しない。校正対象の応答に合う低次元モデルであり、空間分布や条件外への一般化は別途確認が必要。

出典：dacore/cht_rom.py; dacore/calibrate.py; docs/03_rom_derivation.md

## 109. ROMの熱収支と校正

熱容量の合計が物理的な目安に近いことは整合性の一指標だが、9係数の一意性や各ノード容量の正しさを保証しない。midの熱容量は下限近傍で、識別性の評価余地がある。

出典：config/rom_calibrated.yaml; dacore/calibrate.py

## 110. ROM版のアンサンブルは60行×7列

ROM版では各ノード温度も独立に乱数化し、qとhには前進時に微小なランダムウォークを加える。実ソルバ版と初期分布・推定パラメータが異なるので性能の直接比較には使わない。

出典：dacore/twin_fem.py; ensemble.py; cht_rom.py; config/da_config.yaml

## 111. フォルダ構成（ROM版）

ROM版=dacore/だけで完結（純NumPy、実ケース不要）。EnKFは共通。

出典：104のdacore

## 112. 計算Bの条件：温度2点だけで同化し、同化なしと比較

ROM章では計算A〜Dを区別して説明する。温度のみの既定計算と、別ドライバの変位追加比較では観測ノイズの生成順序が異なる。同じ乱数シードでも同じノイズ系列にはならない。

出典：run/run_rom_fem.py → run_rom_fem_twin(...), assim_disp=False

## 113. 計算Bの結果：観測していない代表温度も補正される

図はROM自身を真値に使った双子実験。青が同化、灰が同化なし、破線が真値。OpenFOAMによる600秒の全メンバー計算結果ではない。

出典：run/run_rom_fem.py; dacore/twin_fem.py

## 114. 計算Bの結果：同化なしと温度RMSEを比較

2026-09-10にrun_rom_fem_twinを現設定で再実行し数値を確認。最終core誤差は+0.04039K。実ソルバの温度のみ同化との比較ではない。

出典：dacore/twin_fem.py: run_rom_fem_twin(); config/da_config.yaml

## 115. ROM変位の理論式：温度差の線形写像

理論式を省略せず、次の数値例に対応させる。

出典：dacore/enkf.py: enkf_update(); 理論章のEnKF導出

## 116. ROM変位①：実際のDと、説明用の温度差を用意する

座標はdacore/node_locations.py、Dはconfig/displacement_operator.yamlの実値。入力温度と補間先は説明用で、保存された特定時刻の解析結果ではない。Wは可視化用の温度補間、Dは上面2成分の変位近似であり、全節点変位を与える行列ではない。

出典：実座標・校正D＋説明用温度：presentation/rom_numeric_example.py

## 117. ROM変位②：ヒータ側はDの1行目と温度差の内積

座標はdacore/node_locations.py、Dはconfig/displacement_operator.yamlの実値。入力温度と補間先は説明用で、保存された特定時刻の解析結果ではない。Wは可視化用の温度補間、Dは上面2成分の変位近似であり、全節点変位を与える行列ではない。

出典：実座標・校正D＋説明用温度：presentation/rom_numeric_example.py

## 118. ROM変位③：2行まとめて計算すると、上面2成分が得られる

座標はdacore/node_locations.py、Dはconfig/displacement_operator.yamlの実値。入力温度と補間先は説明用で、保存された特定時刻の解析結果ではない。Wは可視化用の温度補間、Dは上面2成分の変位近似であり、全節点変位を与える行列ではない。

出典：実座標・校正D＋説明用温度：presentation/rom_numeric_example.py

## 119. ROM変位④：同化との関係を数値例で確認する

座標はdacore/node_locations.py、Dはconfig/displacement_operator.yamlの実値。入力温度と補間先は説明用で、保存された特定時刻の解析結果ではない。Wは可視化用の温度補間、Dは上面2成分の変位近似であり、全節点変位を与える行列ではない。

出典：実座標・校正D＋説明用温度：presentation/rom_numeric_example.py

## 120. 計算Cの目的：温度観測に変位を足すとどう変わるか

ROM章では計算A〜Dを区別して説明する。温度のみの既定計算と、別ドライバの変位追加比較では観測ノイズの生成順序が異なる。同じ乱数シードでも同じノイズ系列にはならない。

出典：run/compare_rom_disp.py: main(), run()

## 121. 計算Cの前処理：Dではなく、単位応答からMを作る

ROM章では計算A〜Dを区別して説明する。温度のみの既定計算と、別ドライバの変位追加比較では観測ノイズの生成順序が異なる。同じ乱数シードでも同じノイズ系列にはならない。

出典：run/compute_M_operator.py; config/M_frontistr_operator.npy（m/K）

## 122. 計算Cの読み方：同じRMSEでも、計算Bとは別系列

ROM章では計算A〜Dを区別して説明する。温度のみの既定計算と、別ドライバの変位追加比較では観測ノイズの生成順序が異なる。同じ乱数シードでも同じノイズ系列にはならない。

出典：dacore/twin_fem.py; run/compare_rom_disp.py の乱数生成・呼出順序

## 123. 計算Cの結果：4つの観測条件で温度RMSEを比較

正しい演算子M＋アフィン評価なら、ROMでも変位は効く（特に温度点が少ないほど）。

出典：run/compare_rom_disp.py; docs/11_rom_concept.md（保存図の丸め値）

## 124. IDWの理論式：距離の重みから温度場へ

理論式を省略せず、次の数値例に対応させる。

出典：dacore/enkf.py: enkf_update(); 理論章のEnKF導出

## 125. IDW①：5点の座標と、補間したい位置を決める

座標はdacore/node_locations.py、Dはconfig/displacement_operator.yamlの実値。入力温度と補間先は説明用で、保存された特定時刻の解析結果ではない。Wは可視化用の温度補間、Dは上面2成分の変位近似であり、全節点変位を与える行列ではない。

出典：実座標・校正D＋説明用温度：presentation/rom_numeric_example.py

## 126. IDW②：距離の二乗の逆数を、合計1に正規化する

座標はdacore/node_locations.py、Dはconfig/displacement_operator.yamlの実値。入力温度と補間先は説明用で、保存された特定時刻の解析結果ではない。Wは可視化用の温度補間、Dは上面2成分の変位近似であり、全節点変位を与える行列ではない。

出典：実座標・校正D＋説明用温度：presentation/rom_numeric_example.py

## 127. IDW③：Wの1行と5温度の内積で、1点の温度を求める

座標はdacore/node_locations.py、Dはconfig/displacement_operator.yamlの実値。入力温度と補間先は説明用で、保存された特定時刻の解析結果ではない。Wは可視化用の温度補間、Dは上面2成分の変位近似であり、全節点変位を与える行列ではない。

出典：実座標・校正D＋説明用温度：presentation/rom_numeric_example.py

## 128. IDW④：補間先を並べれば、温度分布は行列積になる

座標はdacore/node_locations.py、Dはconfig/displacement_operator.yamlの実値。入力温度と補間先は説明用で、保存された特定時刻の解析結果ではない。Wは可視化用の温度補間、Dは上面2成分の変位近似であり、全節点変位を与える行列ではない。

出典：実座標・校正D＋説明用温度：presentation/rom_numeric_example.py

## 129. 温度GIFの生成過程①：各時刻の5温度をメッシュへ補間

現行GIF生成設定assim_disp=Trueで再計算したROMの20秒例。既存GIFの作成設定を遡って検証したものではない。 変位GIFでは各フレーム両側のFrontISTR静解析を実行する設計。今回の資料修正でGIFやFrontISTRを再実行したわけではない。

出典：run/make_da_gifs.py: main, temperature_gif, displacement_gif; rom_animation_numeric_example.json

## 130. 温度GIFの生成過程②：20秒の温度で1フレームを追う

現行GIF生成設定assim_disp=Trueで再計算したROMの20秒例。既存GIFの作成設定を遡って検証したものではない。 変位GIFでは各フレーム両側のFrontISTR静解析を実行する設計。今回の資料修正でGIFやFrontISTRを再実行したわけではない。

出典：run/make_da_gifs.py: main, temperature_gif, displacement_gif; rom_animation_numeric_example.json

## 131. ROM：温度分布の時刻歴（データ同化 vs 真値, 0→600s）

でたらめ初期→同化で温度分布が真値に一致（全時間）。

出典：104のmake_da_gifs

## 132. 変位GIFの生成過程①：補間した節点温度をFrontISTRへ渡す

現行GIF生成設定assim_disp=Trueで再計算したROMの20秒例。既存GIFの作成設定を遡って検証したものではない。 変位GIFでは各フレーム両側のFrontISTR静解析を実行する設計。今回の資料修正でGIFやFrontISTRを再実行したわけではない。

出典：run/make_da_gifs.py: main, temperature_gif, displacement_gif; rom_animation_numeric_example.json

## 133. 変位GIFの生成過程②：熱ひずみ・構造解析・描画

現行GIF生成設定assim_disp=Trueで再計算したROMの20秒例。既存GIFの作成設定を遡って検証したものではない。 変位GIFでは各フレーム両側のFrontISTR静解析を実行する設計。今回の資料修正でGIFやFrontISTRを再実行したわけではない。

出典：run/make_da_gifs.py: main, temperature_gif, displacement_gif; rom_animation_numeric_example.json

## 134. ROM：変位分布の時刻歴（データ同化 vs 真値, 0→600s）

温度が直れば、そこから決まる変位分布も真値に一致。

出典：104のmake_da_gifs

## 135. 第2部まとめ：軽くて全時間を可視化できる

ROMは軽く全時間を見られる。分布は粗いが仕組み理解に最適。

出典：results/rom

## 136. 全体のまとめ

章区切り

出典：104ケースの設定・実装を基に作成

## 137. 全体のまとめ

少数観測＋EnKFで温度場とQを推定。実ソルバとROMを使い分ける。

出典：104全体

## 138. 目的への答え：この条件では温度場とQを推定できた

まとめではROMの結果を実ソルバの結果に混ぜない。最も重要な数値と、研究上の次の問いに絞って説明する。

出典：results/openfoam_fem_enkf_summary.yaml; 本資料の検証範囲

## 139. 現在の結果から言えること

変位を追加する実装の成立と、その追加がどれだけ有利かという科学的評価は別の問い。結論を現時点のエビデンスに合わせる。

出典：results/openfoam_fem_enkf_summary.yaml; ROM再計算結果

## 140. 変位観測の追加効果は、統制比較で検証する

これは今後の計画であり実施結果ではない。103の既存結果をそのまま統制比較とみなすには設定と乱数、時間条件をそろえる確認が必要。

出典：今後の検証計画（本発表で未実施）

## 141. 何がわかり、何がまだわからないか

聞き手が成果を過大評価しないよう、研究の到達点を明示する。このあと対象・計算・結果の根拠を示す。ROMの結果は実ソルバ版と切り分け、補足へ移動する。

出典：104の保存結果と実装に基づく検証範囲

## 142. 補足資料

章区切り

出典：104ケースの設定・実装を基に作成

## 143. 何をしたか①：正解を作り、4成分だけを観測した

真値側と推定側を分けて説明する。正解の全温度をフィルタへ直接入力するわけではない。フィルタが受け取るのは4成分の合成観測。実測による検証とは異なる。

出典：daof/of_fem_twin.py: truth loop, obs_values

## 144. 得られた結果：温度分布と発熱量が正解に近づいた

結果の詳細グラフの前に、研究目的への答えを示す。単一シード・5メンバー・20秒の数値結果である。初期値との比較であり、同化なし20秒との統制比較ではない。

出典：results/openfoam_fem_enkf_history.csv; openfoam_fem_enkf_summary.yaml

## 145. 何ケース、どこで計算したか

base_caseから各ケースを複製する。ケースを再生成する実行は既存ケースを削除して作り直すため、本資料作成では重いソルバは再実行していない。

出典：run/run_openfoam_fem_enkf.py: WORKDIR; daof/of_case.py: prepare_member()

## 146. 初期アンサンブルの具体値

初期温度は280〜320K、Qは5〜30Wの一様乱数。設定seed=20260908、初期アンサンブルの乱数生成器はseed+1=20260909。全セル独立なランダム初期場ではない。

出典：daof/of_fem_twin.py: rng.uniform(); docs/06_ensemble_calculation.md

## 147. 熱解析と構造解析の役割

流体構造の双方向連成と混同しない。EnKFが戻すのは温度とQであり、構造変位そのものを熱ソルバへ書き戻さない。

出典：daof/of_case.py; fem/fem_obs.py; 102_1/python/run_thermal_expansion.py

## 148. 結果表示用の標準偏差は、分母がN

NumPy既定のddof=0で計算した標準偏差を使用。初期Qの具体値はインフレーション前の例。信頼区間や標準誤差と混同しない。

出典：daof/of_fem_twin.py: Q.std(), Yf.std(); docs/06_ensemble_calculation.md

## 149. 評価指標：平均温度場と真値のRMSE

体積加重RMSEではなく、コード上はセル値の単純な二乗平均。保存CSVのrmse_obs_Kはrmse_fieldのコピーなので別指標として使わない。

出典：daof/of_fem_twin.py: rmse(); run/run_openfoam_fem_enkf.py: hist_c

## 150. 5メンバーで20,696セルを扱う意味

20,696自由度を4測定から完全同定といった表現を避ける。モデルが結ぶ低次元の変動方向を用いて推定している。今後、空間的な初期温度誤差を増やすと難易度が上がる。

出典：dacore/enkf.py; daof/of_fem_twin.py: initial ensemble

## 151. 実験に進む前に評価すべき誤差

現状で誤差の全てをRの対角に含めたと解釈しない。物性や拘束の誤差は温度推定の偏りとして現れる可能性がある。実験導入時の研究計画として示す。

出典：今後の検証計画（本発表で未実施）

## 152. 再現性：保存される情報と不足する情報

現時点では共分散の計算法はコードで追えるが、保存ケースだけから当時の全予報共分散を厳密に復元できるわけではない。将来はnpzなどにZf/Yf/y/R/Zaを保存する。

出典：docs/06_ensemble_calculation.md; daof/of_fem_twin.py; dacore/enkf.py

## 153. プログラムを追う順番

入口から全体ループへ進み、ソルバの入出力、観測演算子、更新本体の順に読む。ファイル名は104ケースを基準とする。

出典：docs/06_ensemble_calculation.md

## 154. 観測ノイズとEnKFの観測摂動は別

同じRを使っていても乱数の役割が異なる。再現時は乱数の呼び出し順序も結果に影響する。温度のみ比較を作る際は共通温度ノイズが変わらないよう設計する。

出典：daof/of_fem_twin.py: rng_obs, rng; dacore/enkf.py: eps

## 155. どの「平均」を見ているか

質問で平均の意味が曖昧になったときに使用する。予報変位の平均と同化後温度の平均も、集計する時点が異なる。

出典：102_1/run_one_time(); daof/of_fem_twin.py

## 156. 「観測演算子は非線形か」という問い

docs/00と01の非線形という表現をそのまま学術的な主張に使わない。固定線形構造系であればK_u u=f_th(T)なので、温度差に対する線形作用として解釈できる。ヤコビアンを計算せず使える実装の利点は残る。

出典：fem/fem_obs.py; 102_1/python/fistr_case.py: TYPE=STATIC; dacore/enkf.py

## 157. 可視化の由来を混同しない

静止図の左は初期状態であり、同化なしで時間発展した結果ではない点も補足する。発表スライドの各図には由来と時刻条件を記載した。

出典：run/export_paraview.py; plot_da_temperature_field3d.py; plot_da_displacement_field.py

## 158. 文献と結果の一次資料

Evensen文献はEnKFの背景として引用し、この円筒ケースの精度を裏付ける文献とはしない。文献の書誌情報はオンライン一次掲載元で確認した。summary.yamlの正式ファイル名はopenfoam_fem_enkf_summary.yaml。

出典：文献掲載元を2026-09-10確認; 104ケースのローカル一次資料
