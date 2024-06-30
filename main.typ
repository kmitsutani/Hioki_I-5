// MIT No Attribution
// Copyright 2024 Shunsuke Kimura

#import "jdocuments/jnote.typ": main, appendix, thebibliography
#import "@preview/physica:0.9.3"
#show: main.with(
  title: [日置『場の量子論』 I-5], 
  authors: [三ツ谷和也],
)

= Boson と Fermion

通常の散乱によって変化することのない量子数がすべて同じ粒子からなる系を考える．
この系の$N$粒子状態を$Phi^((N))$とすると，粒子$i$と粒子$j$の入れ替えは区別ができず
物理を変えないのでそのような入れ替えの結果量子状態が
$Phi^((N))_((i arrow.l.r j))$になったとするとこれは元の量子状態と位相しか違わず
$
Phi^((N))_((i arrow.l.r j)) = C_(i j) Phi^((N))
$
となる．ただしここに$C_(i j)$は$abs(C_(i j))^2=1$を満たす複素係数である．
$i arrow.l.r j$ の交換をもう一度 $Phi^((N))_((i arrow.l.r j))$ に施すと，
位相も含めて完全に元の状態に戻る．
したがって$C_(i j)^2 = 1 arrow.double.l.r C_(i j) = plus.minus 1$となる．
$C_(i j) = 1$となる粒子種をBoson, $C_(i j) = -1$ となる粒子種を Fermion という．

= 多粒子の波動関数

粒子間の相互作用がない非摂動（あるいは摂動の零次）的状況を考えると，
$Phi^((N))(bold(r_1),dots,bold(r_N))$ は
$1$粒子の正規完全直交系$phi_alpha^((1))$の積の線形結合で
$
Phi^((N))(bold(r_1),dots,bold(r_2)) = 
  sum_({s_i}) C({s_i}) phi_(s_i)^((1))(bold(r_1))dots.c phi_(s_N)^((1))(bold(r_N))
$
とかける（独立粒子近似）．
ここで$s_i$は$i$番目の粒子を$alpha=1,2,dots$のどれに割り当てるかの組み合わせであり，
$c({s_i})$はそれぞれの割り当てに対する結合係数である．

ここで状態$phi_alpha$に粒子が$n_alpha$個ある状態を
$
physica.ket(n_1\,n_2\, dots\, n_alpha\, dots)
$
と書くものとする．

以下では簡単のためにすべてがボソンの場合を考える
#footnote[本テキストでは Fermion の取り扱いについて一切書かれていないが
それに関しては例えば砂川 @SunakawaQM の第 $6$ 章を見よ．]
．

= 多粒子系における行列要素に関する物理的考察

以下では$1$粒子状態の正規直行完全系を${physica.ket(alpha)}_(alpha=1,2,dots)$とする．
すなわち
$
physica.braket(alpha, beta) = delta_(alpha,beta), 
sum_alpha physica.ketbra(alpha, alpha) = hat(cal(1))
$

本節では簡単のために二準位系${physica.ket(alpha)} = {physica.ket(1), physica.ket(2)}$
を考えオブザーバブルの行列要素を考える．
ここで独立粒子近似の範囲で扱える粒子間相互作用を含まないオブザーバブル$hat(F)$の行列要素を考える．
粒子間相互作用を含まないという条件から
$hat(F)$は$1$粒子状態に作用するオブザーバブル$hat(F)^((1))$を用いて
$
hat(hat(F)) = &hat(F)^((1)) times.circle hat(1)_2 times.circle dots.c times.circle hat(1)_N \
  &+ hat(1)_1 times.circle hat(F)^((1)) times.circle dots.c times.circle hat(1)_N \
  &dots \
  &+ hat(1)_1 times.circle dots.c times.circle hat(1)_(N-1) times.circle hat(F)^((1)) \
$ <decomposition_of_F>
と展開される．

また状態は
$
physica.ket(n_1\, n_2) = c_n (n_1, n_2) sum_P physica.ket(P(1)) physica.ket(P(2))
                                                  dots physica.ket(P(N))
$
ここで$P(i)$は$i$番目の粒子をどの状態に割り当てるかを表す関数である
（すなわち今の場合$forall i(P(i) in {1, 2})$）．
$c_n(n_1, n_2)$は規格化因子であるがこれは
$
physica.braket(n_1\, n_2, n_1\, n_2) &= abs(c_n (n_1, n_2))^2 sum_(P,P^prime)
physica.braket(P(1),P^prime(1))times dots times physica.braket(P(N),P^prime(N)) \
&= abs(c_n (n_1, n_2))^2 sum_P = abs(c_n (n_1, n_2))^2 binom(N, n_1)
$
となるので，オーバーオールの位相を気にせず$c_n (n_1, n_2) in bb(R)$とすれば
$c_n (n_1, n_2) = sqrt((n_1 ! n_2 ! )\/ N !)$となる．

式 @decomposition_of_F から $hat(F)physica.ket(n_1\, n_2)$ には
$hat(F)^((1))physica.ket(P(i))$が繰り返し登場するのがわかるので$hat(F) physica.ket(P(i))$
を以下に書き下すと
$
hat(F)^((1))physica.ket(P(i)) &= hat(bb(1))hat(F)^((1))physica.ket(P(i)) \
  &= sum_alpha physica.ketbra(alpha, alpha)hat(F)^((1))physica.ket(P(i)) \
  &= physica.bra(1) hat(F)^((1))physica.ket(P(i))physica.ket(1)
  + physica.bra(2) hat(F)^((2))physica.ket(P(i))physica.ket(2)
$
となる．この式を物理的に解釈すると$hat(F)^((1))$の作用によって
状態$physica.ket(P(i))$から$physica.ket(1)$および$physica.ket(2)$への遷移が起きていることになる．
したがって遷移後の状態の占有数を$(n_1^prime, n_2^prime)$とすると
$(n_1^prime, n_2^prime) in {(n_1 - 1, n_2 + 1), (n_1, n_2), (n_1 + 1 , n_2 - 1)}$となり，
この占有数変化に対応する行列要素だけが値を持つ．

= 多粒子系における行列要素の具体的計算
前節での考察から非零になる行列要素が絞り込まれたのでその一つ
$physica.bra(n_1 - 1\, n_2 + 1) hat(F) physica.ket(n_1\, n_2)$を計算する．
$
physica.bra(n_1 - 1\, n_2 + 1)hat(F)physica.ket(n_1\, n_2) 
= c_n (n_1, n_2) physica.bra(n_1 - 1\, n_2 + 1)sum_P &
  {
     [hat(F)^((1))physica.ket(P(1))]physica.ket(P(2))dots physica.ket(P(N)) \
    &+ physica.ket(P(1))[hat(F)^((1)) physica.ket(P(2))] dots physica.ket(P(N)) \
    &dots \
    &+ physica.ket(P(1)) dots physica.ket(P(N-1))[hat(F)^((1)) physica.ket(P(N))]
  }
$
占有数に関する考察から$hat(F)^((1))physica.ket(P(i))$は実質的に
$delta_(P(i),1) f_(2, P(i)) physica.ket(2)$としてふるまう．
ただしここで$hat(F)^((1))$の行列要素を
$physica.bra(i)F^((1))physica.ket(j) = f_(i,j)$と表記した．
したがってうえの行列要素は
$
physica.bra(n_1 - 1\, n_2 + 1)hat(F)physica.ket(n_1\, n_2)
= c_n (n_1, n_2) f_(2, 1) physica.bra(n_1-1\,n_2+1){
  &sum_(P(1)=1)physica.ket(2)physica.ket(P(2))dots.c physica.ket(P(N))\
 &+sum_(P(2)=1)physica.ket(P(1))physica.ket(2)dots.c physica.ket(P(N))\
 &dots\
 &+sum_(P(N)=1)physica.ket(1)dots physica.ket(P(N-1)) physica.ket(2)
}
$
と書ける．
ここで $physica.bra(n_1-1\,n_2+1)$ も$1$粒子状態であらわに書くと
$
& physica.bra(n_1-1\,n_2+1)hat(F)physica.ket(n_1\,n_2)\
& =c_n (n_1 - 1, n_2 + 1) c_n (n_1, n_2)f_(1,2)sum_(i=1)^N sum_(P^prime, P\ P(i) = 1)
  physica.braket(P^prime (1), P(1)) dots physica.braket(P^prime (i),2)
  dots physica.braket(P^prime (N), P(N))
$
この表式の和の部分は$i$に関して対象になっており和の対象はどの$i$ でもおなじである．
そこで$i$を固定して考えると$P,P^prime$に関する和は
$P^prime (i) = 2, P(i) = 1$
という条件下で $forall j eq.not i (P^prime (j) = P(j) )$ を満たす$P,P^prime$ 
の組の数に等しい．
これは$N-1$個の粒子を準位$physica.ket(1)$に$n_1-1$個，
準位$physica.ket(2)$に$n_2$個それぞれ割り当てる場合の数に等しい．
したがって
$
&physica.bra(n_1 - 1\,n_2 + 1)hat(F)physica.ket(n_1\,n_2) \
&= [((n_1-1)!(n_2+1)!)/N!]^(1\/2)[((n_1 ! n_2 !)/N!)]^(1\/2) N binom(N-1, n_1 - 1) f_(2,1)\
&= [((n_1-1)!(n_2+1)!)/N!]^(1\/2)[((n_1 ! n_2 !)/N!)]^(1\/2) [N! / ((n_1-1)! n_2 !)] f_(2,1)\ 
&= sqrt(n_1(n_2+1)) f_(2,1)
$
と計算できる（$n_1 !\/(n_1-1)! = n_1$などに注意）．

同様の計算により $physica.bra(n_1+1\,n_2-1)hat(F)physica.ket(n_1\,n_2)=sqrt(n_2(n_1+1))f_(1,2)$

== $mu$準位系への拡張
$1$粒子状態として${physica.ket(1),dots, physica.ket(mu)}$を考える．
この場合の占有数表示は$physica.ket(n_1\,dots\,n_mu)$というものになる．
いま二準位系の場合と同様に$1$粒子演算子を多粒子状態に作用するように細工した
$hat(F)$の行列要素
$physica.bra(n_1+d_1\,n_2+d_2\,dots\,n_mu +d_mu)
hat(F)
physica.ket(n_1\,n_2\,dots\,n_mu)$を考える．
ここで$d_i$は$sum d_i = 0$を満たす整数の組である．
行列要素のどの部分が残るかについて考察するために
二準位系の場合と同様に$hat(F) physica.ket(n_1\,dots\,n_mu)$について考察することから始める．
そのために占有数表示された状態を$1$粒子状態で展開すると
$
hat(F) physica.ket(n_1\,n_2\,dots \,n_mu)=c_N(n_1, n_2, dots, n_mu)
  hat(F) sum_(P)physica.ket(P(1)) physica.ket(P(2)) dots physica.ket(P(N))
$
となる．ここで$c_N(n_1, n-2, dots, n_mu) = sqrt((product n_alpha !)\/N!)$ は規格化因子．
この形から二準位系の時と同様に$hat(F)^((1)) physica.ket(P(i))$ が出てくるのがわかるのでこの因子
をあらわに書くと
$
hat(F)^((1))physica.ket(P(i)) 
  = sum_(alpha=1)^mu physica.ketbra(alpha)hat(F)^((1))physica.ket(P(i))
  = sum_(alpha=1)^mu f_(alpha, P(i)) physica.ket(alpha)
$
となり，$physica.ket((P(i)))arrow.r physica.ket(alpha) "for" alpha = 1dots mu$
という遷移に対応する行列要素のみが残ることがわかる．
したがって二準位系の場合と同様に行列要素が非零になるのは対角線とその隣のみである．

有限の値を持つ行列要素
$physica.bra(n_alpha - 1\, n_beta + 1)hat(F)physica.ket(n_alpha\, n_beta)$
も二準位系の場合とほぼ同様に以下のように計算できる．
$
physica.bra(n_alpha-1\,n_beta+1)sum_P &
{
[sum_(sigma=1)^mu f_(sigma, P(1))physica.ket(sigma)] physica.ket(2) dots physica.ket(P(N))\ &
+ physica.ket(P(1)) [sum_(sigma=1)^mu f_(sigma, P(2))physica.ket(sigma)] dots physica.ket(P(N))\ &
dots \ &
+ physica.ket(P(1)) dots physica.ket(P(N-1)) [sum_(sigma=1)^mu f_(sigma, P(N))physica.ket(sigma)] 
}
$
今回の場合$hat(F)$の作用で引き起こされる遷移が$physica.ket(alpha) arrow.r physica.ket(beta)$
であると解釈できる要素のみが残るため$[~]$の部分は
$f_(beta, alpha) physica.ket(beta)$に置き換えられる．

$
physica.bra(n_alpha-1\,n_beta+1)sum_P &
{
 f_(beta, alpha) physica.ket(beta) physica.ket(2) dots physica.ket(P(N))\ &
+ physica.ket(P(1)) f_(beta, alpha) physica.ket(beta) dots physica.ket(P(N))\ &
dots \ &
+ physica.ket(P(1)) dots physica.ket(P(N-1)) f_(beta, alpha) physica.ket(beta) 
}
$
$physica.bra(n_alpha -1\, n_beta + 1)$も$1$粒子状態であらわに書き表すと
$
&physica.bra(n_alpha -1\,n_beta+1)hat(F)physica.ket(n_alpha\,n_beta) \
&=c_N(n_alpha - 1\,n_beta + 1)c_N(n_alpha, n_beta) f_(beta, alpha)
  sum_i sum_(P, P^prime\ P(i)=alpha\ P^prime (i) =beta) 
  physica.braket(P^prime (1), P(1)) dots.c physica.braket(P^prime (N), P(N))
$
和に関して二準位系の場合と同じロジックで処理すると
$
physica.bra(n_alpha -1 \, n_beta + 1)hat(F)physica.ket(n_alpha\, n_beta) 
= f_(beta, alpha) sqrt((n_beta+1)/n_alpha) [(product_sigma n_sigma !) / N!] 
  [N! / (product_sigma n_sigma !)] n_alpha = sqrt(n_alpha (n_beta+1)) f_(beta, alpha)
$
となる．

== 期待値の計算

演算子$hat(F)^((1))_i := hat(bb(1))_1
  times.circle dots hat(bb(1))_(i-1) times.circle hat(F)^((1))
  times.circle hat(bb(1))_(i+1) times.circle dots hat(bb(1))_N$ 
#let expval=physica.expectationvalue
#let ket=physica.ket
#let bra=physica.bra
#let braket=physica.braket
#let ketbra=physica.ketbra
を導入すると$hat(F)$の期待値は
$
expval(F) = c_N^2 sum_(P, P^prime) sum_(i=1)^N 
  bra(P^prime (1))bra(P^prime (2)) dots.c bra(P^prime (N))
  hat(F)^((1))_i
  ket(P^prime (1))ket(P^prime (2)) dots.c ket(P^prime (N))

$
と書ける．
占有数状態に関する対角和なので占有数分布を変えない項が残る．
したがって
$
[hat(F)^((1))physica.ket(P(i))] 
  arrow f_(P(i), P(i)) ket(P(i))
$
とすればよい．したがって
$
expval(F) 
  &= c_N^2 sum_i^N sum_(P, P^prime) 
    f_(P(i), P(i)) product_i delta_(P^prime (i), P(i)) \
  &= c_N^2 sum_i^N sum_(sigma=1)^mu sum_(P, P^prime\ P(i) = sigma) f_(sigma, sigma) product delta_(P^prime (i), P(i))
$ <eq21>
となる．ここで$i$および$P,P^prime$に関する和を取り出した部分
$
sum_(i=1)^N sum_(P,P^prime\ P(i) = sigma) 
  f_(sigma, sigma) product delta_(P^prime (i), P(i)) equiv S_sigma
$
を見積もる．ある$i$に対して$sum_(P,P^prime)$ 以下は 
「$P(i)$ を $sigma$ に固定する条件の元で $P(dot.c)$ と $P^prime (dot.c)$ が
写像として一致する場合の数に$f_(sigma, sigma)$を乗じたもの」であるので
$f_(sigma, sigma) binom(N-1, dots\, n_sigma - 1\, dots)$ となるしたがって 
$
S_sigma 
  &= sum_(i=1)^N f_(sigma, sigma) (N-1)!/(product_alpha n_alpha !) n_sigma \
  &= f_(sigma, sigma) N! / (product_alpha n_alpha !) n_sigma

$
ここで式 <eq21> に $S_sigma$ を代入すると
$
expval(F) &= (product_alpha n_alpha ! ) / N! N! / (product_alpha n_alpha ! )
  sum_sigma f_(sigma, sigma) n_sigma
  &= sum_sigma f_(sigma, sigma) n_sigma
$
となり$F$の期待値が求まる．


= Fock空間の導入
#let nket(n, diff: "")=$ket(n_1\,n_2\,dots\,n_#n #diff\,dots)$
#let nbra(n, diff: "")=$bra(n_1\,n_2\,dots\,n_#n #diff\,dots)$
#let nbraket(n)=$braket(n^prime_1\,n^prime_2\,dots\,n^prime_#n\,dots, n_1\,n_2\,dots\,n_#n\,dots)$

$1$粒子ヒルベルト空間の正規直行基底を${ket(alpha)}_(alpha=1,2,dots)$とする．
そしてこの$1$粒子状態$ket(sigma)$を占有する粒子の数が$n_sigma$であるような
正規直交化された$N$粒子状態を
$
nket(sigma)
$
とかく．正規直交性は
$
nbraket(sigma) = product_sigma delta_(n_sigma^prime, n_sigma)
$
とあらわされる．この正規直交性より少し考えると完全性も従う．
これらのベクトルで張られるヒルベルト空間をFock空間といい上記の占有数表示されたベクトルによる
Fock空間の元の表示をFock表示という．

== 生成消滅演算子

この表示の下で$1$粒子状態$sigma$の消滅演算子を調和振動子の場合と同様に

$
hat(a)_sigma nket(sigma) = sqrt(n_sigma) nket(sigma, diff: -1)
$
と作用する演算子であると定義する．左から$nbra(sigma, diff: -1)$ を作用させると
$
nbra(sigma, diff: -1)hat(a)_sigma nket(sigma) = sqrt(n_sigma)
$
となる．これのエルミート共役をとると$sqrt(n_sigma)^dagger = sqrt(n_sigma)$なので
$
nbra(sigma) hat(a)^dagger_sigma nket(sigma, diff: -1) = sqrt(n_sigma)
$
となる．これは任意の$n_sigma$について成立するので$n_sigma$の値を一つずらしても成立する．
したがって
$
hat(a)^dagger_sigma nket(sigma) = sqrt(n_sigma+1) nket(sigma, diff: +1)
$

== 生成消滅演算子の正準交換関係

#let nket2(s1, s2, d1:"", d2: "")=$ket(n_1\, n_2\, dots\, n_#s1 #d1\,dots n_#s2 #d2\,dots)$

以下では前節で導入した生成消滅演算子の正準交換関係について計算していく．
まず $[hat(a)_sigma, hat(a)_rho]$ について考える．$sigma=rho$ の時は自明に零なので
$sigma eq.not rho$ の場合を検討する．
$
[hat(a)_sigma, hat(a)_rho] 
  &= (hat(a)_sigma hat(a)_rho - hat(a)_rho hat(a)_sigma) nket2(sigma, rho)\
  &= hat(a)_sigma sqrt(n_rho) nket2(sigma, rho, d2: -1) 
      - hat(a)_rho sqrt(n_sigma) nket2(sigma, rho, d1: -1)\
  &= (sqrt(n_sigma)sqrt(n_rho) - sqrt(n_rho)sqrt(n_sigma)) nket2(sigma, rho, d1: -1, d2: -1)\
  &=0
$
生成演算子どうしの交換関係も同様の計算により$0$となる．
最後に$[hat(a)_sigma, hat(a)^dagger_rho]$を計算する．
$sigma eq.not rho$の場合と$sigma eq rho$の場合に分けて計算する．
まず $sigma eq.not rho$ の場合は以下の通り
$
[hat(a)_sigma, hat(a)^dagger_rho] nket2(sigma, rho)
  &= (hat(a)_sigma hat(a)^dagger_rho - hat(a)^dagger_rho hat(a)_sigma) nket2(sigma, rho)\
  &= hat(a)_sigma sqrt(n_rho + 1) nket2(sigma, rho, d2: +1) 
      - hat(a)_rho^dagger sqrt(n_sigma) nket2(sigma, rho, d1: -1) \
  &= sqrt(n_sigma (n_rho +1)) [nket2(sigma, rho, d1:-1, d2:+1) 
        - nket2(sigma, rho, d1:-1, d2:+1)] \
  &=0
$
次に $sigma = rho$の場合を計算する．
$
[hat(a)_sigma, hat(a)^dagger_sigma] nket(sigma) 
  &= (hat(a)_sigma hat(a)_sigma^dagger -hat(a)^dagger_sigma hat(a)_sigma) nket(sigma)\
  &= sqrt(n_sigma + 1) hat(a)_sigma nket(sigma, diff: +1) 
    - sqrt((n_sigma)) hat(a)^dagger_sigma nket(sigma, diff: -1) \
  &= (sqrt(n_sigma+1)sqrt(n_sigma+1) - sqrt(n_sigma)sqrt(n_sigma))nket(sigma) \
  &= nket(sigma)
$
となる．これらをまとめると
$
[hat(a)_sigma,hat(a)_rho] = [hat(a)^dagger_sigma, hat(a)^dagger_rho] 
  = 0 times hat(bb(1))_(cal(H)_"Fock")\
[hat(a)_sigma,hat(a)^dagger_rho] = delta_(sigma, rho) times hat(bb(1))_(cal(H)_"Fock")\
$ <CaCR>
となる．ここで$hat(bb(1))_(cal(H)_"Fock")$は今考えているFock空間における恒等変換である．

以下では誤解を招く可能性がない限り $hat(bb(1))_(cal(H)_"Fock")$ を明示的には書かない．

= 演算子の「第二量子化」
#let nbra2(s1, s2, d1:"", d2: "")=$bra(n_1\, n_2\, dots\, n_#s1 #d1\,dots n_#s2 #d2\,dots)$
少し前に導いた演算子$hat(F)$の行列要素の表式
$
nbra2(n_sigma,n_rho,d1:-1,d2:+1) hat(F) nket2(n_sigma, n_rho) 
  = sqrt(n_sigma (n_rho + 1)) bra(rho) hat(F)^(1) ket(sigma)
$
生成演算子を$nket(mu)$に作用させたときに出る係数と比較すると
$
hat(F) = sum_(sigma,rho) hat(a)^dagger_rho bra(rho) F^((1)) ket(sigma) hat(a)_sigma
$
と置けばこれがもともとのの$hat(F)$の行列要素を再現することは容易に確認できる．
初め$hat(F)$を
$
hat(F) = &hat(F)^((1)) times.circle hat(1)_2 times.circle dots.c times.circle hat(1)_N \
  &+ hat(1)_1 times.circle hat(F)^((1)) times.circle dots.c times.circle hat(1)_N \
  &dots \
  &+ hat(1)_1 times.circle dots.c times.circle hat(1)_(N-1) times.circle hat(F)^((1)) \
$
と定義したがこのような煩雑な記述をせずとも生成消滅演算子を用いれば各粒子に独立に
$1$粒子演算子を作用させる $cal(H)_"Fock"$ 上の演算子を簡潔に求めることができる．

一般に$1$粒子演算子 $hat(o):cal(H)_1 arrow cal(H)_1$ に対して
$
hat(cal(B))(hat(o)) equiv sum_(alpha, beta)
hat(a)^dagger_alpha
bra(alpha)hat(o)ket(beta)
hat(a)_beta
$
で定義される操作 $hat(cal(B))(hat(o))$ を $hat(o)$ の第二量子化という．

演算子の第二量子化でもちいる$1$粒子の正規完全直交系は正規完全直交系であればなんでもよい
#footnote[
これについては田崎 @TasakiFock に証明がみられる．
生成消滅演算子を一般の状態を系に生成したり消滅させたりするというものと定義するところが
証明の肝に見えるので普通の教科書では証明を見つけるのは難しいかもしれない．
]．
そこで基底として$1$粒子ハミルトニアン$hat(H)^((1))$の固有状態を選び
$1$粒子ハミルトニアンの第二量子化を行うと，
$
hat(cal(B))(hat(H)^((1))) &= sum_(alpha, beta) 
  hat(a)_alpha^dagger bra(alpha) hat(H)^((1)) ket(alpha) hat(a)_alpha \
  &=sum_(alpha, beta) epsilon_beta delta_(alpha, beta) hat(a)_alpha^dagger hat(a)_beta \
  &= sum_alpha epsilon_alpha hat(n)_alpha
$
となる（ここで $epsilon_alpha$ は準位$ket(alpha)$ のエネルギー固有値）．
これはI-$4$節で導いたスカラー場のハミルトニアン(I-39)と同様の形となる．

= 場の「第二量子化」

古典論から量子論の基礎方程式を導く際に実験結果を見ながらの論理的跳躍が必要であったように，
粒子系の量子論から相対論的場の量子論に演繹的にたどり着けなくても何の問題もない．
筆者は第二量子化は場の理論黎明期の物理学者を刺激した足がかりに過ぎず
場の量子論の現代的な導入には不要であるという立場である．
しかしながら一応テキストをなぞる．

#let Schrodinger=[Schro\u{308}dinger]
== "#Schrodinger 方程式"の導出
$1$粒子波動関数$Phi (bold(x), t)$ を $t=0$で$1$粒子ハミルトニアン$H^((1))$
の固有状態 $phi_alpha (bold(x))$ で $Phi(bold(x), t=0) = sum_alpha a_alpha phi_alpha (bold(x))$
と展開する．ここで規格化条件より $sum_alpha abs(a_alpha)^2 = 1$ となる．
各 $phi_alpha$ は $hat(H)^((1))$ の下で $exp(-i epsilon_alpha t)$ で時間発展する．
これらを合わせると
$
Phi(bold(x), t) = sum_alpha a_alpha phi_alpha (bold(x)) e^(-i epsilon_alpha t)
$
となる（ここまでは波動関数として扱っている．そうでないと時間発展が合わない．）．
ここで $a_alpha arrow hat(a)_alpha$ の読みかえを行って
波動関数 $Phi$ を場の演算子 $hat(Phi)$ に置き換える．すなわち
$
hat(Phi)(bold(x), t) = sum_alpha hat(a)_alpha phi_alpha (bold(x)) e^(-i epsilon_alpha t)
$
とする．そうすると規格化条件は$sum_alpha hat(n)_alpha = hat(bb(1))$ と
置き換えるのが妥当であろう．すると
$
nket(sigma) in {
  ket(1\, 0\, dots\, 0\, dots),
  ket(0\, 1\, dots\, 0\, dots), dots }

$
すなわち状態は$1$粒子ヒルベルト空間$cal(H)_1$ に住んでいることになる．
$
i partial / (partial t) hat(Phi) (bold(x), t) 
    = sum_alpha hat(a)_alpha phi_alpha (bold(x)) epsilon_alpha e^(-i epsilon_alpha t)
$
$
hat(H)^((1)) hat(Phi)(bold(x), t) 
  = sum_alpha hat(a)_alpha epsilon_alpha phi_alpha(bold(x)) e^(-i epsilon_alpha t)
$
なので
$
  i partial / (partial t) hat(Phi)(bold(x), t) = hat(H)^((1)) hat(Phi)(bold(x), t)
$ <SecondaryQuantizedSchroedingerEq>
となる．ただし$1$粒子ハミルトニアンを作用させた時に
$H^((1)) phi_alpha(bold(x)) = epsilon_alpha phi_alpha(bold(x))$ および
$[hat(H)^((1)), hat(a)_alpha] = 0$ を用いた（後者が成立するかはよくわからない）．

場の演算子 $hat(Phi)$ とそのエルミート共役 $hat(Phi)^dagger$ 
は以下の同時刻交換関係を満たす．
$
[hat(Phi)(bold(x), t), hat(Phi)^dagger (bold(x)^prime, t)]
  &= sum_(alpha, beta) [hat(a)_alpha, hat(a)_beta^dagger]
    phi_alpha (bold(x)) phi_beta^* (bold(x^prime))
    exp(i (epsilon_beta - epsilon_alpha) t) \
  &= sum_(alpha, beta) delta_(alpha, beta) 
    phi_alpha (bold(x)) phi_beta^* (bold(x^prime))
    exp(i (epsilon_beta - epsilon_alpha) t) \
  &= sum_alpha phi_alpha (bold(x)) phi_alpha^*(bold(x^prime)) \
  &= delta(x - x^prime)
$
ここで$phi_alpha (bold(x))$の完全性を用いた#footnote[
ブラケット記法を経由するとこれは
$sum_alpha phi_alpha (bold(x)) phi_alpha^* (bold(x^prime)) 
= sum_alpha braket(bold(x)^prime, phi_alpha) braket(phi_alpha, bold(x))
=braket(bold(x)^prime, bold(x)) = delta(bold(x) - bold(x)^prime)$
と書ける．
ブラケットを用いない場合完全性条件
$forall xi in cal(H)_1 (sum_alpha phi_alpha (phi_alpha, xi) = xi)$
に位置表示の内積 $(psi, xi) = integral psi(bold(x))^* xi(bold(x)) d^3x$
を代入して$ forall xi (integral d^3 bold(x)^prime
[sum_alpha phi_alpha (bold(x))phi_alpha^* (bold(x)^prime)]
 xi(bold(x)^prime) ) = xi(bold(x))$ から
 $sum_alpha phi_alpha (bold(x)) phi_alpha^* (bold(x)^prime) = delta(bold(x) - bold(x)^prime)$

]
．
ここで一度古典場に戻ると式 @SecondaryQuantizedSchroedingerEq 
と同じ形の "#Schrodinger 方程式" は
ラグランジアン密度
$
cal(L)(bold(x), t) 
  = Phi^*(bold(x)) [i partial / (partial t) + 1/(2m) laplace - V(bold(x))] Phi(bold(x),t)
$
から導かれる．実際
$
(delta cal(L))/(delta dot(Phi)) &= i Phi^*(bold(x), t)\
(delta cal(L))/(delta nabla Phi) &= -(nabla Phi^*(bold(x), t))/(2m)\
(delta cal(L))/(delta Phi) &= -V(x) Phi^*(bold(x), t)
$
二行目を導くにあたってラグランジアン密度が積分の中に入っていることを利用して
部分積分により
$Phi^* laplace Phi arrow - nabla Phi^* dot.c nabla Phi$
とした#footnote[
何らかのベクトル解析の公式を用いたりしたのではなくラプラシアンの各成分を分けて書いたのち
各成分に対して部分積分を適用してまとめなおしたものである．
]．これらより E-L 方程式は
$
-V(bold(x)) Phi^* - i partial_t Phi^* + (nabla dot.c nabla Phi^*) / (2m)
arrow.double.r.l i partial_t Phi(bold(x),t)= [- laplace / (2m) + V(bold(x))] Phi(bold(x), t)
$
となり#Schrodinger 方程式に一致する．
このラグランジアン密度に対して$Phi$に共役な運動量密度$Pi$を導入すると
$ Pi(bold(x), t) = (delta cal(L))/(delta dot(Phi)) &= i Phi^*(bold(x), t) $
となる．
もしこの共役運動量密度を上で導入した場の演算子$hat(Phi)$と同様に演算子化（量子化）すると
$hat(Pi)$ と $hat(Phi)$ の同時刻交換関係は
$
[hat(Phi)(bold(x),t), hat(Pi)(bold(x)^prime, t)]
=  [hat(Phi)(bold(x), t), i Phi^dagger (bold(x), t)] = i delta(bold(x)-bold(x)^prime)
$
となりこれはI-4で登場したスカラー場の同時刻交換関係と一致する．

= 生成消滅演算子の交換関係と場の演算子の同時刻交換関係の同値性について
$hat(a)^((1))_alpha = hat(a)_alpha, hat(a)^((2))_beta = hat(a)^dagger_beta$ の表記法および，
完全反対称テンソル $epsilon.alt^(i,j)$ (
$epsilon.alt^(1,2) = - epsilon.alt^(2,1) = 1$
それ以外の要素は $0$)
を導入すれば式 @CaCR の形に（複数の等式で）書かれた生成消滅演算子の交換関係は
$
[hat(a)_alpha^((i)), hat(a)_beta^((j))] = epsilon.alt^(i,j) delta_(alpha, beta)
$
と一つの式で書ける．この式と同様に
$Phi^((1))(bold(x),t) = Phi(bold(x),t)$,
$Phi^((2))(bold(x), t) = Pi(bold(x), t) = i Phi^*(bold(x), t)$
と定義すれば先ほど導出した場の演算子の同時刻交換関係は
#let _x=$bold(x),t$
#let _y=$bold(x^prime),t$
$
[hat(Phi)^((i)) (#_x), hat(Phi)^((j)) (#_y)] 
  = i epsilon.alt^(i,j) delta(bold(x) - bold(x)^prime)
$
と書ける．
これらの表式を複数の等式であらわされる正交換関係
や同時刻交換関係を簡潔に記述するために用いる．
少なくとも今回は実際の計算に用いても特に一般化のメリットは得られないので計算には用いない．
上記の同時刻交換関係には明示的に導出していない成分も含まれるが$Phi$どうしや$Pi$どうしの交換関係は
生成演算子同士や消滅演算子同士の交換関係に帰着されるので$0$になるのはほぼ自明であろう．
したがって前節の議論で
$[hat(a)_alpha^((i)), hat(a)_beta^((j))] = epsilon.alt^(i,j) delta_(alpha, beta)$
から $ [hat(Phi)^((i)) (#_x), hat(Phi)^((j)) (#_y)] 
= i epsilon.alt^(i,j) delta(bold(x) - bold(x)^prime)$
が導出されることは示したといえる．

ここではその逆を示し両者が同値であることの証明を与える．まず $[hat(Phi)(#_x), hat(Phi)(#_y)] = 0$ の帰結を考える．
交換関係の $hat(Phi) (#_x)$ に $sum_alpha phi_alpha (bold(x)) exp(-i epsilon_alpha t)$ を代入すると
$
sum_(alpha, beta) [hat(a)_alpha, hat(a)_beta]
      phi_alpha (bold(x)) phi_beta (bold(x)^prime) exp(-i(epsilon_alpha + epsilon_beta) t)
      = 0
$
エネルギー固有値は非負なので$alpha != 0 or beta != 0$の時は$[hat(a)_alpha, hat(a)_beta] = 0$．
$alpha=beta=0$ の時は交換関係の性質より $0$．
結局 $forall alpha, beta ([hat(a)_alpha, hat(a)_beta] = 0)$
同様にして $forall alpha, beta ([hat(a)^dagger_alpha, hat(a)^dagger_beta]=0)$も導ける．


最後に $[hat(Phi) (#_x), hat(Pi) (#_y)] = i delta (bold(x) - bold(x)^prime)$
について考える．
同じように$Phi$および$Pi$を生成消滅演算子で展開すると．
$
sum_(alpha, beta) [hat(a)_alpha, hat(a)^dagger_beta]
      phi_alpha (bold(x)) phi^*_beta (bold(x)^prime) exp(-i(epsilon_alpha - epsilon_beta) t)
      = delta(bold(x) - bold(x)^prime)
$
右辺は時間依存性がないので左辺の時間依存性を打ち消す必要がある．
ここで「エネルギー固有値に縮退がない」という仮定を置くと
$alpha=beta$の項だけ拾えばその目的が達成される．
そこで適当な係数 $C$ を用いて
$
[hat(a)_alpha, hat(a)_beta^dagger] = C delta_(alpha, beta)
$
であればよいということになる．この係数は
$
sum_(alpha,beta)[ ~ ] = C sum_alpha
  phi_alpha (bold(x)^prime) phi^*_alpha (bold(x)) 
= C delta(bold(x) - bold(x)^prime)
$
となる．よって $C=1, [hat(a)_alpha, hat(a)^dagger_beta] = delta_(alpha, beta)$
が導かれた．
以上より エネルギー固有値に縮退がない場合に
$Phi$ と $Pi$ の間の同時刻交換関係と ${a_alpha}$ と ${a^dagger_beta}$
の間の正準交換関係が同値であることが示された．

= 「第二量子化」

電磁場の量子論のために近似モデルとして考案されたスカラー場の量子論は便利なものであった．
そのため多くの物理学者が粒子系にもこの理論を適用できないかと考えた
．
とはいえ初めから場の量がある電磁場の系とは異なり粒子系には場の自由度は存在しないため
それは困難であるように見えた．
結局本節で見たような試行錯誤により
波動関数を演算子化し粒子系の基礎方程式を導くラグランジアン密度を用いて
場の演算子に共役な運動量場との間にスカラー場と同様の
同時刻交換関係を課せば通常の正準交換関係を課したのと
同等な物理を場の形式で扱えることがわかった
#footnote[このあたりの歴史はワインバーグ @WeinbergJa1 １章に詳しい．]．
この手続きはあたかも量子論に移行して初めて現れる
波動関数をさらに量子化したように見えるので第二量子化と呼ばれている．
しかし交換関係を要請したのは一回であり量子化も一回である．
さらに要請する交換関係は同値なので物理としても同じである．

// Appendix

#show: thebibliography.with(bibliography("refs.yml"))
