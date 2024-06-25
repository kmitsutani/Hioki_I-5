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

以下では簡単のためにすべてがボソンの場合を考える．

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
ただしここで$hat(F)^((1))$の行列要素を$physica.bra(i)F^((1))physica.ket(j) = f_(i,j)$と表記した．
$
physica.bra(n_1 - 1\, n_2 + 1)hat(F)physica.ket(n_1\, n_2)
= c_n (n_1, n_2) f_(2, 1) physica.bra(n_1-1\,n_2+1){
  &sum_(P(1)=1)physica.ket(2)physica.ket(P(2))dots.c physica.ket(P(N))\
 &+sum_(P(2)=1)physica.ket(P(1))physica.ket(2)dots.c physica.ket(P(N))\
 &dots\
 &+sum_(P(N)=1)physica.ket(1)dots physica.ket(P(N-1)) physica.ket(2)
}
$
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
#let exp=physica.expectationvalue
#let ket=physica.ket
#let bra=physica.bra
#let braket=physica.braket
#let ketbra=physica.ketbra
を導入すると$hat(F)$の期待値は
$
exp(F) = c_N^2 sum_(P, P^prime) sum_(i=1)^N 
  bra(P^prime (1))bra(P^prime (2)) dots.c bra(P^prime (N))
  hat(F)^((1))_i
  ket(P^prime (1))ket(P^prime (2)) dots.c ket(P^prime (N))

$
占有数状態に関する対角和なので占有数分布を変えない項が残る．
したがって
$[hat(F)^((1))physica.ket(P(i))] arrow sum_sigma^mu f_(sigma, sigma)physica.ket(sigma)$
とすればよい．
今上式を$exp(F) = sum_i F_i$と書いたとすると式の対称性より$f_i$は$i$に依らない．
そこで$i=N$ととると
$
F_n &= c_N^2 sum_(P, P^prime) delta_(P^prime (1), P(1))dots delta_(P^prime (N-1), P(N-1))
  [sum_(alpha=1)^mu f_(alpha, alpha) braket(P^prime (n), alpha)] \
  &= c_N^2 sum_(alpha=1)^mu f_(alpha, alpha) sum_(P, P^prime)
  delta_(P^prime (1), P(1)) delta_(P^prime (2), P(2)) dots delta_(P^prime (N), alpha)
$
この式の$P, P^prime$に関する和は$P^prime (N) = alpha$の制約下で
$1dots n_1$に$P$および$P^prime$を作用させた結果が一致するような状態割り当ての数なので
$(N-1)!\/(n_1 ! n_2 ! dots (n_alpha - 1)! dots n_N !)$となる．
したがって
$
F_N = sum_(alpha=1)^mu f_(alpha, alpha) 
  [(product_sigma n_sigma !)/ N!] [N! / (product_sigma n_sigma)] n_alpha / N 
  = f_(alpha, alpha) n_alpha / N 
$
これを$exp(F) = sum_(i=1)^N F_i$ に代入すると
$
exp(F) = F_N times N = sum_(alpha=1)^(mu) bra(alpha) hat(F)^((1)) ket(alpha) n_alpha
$
となる．

= Fock空間の導入

= Fock空間上での生成消滅演算子の正準交換関係

= 「第二量子化」


= 第二量子化という言葉について

大抵の教科書には

- 「誤解をを招く言葉であるが広く広まっているので使用する．」
- 「（偉大な）先人たちの黎明期の混乱を反映する言葉である．」
というようなことが書いてある．

しかし具体的にどのような混乱があったかまでは書いていないのでいまいち釈然としない．
いくつか積んである場の量子論の教科書を読んでみたらワインバーグが詳しそうであった．

// Appendix
#show: appendix

= 量子論における粒子に関する（与太）話

= 速習ブラケット

#show: thebibliography.with(bibliography("refs.yml"))
