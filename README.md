# Algorithms

Notes and homework for Coursera's "Algorithms" online class :books:

## Divide and Conquer, Sorting and Searching, and Randomized Algorithms

### Introduction

#### Why study algorithms?

foobar

#### Integer multiplication

input: 2 n-digit numbers x and y
output: the product x.y

Primitive operation (for now): add or multiply 2 single-digit numbers.

Algo is "correct": the algo will eventually terminate and give the correct answer.

What about the amount of time needed?

upshot: nOperations <= some constant
upshot: can we do better?

#### Karatsuba multiplication

**Example:**
x = 5678, a = 56, b = 78
y = 1234, c = 12, d = 34

step1: a.c = 672
step2: d.b = 2652
step3: (a+b)(c+d) = 6164
step4: (3) - (2) - (1) = 2840
step5: (1) * 10^4 + (2) + (4) * 10^2 = 7006652, the output

**Recursive algorithm:**
write x = 10^(n/2) * a + b, y = 10^(n/2) * c + d where a, b, c, d are n/2-digit numbers.

x.y = ...
    = 10^n * ac + 10^(n/2) * (ad + bc) + bd (¤)

idea: recursively compute ac, ad, bc, bd then compute ¤

step1: recursively compute ac
step2: recursively compute bd
step3: recursively compute (a + b)(c + d) = ac + ad + bc + bd
Gauss's trick: (3) - (1) - (2) = ad + bc

upshot: only need 3 recursive multiplications and some additions

#### About the course

spoilers

#### Merge Sort: Motivation and Example

Improves over Selection, Insertion, Bubble sorts
Good example for worst case behavior and asymptotic analysis

Recursion tree methods --> "Master method"

input: array of n numbers
output: array of n numbers, sorted from min to max

Assumption: numbers are distincts (does not change much)

#### Merge Sort: Pseudocode / Merge Sort: Analysis

ignore base cases

```
C = output[length = n]
A = first sorted array [n/2]
B = second sorted array [n/2]
i = 1
j = 1

for k = 1 to n
  if A(i) < B(j)
    C(k) =  A(i)
    i++
  else
    C(k) = B(j)
    j++
end
```

Running time ? n numbers to sort

Upshot: <= 4n + 2
        <= 6n (since n >= 1)

Claim: merge sort requires <= 6n * (log(n) + 1) operations to sort n

At each level j=0,1,...,log(n) there are 2^j subproblems each of size n/(2^j)

#### Guiding Principles for Analysis of Algorithms

- worst case analysis
- don't pay attention to constant factors and lower order terms
- asymptotic analysis (focus on large input sizes)

Fast algorithm ~= worst-case running time grows slowly with input size

Holy grail: linear runnning time

### Asymptotic analysis

#### The gist

Importance: vocabulary, "sweet spot", sharp enough to differentiate algo

High-level idea: suppress constant factors and lower-order terms.

constant factors: too system-dependant
lower-order terms: irrelevant for large input

--> 6n(log(n) + 1) ==> nlog(n)

**examples:**

- One loop: does an array contains a given integer  --> O(n)
- Two loops: does A or B contain t? --> O(n)
- Two nested loops: do A and B have a number in common --> O(n²)
- Two nested loops (II): does A have duplicate entries --> O(n²)

#### Big-Oh Notation

T(n) = function on N

`T(n) = O(f(n)) <--> E c, n0 > 0 / T(n) <= c * f(n) A n >= n0`

#### Basic Examples

##### Example 1

Claim:  if T(n) = a_k * n^k + ... + a_0 * n^0
        then T(n) = O(n^k)

Proof: choose n0 = 1 and c = |a_k|+...+|a_0|

##### Examples 2

Claim: A k >= 1, n^k != O(n^(k-1))

Proof: by contradiction. Suppose n^k = O(n^(k-1))
`E c, n0 > 0 / n^k <= c * n^(k -1), A n >= 0`
ie `n <= c, A n >= 0`

#### Big Omega and Theta

Omega notation:
`T(n) = Om(f(n)) <--> E c, n0 / T(n) >= c * f(n) A n >= n0`

Theta notation: O(f(n)) AND Om(f(n))
`T(n) = th(f(n)) <--> E c1, c2, n0 / c1 * f(n) <= T(n) <= c2 * f(n) A n >= n0`

Little Oh notation:
`T(n) = o(f(n)) <--> A c > 0, E n0 / T(n) <= c * f(n) A n >= n0`

### Divide and conquer algorithms

#### O(n log n) Algorithm for Counting Inversions I

The divide and conquer paradigm:
- Divide the problem into smaller subproblems
- Conquer the subproblems via recursive calls
- Combine solutions of subproblems into the original solution

Input: array A containing the numbers 1, 2, 3, ... in some arbitrary order
Output: number of inversions = number of pairs (i,j) of indices with i < j and A[i] > A[j]
Motivation: numerical similarity measure between two ranked lists
Largest possible number of inversions: (n C 2) = n(n - 1) / 2

Note: `(n choose k) = n! / (k!(n - k)!)`

Brute force algorithm: O(n²) time

3 types of inversions, (i,j) <= n / 2; (i,j) > n / 2; i <= n / 2 < j (split inversion)

High-level algorithm:
```
CountInversions(array A of length n)
  if (n == 1) return 0
  else
    x = CountInversions(first half of A)
    y = CountInversions(second half of A)
    z = CountSplitInversions(A, n) // TODO
    return x + y + z
```

#### O(n log n) Algorithm for Counting Inversions II

Idea: piggyback on merge sort

```
SortAndCountInversions(array A of length n)
  if (n == 1) return (A, 0)
  else
    (B, x) = SortAndCountInversions(first half of A)
    (C, y) = SortAndCountInversions(second half of A)
    (D, z) = MergeAndCountSlipInversions(B, C, n) // TODO
    return (D, x + y + z)
```

MergeAndCountSlipInversions: when an element of 2nd array C gets copied into output D, increment total by number of elements remaining in 1st array B

Run time: O(nlog(n))

#### Strassen's Subcubic Matrix Multiplication Algorithm

n by b matrices

Classic multiplication: O(n^3)

Recursive algorithm n#1:

Idea: write X = (A, B; C, D) and Y = (E, F; G, H) where A...H are n/2 by n/2 matrices

Then: XY = (AE + BG, AF + BH; CE + DG, CF + DH)

--> O(n^3)

Strassen's algorithm (1969)

- Step 1: recursiely compute only 7 (cleverly chosen) products
- Step 2: do the necessary (clever) additions and subtractions (still O(n²) time)
- Fact: better than cubic time! (see next lecture)

7 products:
- P1 = A(F - H)
- P2 = (A + B)H
- P3 = (C + D)E
- P4 = D(G - E)
- P5 = (A + D)(E + H)
- P6 = (B - D)(G + H)
- P7 = (A - C)(E + F)

Claim: XY = (P5 + P4 - P2 + P6, P1 + P2; P3 + P4, P1 + P5 - P3 - P7)

#### O(n log n) Algorithm for Closest Pair

input: n points (€ R²) in the plane
output: the closest pair in the plane (euclidian distance)
assumption: all points have distinct x and y coordinates

Brute force: O(n²)

1-D version of closest pair:
- sort points O(nlogn)
- return closest pair of adjacent points O(n)

Goal: O(nlogn) for 2-D version

High level approach:

- make copies of points sorted by x-coordinate (Px) and y-coordinates (Py) O(nlogn)
- use divide-conquer

```
ClosestPair(Px, Py)
  # base case omitted
  let Q = leftmost half part of P
      R = rightmost
  form Qx, Qy, Rx, Ry (takes O(n) time)

  (p1, q1) = ClosestPair(Qx, Qy)
  (p2, q2) = ClosestPair(Rx, Ry)

  Either the closest pair is entirely in Q, either in R
  unlucky case: the pair is splitted in Q and R

  (p3, q3) = ClosestSplitPair(Px, Py)

  return best of 3 pairs
```

If ClosestSplitPair is O(n) then ClosestPair is O(nlogn)

ClosestSplitPair:
key idea: only need to compute the unlucky case. Because the recursions will eventually make all cases unlucky

```
ClosestPair(Px, Py)
  let Q = leftmost half part of P
      R = rightmost
  form Qx, Qy, Rx, Ry (takes O(n) time)

  (p1, q1) = ClosestPair(Qx, Qy)
  (p2, q2) = ClosestPair(Rx, Ry)

  let D = min(d(p1, q2), d(p2, q2))

  (p3, q3) = ClosestSplitPair(Px, Py, D)

  return best of 3 pairs

ClosestSplitPair(Px, Py, D)
  let x = biggest x-coordinate in left of P (O(n))
  let Sy = points of P with x-coordinates in [x - D, x + D], sorted by y-coordinates (O(n))
  let best = D
  let best_pair = null

  for i = 1 to len(Sy) - 1
    for j = 1 to min(7, len(Sy) - i)
      let p, q = i-th, (i + j)-th points of Sy

      if d(p, d) < best
        best = d(p, q)
        best_pair = (p, q)
    # O(1) time
  # O(n) time

  return best_pair
```

Satisfies the correctness requirements u_u
Claim:
let p € Q, q € R be a split pair with d(p, q) < D
Then:
(A) p and q are members of Sy
(B) p and q are at most 7 positions apart in Sy

### The master method (aka master theorem)

#### Motivation

Re: Karatsuba

T(n) = maximum number of operations this algorithm needs to multiply two n-digit numbers

Recurrence: express T(n) in term of running time of recursive calls

- Base case: T(1) <= a constant
- A n > 1: T(n) <= 3T(n / 2) + O(n)

#### Formal Statement

assumption: all subproblems have equal size (previously: n / 2)

- A n > n_base_case: T(n) <= a * T(n / b) + O(n^d)
- a: number of recursive calls (>= 1)
- b: input size shrinkage factor (> 1)
- d: exponent in running time of "comble step" (>= 0)
- (a, b, d) independent of n

```
if T(n) <= a * T(n / b) + O(n^d)
then T(n) = {
  if a = b^d: O(n^d * log(n))
  if a < b^d: O(n^d)
  if a > b^d: O(n^log_b(a))
```

Only gives upper bond, but can use theta

log: no base ? base does not matter (constant factor change) UNLESS in exponent, then the constant factor matters, a lot.

#### Examples

Merge sort:
a = 2
b = 2
d = 1
--> a = b^d --> O(nlogn)

Binary search:
a = 1
b = 2
d = 0
--> a = b^d --> O(logn)

Karatsuba without Gauss's trick:
a = 4
b = 2
d = 1
--> a > b^d --> O(n^(log2(4))) = O(n²)

Karatsuba:
a = 3
b = 2
d = 1
--> a > b^d --> O(n^(log2(3))) = O(n^1.59)

Strassen's matrix multiplication:
a = 7
b = 2
d = 2
--> a > b^d --> O(n^(log2(7))) = O(n^2.81)

#### Proof

Assume: recurrence is
(i)  T(1) < c
(ii) T(n) < a * T(n / b) + c * n^d
and n is a power of b

Idea: generalize merge sort analysis ie use a recursion tree

at each level j=0,1,2,...,log_b(n) there are a^j subproblems each of size n / b^j

The recursion tree

Level 0: n
Level 1: n/b, n/b, ... (a branches)
...
Level log_b(n): base cases

Work at a single level (ignoring work in recursive calls): <= a^j * c * (n / b^j)^d
where:
- a^j is the # of level-j subproblems
- n / b^j is the size of each level-j subproblem
- c * (n / b^j)^d is the work per level-j subproblem

Total work: summing over all levels j:
total work <=  c * n^d * sum(j=0...log_b(n), (a / b^d)^j) (¤)

How to think about (¤)

- a = rate of subproblem proliferation (RSP)
- b^d = rate of work shrinkage (RWS) (per subproblem)

- If RSP < RWP then the amount of work is decreasing with the recursion level j (expect O(n^d * logn))
- If RSP > RWP then the amount of work is increasing with the recursion level j (extect O(n^d))
- If RSP = RWP then the amount of work is the same at every recursion level j (expect O(# leaves))

If a = b^d then (¤) = c * n^d * (log_b(n) + 1)
                    = O(n^d * log(n)) We can suppress constant factors

A r != 1, 1 + r + r² + ... + r^k = (r^(k + 1) - 1) / (r - 1)
upshot:
- If r < 1, <= 1 / (r - 1) = a constant (independent of k)
- If r > 1, <= r^k * (1 + 1 / (r - 1))

### Quicksort algorithm

#### Quicksort: Overview

- Greatest hit algorithm
- Prevalent in practice
- O(nlogn), works in place

input: array of n unsorted numbers
output: array of n sorted numbers
assume: all array entries distinct

key idea: partition array arounda pivot element
- pick element of array
- rearrange array so that left to pivot === less than pivot, right === greater

Note: puts pivot in its "rightful position"

Partition:
- linear O(n) time, no extra memory
- reduces problem size

High-level description:
```
QuickSort(A of length n)
  if n == 1 return
  p = ChoosePivot(A)
  partition A around p
  recursively sort first and second part (not including p)
```

#### Partitioning Around a Pivot

The easy way: using O(n) extra memory. bouh!

In-place implementation:
Assume: pivot is the first element of array
High-level idea:
- single scan through array
- invariant: everything looked at so far is partitioned

```
Partition(A, l, r) // input = A[l...r]
  p = A[l]
  i = l + 1
  for j = l + 1 to r
    if A[j] < p // if not, do nothing
      swap A[j] and A[i]
      i++
  swap A[l] and A[i - 1]
```

Running time: O(n) where n = r - l + 1 (O(1) work per array entry)

Correctness: the for loop maintains the invariants (A[l + 1]...A[i - 1] < p and A[i]...A[j - 1] > p)

#### Correctness of Quicksort

By induction. Think about the partitioned array. The "parts" are smaller than n, hence the proof.

#### Choosing a Good Pivot

The running time of QuickSort depends on the quality of the pivot

If pivot is always the first element, on a sorted array QuickSort performs theta(n²)
Because n + (n - 1) + ... + 1 < n²

If pivot is always the median element of te array, on any array QuickSort performs theta(nlogn)

How to choose pivots?
Key idea: random pivots!

QuickSort theorem:
for every input array of length n, the average running time of QuickSort with random pivots is O(nlogn)

### QuickSort analysis

#### Analysis: A Decomposition Principle

Fixed input array A of length n

sample space Om = all possible outcomes of random choices in QuickSort (i.e. pivot sequences)

random variable A s € Om, E C(s) the # of comparisons between two input elements

lemma: running time of QuickSort dominated by companrisons

notation:
z_i = i-th smallest item in A

X_ij(s) = # of times z_i and z_j get compared in pivot sequence s (indicator random variable: ie can take only 0 or 1)

A s, C(s) = sum(i=1...n - 1, sum(j=i + 1...n, X_ij(s)))

By linearity of ... this is cumbersome to type

A general decomposition principle:
- Identify random variable Y that you really care about
- Express Y as a sum of indicator random variable: Y = sum(X)
- Apply linearity of expectation: E[Y] = sum(E[x]) = sum(P(x = 1))

### Linear-time Selection

#### Randomized Selection - Algorithm

input:  array A of n disctinct numbers and a number i
output: i-th order statistic ie the i-th smallest element of A
example: median

The 1st order statistic is the min --> trivial with linear scan
The nth order statistic is the max --> idem

Easy to do with sorting, but nlogn

```
RandomizedSelection(A of length n, order statistic i)
  if n == 1 return A[0]
  choose pivot p uniformly at random
  partition A around p

  let j = new index of p

  if i == j return p // lucky

  if j > i return RandomizedSelection(first part of A of length j - 1, i)
  if j < i return RandomizedSelection(second part of A of length n - j, i - j)
```

Claim: RandomizedSelection is correct (proof: by induction)

Running time ? Depends on quality of the pivot
worst pivot: theta(n²)
key: find a pivot giving a balanced split
best pivot: the median (lol)
--> would get recurrence T(n) <= T(n / 2) + O(n) --master method--> T(n) = O(n)
Hope: random pivot is "pretty good" "often enough"

RandomizedSelection theorem: for every input array of length n, the average running time is O(n)

#### Randomized Selection - Analysis

Notation: RS is in phase j if the current array size is between (3/4)^(j+1) * n and (3/4)^j * n
X_j = number of recursion calls during phase j

...

coin flip: E[heads] = 1 + 0.5 * E[heads]
- 1 is the flip
- 0.5 is the probability that the flip fails, in which case we would try again

#### Deterministic Selection - Algorithm

key idea: use the "median of medians" as the pivot

```
ChoosePivot(A of length n)
  break A into n / 5 groups of size 5 each
  sort each group
  copy the n / 5 medians (ie middle elements of each sorted group) into new array C
  recursively compute median of C
  return it as the pivot
  partition around it

i.e.:

DSelect(A of length n, order statistic i)
  Break A into n / 5 groups of size 5, sort the groups
  C = the n / 5 "middle elements"
  p = DSelect( of length n / 5, n / 10)
  Partition A around P
  if j = i return p
  if j < i return DSelect(first part of A of length j - 1, i)
  if j > i return DSelect(second part of A of length n - j, i - j)
```

Not as good as RSelect in practice
- worse constants
- not in-place
