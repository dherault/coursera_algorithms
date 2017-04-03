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

Holy grail: linear running time

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

#### Deterministic Selection - Analysis

#### Omega(n log n) Lower Bound for Comparison-Based Sorting

### Graphs and the contraction algorithm

#### Graphs and Minimum Cuts

"Only" 20 yo algorithm

Cuts of a graph (V, E) a partition of V into two non-empty sets A and B

A graph with n vertices has 2^n possible cuts

Input: an undirected graph G = (V, E)
Goal: compute a cut with fewest number of crossing edges (i.e. a min-cut)

#### Graph Representations

Let an undirected graph of n vertices
Min # of edges: n - 1
Max # of edges: n(n - 1) / 2

n = # of vertices, m = # of edges
In most applications, m is Omega(n) and O(n²)

- Sparce graph: m is O(n) or close to it
- Dense graph: m is closer to O(n²)

Adjacency matrix: represent G using a n by n matrix

A_ij = 1 <=> G has an i-j edge
(requires n² space)

Adjacency list:
- array or list of vertices (theta(n))
- array or list of edges (theta(m))
- each edge points to its endpoints (theta(m))
- each vertex points to edges incident to it (theta(m)) (same pointers as the 3rd item)
--> requires theta(n + m) space (or theta(max(n, m)))

Which is better: it depends on graph density and operations needed.

#### Random Contraction Algorithm

```
while there are more than 2 vertices left:
  - pick a remaining edge (u, v) at random
  - merge (or "contract") u and v into a single vertex
  - remove self loops
return cut represented by final 2 vertices
```

#### Analysis of Contraction Algorithm

What is the probability of success?

Fix G = (V, E), n and m
Fix A and B the minimum cut
Let k be the # of edges crossing (A, B), those edge form F

if an edge of F is contrated, the algo fails
if not, succeeds

First iteration, probability to fail: k/m

Key observation: The degree of each vertex is at least k

We have sum(v, degree(v)) = 2m and sum(v, degree(v)) >= kn
--> m >= kn / 2

P(fail1) = k/m <= 2 / n

Second iteration:
P(success1, success2) = P(success2|succes1) * P(success1)

P(success1) >= 1 - 2/n
P(success2|success1) = 1 - k/remaining edges

remaining edges >= m - 1 >= k(n - 1) / 2

Total probability: >= (1 - 2/n)(1 - 2/(n - 1))(1 - 2/(n - 2))...(1 - 2/(n - (n - 3)))
i.e. >= 2 / (n(n - 1 )) >= 2 / n²

Problem: low success proba... but non trivial! (# of cuts is exp, proba is ²)

Solution: run the algo a large number N of times and remember the smallest cut found

How many trials needed?

P(all trials fail) <= (1 - 1 / n²)^N

We have 1 + x <= e^x

If we take N = n², P(all trials fail) <= (e^(-1/n²))^n² = 1 / e

If we take N = n²ln(n), P(all trials fail) <= 1 / n

Running time ?
polynomial in n and m but slow
Omega(n²m)

But: better algo exist

#### Counting Minimum Cuts

A graph may ave more than one min cut

Eg: a tree has n - 1 min cut

Question: what is the largest number of min cut a graph can have?
Answer: n choose 2 = n(n - 1) / 2

## Graph Search, Shortest Paths, and Data Structures

### Graph search and connectivity

#### Graph Search - Overview

Generic graph search
Goals:
- Find everything findable from a given vertex
- Don't explore anything twice
- in O(n + m) time

Generic algo: given a grap G and a vertex s
- initially s explored, all other vertices unexplored
- while possible:
  - choose an edge (u, v) with u explored and v unexplored
  - mark v explored

Claim: at the end of the algo, v explored <=> G has a path from s to v

Breath-First Search (BFS)
- explore nodes in layer
- can compute shortest paths
- can compute connected components of an undirected graph
--> O(n + m) time using a queue (FIFO)

Depth-First Search (DFS)
- explore agressively like a maze, backtrack only when necessary
- compute topological ordering of a DAG
- compute connected components in directed graphs

#### Breadth-First Search (BFS): The Basics

```
BFS(graph G, start vertex s)
  mark s as explored
  let Q = queue FIFO initialized with s
  while Q != 0:
    v = Q.shift() // The first
    for each edge (v, w):
      if w unexplored:
        mark w as explored
        Q.append(w) // At the end
```

Running time: O(ns + ms), where ns is the # of nodes reachable from s

#### BFS and Shortest Paths

Goal: compute dist(v), te fewest # of edges between s and v

extra code:
- initialize dist(v) = { 0 is s == v, +infinity otherwise }
- when considering edge (v, w):
  - if w inexplored: dist(w) = dist(v) + 1

#### BFS and Undirected Connectivity

Connected components: "the pieces" of G

Formal definition: equivalence classes of the relation u~v <=> E u-v path in G

Goal: compute all connected components
Why: check if network is disconnected

```
ComputeAllComponents(G)
  mark all nodes as unexplored
  for i=1 to n
    if i not explored
      BFS(G, i)
# note: maybe some code to describe the found components
```

Running time: O(n + m)

#### Depth-First Search (DFS): The Basics

```
DFS(G, starting vertex s)
  mark s as explored
  for every edge (s, v):
    if v is unexplored:
      DFS(G, v)
```

Could also mimick BFS with a stack instead of a queue.

Running time: O(ns + ms)

#### Topological Sort

Definition: a topological ordering of a directed graph G is a labelling F of G's nodes such that:
- the f(v)'s are the set {1,2,...,n}
- (u,v) € G => f(u) < f(v)

Motivation: Sequence tasks while respecting precedence constraints

Note: cycle => no topo ordering

Note: every DAG has a sink vertex.

Straightforward solution:
- let v be a sink vertex of G
- set f(v) = n
- recurse on G - {v}

Topological sort via DFS (Slick)

```
TopologicalSort(G)
  mark all vertices unexplored
  current_label = n
  for each vertex v:
    if v not yet explored:
      DFSForTopologicalSort(G, v)

DFSForTopologicalSort(G, s)
  mark s as explored
  for every edge (s, v):
    if v not yet explored:
      DFS(G, v)
  set f(s) = current_label
  current_label--
```

Running time: O(m + n)
Reason: O(1) time per node, O(1) time per edge
Correctness: need to show that if (u,v) is an edge, then f(u)<f(v)

#### Computing Strong Components: The Algorithm

Definition: the Strongly Connected Components (SCC) of a directed graph are the equivalence classes of the relation: u~v <=> E path u->v and a path v->u in G

Kosaraju's two-pass algorithm O(m+n)

```
KosarajuSCC(G)
  let Grev = G with all arcs reversed (could also run DFS going backward)
  run DFS-loop on Grev (goal: compute "magical ordering" of nodes)
    let f(v) = "finising time" of each vertex v
  run DFS-loop on G (goal: discover the SCC one by one)
    processing nodes in decreasing order of finising time
    save the leaders

global variable t = 0 (# of nodes processed so far)
global variable s = Null (current source vertex)

DFS-loop(G)
  Assume nodes labeled 1 to n
  for i=n to 1:
    if i not yet explored:
      s = i
      DFS(G, i)

DFS(G, i)
  mark i as explored (in a given DFS-loop)
  set leader(i) = node s
  for each edge (i,j):
    if j not yet explored:
      DFS(G, j)
  t++
  set f(i) = t
```

#### Computing Strong Components: The Analysis

2nd pass of DFS-loop begings somewhere in a sink SCC C*
First pass of DFS-loop discovers C* and nothing else

### Dijkstra's shortest-path algorithm

#### Dijkstra's Shortest-Path Algorithm

Input: a (directed) Graph G = (V, E) (n, m)
each edge has a nonnegative length "le"
source vertex s

Output: for each v€V, compute L(v) = length of the shortest s-v path in G

Assumption:
- for convinience: A v€V, E s->v path
- important: no negative edge length: le >= 0 A e € E

BFS computes the shortest path iff le=1 A e € E

```
DijkstraShortestPath(G, s)
  X = {s}             # Vertices processed so far
  A[s] = 0            # Computed shortest path distances
  B[s] = empty path   # Computed shortest paths (do not include in actual implementation)

  while X != V:
    among all edges (v, w) € E with v € X, w not € X, pick the one that minimizes A[v] + l_vw (Dijkstra's greedy criterion), call it (v*, w*)
    add w* to X
    set A[w*] = A[v*] + l_v*w*
    set B[w*] = B[v*] + (v*, w*)
```

#### Dijkstra's Algorithm: Implementation and Running Time

As such, running time is theta(nm). But we can do better.

Heap: perform insert, extract min in O(logn) time
- perfectly balanced binary tree (height: log_2(n))
- extract  min by swapping up last leaf, bubbling down
- insert bia bubling up
Also: will need the ability to delete from the middle of the heap (bubble up or down, as needed)

Two invariants:
- elements in heap = vertices of V - X
- for v € V - X, key[v] = smallest Dijkstra greedy score of an edge (u, v) € E with U € X
(key means value somehow) or +infinity if no such edge

With those invariants, extract-min yields correct vertex w* to add to X next

Maintaining the 2nd invariant:
When w is extracted from heap, (ie added to X)
```
for each edge (w, v) € E:
  if v € V - X (ie in heap):
    remove v from heap
    recompute key[v] = min(key[v], A[w] + l_wv)
    reinsert v into heap
```

Running time analysis:
- Dominated by heap operations, O(logn) each
- (n - 1) extract mins
- each edge (v, w) triggers at most one Delete/Insert combo (if v added to X first)
--> # of heap operations is O(n + m) = O(m) since the graph is weakly connected
--> Running time: O(mlogn)

### Heaps

#### Data Structures: Overview

Point: organize data so that it can be accessed quickly and usefully

Ex: lists, stacks, queues, heaps, search trees, hashtables, bloom filters, union-find, etc...

#### Heaps: Operations and Applications

- insert O(logn)
- extract min O(logn)
- heapify O(n)
- delete O(logn)

Application: sorting: HeapSort

Application: Median maintenance

Input: a sequence of n numbers
Output: for each step, the median of the numbers
Constraint: O(log i)

Solution: 1 extractMin (Hhigh) and 1 extractMax (Hlow) heap
Maintain invariant: i/2 smallest (largest) elements in Hlow (Hhigh)

### Balanced binary search trees

#### Balanced Search Trees: Operations and Applications

Things to do with a sorted array:
- Search (binary search) O(logn)
- Select (a given order statistic, ex: min/max) O(1)
- Predecessor / successor O(1)
- Rank (know # of keys less than a given value) O(logn)
- Output in sorted order O(n)

--> Apply on a static dataset

BST: like sorted array + fast log insertion and deletion

- Search O(logn)
- Select O(logn)
- Min/Max O(logn)
- Predecessor / successor O(logn)
- Rank O(logn)
- Output O(n)
- Insert / Delete O(logn)

#### Binary Search Tree Basics

Exactly one node per key
most basic version each child has left and right child + parent pointers

Search tree property: left children < current key, right children > current key

Height of a BST:
Number or "levels" (arcs)
Note: many possible search tree for a set of keys
--> height or depth: log2n <= h <= n

Search:
- Start at the root
- Traverse left/right child pointrs as needed

Insert:
- Search for k (unseccessfully)
- Rewire final Null pointer to point to new node with key k

--> theta(height)

Min/Max:
- Start at the root
- Follow left/right child pointers, return last key found

Predecessor of k:
- Easy case: if k's left subtree nonempty, return max key in left
- Otherwise: follow parent pointers until you get to a key less than k

In-order traversal
- Let r = root of search tree, with subtrees Tl and Tr
- Recurse on Tl
- print out r's key
- Recurse on Tr

Deletion (of a key k)
- Search for key
- easy case (k has no children): just delete k
- medium case (k has only one child): just "splice out" k, replace with child
- difficult case (k has 2 children):
- compute k's predecessor l
- swap k and l (note: k does not have a right child anymore)

--> theta(height)

Select and rank:
Idea: store a little bit of extra info at each tree node about the tree itself
Example augmentation: size(x) = # of tree nodes in subtree rooted at x = size(y) + size(z) + 1 (children + 1)

How to select ith order statistic (with such a tree)
- start at root x, with children y and z
- let a = size(y) (0 if no left child)
- if a = i - 1 return x's key
- if a >= i recursively compute ith order statistic of search tree rooted at y
- if a < i - i recursively compute (i - a - 1)th order statistic of tree rooted at z

--> theta(height)

#### Red-Black Trees

Balanced, logn garanteed

Idea: ensure that height stays logarithmic
--> Ops run at O(logn)

Example: red-black trees (see also AVL trees, splaytrees, b-trees, b+-trees)

Red-black invariants:

- Same as BST, plus more
- Each node red or black
- Root is black
- No 2 reds in a row (red node ==> children must be black)
- Every root-Null path (unsuccessful search) path has same number of black nodes

Example: a chain of length 3 cannot be a red-black tree. (4th invariant violated)

Example: perfectly balanced search tree

Height garantee: every red-black tree with n nodes has height O(logn) [h <= 2log2(n+1)]

#### Rotations

Left rotations: locally rebalance subtrees at a node in O(1) time

Of a parent x and right child y (x < y).

P - x - A (left st of x), (y - B (left st of y), C (right sb of y))
-->
P - y - (x - A, B), C

### Hashing; bloom filters

#### Hash Tables: Operations and Applications

purpose: maintain a possibly evolving set of stuff. (gné?)

Supported operation: insert, delete, lookup (using a "key")
--> O(1)

#### Hash Tables: Implementation Details

Setup: know the universe U that we ant to store (all IP address, all chessboard position, ...)

Goal: want to maintain a set S € U

Naive solutions:
- array-based, indexed by u, O(1) time but theta(|U|) space
- list-based, theta(|S|) space but theta(|S|) memory

Solution:
- pick n = # of "buckets" with n ~= |S| (assume |S| does not vary too much)
- choose a hash function h: U --> {0, 1, 2, ..., n - 1}
- Use array A of length n, store x in A[h(x)]

Birthday "paradox" => colisions will happen

Resolving colisions:

solution #1: separate chaining
- keep linked list in each bucket
- given a key/object x, perform insert/delete/lookup in the list in A[h(x)]

solution #2: open addressing (only o object per bucket)
- hash function now specifies probe sequence h1(x), h2(x), ...
- keep trying until find an open slot

What makes a good hash function ?

Note: in hash table with chaining, insert is O(1), O(list length) for insert/delete
Point: performance depends on a good hash function

Properties of a good has function:
- should "spread the data" as much as possible --> performance
- does not work too hard, fast to evaluate

Bad hash function, keys = phone numbers (10 digits)
|u| = 10^10, n = 10^3
- first 3 digits (terrible, some buckets will be filled, other empty)
- last 3 digits (not uniformly distributed)

next example: keys = memory locations (will be multiples of a number of 2)
- bad hash function: x mod 1000 (some buckets will be empty (odd numbers))

Quick and dirty hash functions:
- turn the object into a (preferably big) number (eg: subroutine to convert strings into intergers)
- compress that number (eg: the mod function)

How to choose the number of buckets:
- choose n to be a prime number (no multiple in common with data)
- not too close to a power of 2 or 10

#### Pathological Data Sets and Universal Hashing Motivation

The load of a hash table:
alpha = # of things that have been inserted / # of buckets of hash table

For good hash table performance:
- alpha = O(1) is necessary condition for operations to run in constant time
- ith open addressing, need alpha << 1

--> Need to control the load, make sure the # of buket increases with the # of items in the HT

Pathological data sets

--> Need a good hash function (spreads the data evenly)
--> super good hash function does not exist

Hacking prevention:
- use a cryptografic hash function (ex: SHA-2), infeasible to reverse engineer
- use randomization: design a family of hash functions such as for any data set, on average, the hash function will perform well

#### Universal Hashing: Definition and Example

What is a good random hash function ?

Def: let H be a set of hash functions U --> {0, ..., n - 1}

H is universal <=> A (x, y) € U², x != y, P(h(x) = h(y), h € H) <= 1 / n
where h is chosen uniformly at random

Example: IP address (form: (x1, x2, x3, x4) with x € {0, ..., 255})

Let n = a prime

Construction:   

### Bloom filters

#### Bloom Filters: the basics

Raison d'être lookups

Comparaison with HT:
- Pro: more space efficient
- Con: can't store objects, no deletions (with basic implementation), small fasle negative probability

Applications:
- Early spellchecker (is this a legitimate word)
- Canonical: list of forbiden passwords
- Modern: network routers

Under the hood:
- Array of n bits (0 or 1)
- k hash functions h1, ..., hk (k small constant)

Insert: for i = 1...k, set A[hi(x)] = 1
Lookup: return true if A[hi(x)] = 1 A i = 1...k

Note: no false negative, but posibly false positive

#### Heuristic Analysis

Intuition: sould be a trade-off between space and error probability

Proba that a bit (in BF) has been set to 1 = 1 - (1 - 1/n)^(k * |S|)

ie P <~= 1 - e^(-k|S|/n) = 1 - e^(-k/b) where b is # bits per object

So under assumptions, for x !€ S, false positive probalbility is <= (1 - e^-k/b)^k

How to set k ? k ~= ln2 * b

## Greedy Algorithms, Minimum Spanning Trees, and Dynamic Programming

### Greedy Algorithms

#### Introduction to Greedy Algorithms

Definition: iteratively make "myopic" decisions, hope everything works at the end

Example: Dijkstra's shortest path

#### Application: Optimal Caching

The caching problem:
- small fast memory (the cache)
- big slow memory
- process sequences of page request
- on a fault, ie cache miss, need to evict something from the cache

Theorem: the "furthest-in-future" algo is optimal (minimizes the number of cache miss)

#### Scheduling problem: Definition

Setup:
- one shared resource
- many "jobs" to do

Question: In what order to do them ?

Assume:
- Each job as a weight wj (priority)
- and a length lj

Completion time cj of job j is the sum of job length up to and including j

Goal (objective function): minimizes the weighted sum of completion times
ie min sum(j=1..n, wj*cj)

#### A Greedy Algorithm

Good idea: look at special cases, then move to general case

Edge case ideas: larger weights first, smaller length first

what if wi > wi and li > lj ?

Idea: assign a score: w - l or w/l

Breaking a greedy algorithm:
- Find examples where the two algorithms produce different input
- ex: algo 1 ot always correct, algo 2 always


#### Minimum spanning trees

Informal definition: connect some points as cheaply as possible
Applications: clustering, networking

Fast greedy algorithm:
- Prim's algorithm, 1957
- Kruskal's algorithm, 1956
=> O(mlogn) time unsing suitable data structures

Input: undirected graph G = (V, E), n,m, and a cost c for each edge

Assume adjacency list representation

Output: minimum cost tree T <= E that spans all vertices
- Has no cycles
- The subgraph (V, T) is connected (ie E a path for any 2 vertices)

Assumptions:
- G is connected
- Edge costs are distinct (algos still correct without this assumption)

#### Prim's MST Algorithm

```
X = [s], chosen arbitrarily
T = empty [invariant: X = vertices spanned by tree-so-far T]
while X != V
  let e = (u, v) be the cheapest edge of G with u € X, v !€ X
  add e to T
  add v to X
```

O(n) iterations, O(m) per iteration => O(nm) running time for naive implementation

Empty cut lemma: A graph is not connected <=> E cut (A, B) with no crossing edges

Double crossing lemma: Suppose the cycle C € E has an edge crossing the cut (A, B), then so does some other edge of C

Lonely cut corollary: If e is the only edge crossing some cut (A, B), then it is not in any cycle

The cut property:
Consider an edge e of G. Suppose there is a cut (A, B) such that e is the cheapest edge of G that crosses it. Then e belongs to the MST of G.

#### Fast Implementation

Using a Heap

Using heap to store edges with keys = cost => O(mlogn)

But there is better.

Maintain:
- Invariant 1: elements in heap = vertices of V - X
- Invariant 2: key[v] = cheapest edge (u, v) cost with u € X

Heap population:
- Scan edges, find cheapest edge for each vertice of V - X O(m), populate heap O(nlogn) --> O(m + nlogn) = O(mlogn) since m >= n-1 because G is connected

How to maintain invariant 2:
When adding a new vertex, we need to recompute keys of the neighbours

```
when V added to X:
  For each edge (v, w) € E:
    If w € V - X
      Delete w from heap
      recompute key[w] = min(key[w], c_vw)
      re insert w into heap
```

Running time with heaps:
- n - 1 inserts durring preprocessing
- n - 1 extract min per iteration
- each edge (v, w) triggers one delete/insert op
--> O(m) heap operations (since m >= n - 1)
--> O(mlogn)

### Kruskal's minimum spanning tree algorithm

#### Kruskal's MST Algorithm

Assumption: G is connected, costs are distincts (to ease proofs)

Cut property: If e is the cheapest edge crossing some cut (A, B) then e belongs to the MST

Kruskal's MST algo:
```
sort edges in order of increasing cost
T = empty
for i = 1 to m:
  if T u {i} has no cycle:
    add edge at sorted index i to T
```

#### Implementing Kruskal's Algorithm via union-find

How to check for cycles: use BFS or DFS in the graph (V, T), which contains n - 1 edges

algo --> O(nm)

Union-find data structure: maitain a partition of a set of objects

Find: return name of the group x belongs to
Union: fuse two groups into a single one

In Kruskal's: objects = vertices, groups = connected components, wen adding a new edge, check its groups and fuse them if possible

Motivation: O(1) cycle check

Idea:
- Maintain one linked structure per connected component
- each component has an arbitrary leader

Invariant:
- Each vertex points to its leader of its component ("name" of component is inherited from the leader)

key point:
- Cycle <=> edge (u,v), u and v have the same leader

Maintaining the invariant:
When new edge (u, v) added to T, connectected components of u and v merge.
When two components merge, have smaller one  inherit the leader of the larger one. (maintain a count of each component)

Running time:
- O(mlogn) for sorting
- O(m) for cycle checks (O(1) per iteration)
- O(nlogn) for leader pointe updates
- --> O(mlogn) total

#### Application to Clustering

Goal: given n points, classify them into "coherent groups"

Assumptions:
- As input, given a (dis)similarity measure - a distance d(p, q) between each pair of points
- d is symmetric
- we know k = # of clusters desired

Call points p, q separated if they're assigned to different clusters

```
assign each point to a separate clusters

repeat until only k clusters:
  let p, q be the closest pair of separated points
  merge the 2 clusters of p and q into one
```

--> Called single-link clustering

### Huffman codes

#### Introduction: prefix-free variable encoding

Binary code: map every char of an alphabet (S: say 32 chars) to a binary string

Obvious encoding: 32 5-bit binary strings (eg: ASCII)

Can we do better ? Yes if some chars are more frequent than others, using a variable length code

Suppose { A, B, C, D } = { 00, 01, 10, 11 }. Assume we use { 0 01 10 1 } instead

Then 001 could be AAD or AB

Prefix-free codes: make sure that for every pair i, j € S, neither the encoding F(i), F(j) is a prefix of the other.

Example: { 0, 10, 110, 111 } prefix-free

Example: A: 60% 0, B: 25% 10, C: 10% 110, D: 5% 111 --> 1.55 bits needed of average per char

#### Problem definition

Codes as tree:

Goal: find the best prefix-free encoding for a given set of character frequencies

Fact: binary codes <--> binary trees

Left child: 0, right child: 1

Prefix-free <--> labelled node are leaves

To decode: repeatedly follow path from root until you hit a leaf

Input: probability pi for each character i € S

Notation: if T = tree with leaves <--> symbols of S, then L(T) (average encoding length) = sum(on S, pi * depth of i in T)

#### A greedy algorithm

Natural but sub-optimal idea: top-down, divide and conquer
- Partition S into S1 and S2 each with 50% of total frequency
- Recursively compute T1 for S1 and T2 for S2, return

Huffman's (optimal) idea:
- build tree bottom up

```
HuffmanCodes(S) {
  if |S|==2 return A --0-- x --1-- B
  Let a, b € S ave the smallest frequencies
  Let S' = S with a, b replaced wih new symbol ab
  define pab = pa + pb
  T' = HuffmanCodes(S')
  Extend T' (with leaves = S') to a tree T (with leaves = S) by splitting leaf ab into a --0-- x --1-- b
  return T
}
```

Running time: O(n²) where n = |S|
Speed up: heap
Even better: sorting + 2 queues

### Introduction to dynamic programming

#### Weighted Independent Sets in Path Graphs

Input: a path graph G = (V, E) with nonnegative weights on vertices

Output: subset of non-adjacent vertices (an independant set) of maximum total weigth

#### WIS in Path Graphs: Optimal Substructure

Optimal substructure:
critical step: reason about structure of an optimal solution

Upshot: a max-weight IS must be either
- a max-weight IS of G', or
- vn + a max-weight IS of G''

If we knew wether or not vn is in the max-weight IS, then  could recursively compute the max-weight IS of G' or G''.

Let Gi = first i vertices of G

Plan: populate array A left to right with A[i] = Value of max-weight IS of Gi

```
A[0] = 0
A[1] = w1

For i=1:n {
  A[i] = max(A[i - 1], A[i - 2] + wi)
}
```
#### WIS in Path Graphs: A Reconstruction Algorithm

Reconstruction algorithm
```
Let A = filled in array A

Let S = empty
While i >= 1 {
  If A[i - 1] >= A[i - 2] + wi {
    i -= 1
  }
  else {
    add vi to S
    i -= 2
  }
}

return S
```

#### Principles of Dynamic Programming

- Identify a small number of subproblems
- Can quickly + correctly solve "larger" subproblems given the solutions to the smaller subproblems
- After solving all subproblems, can quickly compute the final solution
