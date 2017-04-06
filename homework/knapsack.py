#!/usr/bin/python3.5
import time
start_time = time.time()

# from math import gcd

# Let A = 2D array
# initialize A[0, x] = 0 for x = 0, 1, ..., W
# For i = 1:n {
#   For x = 0:W {
#     A[i, x] = max(A[i - 1, x], A[i - 1, x - wi] + vi or 0 if wi > x)
#   }
# }
# Return A[n, W]

def knapsack(items, W):
  n = len(items)
  W_range = range(0, W + 1)

  last_row = [0 for x in W_range]

  for i in range(1, n + 1):
    print(100 * i / n)

    current_row = []
    (v, w) = items[i - 1]

    for x in W_range:
      if x < w:
        current_row.append(last_row[x])
      else:
        current_row.append(max(last_row[x], last_row[x - w] + v))

    last_row = current_row

  return last_row[W]

if __name__ == '__main__':
  input_file_name = 'knapsack_big_input.txt'

  W = 0
  items = []

  with open(input_file_name) as f:
    W = int(f.readline().rstrip().split(' ')[0])

    for line in f:
      [value, weight] = [int(x) for x in line.rstrip().split(' ')]

      items.append((value, weight))

  # print(W)
  # print(len(items))

  # d = W
  #
  # for (v, w) in items:
  #   d = gcd(d, w)
  #
  # print(W)
  # print(d)
  # print(W / d)

  optimal_value = knapsack(items, W)

  print('---')
  print(optimal_value) # small input: 2493893
  print('---')



print("--- %s seconds ---" % (time.time() - start_time))
