#!/usr/bin/python3.5
import time
start_time = time.time()

def mwis(weights):
  n = len(weights) - 1
  A = [0, weights[1]]

  for i in range(2, n + 1):
    A.append(max(A[i - 1], A[i - 2] + weights[i]))

  i = n
  S = set()

  while i >= 1:
    if A[i - 1] >= A[i - 2] + weights[i]:
      i -= 1
    else:
      S.add(i)
      i -= 2

  return S



if __name__ == '__main__':
  input_file_name = 'mwis_input.txt'

  weights = [None]

  with open(input_file_name) as f:
    f.readline() # 1000 weights

    for line in f:
      weights.append(int(line.rstrip()))

  S = mwis(weights)

  v = [1, 2, 3, 4, 17, 117, 517, 997]

  print("".join(["1" if i in S else "0" for i in v]))

print("--- %s seconds ---" % (time.time() - start_time))
