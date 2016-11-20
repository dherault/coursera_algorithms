#!/usr/bin/python3.5
# import sys

def choose_pivot_l(a, l, r):
  return a[l]

def choose_pivot_r(a, l, r):
  return a[r]

def quick_sort(a, choose_pivot, l = 0, r = None):
  if r == None: r = len(a)

  # print('quick_sort:', l, r, a, a[l:r])

  if r - l < 2: return a

  p = choose_pivot(a, l, r)
  i = l + 1

  # print('pivot:', p)

  for j in range(l + 1, r):
    if a[j] <= p:
      a[i], a[j] = a[j], a[i]
      i += 1

  a[l], a[i - 1] = a[i - 1], a [l]

  # print('i:', i)
  # input()

  quick_sort(a, choose_pivot, l, i - 1)
  quick_sort(a, choose_pivot, i, r)

  return a

def main():
  # input = sys.argv[1:]

  input1 = [5, 1, 1, 12, 12, 2, 3, 4, 11, 8, 9, 10, 7, 6]
  input2 = [1, 1, 1]
  # input3 = [1, 2, 3]
  input3 = [10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10]

  print('input', input1)
  print('output', quick_sort(input1, choose_pivot_l))
  print('input', input2)
  print('output', quick_sort(input2, choose_pivot_l))
  print('input', input3)
  print('output', quick_sort(input3, choose_pivot_l))

if __name__ == '__main__': main()
