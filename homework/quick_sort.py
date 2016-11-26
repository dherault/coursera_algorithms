#!/usr/bin/python3.5
# import sys

def swap(a, x, y):
  a[x], a[y] = a[y], a[x]

def choose_pivot_l(a, l, r):
  return l

def choose_pivot_r(a, l, r):
  return r - 1

def choose_pivot_median_of_three(a, l, r):
  n = r - l
  i1 = l
  i2 = l + (int(n / 2) if n % 2 else (int(n / 2) - 1))
  i3 = r - 1

  m1 = a[i1]
  m2 = a[i2]
  m3 = a[i3]

  if (m2 >= m1 and m1 >= m3) or (m3 >= m1 and m1 >= m2): return i1
  if (m1 >= m2 and m2 >= m3) or (m3 >= m2 and m2 >= m1): return i2
  return i3

def quick_sort(a, choose_pivot, l = 0, r = None, m = [0]): # /!\ m
  if r == None: r = len(a)

  if r - l < 2: return a

  p_index = choose_pivot(a, l, r)
  p = a[p_index]

  swap(a, p_index, l) # Put the pivot where it does not bother anyone

  i = l + 1
  m[0] += r - l - 1

  for j in range(l + 1, r):
    if a[j] <= p:
      swap(a, i, j)
      i += 1

  swap(a, l, i - 1) # Put the pivot we it should be

  quick_sort(a, choose_pivot, l, i - 1, m)
  quick_sort(a, choose_pivot, i, r, m)

  return m[0]

def main():
  # input = sys.argv[1:]

  # input1 = [5, 1, 1, 12, 7, 12, 10, 6]
  # input2 = [1, 1, 1]
  # input3 = [10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10]

  # print(input1, choose_pivot_median_of_three(input1, 0, len(input1)))
  # print(input2, choose_pivot_median_of_three(input2, 0, len(input2)))
  # print(input3, choose_pivot_median_of_three(input3, 0, len(input3)))

  # print('input1', input1)
  # print('input2', input2)
  # print('input3', input3)
  #
  # print('output1', quick_sort(input1, choose_pivot_median_of_three))
  # print('output2', quick_sort(input2, choose_pivot_median_of_three))
  # print('output3', quick_sort(input3, choose_pivot_median_of_three))
  #
  # print(input1)
  # print(input2)
  # print(input3)

  # print(choose_pivot_median_of_three(input1, 0, len(input1)))

  with open('quick_sort_input.txt') as f:
  # with open('quick_sort_test_1000.txt') as f:
    input = [int(x) for x in f.readlines()]
    # print('choose_pivot_l:', quick_sort(input, choose_pivot_l))
    # print('choose_pivot_r:', quick_sort(input, choose_pivot_r))
    print('choose_pivot_median_of_three:', quick_sort(input, choose_pivot_median_of_three))

if __name__ == '__main__': main()
