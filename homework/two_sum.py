#!/usr/bin/python3.5
# import time
# import math

def two_sum():
  # input_file_name = 'two_sum_test_range_is_3_10_result_is_8.txt'
  # min_sum = 3
  # max_sum = 10

  input_file_name = 'two_sum_input.txt'
  min_sum = -10000
  max_sum = 10000

  hash_table = set() # I have to be fast today
  target_values = set(range(min_sum, max_sum + 1))

  with open(input_file_name) as f:
    for line in f:
      hash_table.add(int(line.rstrip()))

  found_values = set()

  for s in target_values:
    print(s)
    for n in hash_table:
      m = s - n
      if m != n and m in hash_table:
        found_values.add(s)
        break

  print('result:')
  print(len(found_values))



if __name__ == '__main__': two_sum()
