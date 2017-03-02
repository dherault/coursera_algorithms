#!/usr/bin/python3.5
# import time
# import math

class MinHeap:
  def __init__(self):
    self.list = [0]
    self.size = 0

  def getRoot(self):
    if self.size == 0:
      return None
    else:
      return self.list[1]

  def insert(self, value):
    self.list.append(value)
    self.size += 1
    self.percolateUp(self.size)

  def percolateUp(self, i):
    half_i = i // 2

    while half_i:

      if self.list[i] < self.list[half_i]:
        self.list[i], self.list[half_i] = self.list[half_i], self.list[i]

      i = half_i
      half_i //= 2

  def extractMin(self):
    self.list[self.size], self.list[1] = self.list[1], self.list[self.size]
    self.size -= 1

    retval = self.list.pop()

    self.percolateDown(1)

    return retval

  def percolateDown(self, i):
    while 2 * i <= self.size:
      minChildIndex = self.minChildIndex(i)

      if self.list[i] > self.list[minChildIndex]:
        self.list[i], self.list[minChildIndex] = self.list[minChildIndex], self.list[i]

      i = minChildIndex

  def minChildIndex(self, i):
    if 2 * i + 1 > self.size:
      return 2 * i
    elif self.list[2 * i] < self.list[2 * i + 1]:
      return 2 * i
    else:
      return 2 * i + 1

class MaxHeap:
  def __init__(self):
    self.list = [0]
    self.size = 0

  def getRoot(self):
    if self.size == 0:
      return None
    else:
      return self.list[1]

  def insert(self, value):
    self.list.append(value)
    self.size += 1
    self.percolateUp(self.size)

  def percolateUp(self, i):
    half_i = i // 2

    while half_i:

      if self.list[i] > self.list[half_i]:
        self.list[i], self.list[half_i] = self.list[half_i], self.list[i]

      i = half_i
      half_i //= 2

  def extractMax(self):
    self.list[self.size], self.list[1] = self.list[1], self.list[self.size]
    self.size -= 1

    retval = self.list.pop()

    self.percolateDown(1)

    return retval

  def percolateDown(self, i):
    while 2 * i <= self.size:
      maxChildIndex = self.maxChildIndex(i)

      if self.list[i] < self.list[maxChildIndex]:
        self.list[i], self.list[maxChildIndex] = self.list[maxChildIndex], self.list[i]

      i = maxChildIndex

  def maxChildIndex(self, i):
    if 2 * i + 1 > self.size:
      return 2 * i
    elif self.list[2 * i] > self.list[2 * i + 1]:
      return 2 * i
    else:
      return 2 * i + 1

def median_maintenance():
  # input_file_name = 'median_maintenance_test_result_is_142.txt'
  # input_file_name = 'median_maintenance_test_result_is_9335.txt'
  input_file_name = 'median_maintenance_input.txt'

  minh = MinHeap()
  maxh = MaxHeap()

  median = 0
  total = 0

  with open(input_file_name) as f:
    for line in f:
      n = int(line.split('\n')[0])

      print(n)

      if minh.size < maxh.size:
        if n < median:
          minh.insert(maxh.extractMax())
          maxh.insert(n)
        else:
          minh.insert(n)

        median = maxh.getRoot()

      elif minh.size > maxh.size:
        if n >= median:
          maxh.insert(minh.extractMin())
          minh.insert(n)
        else:
          maxh.insert(n)

        median = maxh.getRoot()

      else:
        if n < median:
          maxh.insert(n)
          median = maxh.getRoot()
        else:
          minh.insert(n)
          median = minh.getRoot()

      total += median

      # print(median)
      print(total % 10000)
      # print(maxh.list)
      # print(minh.list)


if __name__ == '__main__': median_maintenance()
