#!/usr/bin/python3.5
import time
start_time = time.time()

import sys
from heapq import *

sys.setrecursionlimit(1500)

V = {}

class BinaryTree:
  """A special BT for Huffman's algorithm"""

  def __init__(self, left = None, right = None, value = None):
    self.value = value
    self.left = BinaryTree(None, None, left) if left else None
    self.right = BinaryTree(None, None, right) if right else None
    V[self.value] = self

  def getLength(self, fn):
    if not (self.left and self.right):
      return 0

    leftLength = self.left.getLength(fn) if self.left else 0
    rightLength = self.right.getLength(fn) if self.right else 0

    return 1 + fn(leftLength, rightLength)


# if |S|==2 return A --0-- x --1-- B
# Let a, b € S ave the smallest frequencies
# Let S' = S with a, b replaced wih new symbol ab
# define pab = pa + pb
# T' = HuffmanCodes(S')
# Extend T' (with leaves = S') to a tree T (with leaves = S) by splitting leaf ab into a --0-- x --1-- b
# return T

def huffman(alphabet):

  if len(alphabet) == 2:
    return BinaryTree(alphabet[0][1], alphabet[1][1])

  (pa, a) = heappop(alphabet)
  (pb, b) = heappop(alphabet)

  ab = str(a) + '_' + str(b)

  heappush(alphabet, (pa + pb, ab))

  T = huffman(alphabet)

  V[ab].__dict__.update(BinaryTree(a, b).__dict__)

  return T


if __name__ == '__main__':
  input_file_name = 'huffman_input.txt'

  # Build alphabet heap from file
  alphabet = []

  with open(input_file_name) as f:
    f.readline() # 1000 chars

    char_name = 0
    for line in f:
      char_name += 1
      heappush(alphabet, (int(line.rstrip()), char_name))

  tree = huffman(alphabet)

  print(tree.getLength(min)) # 9 bits
  print(tree.getLength(max)) # 19 bits



print("--- %s seconds ---" % (time.time() - start_time))
