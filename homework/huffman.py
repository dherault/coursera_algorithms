#!/usr/bin/python3.5
import sys
from heapq import *

# Run you fools
sys.setrecursionlimit(1500)

class BinaryTree:
  """A special BT for Huffman's algorithm"""

  def __init__(self, left = None, right = None, value = None):
    self.value = value
    self.left = BinaryTree(None, None, left) if left else None
    self.right = BinaryTree(None, None, right) if right else None

  def replaceLeaf(self, leaf_value, tree):
    if self.value == leaf_value:
      self.value = tree.value
      self.left = tree.left
      self.right = tree.right
      return

    if self.left:
      self.left.replaceLeaf(leaf_value, tree)

    if self.right:
      self.right.replaceLeaf(leaf_value, tree)

  # The assignment asks for the max and min length of the tree
  def getLength(self, fn):
    if not (self.left and self.right):
      return 0

    leftLength = self.left.getLength(fn) if self.left else 0
    rightLength = self.right.getLength(fn) if self.right else 0

    return 1 + fn(leftLength, rightLength)


# if |S|==2 return A --0-- root --1-- B
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

  T.replaceLeaf(ab, BinaryTree(a, b))

  return T


if __name__ == '__main__':
  input_file_name = 'huffman_input.txt'

  # alphabet is a heap
  alphabet = []

  with open(input_file_name) as f:
    f.readline()
    char_name = 0
    for line in f:
      char_name += 1
      char_weight = int(line.rstrip())

      heappush(alphabet, (char_weight, char_name))

  tree = huffman(alphabet)

  print(tree.getLength(max))
  print(tree.getLength(min))
