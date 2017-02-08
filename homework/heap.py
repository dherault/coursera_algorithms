class Heap:
  """Min-Heap simple implementation"""
  def __init__(self): # Maybe add a list arg
    self.list = [0]
    self.size = 0

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

  def delMin(self):
    self.list[self.size], self.list[1] = self.list[1], self.list[self.size]
    self.size -= 1

    retval = self.list.pop()

    self.percolateDown(1)

    return retval

  def percolateDown(self, i):
    while 2 * i <= self.size:
      minChild = self.minChild(i)

      if self.list[i] > self.list[minChild]:
        self.list[i], self.list[minChild] = self.list[minChild], self.list[i]

      i = minChild

  def minChild(self, i):
    if 2 * i + 1 > self.size:
      return 2 * i
    elif self.list[2 * i] < self.list[2 * i + 1]:
      return 2 * i
    else:
      return 2 * i + 1


# h = Heap()
# h.insert(20)
# h.insert(7)
# h.insert(8)
# h.insert(11)
# h.insert(9)
# h.insert(5)
# h.insert(12)
# h.insert(1)
# h.insert(24)
# h.insert(3)
# h.insert(4)
# h.insert(111)
# h.insert(6)
# h.insert(2)
#
# print(h.list)
#
# print(h.delMin())
# print(h.delMin())
# print(h.delMin())
# print(h.delMin())
# print(h.delMin())
# print(h.delMin())
# print(h.delMin())
