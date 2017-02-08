#!/usr/bin/python3.5
# import time
import math

class DataHeap:
  """Min-Heap with additionnal data field"""
  def __init__(self):
    self.list = [(0, None)]
    self.size = 0

  def insert(self, value, data = None):
    self.list.append((value, data))
    self.size += 1
    self.percolateUp(self.size)

  def percolateUp(self, i):
    half_i = i // 2

    while half_i:

      if self.list[i][0] < self.list[half_i][0]:
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

      if self.list[i][0] > self.list[minChildIndex][0]:
        self.list[i], self.list[minChildIndex] = self.list[minChildIndex], self.list[i]

      i = minChildIndex

  def minChildIndex(self, i):
    if 2 * i + 1 > self.size:
      return 2 * i
    elif self.list[2 * i][0] < self.list[2 * i + 1][0]:
      return 2 * i
    else:
      return 2 * i + 1

  def extractByData(self, data):
    for index, x in enumerate(self.list): # O(n) ...
      if x[1] == data:
        self.list[index] = (0, data)

        self.percolateUp(index)
        self.extractMin()

        return x

    return None

def dijkstra_shortest_path():
  # input_file_name = 'dijkstra_shortest_path_test_result_is_0_1_2_3_4_4_3_2.txt'
  # input_file_name = 'dijkstra_shortest_path_test_result_is_0_3_5_8_5_7_11_4_6_10_10.txt'
  input_file_name = 'dijkstra_shortest_path_input.txt'

  vertices = []
  edges = []

  with open(input_file_name) as f:
    for line in f:
      line_data = line.split('\n')[0].split('\t')

      vertex = int(line_data.pop(0))

      vertices.append(vertex)

      for edge_data in line_data:
        other_vertex, edge_length = edge_data.split(',')
        edges.append((vertex, int(other_vertex), int(edge_length)))

  # print(len(vertices))
  # print(len(edges))

  n_vertices = len(vertices)
  shortest_paths = {}

  # 1 is the start vertex
  start_vertex = vertices.pop(0)
  shortest_paths[start_vertex] = 0

  heap = DataHeap()

  for v in vertices:
    heap.insert(math.inf, v)

  def maintain_second_invariant(removed_vertex):
    # print('maintain_second_invariant', removed_vertex)
    # print(heap.size)

    shortest_path = shortest_paths[removed_vertex]

    for edge in edges:
      if edge[0] == removed_vertex:
        extracted_data = heap.extractByData(edge[1]) # O(n)...

        # print('extracted_data', extracted_data, shortest_path + edge[2])
        # print(heap.list)
        if extracted_data:
          heap.insert(min(extracted_data[0], shortest_path + edge[2]), extracted_data[1])


  maintain_second_invariant(start_vertex)

  # print(heap.list)
  # print()
  # print(heap.extractMin())

  while heap.size > 0:
    score, vertex = heap.extractMin()

    # print(vertex)

    shortest_paths[vertex] = score

    maintain_second_invariant(vertex)

  # print(shortest_paths)
  # print('')
  print(','.join([str(shortest_paths[x]) for x in [7,37,59,82,99,115,133,165,188,197]]))
  # print('')
  # print(shortest_paths[7])
  # print(shortest_paths[37])
  # print(shortest_paths[59])
  # print(shortest_paths[82])
  # print(shortest_paths[99])
  # print(shortest_paths[115])
  # print(shortest_paths[133])
  # print(shortest_paths[165])
  # print(shortest_paths[188])
  # print(shortest_paths[197])



  #while foobar

  # (score, nearest_vertex) = heap.extractMin()
  #
  # if shortest_paths.has_key(nearest_vertex):
  #   shortest_paths[nearest_vertex] +=
  # start_vertex = vertices.pop(0)
  # shortest_paths[start_vertex] = 0
  #
  #
  # outgoing_edges_from_start_vertex = [e for e in edges if e[0] == start_vertex]
  #
  # print(outgoing_edges_from_start_vertex)
  #
  # heap = DataHeap()
  #
  # for v in vertices:
  #   existing_edge_length_from_start_vertex = [e[2] for e in outgoing_edges_from_start_vertex if e[1] == v]
  #
  #   if len(existing_edge_length_from_start_vertex):
  #     heap.insert(min(existing_edge_length_from_start_vertex), v)
  #   else:
  #     heap.insert(math.inf, v)

if __name__ == '__main__': dijkstra_shortest_path()
