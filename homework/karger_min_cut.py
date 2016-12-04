#!/usr/bin/python3.5
# import sys
import random
# import copy
import math

def build_graph(input):
  vertices = []
  edges = []

  for line in input.split('\n'):
    if not line: continue

    row = line.split()
    vertex = row.pop(0)

    vertices.append(vertex)

    for other_vertex in row:
      if other_vertex in vertices: continue
      edges.append([vertex, other_vertex])

  # print(vertices)
  # print(edges)

  return vertices, edges

def min_cut(g):
  vertices, edges = g

  # save original edges
  # for edge in edges:
  #   edge.append([edge[0], edge[1]])

  # or not, to speed things up

  while len(vertices) > 2:
    # pick one edge uniformly at random and remove it
    popped_edge = edges.pop(random.randrange(len(edges)))

    # u, v, save = popped_edge
    u, v = popped_edge

    # print(v, 'is contracted into', u)

    # remove the corresponding vertex
    vertices.pop(vertices.index(v))

    # update the edges
    for edge in edges:
      if edge[0] == v: edge[0] = u
      if edge[1] == v: edge[1] = u

    # remove loops
    edges[:] = [edge for edge in edges if edge[0] != edge[1]]

  # return original min-cut
  # return [edge[2] for edge in edges]

  # or not, to speed things up
  return edges

def loop_min_cut(g):
  best = None
  n_best = len(g[1])

  n = len(g[0])

  print('g:', n, n_best)
  # iterations = int(n * n * math.log(n))
  iterations = n * n

  print('will perform', iterations, 'iterations')

  for i in range(iterations):
    print(100 * i / iterations, '%')

    g_copy = (g[0][:], [edge[:] for edge in g[1]])

    # result = min_cut(copy.deepcopy(g))
    result = min_cut(g_copy)
    n_result = len(result)

    if n_result < n_best:
      best = result
      n_best = n_result

  return best

def main():
  # input = sys.argv[1:]

  # with open('karger_min_cut_test_result_is_3.txt') as f:
  with open('karger_min_cut_input.txt') as f:
    g = build_graph(f.read())

    print(len(loop_min_cut(g)))




if __name__ == '__main__': main()
