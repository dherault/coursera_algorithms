#!/usr/bin/python3.5

def pick_cost(edge):
  return edge[2]

def prim_MST():
  # input_file_name = 'prim_MST_test_result_is_7.txt'
  input_file_name = 'prim_MST_input.txt'

  T = []
  E = []
  V = set()

  with open(input_file_name) as f:
    for line in f:
      edge = [int(x) for x in line.rstrip().split(' ')]

      if len(edge) == 3:
        V.add(edge[0])
        V.add(edge[1])
        E.append(edge)

  X = { V.pop() }

  while len(V):
    candidate_edges = []

    for edge in E:
      if (edge[0] in X and edge[1] in V) or (edge[1] in X and edge[0] in V):
        candidate_edges.append(edge)

    best_edge = min(candidate_edges, key=pick_cost)

    T.append(best_edge)

    if best_edge[0] in X:
      V.remove(best_edge[1])
      X.add(best_edge[1])
    else:
      V.remove(best_edge[0])
      X.add(best_edge[0])


  answer = sum(pick_cost(edge) for edge in T)

  print('Done, ' + str(answer))

if __name__ == '__main__':
  prim_MST()
