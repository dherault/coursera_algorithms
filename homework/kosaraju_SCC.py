#!/usr/bin/python3.5
# import time

def kosaraju_strongly_connected_components():
  # input_file_name = 'kosaraju_SCC_test_result_is_3_3_3_0_0.txt'
  # input_file_name = 'kosaraju_SCC_test_result_is_6_3_2_1_0.txt'
  input_file_name = 'kosaraju_SCC_input.txt'

  print('Creating adjacency lists')

  adjacency_list = []
  reversed_adjacency_list = []

  n = 0
  rn = 0

  with open(input_file_name) as f:
    for line in f:

      uv = line[:-1].strip().split(' ')
      u = int(uv[0]) - 1 # - 1 is problem-specific (teacher's vertices id start from one)
      v = int(uv[1]) - 1 # also, we assume that if vertices 1 and 3 exist, vertex 2 exists too
      max_uv = max(u, v) # even if no edge refers to vertex 2

      while n <= max_uv:
        adjacency_list.append([])
        n += 1

      while rn <= max_uv:
        reversed_adjacency_list.append([])
        rn += 1

      adjacency_list[u].append(v)
      reversed_adjacency_list[v].append(u)

  # print(adjacency_list)
  # print(reversed_adjacency_list)

  print('Adjacency lists created:', n, 'vertices')
  print('Running first depth-first search')

  stack = [] # A stack is used to perform DFS without recursion and avoid stack overflow
  orderings = []
  explored = [False] * n

  for i in range(n - 1, -1, -1):
    if not explored[i]:
      explored[i] = True
      stack.append(i)

      # print('considering', i)
      # print('starting explored_vertices:', explored_vertices)

      while len(stack):
        current_vertex = stack[-1]
        all_neighbours_explored = True

        for j in reversed_adjacency_list[current_vertex]:
          if not explored[j]:
            all_neighbours_explored = False
            explored[j] = True
            stack.append(j)
            break

        if all_neighbours_explored:
          stack.pop()
          orderings.append(current_vertex)

          # print('All neighbours of', current_vertex, 'explored. order:', t)

        # print('explored_vertices:', explored_vertices)
        # print('stack:', stack)
        # print('stack size:', len(stack), 't:', len(orderings))

  # print('orderings:', orderings)
  print('Running second depth-first search')

  orderings.reverse()
  explored = [False] * n
  leaders_dict = {}

  for i in orderings:
    if not explored[i]:
      explored[i] = True
      stack.append(i)
      # We use a set to prevent member testing: (if current_vertex not in leaders_dict[i + 1])
      # Idk what is better: python noob
      # We add + 1 to retrieve the correct vertex id
      leaders_dict[i + 1] = set()

      while len(stack):
        current_vertex = stack[-1]
        all_neighbours_explored = True

        leaders_dict[i + 1].add(current_vertex + 1)

        for j in adjacency_list[current_vertex]:
          if not explored[j]:
            all_neighbours_explored = False
            explored[j] = True
            stack.append(j)
            break

        if all_neighbours_explored:
          stack.pop()

  # print(leaders_dict)
  # print(len(leaders_dict))

  # A conventional implementation would return here
  # return leaders_dict

  print('Computing assigment result')

  results = []
  items = leaders_dict.items()

  while len(results) < 5:
    best_leader = None
    best_n_followers = 0

    for (leaders, followers) in items:
      n_followers = len(followers)

      if n_followers > best_n_followers:
        best_n_followers = n_followers
        best_leader = leaders

    results.append(str(best_n_followers))

    if best_leader is not None:
      del leaders_dict[best_leader]

  print('Final result:')
  print(','.join(results))

if __name__ == '__main__': kosaraju_strongly_connected_components()
