#!/usr/bin/python3.5
# I stole this code, it taught me a lesson about stealing code, I'll try not to do that anymore to improve my game
# https://github.com/dbarabanov/coursera/blob/master/algorithms_2/assignment_2/question_2.py

def invert_bit_string(bit_string):
  if bit_string == '0':
    return '1'

  return '0'

def get_cost_one_and_two_vertex(v):
  vertices = []

  for i in range(len(v)):
    vertices.append(v[:i] + invert_bit_string(v[i]) + v[i + 1:])

    for j in range(i + 1, len(v)):
      vertices.append(v[:i] + invert_bit_string(v[i]) + v[i + 1:j] + invert_bit_string(v[j]) + v[j + 1:])

  return vertices

def clustering_big():
  input_file_name = 'clustering_big_input.txt'

  vertices = []
  leaders = {}

  def find_leader(vertex):
    leader = leaders[vertex]

    while leader != leaders[leader]:
      leader = leaders[leader]

    return leader

  with open(input_file_name) as f:
    f.readline() # pop first line

    for line in f:
      vertex = "".join(line.rstrip().split(' '))
      vertices.append(vertex)
      leaders[vertex] = vertex

  num_clusters = len(leaders)

  # Merge 1 and 2-cost edges
  for vertex in vertices:
    # Find true leader
    leader = find_leader(vertex)

    for similar_vertex in get_cost_one_and_two_vertex(vertex):
      if leaders.get(similar_vertex):

        similar_leader = find_leader(similar_vertex)

        if leader != similar_leader:
          leaders[similar_leader] = leader
          num_clusters -= 1

  print('Minimun number of clusters: ' + str(num_clusters))



if __name__ == '__main__':
  clustering_big()
