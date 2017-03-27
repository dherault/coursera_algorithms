#!/usr/bin/python3.5

def clustering():
  # input_file_name = 'clustering_test_k_is_2_result_is_100.txt'
  input_file_name = 'clustering_input.txt'

  k = 4
  edges = []
  vertex_to_cluster = {}
  cluster_to_vertices = {}

  with open(input_file_name) as f:
    for line in f:
      edge = [int(x) for x in line.rstrip().split(' ')]

      if len(edge) == 3:
        edges.append(edge)
        cluster_to_vertices[edge[0]] = [edge[0]]
        cluster_to_vertices[edge[1]] = [edge[1]]
        vertex_to_cluster[edge[0]] = edge[0]
        vertex_to_cluster[edge[1]] = edge[1]

  edges.sort(key=lambda e: e[2])
  # print(edges)
  # print(vertex_to_cluster)
  # print(cluster_to_vertices)

  while True:
    a = None
    b = None

    for index, edge in enumerate(edges):
      if vertex_to_cluster[edge[0]] != vertex_to_cluster[edge[1]]:
        a = edge[0]
        b = edge[1]
        break

    cluster_of_a = vertex_to_cluster[a]
    cluster_of_b = vertex_to_cluster[b]

    for vertex in cluster_to_vertices[cluster_of_b]:
      vertex_to_cluster[vertex] = cluster_of_a
      cluster_to_vertices[cluster_of_a].append(vertex)

    del cluster_to_vertices[cluster_of_b]

    if len(cluster_to_vertices) == k:
      break

  max_spacing = None

  for edge in edges:
    if vertex_to_cluster[edge[0]] != vertex_to_cluster[edge[1]]:
      max_spacing = edge[2]
      break

  print('Max spacing with ' + str(k) + ' clusters: ' + str(max_spacing))

if __name__ == '__main__':
  clustering()
