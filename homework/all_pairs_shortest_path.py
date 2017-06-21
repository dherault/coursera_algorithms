#!/usr/bin/python3.5
import math
# Floyd-Warshall Algorithm
def floyd_warshall(n, edges):
    range_n = range(n)
    A = [[[0 if i == j else math.inf] for j in range_n] for i in range_n]

    for (s, t, c) in edges:
        A[s - 1][t - 1][0] = c

    print(A)

if __name__ == '__main__':
    input_file_name = 'all_pairs_shortest_path_input_1.txt'

    with open(input_file_name) as f:
        n = int(f.readline().split()[0])
        input = [tuple(map(int, x.split())) for x in f.readlines()]


    floyd_warshall(n, input)
