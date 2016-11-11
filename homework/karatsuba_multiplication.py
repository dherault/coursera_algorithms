#!/usr/bin/python
import sys


def len_int(num):
  return len(str(num))


def split_at(num, index):
  return (num / 10 ** index, num % 10 ** index)


def karatsuba(x, y):
  if x < 10 or y < 10: return x * y

  n_by_2 = max(len_int(x), len_int(y)) / 2

  (a, b) = split_at(x, n_by_2)
  (c, d) = split_at(y, n_by_2)

  ac = karatsuba(a, c)
  bd = karatsuba(b, d)

  return ac * 10 ** (2 * n_by_2) + (karatsuba(a + b, c + d) - ac - bd) * 10 ** n_by_2 + bd


def main():
  args = sys.argv[1:]

  if not args or len(args) < 2:
    print 'Bad arguments'
    return

  x = int(args[0])
  y = int(args[1])

  print 'Multiplying %d and %d\n' % (x, y)

  print 'Classic multiplication:'
  print x * y
  print 'Karatsuba:'
  print karatsuba(x, y)

if __name__ == '__main__': main()
