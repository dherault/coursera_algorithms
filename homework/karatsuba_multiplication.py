#!/usr/bin/python
import sys


def len_int(num):
  return len(str(num))


def split_at(num, index):
  string = str(num)

  return (int(string[0:index]), int(string[index:]))

# Not optimized
def karatsuba(x, y):
  print 'karatsuba x: %d _ y: %d' % (x, y)

  if x < 10 or y < 10:
    result = x * y
    print 'low numbers, result: %d' % result
    return result

  n = max(len_int(x), len_int(y))
  n_by_2 = int(n / 2)

  (a, b) = split_at(x, n_by_2)
  (c, d) = split_at(y, n_by_2)

  print 'a: %d _ b: %d _ c: %d _ d: %d' % (a, b, c, d)

  ac = karatsuba(a, c)
  bd = karatsuba(b, d)
  ad_plus_bc = karatsuba(a + b, c + d) - ac - bd

  result = ac * (10 ^ n) + ad_plus_bc * (10 ^ n_by_2) + bd

  print 'result: ' + str(result)

  return result


def main():
  args = sys.argv[1:]

  if not args or len(args) < 2:
    print 'Bad arguments'
    return

  x = int(float(args[0]))
  y = int(float(args[1]))

  print 'Multiplying ' + str(x) + ' with ' + str(y)

  print karatsuba(x, y)

if __name__ == '__main__': main()
