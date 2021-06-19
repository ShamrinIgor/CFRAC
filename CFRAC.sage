debugMode = false

def CFRAC(N):
  # Шаг 1...
  lengthN = len(str(N))
  fb = 200
  ub = 14400
  fBase = []
  count = 1
  while len(fBase) < fb:
    fBase = [-1] + list(filter(lambda p: kronecker(N,p) == 1, primes_first_n(fb*2*count)))
    count += 1
  fBase = fBase[0:fb]
  if debugMode:
    print(fBase)

  # Чекаем k 😳
  print("Чекаем k 😳: ",getK(N, fBase))

  # Шаг 2...
  cf = continued_fraction(sqrt(N))
  ui = []
  i = 0
  F = IntegerModRing(N)
  while len(ui) < fb + 2:
    numerator = F(cf.numerator(i)).lift()
    q = F((-1)^(i+1)*numerator^2).lift()*(-1)^(i+1)
    #minRoot = min(F(numerator), n - F(numerator))
    #if abs(numerator-q) >= ub:
    #  break
    if debugMode:
      print(numerator, q)

    parFact = product(list(map(lambda x: x[0]^x[1],list(zip(fBase, BFac(q, fBase))))))
    if q == parFact or isLpUB(q//parFact, ub, fBase[-1]):
    # if q == parFact:
      ui.append((numerator, q))
      print(numerical_approx(len(ui)/(fb + 1))*100)
    i += 1
  if debugMode:
    print(ui)

  # Шаг 3...
  powers = list(map(lambda x: BFac(x[1], fBase), ui))
  gf2 = IntegerModRing(2)
  A = matrix(gf2, powers)
  if debugMode:
    print(A)
  solIndex = 0
  for sol in A.left_kernel().basis():
    qf = product(list(map(lambda x: x[1][1],list(filter(lambda x: x[0] == 1,list(zip(sol, ui)))))))
    af = product(list(map(lambda x: x[1][0],list(filter(lambda x: x[0] == 1,list(zip(sol, ui)))))))
    x = F(af)
    y = F(sqrt(qf))
    if x != y and x != F(-y):
      factor = gcd(x-y, N)
      #print(n, factor)
      #print(n/factor)
      #print('n = ', factor, ' * ', n/factor)
      print('factor: ', factor)
      return
  print('Cringe')

def BFac(N, fBase):
  temp = N
  pwsVec = [0 for i in fBase]
  if temp < 0:
    pwsVec[0] = 1
  num = 1
  for f in fBase[1:]:
    i = 0
    while temp % f == 0 and temp != 0:
      temp = temp // f
      i += 1
    pwsVec[num] = i
    temp = N
    num += 1
  return pwsVec


def isLpUB(x, UB, largestPrimeInFB):
  # print(x)
  if not is_prime_power(x):
    return False
  base = next_prime(largestPrimeInFB)
  while not x.is_power_of(base) and base <= UB:
    base = next_prime(base)
  return x.is_power_of(base)


def getK(N, fBase):
# 🚨🚨🚨 Норм тут кстатить в инт?
  return int(max(list(map(lambda k: knuthFunc(k, N, fBase[1:]), list(range(1,100))))))

def knuthFunc(k, N, primes):
  # Наверное тут должны быть числа из fBase
  return sum(list(map(lambda x: f(x, k, N)*log(x), primes))) - log(k)/2

def f(p, k, N):
  if p%2 != 0:
    if p%k == 0:
      # print("Cringe 1")
      return 1/(p + 1)
    else:
      # print("Cringe 2")
      return 2*p/(p^2 - 1)
  else:
    if k%2 == 0 or (k*N)%4 == 3:
      # print("Cringe 3")
      return 1/3
    if (k*N)%8 == 5:
      # print("Cringe 4")
      return 2/3
    if (k*N)%8 == 1:
      # print("Cringe 5")
      return 4/3
    # print("CRINGE", p, k, N)
