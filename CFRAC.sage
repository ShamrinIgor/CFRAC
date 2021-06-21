import time

debugMode = false

def CFRAC(N, showProgress = False):
  # Шаг 1...
  FBandUB = getFBandUB(N)
  fb = FBandUB[0]
  ub = FBandUB[1]
  fBase = getFBase(N, fb)

  if debugMode:
    print(fBase)

  k = getK(N, fBase)

  kN = k*N
  fBase = getFBase(kN, fb)

  # Шаг 2...
  cf = continued_fraction(sqrt(kN))
  ui = []
  i = 0
  F = IntegerModRing(kN)
  while len(ui) < fb + 2:
    numerator = F(cf.numerator(i)).lift()
    q = F((-1)^(i+1)*numerator^2).lift()*(-1)^(i+1)
    if debugMode:
      print(numerator, q)

    parFact = product(list(map(lambda x: x[0]^x[1], list(zip(fBase, BFac(q, fBase))))))
    if q == parFact or isLpUB(q//parFact, ub, fBase[-1]):
      ui.append((numerator, q))
      progress = numerical_approx(len(ui)/(fb + 1))*100
      if showProgress:
        printProgressBar(len(ui), fb + 1, prefix = ' Progress:')
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

  factors = set()
  for sol in A.left_kernel().basis():
    qf = product(list(map(lambda x: x[1][1], list(filter(lambda x: x[0] == 1, list(zip(sol, ui)))))))
    af = product(list(map(lambda x: x[1][0], list(filter(lambda x: x[0] == 1, list(zip(sol, ui)))))))
    x = F(af)
    y = F(qf).nth_root(2)
    if x != y and x != F(-y):
      fac = gcd(x-y, kN).lift()
      if (1 < fac and fac < N) and (N%fac == 0):
        factors.add(fac)
  print("factors = ", factors)

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
  if not is_prime_power(x):
    return False
  base = next_prime(largestPrimeInFB)
  while not x.is_power_of(base) and base <= UB:
    base = next_prime(base)
  return x.is_power_of(base)


def getK(N, fBase):
  return max(list(map(lambda k: (k, knuthFunc(k, N, fBase[1:])), list(range(1,101,2)))), key = lambda x: x[1])[0]

def knuthFunc(k, N, primes):
  return sum(list(map(lambda x: f(x, k, N)*log(x), primes))) - log(k)/2

def f(p, k, N):
  if p%2 != 0:
    if p%k == 0:
      return 1/(p + 1)
    else:
      return 2*p/(p^2 - 1)
  else:
    if k%2 == 0 or (k*N)%4 == 3:
      return 1/3
    if (k*N)%8 == 5:
      return 2/3
    if (k*N)%8 == 1:
      return 4/3

def getFBase(N, fb):
  fBase = []
  count = 1
  while len(fBase) < fb:
    fBase = [-1] + list(filter(lambda p: kronecker(N,p) == 1, primes_first_n(fb*2*count)))
    count += 1
  return  fBase[0:fb]

def getFBandUB(N):
  numberOfDigits = len(str(N))
  if numberOfDigits <= 20:
    return (60, 3000)
  if 21 <= numberOfDigits <= 23:
    return (150, 10000)
  if 24 <= numberOfDigits <= 25:
    return (200, 14400)
  if 26 <= numberOfDigits <= 28:
    return (300, 22500)
  if 29 <= numberOfDigits <= 30:
    return (400, 29000)
  if 31 <= numberOfDigits <= 32:
    return (450, 36000)
  if 33 <= numberOfDigits <= 34:
    return (500, 36000)
  if 35 <= numberOfDigits <= 36:
    return (550, 36000)
  if 37 <= numberOfDigits <= 38:
    return (600, 44000)
  if 39 <= numberOfDigits <= 40:
    return (650, 53000)
  if 41 <= numberOfDigits:
    return (1000, 63000)

def runCFRACtests():
  numbers = [89798901815566097617, #20
             814549160699567904701, #21
             8145491579225959753183, #22
             49166729254556864996653, #23
             567145958951392379680547, #24
             4340856958597307344296589, #25
             45550029959058235964047399, #26
             132264170437550839716317897, #27
             1322641709009802963596389843, #28
             23943529316361761474523364773, #29
             291163170528478955353327561417, #30
             9592861349514945348515183836369, #31
             26998418696705720904591837114989, #32
             232474835216173169682008181183211, #33
             1391728916433545665546310683013867, #34
             15902168834583676234510442528251181, #35
             2^128 + 1] #39

  resultTimes = []
  for num in numbers:
    print("Факторизуем число: {}".format(num))
    start = time.time()
    CFRAC(num, showProgress = True)
    end = time.time()
    print("Time: ", end - start, "s.")
    resultTimes.append(end - start)
  print("Итоговый массив времени: ", resultTimes)


# Вспомогательная функция для генерации тестовых чисел
def generateTestNumbers():
  a = random_prime(1856707808261)
  b = random_prime(1856707808261)
  ab = a*b
  while len(str(ab)) != 20:
   a = random_prime(1856707808261)
   b = random_prime(1856707808261)
   ab = a*b
   print(len(str(ab)))
  print('a = ', a, 'b = ', b)
  print('a*b = ', ab)
  print('len a*b = ', len(str(ab)))

# Вывод прогресса в консоль
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, printEnd = "\r"):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    print(f'\r{prefix} {percent}% {suffix}', end = printEnd)
    if iteration == total:
        print()
