import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt

def in_pairs(pair):
  #dictionary to check if pairs are complementary
  pairs = {'A': 'U', 'U':'A', 'G':'C', 'C':'G'}
  if pair in pairs.items():
    #true
    return 1
  #false
  return 0

def build_matrix(seq):
  m = len(seq)
  D = np.empty((m, m))
  D.fill(np.NAN)
  for i in range(m):
   for j in range(i+1):
    D[i][j] = 0 
  return D


def nussinov(seq):
  D = build_matrix(seq)
  m = len(seq)
  for h in range(1, m):
    for i in range(m - h):
      j = h + i
      if i < j:
        diag = D[i + 1][j - 1] + in_pairs((seq[i], seq[j]))
        down = D[i + 1][j]
        left = D[i][j - 1]
        skip = -np.inf
        for k in range(i + 1, j):
          if (D[i][k] + D[k+1][j] > skip):
            skip = D[i][k] + D[k+1][j]
        D[i][j] = max(diag, down, left, skip)
      else:
        D[i][j] = 0
  temp_record = []
  record = traceback(D, 0, m-1, temp_record)
  ans = reading_record(seq, record)
  return ans

def traceback(D, i, j, record):
  if i >= j:
    return record
  elif D[i+1][j] == D[i][j]:
    traceback(D, i+1, j, record)
  elif D[i][j-1] == D[i][j]:
    traceback(D, i, j-1, record)
  elif D[i+1][j-1] + 1 == D[i][j]:
    record.append((i, j))
    traceback(D, i+1, j-1, record)
  else:
    for k in range(i+1, j):
      if D[i][k] + D[k+1][j] == D[i][j]:
        traceback(D, i, k, record)
        traceback(D, k+1, j, record)
        break
  return record

def reading_record(seq, record):
  output = ["-"]*len(seq)
  for (i, j) in record:
    if i < j:
      output[i] = '('
      output[j] = ')'
    else:
      output[i] = ')'
      output[j] = '('
  return "".join(output)


if __name__ == "__main__":

    path = 'C:\\Users\\19782\\Downloads\\bpRNA_test\\'
    files = os.listdir(path)
    x = []
    y = []
    a = []
    b = []

    for file in files:
        if os.path.isfile(os.path.join(path, file)):
            f = open(os.path.join(path, file), 'r')
            lines = f.readlines()
            rna = lines[3]
            ans = nussinov(rna)
            c = 0
            c = ans.count('(')
            counter = 0
            counter2 = 0
            counter = lines[4].count('(') 
            counter2 = lines[4].count('[')
            if (counter2 > 0):
                a.append(counter + counter2)
                b.append(c)
            else:
                x.append(counter)
                y.append(c)
    np.savetxt('x_data.dat', x)
    np.savetxt('y_data.dat', y)
    np.savetxt('a_data.dat', a)
    np.savetxt('b_data.dat', b)    
    plt.scatter(x, y)
    plt.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)))
    plt.xlim([0, 250])
    plt.ylim([0, 250])
    plt.show()
    plt.savefig('nussinov_scatter_xy.png')
    plt.scatter(a, b)
    plt.plot(np.unique(a), np.poly1d(np.polyfit(a, b, 1))(np.unique(a)))
    plt.xlim([0, 250])
    plt.ylim([0, 250])
    plt.show()
    plt.savefig('nussinov_scatter_ab.png')   

