import numpy as np


def EvenOdd(List):
    ListShape = np.reshape(List, (int(len(List)/2), 2))
    Even = ListShape[:, 0]
    Odd = ListShape[:, 1]
    return Even, Odd

def FFT(P):
    n = int(len(P))
    if n == 1:
        return P
    w = np.cos(2*np.pi/n) + np.sin(2*np.pi/n)*1j
    P_e, P_o = EvenOdd(P)
    y_e, y_o = FFT(P_e), FFT(P_o)
    y = np.zeros(n, dtype = 'complex_')
    for i in range(int(n/2)):
        y[i] = y_e[i] + w**i*y_o[i]
        y[int(i + n/2)] = y_e[i] - w**i*y_o[i]
    return y

    
X  = np.array([5, 3, 2, 1])
FFT = FFT(X)
print(FFT)
