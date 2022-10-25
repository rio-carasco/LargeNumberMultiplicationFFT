import numpy as np
from scipy.fft import fft, ifft
from random import randint

def randomDigits(n):                            # Generates random digits of specific and equal length (e.g between 10^9 = 1000000000 and (10^10 - 1) = 9999999999)
    range_start = 10**(n-1)
    range_end = (10**n)-1
    return randint(range_start, range_end)

def Pwr2(Number):                               # Checks if a number is a power of 2; if not, log2 will return a non-integer value
    if np.log2(Number).is_integer():
        return True
    else:
        return False


def Info(A, B):                                 # Finds a value of L such that L+M-1 is a power of 2, precursor for acyclic -> cyclic/negacyclic convolution method
    M = len(B)
    L = 1
    while not Pwr2(L+M-1):
        L+=1
    OL = len(A) + len(B) - 1 
    N = L + M - 1
    return M, L, N, OL

def Num2Arr(Number):                            # Converts a single number to a string and then separates the string into individual integer values placed in an array 
    a = []
    for i in str(Number):
        a.append(int(i))
    a = np.array(a)
    return a
    

def Split(Arr, L):                              # Iteratively slices an array of polynomial values into sub-arrays of specific length 
    if L == 1:                                  # Takes in L from a prior instance of Info() as a testing metric (Galois group for cyclotomic extension: calculating a splitting field)
        return np.array_split(Arr, len(Arr))
    else:
        R = len(Arr)%L                          # Length of array modulo L; tells us how "long" the array is relative to a pertinent value (relates to calculations in power 2, given the earlier role of L)
        NoR = len(Arr) - R                      # Subtracts the excess length of the array, returning an "ideal" array length for splitting to power 2
        NewArr = Arr[:NoR]                      # Selects the region of the array up to the length of NoR
        if len(NewArr) == L:                    # Check to see if the length of the new array is equivalent to L; if so:
            FinalArr = np.empty(2, object)      # np.empty returns a new array of a given shape without filling with zeros. Passing 2 to np.empty presumably initiates an empty array of size 2
            FinalArr[0] = NewArr                # Stores the cut-down region of a passed array in the first entry of the new empty array 
            if R == 0:                          # Check to see if there is a remainder left over (i.e is the length of Arr identified with 0 in modulus L)
                pass
            else:
                RArr = np.array(Arr[-R:])       # If not, stores the cut-off section of the original array in the second entry of FinalArr
                FinalArr[1] = RArr
            return FinalArr          
        NewArrSplit = np.array_split(NewArr, len(NewArr)/L)
        if R == 0:
            pass
        else:
            RArr = np.array(Arr[-R:])
        if R == 0:
            R_Value = 0                          # This function splits an array into sub-arrays of equal power-of-two lengths
        else:
            R_Value = 1
        FinalArr = np.empty(len(NewArrSplit) + R_Value, object)
        for i in range(len(NewArrSplit)):
            FinalArr[i] = NewArrSplit[i]
        if R == 0:
            pass
        else:
            FinalArr[-1] = RArr
        return FinalArr


def Padding(Arr, N):                              # A pad to ensure power-2 calculations for efficient FT
    return np.pad(Arr, (0, N-len(Arr)), 'constant')

A = randomDigits(500)
B = randomDigits(500)

print("First Number {}".format(A))
print("Second Number {}".format(B))

A = Num2Arr(A)                                      # Converts these numbers to individual elements in an array  
B = Num2Arr(B)

M, L, N, OL = Info(A, B)                            # Extracts pertinent info from A,B for decision making on rearrangement/efficient calculation

ASplit = Split(A, L)                                # Split A into tractable sub-arrays of power of two length

OAList = np.empty((len(ASplit), N))                 # Empty array the length of A, each sub-array N long
for i in range(len(ASplit)):                        # Enters array of N padding zeros into each entry of A that isn't filled; ensures parity in length with B? 
    OAList[i] = Padding(ASplit[i], N) 



h = Padding(B, N)                                   # Pads B to power of two length defined by N


FOAList = []
for i in OAList:
    FOAList.append(ifft(fft(i) * fft(h)))           # Multiply FT of a sub-array from A and FT of padded B array 
FOAList = np.array(FOAList)                         # Classic convolution theorem; Schonhage-Strassen. Multiplying each FTd sub-array with the FT of B and then taking the inversion (bit-shifting)
FOAList = np.reshape(FOAList, (len(OAList), N))
FOAList = np.round(np.real(FOAList))

k = round(len(OAList))
height = k
width = (((k - 1) * L) + N)
m = np.zeros((height, width))



for j in range(np.shape(FOAList)[0]):
    m[j,j*L:N+(j*L)] = FOAList[j]

Final = []
for i in range(np.shape(m)[1]):
    Final.append(np.sum(m[:,i]))
Final = np.array(Final)    
Final = Final[:OL]
Final = np.flip(Final, 0).astype(int)


digitArray = []
for i in Final:
    digits = [int(x) for x in str(i)]
    digitArray.append(digits)




dA = len(digitArray)

Len = 0
for i in range(dA):
    Len_i = len(digitArray[i]) + i
    if Len < Len_i:
        Len = Len_i
    else:
        pass


FinalSum = np.zeros((dA, Len))


for i in range(dA):
    FinalSum[dA-i-1, (Len-1-i)-(len(digitArray[i]))+1 : Len-i] = digitArray[i]

Carry = 0
Result = []
for i in range(np.shape(FinalSum)[1]):
    ResultCol = sum(FinalSum[:,-(i+1)]) + Carry
    if ResultCol >= 10:
        Sum = [int(x) for x in str(int(ResultCol))]
        Carry = Sum[0]
        Result.append(int(Sum[1]))
    else:
        Result.append(int(ResultCol))
        Carry = 0
                                   
Result = np.flipud(Result)
Result = ''.join(str(i) for i in Result)
print("Result {}".format(Result))


