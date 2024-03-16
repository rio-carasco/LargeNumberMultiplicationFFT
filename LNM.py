#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 11:16:02 2021

@author: riocarasco
"""
import numpy as np
import py_proto_graphec as opt
#from scipy.fft import fft, ifft
from random import randint

def print_info():
    #print("RMSE: {}".format(err))
    print("No. frames so far: {}".format(opt.get_n_frames(0)))
    print("Runtime so far: {}s".format(opt.get_runtime_s(0)))
    print("Energy cost so far: {}J".format(opt.get_energy_cost_j(0)))

# Using OPT SIMULATOR 2; 1x4 in-plane optical FT core @ 1 GHz (Scalable) #######

def iftoptosim(ftsize, arrayentity):
    print()
    ft_size = ftsize
    target_accuracy = 64
    print("Defining the simulator...")
    opt.define(target_accuracy, ft_size)
    print("Number of simulators: {}".format(opt.n_simulators()))
    print("Accuracy : {}".format(opt.get_accuracy(0)))
    print("FT size: {}".format(opt.get_ft_size(0)))
    ift_out = opt.ift(0, arrayentity)

    return ift_out

def ftoptosim(ftsize, arrayentity):
    print()
    ft_size = ftsize
    target_accuracy = 64
    print("Defining the simulator...")
    opt.define(target_accuracy, ft_size)
    print("Number of simulators: {}".format(opt.n_simulators()))
    print("Accuracy : {}".format(opt.get_accuracy(0)))
    print("FT size: {}".format(opt.get_ft_size(0)))
    ft_out = opt.ft(0, arrayentity)

    return ft_out

##################################################################################

def randomDigits(n):                          
    range_start = 10**(n-1)                   
    range_end = (10**n)-1
    return randint(range_start, range_end)

def Pwr2(Number):
    """
    Check if a number is a power of 2.
    Parameters(Number): The number to be checked.
    Returns: True if the number is a power of 2, False otherwise.
    """
    if np.log2(Number).is_integer():
        return True
    else:
        return False


def Info(A, B):
    """
    Defines variables representing infomation about original numbers
    Parameters(A, B): The two initial numbers
    Returns: 
        - M: Length of array B
        - L: Size used for padding array A before convolution
        - N: Size of the output array after convolution
        - OL: Length of the output array after convolution
    """
    M = len(B)
    L = 1
    while not Pwr2(L+M-1):                     
        L+=1
    OL = len(A) + len(B) - 1                   
    N = L + M - 1
    return M, L, N, OL

def Num2Arr(Number): 
    """
    Converts a Number in string format to an array in which each digit is an element
    Parameters(Number): Any int
    Returns(NumArr): Array
    """
    NumArr = []
    for i in str(Number):
        NumArr.append(int(i))
    NumArr = np.array(NumArr)
    return NumArr
    

def Split(Arr, L): 
    """
    Divides an array into segments of a specified length (L), ensuring that each segment is complete and contains no more than L elements.
    If the length of the input array is not evenly divisible by L, the function ensures that any remaining elements are included in the 
    last segment.
    Parameters(Arr, L): (Array(Initial number), max length of arrays
    Returns(FinalArr)
    """
    #If the desired segment length (L) is 1, the function simply splits the input array into individual elements and returns them as 
    #separate arrays.
    if L == 1:                                    
        return np.array_split(Arr, len(Arr))
    #segment length is greater than 1 
    else:
        #The function calculates the remainder (R) when the length of the input array is divided by L
        R = len(Arr)%L 
        #Determines the length of the input array without considering the remainder (NoR)
        NoR = len(Arr) - R 
        #Creates a new array (NewArr) containing only the elements of the input array up to the nearest multiple of L
        NewArr = Arr[:NoR] 
        #If the length of NewArr is equal to L, indicating that it is a complete segment, the function checks if there is a remainder (R). 
        #If there is a remainder, it creates an additional array (RArr) containing the remaining elements.
        if len(NewArr) == L:                    
            FinalArr = np.empty(2, object)          
            FinalArr[0] = NewArr                    
            if R == 0:                              
                pass
            else:
                RArr = np.array(Arr[-R:])           
                FinalArr[1] = RArr
            return FinalArr   
        #If NewArr is longer than L, indicating that it contains multiple complete segments, 
        #the function splits NewArr into segments of length L using np.array_split
        NewArrSplit = np.array_split(NewArr, len(NewArr)/L)
        if R == 0:
            pass
        else:
            RArr = np.array(Arr[-R:])
        if R == 0:
            R_Value = 0                                                  
        else:
            R_Value = 1
        #If there is a remainder (R), it appends the remaining elements to the last segment
        FinalArr = np.empty(len(NewArrSplit) + R_Value, object)
        for i in range(len(NewArrSplit)):
            FinalArr[i] = NewArrSplit[i]
        if R == 0:
            pass
        else:
            FinalArr[-1] = RArr
        #the function returns an array containing the split segments, with any remaining elements placed in the last segment
        return FinalArr

def Padding(Arr, N):
    """
    Adds a padding of zeros to the arrays to allow the arrays to fit a desired length N.
    Parameters(Arr, N): (Array, Desired length of Array)
    Returns: Array with padding
    """
    return np.pad(Arr, (0, N-len(Arr)), 'constant')          

############################################################################################################################################

#A and B are our two 'Large' numbers. Use randomDigits() to produce a random number with a specified number of digits
A = randomDigits(5000)          
B = randomDigits(5000)

#Convert to Array
A = Num2Arr(A)                 
B = Num2Arr(B)

M, L, N, OL = Info(A, B)       

ASplit = Split(A, L)           

OAList = np.empty((len(ASplit), N)) 
for i in range(len(ASplit)):        
    OAList[i] = Padding(ASplit[i], N) 



h = Padding(B, N)    
h_complex = h.astype(complex)

FOAList = []
for i in OAList:
    i_complex = i.astype(complex)
    ft_i = ftoptosim(len(i), i_complex)                      
    ft_h = ftoptosim(len(h), h_complex)                      
    ifft_prod = iftoptosim(len(h), (ft_i * ft_h))              
    FOAList.append(ifft_prod)                       
    #FOAList.append(ifft(fft(i) * fft(h)))
FOAList = np.array(FOAList)                          
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
print(Result)

# Runtime and energy cost:

print_info()

