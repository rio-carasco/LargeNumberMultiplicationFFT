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
    if np.log2(Number).is_integer():
        return True
    else:
        return False


def Info(A, B):                               
    M = len(B)
    L = 1
    while not Pwr2(L+M-1):                     
        L+=1
    OL = len(A) + len(B) - 1                   
    N = L + M - 1
    return M, L, N, OL

def Num2Arr(Number):                           
    a = []
    for i in str(Number):
        a.append(int(i))
    a = np.array(a)
    return a
    

def Split(Arr, L):                               
    if L == 1:                                    
        return np.array_split(Arr, len(Arr))
    else:
        R = len(Arr)%L                          
        NoR = len(Arr) - R                      
        NewArr = Arr[:NoR]                      
        if len(NewArr) == L:                    
            FinalArr = np.empty(2, object)          
            FinalArr[0] = NewArr                    
            if R == 0:                              
                pass
            else:
                RArr = np.array(Arr[-R:])           
                FinalArr[1] = RArr
            return FinalArr                         
        NewArrSplit = np.array_split(NewArr, len(NewArr)/L)
        if R == 0:
            pass
        else:
            RArr = np.array(Arr[-R:])
        if R == 0:
            R_Value = 0                                                  
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

def Padding(Arr, N):                                         
    return np.pad(Arr, (0, N-len(Arr)), 'constant')          



# A = randomDigits(5)
# B = randomDigits(5)

# A = randomDigits(50)
# B = randomDigits(50)

# A = randomDigits(500)
# B = randomDigits(500)

# A = randomDigits(5000)
# B = randomDigits(5000)

A = randomDigits(5000)          
B = randomDigits(5000)

#print(A)
#print(B)

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

