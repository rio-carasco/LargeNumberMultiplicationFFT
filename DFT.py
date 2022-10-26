#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 16:41:56 2021

@author: riocarasco
"""
# def DFT():
import numpy as np

def DFT(sampleYValues):
    N = len(sampleYValues)
    DFTarray = np.array([])
    for k in range(N):
        X_k = np.array([])
        for n in range(N):
            b = (-2*np.pi*k*n)/(N)
            X_i = sampleYValues[n]*(np.cos(b)+1j*np.sin(b))
            X_k = np.append(X_k, X_i)
        X_k = np.sum(X_k)
        DFTarray = np.append(DFTarray, X_k)
    return DFTarray
        

X  = np.array([5, 3, 2, 1])
DFT = DFT(X)
print(DFT)