# LargeNumberMultiplicationFFT
An extensive research paper on the algorithm, implementation and use is prodivded in "Large_Number_Multiplication_using_Fast_Fourier_Transforms.pdf" in 'main'. Below is a brief description of the code.

## Introduction and Purpose of LNM
This code is designed to perform large number multiplication (lnm) using the Fast Fourier Transform (FFT) algorithm. The algorithm is particularly efficient for multiplying 'Large' numbers and is based on the Cooley-Tukey algorithm. The code is written in Python and utilises NumPy for numerical operations. 

The primary purpose of the code is to challenge/test the efficiency of Photonic Chiplets. The code is currently being exercised by Optalysys LTD for that exact purpose.

## LNM Algorithm Overview
0. **Prequisites:** Both starting numbers must be converted into polynomial form P(x)
1. **Divide and Conquer:** The algorithm starts by dividing the given polynomial P(x)<sub>i</sub> into its odd and even coefficients, creating sub-problems P<sub>even</sub>(x) and P<sub>odd</sub>(x). This process is repeated until the P(x)<sub>i</sub> cannot be divided further.
2. **FFT Component:** The Fourier Transform of the sub-arrays is computed, and then the convolution theorem is applied for multiplication.
3. **Inverse FFT:** After multiplication in the frequency domain, the Inverse Fourier Transform is applied to get the result in the time domain.
4. **Overlap-Add Method:** The final step involves adding up the individual digits, carrying remainders when needed.

## Code Snippets
- **Padding Function:** Ensures that the array length is a power of 2 for efficient FFT.
  ```python
  def Padding(Arr, N):
      return np.pad(Arr, (0, N-len(Arr)), 'constant')

    FFT and Inverse FFT: Utilises NumPy's FFT and inverse FFT functions.

- **FFT and Inverse FFT:** Utilises NumPy's FFT and inverse FFT functions.
  ```python
  FOAList.append(ifft(fft(i) * fft(h)))

- **Overlap-Add Method:** Adds up the individual digits, carrying remainders when needed.
  ```python
    ResultCol = sum(FinalSum[:, -(i+1)]) + Carry

## Dependencies
    Python 3.x
    NumPy

## Sub-Projects
FFT.py and DFT.py are additional code that I developed. These algorithms simply perform the Fast Fourier Transform and Discrete Fourier Transform repectively.

## Applications
The algorithm has applications in various fields including cryptography, machine learning, and AI. It proves very useful for large matrix multiplication, a core part of these fields.

## Author and License
Rio Carasco

The LargeNumberMultiplicationFFT algorithm is distributed on an MIT license and developed as part of a research project in collaboration with Optalysys Ltd.
