#!/usr/bin/env python
# coding: utf-8

#imports
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import math
from sympy.combinatorics.graycode import GrayCode

#Neoklis Vaindirlis : 03116191
#AM=2
#Spyridon Pavlatos : 03116113
AM=5

def dB(a):
    return 10*math.log10(a)

def sample_signal(y,fs,start_sampling,end_sampling):
    """A function that samples a signal y with sampling
    frequency fs in the time interval [start_sampling,end_sampling]"""
    Ts=1/fs
    n=np.arange(start_sampling,end_sampling+1e-8,Ts)
    yd=y(n)
    #yd=[]
    #for i in range(len(n)):
        #yd.append(y(n[i]))
    #yd=np.array(yd)
    return yd,n

fm=1000*AM
A=2
T=1/fm

def trig(t):
   return A*signal.sawtooth(2*math.pi*fm*t,0.5)

#part 1a
fs1=20*fm
fs2=80*fm
#sampling triagonal waves with fs1=20fm and fs2=80fm
yd1,n1=sample_signal(trig,fs1,0,4*T)
yd2,n2=sample_signal(trig,fs2,0,4*T)
#plot sampled signals
plt.figure()
plt.plot(n1,yd1,'ro', markersize=4)
plt.grid()
plt.xlabel("Time(seconds)")
plt.ylabel("Amplitude(V)")
plt.title("Sampled triangle wave 2V,{}Hz(fm) with sampling frequency {}Hz(20fm)".format(fm,fs1))
plt.figure()
plt.plot(n2,yd2,'bo', markersize=4)
plt.grid()
plt.xlabel("Time(seconds)")
plt.ylabel("Amplitude(V)")
plt.title("Sampled triangle wave 2V,{}Hz(fm) with sampling frequency {}Hz(80fm)".format(fm,fs2))
# plot combined graph of sampled signals
plt.figure()
plt.plot(n2,yd2,'bo', label="fs2=80fm", markersize=2)
plt.plot(n1,yd1,'ro', label="fs1=20fm", markersize=2)
plt.grid()
plt.legend()
plt.xlabel("Time(seconds)")
plt.ylabel("Amplitude(V)")
plt.title("Compined graphs")

#part 1b
fs3=5*fm
#sampling with fs=5fm
yd3,n3=sample_signal(trig,fs3,0,4*T)
#plot sampled signal
plt.figure()
plt.plot(n3,yd3,'go',markersize=4, linewidth=0.1)
plt.grid()
plt.xlabel("Time(seconds)")
plt.ylabel("Amplitude(V)")
plt.title("Sampled triangle wave 2V,{}Hz(fm) with sampling frequency {}Hz(5fm)".format(fm,fs3))

#part 1.c.i

def sine(t):
   return np.sin(2*math.pi*fm*t)

#sampling sine waves with fs1=20fm and fs2=80fm
yd4,n4=sample_signal(sine,fs1,0,4*T)
yd5,n5=sample_signal(sine,fs2,0,4*T)
#plot sampled signals
plt.figure()
plt.plot(n4,yd4,'ro', markersize=4)
plt.grid()
plt.xlabel("Time(seconds)")
plt.ylabel("Amplitude(V)")
plt.title("Sampled sine wave 1V,{}Hz(fm) with sampling frequency {}Hz(20fm)".format(fm,fs1))
plt.figure()
plt.plot(n5,yd5,'bo', markersize=4)
plt.grid()
plt.xlabel("Time(seconds)")
plt.ylabel("Amplitude(V)")
plt.title("Sampled sine wave 1V,{}Hz(fm) with sampling frequency {}Hz(80fm)".format(fm,fs2))
#plot combined graph of sampled signals
plt.figure()
plt.plot(n5,yd5,'bo', label="fs2=80fm", markersize=2)
plt.plot(n4,yd4,'ro', label="fs1=20fm", markersize=2)
plt.grid()
plt.legend()
plt.xlabel("Time(seconds)")
plt.ylabel("Amplitude(V)")
plt.title("Combined graphs")
#sampling sine wave with fs=5fm
yd6,n6=sample_signal(sine,fs3,0,4*T)
#plot sampled signal
plt.figure()
plt.plot(n6,yd6,'go',markersize=4,linewidth=0.5)
plt.grid()
plt.xlabel("Time(seconds)")
plt.ylabel("Amplitude(V)")
plt.title("Sampled sine wave 1V,{}Hz(fm) with sampling frequency {}Hz(5fm)".format(fm,fs3))

#part 1.c.ii

def signal_q(t):
    return np.sin(2*math.pi*fm*t)+np.sin(2*math.pi*(fm+1000)*t)

T2=1/math.gcd(fm,fm+1000) #Frequency of sum of sines is equal to gcd of their frequencies
#sampling q(t) with fs1=20fm and fs2=80fm
yd7,n7=sample_signal(signal_q,fs1,0,T2)
yd8,n8=sample_signal(signal_q,fs2,0,T2)
#plot sampled signals
plt.figure()
plt.plot(n7,yd7,'ro', markersize=4)
plt.grid()
plt.xlabel("Time(seconds)")
plt.ylabel("Amplitude(V)")
plt.title("Sampled q(t) with sampling frequency {}Hz".format(fs1))
plt.figure()
plt.plot(n8,yd8,'bo', markersize=4)
plt.grid()
plt.xlabel("Time(seconds)")
plt.ylabel("Amplitude(V)")
plt.title("Sampled q(t) with sampling frequency {}Hz".format(fs2))
#plot combined graphs of sampled signals
plt.figure()
plt.plot(n8,yd8,'bo', label="fs2=80fm", markersize=2)
plt.plot(n7,yd7,'ro', label="fs1=20fm", markersize=2)
plt.legend()
plt.grid()
plt.xlabel("Time(seconds)")
plt.ylabel("Amplitude(V)")
plt.title("Combined graphs")
#sampling q(t) with fs=5fm
yd9,n9=sample_signal(signal_q,fs3,0,T2)
#plot sampled signal
plt.figure()
plt.plot(n9,yd9,'go',markersize=4)
plt.grid()
plt.xlabel("Time(seconds)")
plt.ylabel("Amplitude(V)")
plt.title("Sampled q(t) with sampling frequency {}Hz".format(fs3))

#part 2a
"""We are not going to use the following function to generate n bits Gray Code.
We are going to use sympy's function GrayCode() instead."""
def generateGrayCode(n):
    #function to generate gray code with n bits 
    if n<=0:
        return 0;
    gray=[]
    gray.append("0")
    gray.append("1")
    i=2
    while(i<(1<<n)):
        for j in range(i-1,-1,-1):
            gray.append(gray[j])
        for j in range(0,i):
            gray[j]="0"+gray[j]
        for j in range(i,2*i):
            gray[j]="1"+gray[j]
        i=i<<1
    return gray

#number of bits
if AM%2==1:
    Q=6
else:
    Q=5
L=2**Q #number of quantization levels
D=2*(A+1e-8)/L #step of quantization, +1ε-8 in order to take the interval [-A,A] 

def Mid_Riser(x):
    """A function that takes the input x of a Mid Riser quantizer
    and returns the output of the quantizer"""
    return D*(math.floor(x/D)+1/2)

i=Mid_Riser(-2) 							#lowest quantization level
quantize_levels=[]
for j in range(2**Q):                       #create the quantization levels
    quantize_levels.append(round(i,6))
    i+=D
fs=40*fm #sampling frequency
yd,n=sample_signal(trig,fs,0,4*T) #sampling of triagonal wave with fs=40fm sampling frequency
y_quantized=list(map(Mid_Riser,yd)) #quantization of sampled signal
for j in range(len(y_quantized)): #rounding the quantized values to match quantization levels
    y_quantized[j]=round(y_quantized[j],6)
gray_values=GrayCode(Q) #generate the Gray Code for quantization levels
gray_values=list(gray_values.generate_gray())

#plot the quantizer's output
plt.figure()
plt.xlabel("Time(seconds)")
plt.ylabel("Quantizer Output")
plt.title("6-bits Mid-Riser Quantization of triangle wave 2V,{}Hz(fm) ".format(fm))
plt.yticks(quantize_levels,gray_values, fontsize=5)
plt.step(n,y_quantized,"g",linewidth=0.5,where='post') #parameter where="post" in order to have the correct graph
plt.grid()

#part 2b
def calculate_SNRq(y,y_quantized,N):
    """This function calculates quantization error's variance
    and SNRq"""
    q=[]
    for i in range(N):
        q.append(-y[i]+y_quantized[i])
    q=np.array(q)
    var_q=np.var(q,ddof=1) #calculating quantization's error variance
    x_power=0
    for x in y[0:N-1]:   #calculating signal power
        x_power+=x*x/N
    SNRq=dB(x_power/var_q)
    return var_q,SNRq

var_q, SNRq=calculate_SNRq(yd,y_quantized,10)
print("The quantization error's variance for 10 samples is {} and SNRQ is {}dB.".format(var_q,SNRq))
var_q, SNRq=calculate_SNRq(yd,y_quantized,20)
print("The quantization error's variance for 20 samples is {} and SNRQ is {}dB.".format(var_q,SNRq))
yd_power=0
for x in yd:
    yd_power+=x*x/len(yd)
theoretical_SNRq=dB(3*yd_power*2**(2*Q)/(2*A))
print("The theoretical value of SNRQ is {}dB.".format(theoretical_SNRq))

#part 2c
def polar_RZ(x,Q):
    #function that returns the Polar RZ bitstream of 1 binary number
    b=[]
    for i in range(Q):
        if x[i]=="0":
            b.append((-1)*AM)
        else:
            b.append(AM)
        b.append(0)
    return b

def polar_RZ_bitstream(quantized_signal,t,Q):
    #function that generates the Polar RZ bitstream of a quantized signal in Gray Code representation
    stream=[]
    stream.append(0) #the bitstream has to start from zero
    for i in range(len(quantized_signal)):
        stream+=polar_RZ(quantized_signal[i],Q)
    n=[]
    for i in range(Q*2*len(quantized_signal)+1):
        n.append(i*(t/2))
    return n,stream

def quantized_to_gray(quantized_signal,gray_values,quantize_levels):
    """This function maps quantized signal values to Gray Code representation"""
    a=[]
    for i in range(len(quantized_signal)): 
        x=quantize_levels.index(quantized_signal[i])
        a.append(gray_values[x])
    return a

a=quantized_to_gray(y_quantized[0:len(y_quantized)//4],gray_values,quantize_levels)
n,stream=polar_RZ_bitstream(a,0.001,Q) 
#plot the bitstream
plt.figure()
plt.step(n,stream)
plt.yticks([(-1)*AM,0,AM])
plt.xlabel("Time(seconds)")
plt.ylabel("Amplitude(V)")
plt.title("Bitstream of quantized triangle wave 2V,{}Hz(fm) ".format(fm))
plt.grid()

#part3a
f_m=35
T_m=1/f_m
ka=0.5
fs2=110*fm
def AMDSB(t):
    return (1+ka*np.sin(2*math.pi*f_m*t))*np.sin(2*math.pi*fm*t)
yd,n=sample_signal(AMDSB,fs2,0,4*T_m)
plt.figure()
plt.plot(n,yd,label="AMDSB signal",linewidth=0.5)
plt.ylabel("Amplitude(V)")
plt.xlabel("Time(seconds)")
plt.title("AM modulated signal and the upper envelope")
plt.grid()

#part 3b
"""We will calculate the envelope of the signal using Hilbert Transformation"""
k=np.abs(signal.hilbert(yd))
plt.plot(n,k,label="Upper Envelope") #plotting the upper envelope in the previus figure
plt.legend()
m_k=np.mean(k) 
k=2*(k-m_k)    #subtracting the DC component of the envelope and fixing gain
plt.figure()
plt.plot(n,k)
plt.ylabel("Amplitude(V)")
plt.xlabel("Time(seconds)")
plt.title("Envelope Detector's output - Message signal m(t)")
plt.grid()

"""Demodulation using product-modulator with sine and lowpass filter"""
f=[]
for x in range(len(yd)):
    f.append(yd[x]*math.sin(2*math.pi*fm*n[x]))
b,a=signal.butter(3,700/(fs2*0.5),'lowpass',analog=False)
m=signal.lfilter(b,a,f)
m=np.array(m)
m_m=np.mean(m)
m=m-m_m
m_max=np.max(m)
m=m/m_max #fixing gain and subtracting DC component
plt.figure()
plt.plot(n,m)
plt.ylabel("Amplitude(V)")
plt.xlabel("Time(seconds)")
plt.title("Lowpass FIR output - Message signal m(t)")
plt.grid()

plt.show()
