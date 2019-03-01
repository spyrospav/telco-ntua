#!/usr/bin/env python
# coding: utf-8


import matplotlib.pyplot as plt
import numpy as np
from scipy import special
import scipy.io.wavfile
import math
from sympy.combinatorics.graycode import GrayCode	
import binascii


def sample_signal(y,fs,start_sampling,end_sampling):
    """A function that samples a signal y with sampling
    frequency fs in the time interval [start_sampling,end_sampling]"""
    Ts=1/fs
    n=np.arange(start_sampling,end_sampling+1e-8,Ts)
    yd=y(n)
    return yd,n

def dB(a):
    return 10*math.log10(a)


#Neoklis Vaindirlis : 03116191
#AM=2
#Spyridon Pavlatos : 03116113
AM=5


"""part 1a"""

def polar_NRZ(x):
    """function that returns the Polar NRZ (or BPAM) bitstream of 1 binary number"""
    b=[]
    if x==0:
        b.append((-1)*AM)
    else:
        b.append(AM)
    return b

def polar_NRZ_bitstream(sequence,Tb):
    """function that generates the Polar NRZ bitstream of a quantized signal in Gray Code representation"""
    stream=[]
    #stream.append(0) #the bitstream has to start from zero
    for i in range(len(sequence)):
        stream+=polar_NRZ(sequence[i])
    n=np.arange(0,(len(sequence)+1)*Tb,Tb)
    return n,stream

Tb=0.5 #symbol duration
N=24 #number of samples
a=np.random.randint(2,size=N) #24bits random sequence
n,stream=polar_NRZ_bitstream(a,Tb)
stream.insert(0,0)
plt.figure()
plt.step(n,stream)
plt.grid()
plt.xlabel("Time(seconds)")
plt.ylabel("Amplitude(V)")
plt.title("B-PAM of random 24-digit binary sequence")

"""part 1b"""

def BPAM_constellation_diagram(stream,E):
    """Function that returns constellation points for BPAM"""
    stream=np.delete(stream,0)
    IQ_stream=(stream/AM)*math.sqrt(E)
    return IQ_stream

E_bpam=(AM**2)*Tb
constellation_BPAM=BPAM_constellation_diagram(stream,E_bpam)
plt.figure()
plt.grid()
plt.scatter(np.real(constellation_BPAM),np.imag(constellation_BPAM))
plt.axvline(color="black",label="Decision Threshold")
plt.legend()
plt.xlabel("I")
plt.ylabel("Q")
plt.title("Constellation diagram of B-PAM of random 24-digit binary sequence")

"""part 1c"""
def AWGN(Eb,ratio,N):
    N0=Eb/(10**(ratio/10))
    return np.random.normal(0,np.sqrt(N0/2),size=N)+1j*np.random.normal(0,np.sqrt(N0/2),size=N)

stream_a=np.delete(stream,0)

r1=stream_a+AWGN(E_bpam,5,N) #adding AWGN with SNR=5dB
r1=np.insert(r1,0,0)
plt.figure()
plt.step(n,np.real(r1))
plt.grid()
plt.xlabel("Time(seconds)")
plt.ylabel("Amplitude(V)")
plt.title("B-PAM of random 24-digit binary sequence with added AWGN (SNR=5dB)")

r2=stream_a+AWGN(E_bpam,15,N) #adding AWGN with SNR=15dB
r2=np.insert(r2,0,0)
plt.figure()
plt.step(n,np.real(r2))
plt.grid()
plt.xlabel("Time(seconds)")
plt.ylabel("Amplitude(V)")
plt.title("B-PAM of random 24-digit binary sequence with added AWGN (SNR=15dB)")

"""part 1d"""

constellation_r1=BPAM_constellation_diagram(r1,E_bpam)
plt.figure()
plt.grid()
plt.scatter(np.real(constellation_r1),np.imag(constellation_r1),label="Received signal")
plt.scatter(np.real(constellation_BPAM),np.imag(constellation_BPAM),label="Trasmitted signal")
plt.axvline(color="black",label="Decision Threshold")
plt.legend()
plt.xlabel("I")
plt.ylabel("Q")
plt.ylim(-2.5,2.5)
plt.xlim((-1)*math.sqrt(E_bpam)-math.sqrt(AM),math.sqrt(E_bpam)+math.sqrt(AM))
plt.title("Constellation diagram of B-PAM of random 24-digit binary sequence with added AWGN (SNR=5dB)")

constellation_r2=BPAM_constellation_diagram(r2,E_bpam)
plt.figure()
plt.grid()
plt.scatter(np.real(constellation_r2),np.imag(constellation_r2),label="Received signal")
plt.scatter(np.real(constellation_BPAM),np.imag(constellation_BPAM),label="Trasmitted signal")
plt.axvline(color="black",label="Decision Threshold")
plt.legend()
plt.xlabel("I")
plt.ylabel("Q")
plt.ylim(-2.5,2.5)
plt.xlim((-1)*math.sqrt(E_bpam)-math.sqrt(AM),math.sqrt(E_bpam)+math.sqrt(AM))
plt.title("Constellation diagram of B-PAM of random 24-digit binary sequence with added AWGN (SNR=15dB)")

N1=10**5
b=np.random.randint(2,size=N1)
n1,stream2=polar_NRZ_bitstream(b,Tb)
b=np.insert(stream2,0,0)
b=BPAM_constellation_diagram(b,E_bpam)
SNR=np.arange(0,16,1)
BEP=[]
BEP_th=[]

def compare_BPAM(s,r,N):
    """A function that compares RX signal with TX signal for 
       BPAM and calculates the number of errors and Bit Error rate"""
    count=0
    for j in range(N):
        if np.real(s[j])*np.real(r[j])<0:
            count+=1
    return count,count/N

for ratio in range(16):
    r3=b+AWGN(E_bpam,ratio,N1)
    k,bep=compare_BPAM(b,r3,N1)
    BEP.append(bep)
    BEP_th.append(0.5*special.erfc(math.sqrt(10**(ratio/10))))
    print("Experimental BEP for BPAM with {} SNR={}dB. {} errors at {} samples".format(ratio,bep,k,N1))    
    print("Theoretical BEP for BPAM with {} SNR={}dB.".format(ratio,0.5*special.erfc(math.sqrt(10**(ratio/10)))))

fig, (ax1,ax2) = plt.subplots(2,1,sharex=True, sharey=True)
ax1.set_title("BEP diagram")
ax1.plot(SNR,BEP,"ro--",linewidth=0.5,label="Experimental BEP")
ax1.set_ylabel("BEP (Bit Error Probability)")
ax1.grid()
ax1.legend()
ax2.plot(SNR,BEP_th,"bx--",linewidth=0.5,label="Theoretical BEP")
plt.xlabel("SNR $E_b/N_0$ (dB)")
plt.ylabel("BEP (Bit Error Probability)")
plt.grid()
plt.legend()
plt.show()


"""part 2a+b"""

fc=1 if AM%2==0 else 2

def PSK(M,sequence,Tb,fc,phase=0):
    """M-Phase Shift Keying Modulation. Returns modulated signal"""
    k=[]
    n=[]
    gray_values=GrayCode(math.log(M,2))
    gray_values=list(gray_values.generate_gray())
    for j in range(len(sequence)):
        n1=np.arange(math.log(M,2)*Tb*j,math.log(M,2)*Tb*(j+1),1/1000)
        n=np.append(n,n1)
        k=np.append(k,AM*np.cos(2*math.pi*fc*n1+2*math.pi*gray_values.index(sequence[j])/M)+phase)
    return n,k
   
def PSK_sequence(a,N,M):
    """Generate symbol sequence for M-PSK from bitstream"""
    symbol_sequence=[]
    k=int(math.log(M,2))
    for j in range(N//k):
        symbol=str(a[k*j])
        for i in range(1,k):
            symbol+=str(a[k*j+i])
        symbol_sequence.append(symbol)
    return symbol_sequence
    
bpsk=PSK_sequence(a,N,2)
print("BPSK symbol sequence : {}\n".format(", ".join(bpsk)))
n2,bpsk_mod=PSK(2,bpsk,Tb,fc)
plt.figure()
plt.plot(n2,bpsk_mod)
plt.xlabel("Time(seconds)")
plt.ylabel("Amplitude(V)")
plt.title("BPSK modulation with carrier frequency {}Hz and symbol duration {}s".format(fc,Tb))

qpsk=PSK_sequence(a,N,4)
print("QPSK symbol sequence : {}\n".format(", ".join(qpsk)))
n3,qpsk_mod=PSK(4,qpsk,Tb,fc)
plt.figure()
plt.plot(n3,qpsk_mod)
plt.xlabel("Time(seconds)")
plt.ylabel("Amplitude(V)")
plt.title("QPSK modulation with carrier frequency {}Hz and symbol duration {}s".format(fc,2*Tb))

psk8=PSK_sequence(a,N,8)
print("8-PSK symbol sequence : {}".format(", ".join(psk8)))
n4,psk8_mod=PSK(8,psk8,Tb,fc)
plt.figure()
plt.plot(n4,psk8_mod)
plt.xlabel("Time(seconds)")
plt.ylabel("Amplitude(V)")
plt.title("8-PSK modulation with carrier frequency {}Hz and symbol duration {}s".format(fc,3*Tb))

plt.show()


"""part 3a"""

Es=(AM**2)*Tb
Eb=Es/2

def pi4_QPSK_constellation_diagram(sequence,E):
    theta=math.pi/4
    IQ_stream=[]
    for x in sequence:
        if x=="00":
            IQ_stream.append(math.sqrt(E)*np.exp(1j*theta))
        elif x=="01":
            IQ_stream.append(math.sqrt(E)*np.exp(1j*3*theta))
        elif x=="11":
            IQ_stream.append(math.sqrt(E)*np.exp(1j*5*theta))
        elif x=="10":
            IQ_stream.append(math.sqrt(E)*np.exp(1j*7*theta))
    return IQ_stream

def pi4_qpsk_annotate(E):
    plt.annotate("00",(math.sqrt(2*E)/2,math.sqrt(2*E)/2))
    plt.annotate("01",((-1)*math.sqrt(2*E)/2,math.sqrt(2*E)/2))
    plt.annotate("11",((-1)*math.sqrt(2*E)/2,(-1)*math.sqrt(2*E)/2))
    plt.annotate("10",(math.sqrt(2*E)/2,(-1)*math.sqrt(2*E)/2))
    
qpsk_constellation=pi4_QPSK_constellation_diagram(qpsk,Es)
plt.figure()
plt.grid()
plt.scatter(np.real(qpsk_constellation),np.imag(qpsk_constellation),label="Transmitted signal")
plt.axvline(color="black",label="Decision Threshold")
plt.axhline(color="black")
plt.legend(loc=7,fontsize="small")
plt.xlabel("I")
plt.ylabel("Q")
plt.ylim((-1)*math.sqrt(Es)-math.sqrt(AM),math.sqrt(Es)+math.sqrt(AM))
plt.xlim((-1)*math.sqrt(Es)-math.sqrt(AM),math.sqrt(Es)+math.sqrt(AM))
plt.title("Constellation diagram of pi4-QPSK of random 24-digit binary sequence")
pi4_qpsk_annotate(Es)

"""part 3b"""

plt.figure()
plt.grid()
r3=qpsk_constellation+AWGN(Es,5,N//2)
plt.scatter(np.real(qpsk_constellation),np.imag(qpsk_constellation),label="Transmitted signal")
pi4_qpsk_annotate(Es)
plt.scatter(np.real(r3),np.imag(r3),label="Received signal")
plt.axvline(color="black",label="Decision Threshold")
plt.axhline(color="black")
plt.legend(loc=7,fontsize="small")
plt.xlabel("I")
plt.ylabel("Q")
plt.ylim((-1)*math.sqrt(Es)-math.sqrt(AM),math.sqrt(Es)+math.sqrt(AM))
plt.xlim((-1)*math.sqrt(Es)-math.sqrt(AM),math.sqrt(Es)+math.sqrt(AM))
plt.title("Constellation diagram of pi4-QPSK of random 24-digit binary sequence with added AWGN (SNR=5dB)")

plt.figure()
plt.grid()
r4=qpsk_constellation+AWGN(Es,15,N//2)
plt.scatter(np.real(qpsk_constellation),np.imag(qpsk_constellation),label="Transmitted signal")
pi4_qpsk_annotate(Es)
plt.scatter(np.real(r4),np.imag(r4),label="Received signal")
plt.axvline(color="black",label="Decision Threshold")
plt.axhline(color="black")
plt.legend(loc=7,fontsize="small")
plt.xlabel("I")
plt.ylabel("Q")
plt.ylim((-1)*math.sqrt(Es)-math.sqrt(AM),math.sqrt(Es)+math.sqrt(AM))
plt.xlim((-1)*math.sqrt(Es)-math.sqrt(AM),math.sqrt(Es)+math.sqrt(AM))
plt.title("Constellation diagram of pi4-QPSK of random 24-digit binary sequence with added AWGN (SNR=5dB)")


"""part 3c"""

def received_pi4_qpsk(seq):
    """Demodulation of pi4-QPSK signal"""
    received=[]
    for j in range(len(seq)):
        if np.real(seq[j])>=0 and np.imag(seq[j])>=0:
            received.append("00")
        elif np.real(seq[j])<0 and np.imag(seq[j])>0:
            received.append("01")
        elif np.real(seq[j])<=0 and np.imag(seq[j])<0:
            received.append("11")
        elif np.real(seq[j])>0 and np.imag(seq[j])<0:
            received.append("10")
    return received

def compare_qpsk(seq,received):
    """A function that compares RX signal with TX signal for pi4-QPSK and calculates the number and rate of Bit Error"""
    count=0
    for j in range(len(seq)):
        if seq[j][0]!=received[j][0]:
            count+=1
        if seq[j][1]!=received[j][1]:
            count+=1
    return count, count/(2*len(seq))

c=np.random.randint(2,size=N1)
qpsk1=PSK_sequence(c,N1,4)
constellation_qpsk1=pi4_QPSK_constellation_diagram(qpsk1,Es)
BEP2=[]
BEP2_th=[]

for ratio in range(16):
    r5=constellation_qpsk1+AWGN(Es,ratio,N1//2)
    received=received_pi4_qpsk(r5)
    er,BEP_=compare_qpsk(qpsk1,received)
    BEP2.append(BEP_)
    BEP2_th.append(0.5*special.erfc(math.sqrt(0.5*10**(ratio/10))))
    print("Experimental BEP for pi-4 QPSK with SNR={}dB: {}. {} errors at {} samples".format(ratio,BEP2[ratio],er,N1))
    print("Theoretical BEP for pi-4 QPSK with SNR={}dB: {}".format(ratio,BEP2_th[ratio]))
    

fig, (ax1,ax2) = plt.subplots(2,1,sharex=True, sharey=True)
ax1.set_title("BEP diagram")
ax1.plot(SNR,BEP2,"ro--",linewidth=0.5,label="Experimental BEP")
ax1.set_ylabel("BEP (Bit Error Probability)")
ax1.grid()
ax1.legend()
ax2.plot(SNR,BEP2_th,"bx--",linewidth=0.5,label="Theoretical BEP")
plt.xlabel("SNR $E_s/N_0$ (dB)")
plt.ylabel("BEP (Bit Error Probability)")
plt.grid()
plt.legend()
plt.show()



"""part 3di"""


def text_to_bits(text, encoding='ascii'):
    bits = bin(int.from_bytes(text.encode(encoding), 'big'))[2:]
    return bits.zfill(8 * ((len(bits) + 7) // 8))

def text_from_bits(bits):
    q="".join(bits)
    bytes_=[q[8*i:8*i+8] for i in range(len(q)//8)] #splitting the bitstream into bytes(8bits)
    return "".join([chr(int(b,2)) for b in bytes_])

filename="christmas_carol_even.txt" if AM%2==0 else "christmas_carol_odd.txt"
f=open(filename,"r") 
text=f.read()
f.close()
k=text_to_bits(text)
k=k.replace("0b","")


"""part 3diii"""
qpsk_txt=PSK_sequence(k,len(k),4)

Es_txt=Tb
Eb_txt=Es_txt/2
def QPSK_constellation_diagram(sequence,E):
    IQ_stream=[]
    for x in sequence:
        if x=="00":
            IQ_stream.append(math.sqrt(E)*np.exp(0))
        elif x=="01":
            IQ_stream.append(math.sqrt(E)*np.exp(1j*math.pi/2))
        elif x=="11":
            IQ_stream.append(math.sqrt(E)*np.exp(1j*math.pi))
        elif x=="10":
            IQ_stream.append(math.sqrt(E)*np.exp(1j*3*math.pi/2))
    return IQ_stream

x=[-5,5]
y1,y2=[-5,5],[5,-5]

def qpsk_annotate(E):
    plt.annotate("00",(math.sqrt(E),0))
    plt.annotate("01",(0,math.sqrt(E)))
    plt.annotate("11",((-1)*math.sqrt(E),0))
    plt.annotate("10",(0,(-1)*math.sqrt(E)))
    
qpsk_txt_constellation=QPSK_constellation_diagram(qpsk_txt,Es_txt)
plt.figure()
qpsk_annotate(Es_txt)
plt.grid()
plt.scatter(np.real(qpsk_txt_constellation),np.imag(qpsk_txt_constellation))
plt.plot(x,y1,color="black",label="Decision Threshold")
plt.plot(x,y2,color="black")
plt.legend(loc=0,fontsize="small")
plt.xlim((-1)*math.sqrt(Es_txt)-1,math.sqrt(Es_txt)+1)
plt.ylim((-1)*math.sqrt(Es_txt)-1,math.sqrt(Es_txt)+1)
plt.xlabel("I")
plt.ylabel("Q")
plt.title("Constellation diagram of QPSK modulation of text file")

r6=qpsk_txt_constellation+AWGN(Es_txt,5,len(qpsk_txt))
plt.figure()
qpsk_annotate(Es_txt)
plt.grid()
plt.scatter(np.real(r6),np.imag(r6),label="Received signal",s=15)
plt.scatter(np.real(qpsk_txt_constellation),np.imag(qpsk_txt_constellation),label="Transmitted signal",s=15)
plt.plot(x,y1,color="black",label="Decision Threshold")
plt.plot(x,y2,color="black")
plt.legend(loc=0,fontsize="small")
plt.xlim((-1)*math.sqrt(Es_txt)-1,math.sqrt(Es_txt)+1)
plt.ylim((-1)*math.sqrt(Es_txt)-1,math.sqrt(Es_txt)+1)
plt.xlabel("I")
plt.ylabel("Q")
plt.title("Constellation diagram of QPSK modulation of text file with added AWGN 5dB")

plt.figure()
r7=qpsk_txt_constellation+AWGN(Es_txt,15,len(qpsk_txt))
qpsk_annotate(Es_txt)
plt.grid()
plt.scatter(np.real(r7),np.imag(r7),label="Received signal",s=15)
plt.scatter(np.real(qpsk_txt_constellation),np.imag(qpsk_txt_constellation),label="Transmitted signal",s=15)
plt.plot(x,y1,color="black",label="Decision Threshold")
plt.plot(x,y2,color="black")
plt.legend(loc=0,fontsize="small")
plt.xlim((-1)*math.sqrt(Es_txt)-1,math.sqrt(Es_txt)+1)
plt.ylim((-1)*math.sqrt(Es_txt)-1,math.sqrt(Es_txt)+1)
plt.xlabel("I")
plt.ylabel("Q")
plt.title("Constellation diagram of QPSK modulation of text file with added AWGN 15dB")

def received_qpsk(seq):
    received=[]
    for j in range(len(seq)):
        if abs(np.real(seq[j]))>abs(np.imag(seq[j])):
            if np.real(seq[j])>0:
               received.append("00")
            else:
               received.append("11")
        else:
            if np.imag(seq[j])>0:
               received.append("01")
            else:
               received.append("10")
    return received

r6_=received_qpsk(r6)
be_r6,bep_r6=compare_qpsk(qpsk_txt,r6_)
print("Experimental BEP for QPSK with SNR={}dB: {}. {} errors at {}bits".format(4,bep_r6,be_r6,len(k)))
print("Theoretical BEP for QPSK with SNR={}dB: {}".format(4,0.5*special.erfc(math.sqrt(0.5*10**(5/10)))))

r7_=received_qpsk(r7)
be_r7,bep_r7=compare_qpsk(qpsk_txt,r7_)
print("Experimental BEP for QPSK with SNR={}dB: {}. {} errors at {} bits".format(14,bep_r7,be_r7,len(k)))
print("Theoretical BEP for QPSK with SNR={}dB: {}".format(14,0.5*special.erfc(math.sqrt(0.5*10**(15/10)))))

def reconstruct_text(received,SNR):
	received_bitstream=[]
	for c in received:
    		received_bitstream+=c
	received_bitstream="".join(received_bitstream)
	reconstructed_text=text_from_bits(received_bitstream)	
	print("Received text after channel with {}dB AWGN:".format(SNR))
	print(reconstructed_text)
	f_w=open(filename.replace(".txt","")+' with noise {}dB.txt'.format(SNR),'w')
	f_w.write(reconstructed_text)
	f_w.close()

reconstruct_text(r6_,5)
reconstruct_text(r7_,15)

"""part 4"""

soundfile="soundfile1_lab2.wav" if AM%2==1 else "soundfile2_lab2.wav"
rate,data=scipy.io.wavfile.read(soundfile) #read the wav file
n=np.linspace(0,len(data),len(data))
n=n*(1/rate)
plt.figure()
plt.grid()
plt.xlabel("Time(seconds)")
plt.ylabel("")
plt.title("Sound signal")
plt.plot(n,data) #plot the data from wav file

Q=8
L=2**Q #number of quantization levels
max1=max(max(data),abs(min(data)))
D=2*max1/L+1e-8 #step of quantization, +1ε-8 in order to take the interval [-A,A] 

def Mid_Riser(x):
    """A function that takes the input x of a Mid Riser quantizer
    and returns the output of the quantizer"""
    return D*(math.floor(x/D)+1/2)

i=Mid_Riser(min(data)) #lowest quantization level
quantize_levels=[]
for j in range(2**Q):                       #create the quantization levels
    quantize_levels.append(round(i,6))
    i+=D

y_quantized=list(map(Mid_Riser,data)) #quantization of sampled signal
for j in range(len(y_quantized)): #rounding the quantized values to match quantization levels
    y_quantized[j]=round(y_quantized[j],6)
gray_values=GrayCode(Q) #generate the Gray Code for quantization levels
gray_values=list(gray_values.generate_gray())

#plot the quantizer's output
plt.figure()
plt.yticks(quantize_levels,gray_values, fontsize=3)
plt.step(n,y_quantized,"g",linewidth=0.5,where='post') #parameter where="post" in order to have the correct graph
plt.xlabel("Time(seconds)")
plt.ylabel("Quantizer Output")
plt.title("Quantized sound signal (8bit Mid Riser quantization)")
plt.grid()

def quantized_to_gray(quantized_signal,gray_values,quantize_levels):
    """This function maps quantized signal values to Gray Code representation"""
    a=[]
    for i in range(len(quantized_signal)): 
        x=quantize_levels.index(quantized_signal[i])
        a.append(gray_values[x])
    return a

a=quantized_to_gray(y_quantized,gray_values,quantize_levels)
w="".join(a)
qpsk_wav=PSK_sequence(w,len(w),4)

Es_wav=Tb
Eb_wav=Es_wav/2
qpsk_wav_constellation=QPSK_constellation_diagram(qpsk_wav,Es_wav)
plt.figure()
plt.grid()
qpsk_annotate(Es_wav)
plt.scatter(np.real(qpsk_wav_constellation),np.imag(qpsk_wav_constellation),label="Transmitted signal")
plt.plot(x,y1,color="black",label="Decision Threshold")
plt.plot(x,y2,color="black")
plt.legend(loc=0,fontsize="small")
plt.xlim((-1)*math.sqrt(Es_wav)-1,math.sqrt(Es_wav)+1)
plt.ylim((-1)*math.sqrt(Es_wav)-1,math.sqrt(Es_wav)+1)
plt.xlabel("I")
plt.ylabel("Q")
plt.title("Constellation diagram of QPSK of quantized sound signal")

r8=qpsk_wav_constellation+AWGN(Es_wav,4,len(qpsk_wav))
plt.figure()
qpsk_annotate(Es_wav)
plt.grid()
plt.scatter(np.real(r8),np.imag(r8),label="Received signal",s=15)
plt.scatter(np.real(qpsk_wav_constellation),np.imag(qpsk_wav_constellation),label="Transmitted signal",color="orange",s=15)
plt.plot(x,y1,color="black",label="Decision Threshold")
plt.plot(x,y2,color="black")
plt.legend(loc=0,fontsize="small")
plt.xlim((-1)*math.sqrt(Es_wav)-2,math.sqrt(Es_wav)+2)
plt.ylim((-1)*math.sqrt(Es_wav)-2,math.sqrt(Es_wav)+2)
plt.xlabel("I")
plt.ylabel("Q")
plt.title("Constellation diagram of QPSK of quantized sound signal with added AWGN (SNR=4dB)")

plt.figure()
r9=qpsk_wav_constellation+AWGN(Es_wav,14,len(qpsk_wav))
qpsk_annotate(Es_wav)
plt.grid()
plt.scatter(np.real(r9),np.imag(r9),label="Received signal",s=15)
plt.scatter(np.real(qpsk_wav_constellation),np.imag(qpsk_wav_constellation),label="Transmitted signal",color="orange",s=15)
plt.plot(x,y1,color="black",label="Decision Threshold")
plt.plot(x,y2,color="black")
plt.legend(loc=0,fontsize="small")
plt.xlim((-1)*math.sqrt(Es_wav)-2,math.sqrt(Es_wav)+2)
plt.ylim((-1)*math.sqrt(Es_wav)-2,math.sqrt(Es_wav)+2)
plt.xlabel("I")
plt.ylabel("Q")
plt.title("Constellation diagram of QPSK of quantized sound signal with added AWGN (SNR=14dB)")

r8_wav=received_qpsk(r8)
be_r8,bep_r8=compare_qpsk(qpsk_wav,r8_wav)
print("Experimental BEP for QPSK with SNR={}dB: {}".format(4,bep_r8))
print("Theoretical BEP for QPSK with SNR={}dB: {}".format(4,0.5*special.erfc(math.sqrt(0.5*10**(4/10)))))

r9_wav=received_qpsk(r9)
be_r9,bep_r9=compare_qpsk(qpsk_wav,r9_wav)
print("Experimental BEP for QPSK with SNR={}dB: {} ".format(14,bep_r9))
print("Theoretical BEP for QPSK with SNR={}dB: {}".format(14,0.5*special.erfc(math.sqrt(0.5*10**(14/10)))))

def write_soundfile(received,SNR):
    """A function that plots the reconstructed signal
    and creates the wav file after the addition of AWGN"""
    received_=[]
    for i in range((len(received)//4)):
        received_.append(received[4*i]+received[4*i+1]+received[4*i+2]+received[4*i+3])
    u=[]
    for w in received_:
        u.append(gray_values.index(w))
    u=np.array(u)
    u=u.astype("uint8")
    plt.figure()
    plt.plot(n,u)
    plt.xlabel("Time(seconds)")
    plt.ylabel("")
    plt.title("Reconstructed 8PCM sound signal after channel with {}dB AWGN".format(SNR))
    scipy.io.wavfile.write("soundfile{}_with_noise_{}dB.wav".format((AM+1)%2+1,SNR),rate,u)

write_soundfile(r8_wav,4)
write_soundfile(r9_wav,14)

plt.show()

