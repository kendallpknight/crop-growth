import numpy as np
import matplotlib.pyplot as plt
from tkinter import *
#crop growth
#Kendall Knight

#all time steps are 1 day
def plantingtime(j):
    global timeafterplanting
    timeafterplanting = j - 100
def calc_Rl(t):
    global Rl
    top = 35* np.e**(-(7*(t-45))/100)
    bottom = 2 * (np.e**(-(7*(t-45))/100) + 1)**2
    Rl = top/bottom
    return Rl
def calc_LAI(t):
    global Lt, L0, LAI, GroundArea
    Lt = L0 * np.e**(calc_Rl(t) * t)
    LAI = Lt/GroundArea
def calc_Rm():
    global Rm,dW,Ad
    if growthRateArray[len(growthRateArray)-1] == 0:
        Rm = Ad * .0001
        return Rm
    else:
        x = np.arange(0,len(growthRateArray),1)
        coefs = np.polynomial.polynomial.polyfit(x,growthRateArray, 3)
        ffit = np.polynomial.polynomial.polyval(x, coefs)
        dWintergral = np.trapz(growthRateArray, x)
        listofplantfractions = [.25,.25,.25,.25]
        combinedweight = 0
        for i in range(4):
            combinedweight += dWintergral * listofplantfractions[i]
        Rm = (combinedweight + Ad) * .0001
        
def calc_Growth_rate():
    global dW,Cf, A, Rm
    
    calc_A_co2()
    dW = Cf*(A - Rm)
    growthRateArray.append(dW)
    return dW
def calc_A_co2():
    global Al, Am, sigma, k, p, I0, LAI, Am, Ad
    
    for i in range(-1,2):
        level = (.5 + i *np.sqrt(.15))*LAI
        Al = Am*(1-np.e**(-sigma * k*(1-p)*I0*np.e**(-k*level))/Am)
        if i == 0:
            Ad += Al*1.6
        else:
            Ad += Al
    Ad = (Ad*LAI)/3.6
    Ad = Ad * 30/44
    return Ad
#declaration of all variables    
growthRateArray = []
Ad = 0
dW = 0
growthRateArray.append(dW)
LAI = 5
Rl = 0
k = .72
DoE = 0
timeline = np.arange(0,18250,1)
A = 0 
L = 0
Daylenght = 0
L0 = .57
p = .0557
I0 = 0
sigma = .45
Am = 40/36000000
Cf = 1/5.82
Rm = 0
Al = 0
Lt = 0
GroundArea = 8.77
timeafterplanting = 0
for i in range(50):
    for j in range(365):
        if j > 100 and j < 220:
            plantingtime(j)
            calc_Rm()
            calc_LAI(j-100)
            calc_Growth_rate()
            print(j+1)
x = np.arange(0,len(growthRateArray),1)
plt.plot(x,growthRateArray)
