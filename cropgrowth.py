import numpy as np
import matplotlib.pyplot as plt
from tkinter import *
from matplotlib.widgets import RadioButtons
#crop growth, SUCROS87
#Kendall Knight
def plantingtime(j):
    global timeafterplanting #creates a global variable that is the time after planting for each season
    timeafterplanting = j - 101
    
def calc_Rl(t):
    global Rl
    
    toptopex = -7*(t-45)
    
    topex = toptopex/100
    top = 35* np.e**topex
    bottom = 2 * (np.e**topex + 1)**2 #This is the extrapolated function for the growth of leaves 
    Rl = (top/bottom)
    
    return Rl
def calc_LAI():
    global Lt, L0, LAI, GroundArea,timeafterplanting #this is the leaf area index (leaf area/area of ground)
    
    
    LAI = (calc_Rl(timeafterplanting)/GroundArea)
    return LAI
def calc_Rm():
    global Rm,dW,Ad, resetClimb
    if resetClimb:
        Rm = calc_A_co2() * .0001
        resetClimb = False
        
        return Rm
    else:
        if len(growthRateArray) == 2:  #fits it to a linear model for the first 2 points
            i = 1
        elif len(growthRateArray) >= 3 and len(growthRateArray) < 50: #fits to a quadratic 
            i = 2
        else:
            i = 3 #fits to a third degree
        
        x = np.arange(0,len(growthRateArray),1)
        
        dWintergral = np.trapz(growthRateArray, x) #takes the intergral of the growth rate
        listofplantfractions = [.25,.25,.25,.25]
        
        combinedweight = 0
        for i in range(4):
            combinedweight += dWintergral * listofplantfractions[i] #solves for the maintenance respiration using the intergral
        Rm = (combinedweight + Ad) * .0001
        
def calc_Growth_rate():
    global dW,Cf, Ad, Rm, timeafterplanting
    
    calc_A_co2()
    
    dW = Cf*(Ad - Rm)
    growthRateArray.append(dW)
    #main growth rate function that uses all the varibales that are solved for
    return dW
def calc_A_co2():
    global Al, Am, sigma, k, p, I0, LAI,Ad
    Ad = 0
    for i in range(-1,2):
        level = (.5 + i *np.sqrt(.15))*calc_LAI()
        
        Al = Am*(1-np.e**(-sigma * k*(1-p)*I0*np.e**(-k*level))/Am)
        
        if i == 0:
            Ad += Al*1.6
        else:
            Ad += Al
    
    Ad = (Ad*calc_LAI())/3.6
    Ad = Ad * 30/44
    #This solves for the absorbption rate using given constants and the solved LAI
    
    
    return Ad
def reset(i):
    global LAI, Rm, dW, growthRateArray, timeafterplanting, Rl, Ad, dec1, dec2, dec3, dec4,dec5
    
    LAI = 5
    Rm = 0
    dW = 0
    Lt = 0
    Ad = 0
    Rl = 0
    Ad = 0
    
    if i == 1:
        dec1 = growthRateArray
        
    elif i < 10 and i >1:
        for j in range(len(growthRateArray)):
            dec1[j] = dec1[j] + growthRateArray[j]
        
            
    elif i == 10:
        dec2 = growthRateArray
        
    elif i > 10 and i < 20:
        for j in range(len(growthRateArray)):
            dec2[j] = dec2[j] + growthRateArray[j]
        
    elif i == 20:
        dec3 = growthRateArray
        
    elif i > 20 and i < 30:
        for j in range(len(growthRateArray)):
            dec3[j] = dec3[j] + growthRateArray[j]
        
    elif i == 30:
        dec4 = growthRateArray
        
    elif i > 30 and i < 40:
        for j in range(len(growthRateArray)):
            dec4[j] = dec4[j] + growthRateArray[j]
        
    elif i == 40:
        dec5 = growthRateArray
        
    elif i>40:
        for j in range(len(growthRateArray)):
            dec5[j] += growthRateArray[j]
        
        
              #this large statement resets the parameters for the next season and created the decade arrays 
    
    
    growthRateArray = []
    timeafterplanting = 0
    growthRateArray.append(dW)
def average():
    global dec1,dec2,dec3,dec4,dec5
    
    #this averages the decade array sums
    for i in range(len(dec1)):
        dec1[i] = dec1[i]/9
        dec2[i] = dec2[i]/10
        dec3[i] = dec3[i]/10
        dec4[i] = dec4[i]/10
        dec5[i] = dec5[i]/10

def lightfunc(label):
#this is the code for the radio buttons that changes all the titles and axes
    ydict = {'Growth Rates W/':noChange, 'Growth Rates With':[dec1,dec2,dec3,dec4,dec5], 'Intergral':trapplot}
    if label == "Growth Rates With":
        d1.set_visible(True)
        d2.set_visible(True)
        d3.set_visible(True)
        d4.set_visible(True)
        d5.set_visible(True)
        ax.axis([0,120,0,180])
        l.set_xdata(x)
        l.set_ydata(dec1)
        ax.set_title("Crop Growth Rate")
        ax.set_xlabel("Time After Planting")
        ax.set_ylabel("KG DM ha^-1 d^-1")
    elif label == 'Intergral':
        
        l.set_ydata(trapplot)
        l.set_xdata(trapplotx)
        ax.axis([1,5,9620,9623 ])
        d1.set_visible(False)
        d2.set_visible(False)
        d3.set_visible(False)
        d4.set_visible(False)
        d5.set_visible(False)
        ax.set_title("Change in Yield")
        ax.set_xlabel("Decade")
        ax.set_ylabel("KG of DM")
        
    else:
        ydata = ydict[label] 
        l.set_ydata(ydata)
        d1.set_visible(False)
        d2.set_visible(False)
        d3.set_visible(False)
        d4.set_visible(False)
        d5.set_visible(False)
        l.set_xdata(x)
        ax.axis([0,120,0,180])
        ax.set_title("Crop Growth Rate")
        ax.set_xlabel("Time After Planting")
        ax.set_ylabel("KG DM ha^-1 d^-1")
        
    plt.draw()
    
#declaration of all variables with the parameters given by the paper
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
I0 = 81
sigma = .45
Am = 40
Cf = 1/5.82
Rm = 0
Al = 0
Lt = 0
GroundArea = .114
timeafterplanting = 0

#decade array creation
dec1 = np.zeros(120)
dec2 = np.zeros(120)
dec3 = np.zeros(120)
dec4 = np.zeros(120)
dec5 = np.zeros(120)
resetClimb = False

#main loop structure of the code moving  through the days of the season and each season
#it also is the code to decrease the light levels
for i in range(50):
    if i > 0:
        reset(i)
    I0 = I0 * .998
    
    for j in range(365):
        if j > 100 and j < 220:
            if j == 101:
                resetClimb = True
            else:
                resetClimb = False
            plantingtime(j)
            calc_Rm()
            
            calc_Growth_rate()
    if i == 0:
             noChange = growthRateArray
            
average()
x = np.arange(0,len(growthRateArray),1)

#takes the intergral of the growth rate to graph and creates an x arrray for it
trapplot = [np.trapz(dec1,x), np.trapz(dec2, x),np.trapz(dec3, x), np.trapz(dec4,x), np.trapz(dec5, x)]
trapplotx = np.arange(1, len(trapplot)+1,1)

#bunch of tkinter crap to make the radio buttons work
fig, ax = plt.subplots()
l, = ax.plot(x, noChange, lw=2, color='red')
d1, = ax.plot(x,dec1,visible = False, lw = 2)
d2, = ax.plot(x,dec2,visible = False, lw = 2)
d3, = ax.plot(x,dec3,visible = False, lw = 2)
d4, = ax.plot(x,dec4,visible = False, lw = 2)
d5, = ax.plot(x,dec5,visible = False, lw = 2)
plt.subplots_adjust(left=0.4)

axcolor = 'lightgoldenrodyellow'
rax = plt.axes([0.05, .7, 0.3, 0.15], axisbg=axcolor)
radio = RadioButtons(rax, ('Growth Rates W/', 'Growth Rates With', 'Intergral'))
radio.on_clicked(lightfunc)
ax.set_title("Crop Growth Rate")
ax.set_xlabel("Time After Planting")
ax.set_ylabel("KG DM ha^-1 d^-1")
plt.show()
