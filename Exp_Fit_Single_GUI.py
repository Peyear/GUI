import os
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import sys
import tkinter as Tk
import tkinter.filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from scipy.optimize import curve_fit
from scipy.signal import butter, lfilter
root = Tk.Tk()
bufferlines=[]
Buffer_delete=[]
quencherlines=[]
Quencher_delete=[]
time=[]
a=[]
b=[]
c=0
buffer=[]
quencher=[]
def Reset1():
    python = sys.executable
    os.execl(python, python, * sys.argv)
def Command1():
    root.directory = Tk.filedialog.askdirectory()
    CD=root.directory
    label3=Tk.Label(root,text=CD)
    label3.grid(row=0,column=1,sticky=Tk.W,columnspan=6)
    label4=Tk.Label(root,text='Select Run - -')
    label4.grid(row=1,column=1,sticky=Tk.W)
    label5=Tk.Label(root,text='Select Run - -')
    label5.grid(row=2,column=1,sticky=Tk.W)
def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a
def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y
def func1(x, t0, beta, Y0, A1):
    return A1*np.exp(-1*np.power(x/t0,beta))+Y0
def rate1(beta,t0):
    return (beta/t0)*np.power((.002/t0),beta-1)
def slope1(t0, beta, Y0, A1):
    return -(A1*beta*np.exp(-np.power(.002/t0,beta))*np.power(.002/t0,beta - 1))/t0
def func2(x, t0, beta, Y0, A1):
    return A1*np.exp(1-np.power(np.add(1,x/t0),beta))+Y0
def rate2(beta,t0):
    return (beta/t0)*np.power((1+(.002/t0)),beta-1)
def slope2(t0, beta, Y0, A1):
    return -(A1*beta*np.exp(1 - np.power(.002/t0 + 1,beta))*np.power(.002/t0 + 1,beta - 1))/t0
def func3(x, t0, beta, Y0, A1):
    return (-Y0/1.5)+(2.5/1.5)*(Y0-A1)/(1+1.5*(1-np.exp(-1*np.power(x/t0,beta))))+(2.5*A1/1.5)
def rate3(beta,t0):
    return (beta*np.exp(-np.power(.002/t0,beta))*np.power(.002/t0,beta - 1))/t0
def slope3(t0, beta, Y0, A1):
    return -(3*beta*np.exp(-np.power(.002/t0,beta))*((5*Y0)/3-(5*A1)/3)*np.power(.002/t0,beta-1))/(2*t0*((3*np.exp(-np.power(.002/t0,beta)))/2-5/2)*((3*np.exp(-np.power(.002/t0,beta)))/2-5/2))
def func4(x, t0, beta, Y0, A1):
    return (-Y0/1.5)+(2.5/1.5)*(Y0-A1)/(1+1.5*(1-np.exp(1-np.power(1+x/t0,beta))))+(2.5*A1/1.5)
def rate4(beta,t0):
    return (beta*np.exp(1 - np.power(.000/t0 + 1,beta))*np.power(.000/t0 + 1,beta - 1))/t0
def slope4(t0, beta, Y0, A1):
    return (func4(.000000000000001, t0, beta, Y0, A1)-func4(0, t0, beta, Y0, A1))/.000000000000001
def Buffer1():
    global bufferlines
    global Buffer_delete
    global time
    global a
    root.filename =  Tk.filedialog.askopenfilename(initialdir = root.directory,title = "Select file",filetypes = (("DSA","*.dsa"),("all files","*.*")))
    bufferlines = open(root.filename,'r').read().split('\n')
    label4=Tk.Label(root,text=root.filename[-12:])
    label4.grid(row=1,column=1,sticky=Tk.E)
def Quencher1():
    global quencherlines
    global Quencher_delete
    global b
    root.filename =  Tk.filedialog.askopenfilename(initialdir = root.directory,title = "Select file",filetypes = (("DSA","*.dsa"),("all files","*.*")))
    quencherlines=open(root.filename,'r').read().split('\n')
    label5=Tk.Label(root,text=root.filename[-12:])
    label5.grid(row=2,column=1,sticky=Tk.E)
def Graph1():
    global bufferlines
    global Buffer_delete
    global time
    global a
    global quencherlines
    global Quencher_delete
    global b
    global buffer
    global quencher
    global c
    order = 1
    fs = np.fromstring(bufferlines[20], dtype=int, sep=",")
    cutoff = (fs/2)*.99
    Buffer_delete = entry1.get()
    Buffer_delete=[int(x) for x in Buffer_delete.split()]
    Buffer_delete=np.subtract(Buffer_delete,1)
    time=bufferlines[21]
    time=np.fromstring(time, dtype=float, sep=",")
    a=bufferlines[28]
    a=np.int(a)
    Quencher_delete = entry2.get()
    Quencher_delete=[int(x) for x in Quencher_delete.split()]
    Quencher_delete=np.subtract(Quencher_delete,1)
    b=quencherlines[28]
    b=np.int(b)
    f = Figure(figsize=(5, 4), dpi=100)
    aa = f.add_subplot(111)
    buffer=[]
    for i in range(fs[0]):
        if time[i]==0.002:
            A=i
        if time[i]==.2:
            B=i
        if time[i]<.2:
            B=i
    for i in range(a):
        if i in Buffer_delete:
            'nothing'
        else:
            n=i*2+33
            j=i+1
            buffer1=bufferlines[n]
            buffer1=np.fromstring(buffer1, dtype=float, sep=",")
            buffer2=np.max(buffer1)
            buffer3=np.subtract(buffer1,buffer2)
            buffer4=butter_lowpass_filter(buffer3, cutoff, fs, order)
            buffer1=np.add(buffer4,buffer2)
            buffer.append(buffer1)
            aa.plot(time,buffer1,label='%s B' % j)
    buffer=np.array(buffer)
    quencher=[]
    quenchernum=[]
    for i in range(b):
        if i in Quencher_delete:
            'nothing'
        else:
            n=i*2+33
            j=i+1
            quencher1=quencherlines[n]
            quencher1=np.fromstring(quencher1, dtype=float, sep=",")
            quencher2=np.max(quencher1)
            quencher3=np.subtract(quencher1,quencher2)
            quencher4=butter_lowpass_filter(quencher3, cutoff, fs, order)
            quencher1=np.add(quencher4,quencher2)
            quencher.append(quencher1)
            quenchernum.append(j)
            aa.plot(time,quencher1,label='%s Q' % j)
    quencher=np.array(quencher)
    a=[]
    a=quencher.shape
    a=a[0]
    aa.legend(loc='center left', bbox_to_anchor=(.9, 0.5))
    aa.set_xlabel('Time (s)')
    aa.set_ylabel('Relative Flourescence intensity (v)')
    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.draw()
    canvas.get_tk_widget().grid(row=5,column=0, columnspan=7)
    buffer=buffer.mean(0)
    label09=Tk.Label(root,text='Repeat')
    label09.grid(row=6,column=0,sticky=Tk.W)
    label10=Tk.Label(root,text='   r\u00b2   ')
    label10.grid(row=6,column=1,sticky=Tk.W)
    label11=Tk.Label(root,text='   t\u2080   ')
    label11.grid(row=6,column=2,sticky=Tk.W)
    label12=Tk.Label(root,text='Beta')
    label12.grid(row=6,column=3,sticky=Tk.W)
    label13=Tk.Label(root,text='   Y\u2080   ')
    label13.grid(row=6,column=4,sticky=Tk.W)
    label14=Tk.Label(root,text='   A\u2081   ')
    label14.grid(row=6,column=5,sticky=Tk.W)
    label15=Tk.Label(root,text='Rate')
    label15.grid(row=6,column=6,sticky=Tk.W)
    label16=Tk.Label(root,text='Slope')
    label16.grid(row=6,column=7,sticky=Tk.W)
    for i in range(a):
        xdata = time
        xdata = xdata[A:B]
        ydata = np.divide(quencher[i],buffer.mean())
        ydata = ydata[A:B]
        popt, pcov = curve_fit(func1, xdata, ydata, p0=[.001, .1, .7, .3],bounds=([0, .1, 0, 0], [100, 1, 1, 1]))
        residuals = ydata- func1(xdata, *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((ydata-np.mean(ydata))**2)
        label20=Tk.Label(root,text='%.0f' % quenchernum[i])
        label20.grid(row=7+i,column=0,sticky=Tk.W)
        r_squared = 1 - (ss_res / ss_tot)
        label21=Tk.Label(root,text='%.4f' % r_squared)
        label21.grid(row=7+i,column=1,sticky=Tk.W)
        t0=popt[0]
        label22=Tk.Label(root,text='%.4f' % t0)
        label22.grid(row=7+i,column=2,sticky=Tk.W)
        beta=popt[1]
        label23=Tk.Label(root,text='%.4f' % beta)
        label23.grid(row=7+i,column=3,sticky=Tk.W)
        Y0=popt[2]
        label24=Tk.Label(root,text='%.4f' % Y0)
        label24.grid(row=7+i,column=4,sticky=Tk.W)
        A1=popt[3]
        label25=Tk.Label(root,text='%.4f' % A1)
        label25.grid(row=7+i,column=5,sticky=Tk.W)
        rate=rate1(beta,t0)
        label26=Tk.Label(root,text='%.4f' % rate)
        label26.grid(row=7+i,column=6,sticky=Tk.W)
        slope=slope1(t0, beta, Y0, A1)
        label27=Tk.Label(root,text='%.4f' % slope)
        label27.grid(row=7+i,column=7,sticky=Tk.W)
    if c<a:
        c=a
    elif c>a:
        for i in range(c-a):
            label20=Tk.Label(root,text='- - - - -')
            label20.grid(row=7+a+i,column=0,sticky=Tk.W)
            label21=Tk.Label(root,text='- - - - -')
            label21.grid(row=7+a+i,column=1,sticky=Tk.W)
            label22=Tk.Label(root,text='- - - - -')
            label22.grid(row=7+a+i,column=2,sticky=Tk.W)
            label23=Tk.Label(root,text='- - - - -')
            label23.grid(row=7+a+i,column=3,sticky=Tk.W)
            label24=Tk.Label(root,text='- - - - -')
            label24.grid(row=7+a+i,column=4,sticky=Tk.W)
            label25=Tk.Label(root,text='- - - - -')
            label25.grid(row=7+a+i,column=5,sticky=Tk.W)
            label26=Tk.Label(root,text='- - - - -')
            label26.grid(row=7+a+i,column=6,sticky=Tk.W)
            label27=Tk.Label(root,text='- - - - -')
            label27.grid(row=7+a+i,column=7,sticky=Tk.W)
def Graph2():
    global bufferlines
    global Buffer_delete
    global time
    global a
    global quencherlines
    global Quencher_delete
    global b
    global buffer
    global quencher
    global c
    order = 1
    fs = np.fromstring(bufferlines[20], dtype=int, sep=",")
    cutoff = (fs/2)*.99
    Buffer_delete = entry1.get()
    Buffer_delete=[int(x) for x in Buffer_delete.split()]
    Buffer_delete=np.subtract(Buffer_delete,1)
    time=bufferlines[21]
    time=np.fromstring(time, dtype=float, sep=",")
    a=bufferlines[28]
    a=np.int(a)
    Quencher_delete = entry2.get()
    Quencher_delete=[int(x) for x in Quencher_delete.split()]
    Quencher_delete=np.subtract(Quencher_delete,1)
    b=quencherlines[28]
    b=np.int(b)
    f = Figure(figsize=(5, 4), dpi=100)
    aa = f.add_subplot(111)
    buffer=[]
    for i in range(fs[0]):
        if time[i]==0.002:
            A=i
        if time[i]==.2:
            B=i
        if time[i]<.2:
            B=i
    for i in range(a):
        if i in Buffer_delete:
            'nothing'
        else:
            n=i*2+33
            j=i+1
            buffer1=bufferlines[n]
            buffer1=np.fromstring(buffer1, dtype=float, sep=",")
            buffer2=np.max(buffer1)
            buffer3=np.subtract(buffer1,buffer2)
            buffer4=butter_lowpass_filter(buffer3, cutoff, fs, order)
            buffer1=np.add(buffer4,buffer2)
            buffer.append(buffer1)
            aa.plot(time,buffer1,label='%s B' % j)
    buffer=np.array(buffer)
    quencher=[]
    quenchernum=[]
    for i in range(b):
        if i in Quencher_delete:
            'nothing'
        else:
            n=i*2+33
            j=i+1
            quencher1=quencherlines[n]
            quencher1=np.fromstring(quencher1, dtype=float, sep=",")
            quencher2=np.max(quencher1)
            quencher3=np.subtract(quencher1,quencher2)
            quencher4=butter_lowpass_filter(quencher3, cutoff, fs, order)
            quencher1=np.add(quencher4,quencher2)
            quencher.append(quencher1)
            quenchernum.append(j)
            aa.plot(time,quencher1,label='%s Q' % j)
    quencher=np.array(quencher)
    a=[]
    a=quencher.shape
    a=a[0]
    aa.legend(loc='center left', bbox_to_anchor=(.9, 0.5))
    aa.set_xlabel('Time (s)')
    aa.set_ylabel('Relative Flourescence intensity (v)')
    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.draw()
    canvas.get_tk_widget().grid(row=5,column=0, columnspan=7)
    buffer=buffer.mean(0)
    label09=Tk.Label(root,text='Repeat')
    label09.grid(row=6,column=0,sticky=Tk.W)
    label10=Tk.Label(root,text='   r\u00b2   ')
    label10.grid(row=6,column=1,sticky=Tk.W)
    label11=Tk.Label(root,text='   t\u2080   ')
    label11.grid(row=6,column=2,sticky=Tk.W)
    label12=Tk.Label(root,text='Beta')
    label12.grid(row=6,column=3,sticky=Tk.W)
    label13=Tk.Label(root,text='   Y\u2080   ')
    label13.grid(row=6,column=4,sticky=Tk.W)
    label14=Tk.Label(root,text='   A\u2081   ')
    label14.grid(row=6,column=5,sticky=Tk.W)
    label15=Tk.Label(root,text='Rate')
    label15.grid(row=6,column=6,sticky=Tk.W)
    label16=Tk.Label(root,text='Slope')
    label16.grid(row=6,column=7,sticky=Tk.W)
    for i in range(a):
        xdata = time
        xdata = xdata[A:B]
        ydata = np.divide(quencher[i],buffer.mean())
        ydata = ydata[A:B]
        popt, pcov = curve_fit(func2, xdata, ydata, p0=[.001, .1, .7, .3],bounds=([0, .1, 0, 0], [100, 1, 1, 1]))
        residuals = ydata- func2(xdata, *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((ydata-np.mean(ydata))**2)
        label20=Tk.Label(root,text='%.0f' % quenchernum[i])
        label20.grid(row=7+i,column=0,sticky=Tk.W)
        r_squared = 1 - (ss_res / ss_tot)
        label21=Tk.Label(root,text='%.4f' % r_squared)
        label21.grid(row=7+i,column=1,sticky=Tk.W)
        t0=popt[0]
        label22=Tk.Label(root,text='%.4f' % t0)
        label22.grid(row=7+i,column=2,sticky=Tk.W)
        beta=popt[1]
        label23=Tk.Label(root,text='%.4f' % beta)
        label23.grid(row=7+i,column=3,sticky=Tk.W)
        Y0=popt[2]
        label24=Tk.Label(root,text='%.4f' % Y0)
        label24.grid(row=7+i,column=4,sticky=Tk.W)
        A1=popt[3]
        label25=Tk.Label(root,text='%.4f' % A1)
        label25.grid(row=7+i,column=5,sticky=Tk.W)
        rate=rate2(beta,t0)
        label26=Tk.Label(root,text='%.4f' % rate)
        label26.grid(row=7+i,column=6,sticky=Tk.W)
        slope=slope2(t0, beta, Y0, A1)
        label27=Tk.Label(root,text='%.4f' % slope)
        label27.grid(row=7+i,column=7,sticky=Tk.W)
    if c<a:
        c=a
    elif c>a:
        for i in range(c-a):
            label20=Tk.Label(root,text='- - - - -')
            label20.grid(row=7+a+i,column=0,sticky=Tk.W)
            label21=Tk.Label(root,text='- - - - -')
            label21.grid(row=7+a+i,column=1,sticky=Tk.W)
            label22=Tk.Label(root,text='- - - - -')
            label22.grid(row=7+a+i,column=2,sticky=Tk.W)
            label23=Tk.Label(root,text='- - - - -')
            label23.grid(row=7+a+i,column=3,sticky=Tk.W)
            label24=Tk.Label(root,text='- - - - -')
            label24.grid(row=7+a+i,column=4,sticky=Tk.W)
            label25=Tk.Label(root,text='- - - - -')
            label25.grid(row=7+a+i,column=5,sticky=Tk.W)
            label26=Tk.Label(root,text='- - - - -')
            label26.grid(row=7+a+i,column=6,sticky=Tk.W)
            label27=Tk.Label(root,text='- - - - -')
            label27.grid(row=7+a+i,column=7,sticky=Tk.W)
def Graph3():
    global bufferlines
    global Buffer_delete
    global time
    global a
    global quencherlines
    global Quencher_delete
    global b
    global buffer
    global quencher
    global c
    order = 1
    fs = np.fromstring(bufferlines[20], dtype=int, sep=",")
    cutoff = (fs/2)*.99
    Buffer_delete = entry1.get()
    Buffer_delete=[int(x) for x in Buffer_delete.split()]
    Buffer_delete=np.subtract(Buffer_delete,1)
    time=bufferlines[21]
    time=np.fromstring(time, dtype=float, sep=",")
    a=bufferlines[28]
    a=np.int(a)
    Quencher_delete = entry2.get()
    Quencher_delete=[int(x) for x in Quencher_delete.split()]
    Quencher_delete=np.subtract(Quencher_delete,1)
    b=quencherlines[28]
    b=np.int(b)
    f = Figure(figsize=(5, 4), dpi=100)
    aa = f.add_subplot(111)
    buffer=[]
    for i in range(fs[0]):
        if time[i]==0.002:
            A=i
        if time[i]<10000:
            B=i
    for i in range(a):
        if i in Buffer_delete:
            'nothing'
        else:
            n=i*2+33
            j=i+1
            buffer1=bufferlines[n]
            buffer1=np.fromstring(buffer1, dtype=float, sep=",")
            buffer2=np.max(buffer1)
            buffer3=np.subtract(buffer1,buffer2)
            buffer4=butter_lowpass_filter(buffer3, cutoff, fs, order)
            buffer1=np.add(buffer4,buffer2)
            buffer.append(buffer1)
            aa.plot(time,buffer1,label='%s B' % j)
    buffer=np.array(buffer)
    quencher=[]
    quenchernum=[]
    for i in range(b):
        if i in Quencher_delete:
            'nothing'
        else:
            n=i*2+33
            j=i+1
            quencher1=quencherlines[n]
            quencher1=np.fromstring(quencher1, dtype=float, sep=",")
            quencher2=np.max(quencher1)
            quencher3=np.subtract(quencher1,quencher2)
            quencher4=butter_lowpass_filter(quencher3, cutoff, fs, order)
            quencher1=np.add(quencher4,quencher2)
            quencher.append(quencher1)
            quenchernum.append(j)
            aa.plot(time,quencher1,label='%s Q' % j)
    quencher=np.array(quencher)
    a=[]
    a=quencher.shape
    a=a[0]
    aa.legend(loc='center left', bbox_to_anchor=(.9, 0.5))
    aa.set_xlabel('Time (s)')
    aa.set_ylabel('Relative Flourescence intensity (v)')
    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.draw()
    canvas.get_tk_widget().grid(row=5,column=0, columnspan=7)
    buffer=buffer.mean(0)
    label09=Tk.Label(root,text='Repeat')
    label09.grid(row=6,column=0,sticky=Tk.W)
    label10=Tk.Label(root,text='   r\u00b2   ')
    label10.grid(row=6,column=1,sticky=Tk.W)
    label11=Tk.Label(root,text='   t\u2080   ')
    label11.grid(row=6,column=2,sticky=Tk.W)
    label12=Tk.Label(root,text='Beta')
    label12.grid(row=6,column=3,sticky=Tk.W)
    label13=Tk.Label(root,text='   Y\u2080   ')
    label13.grid(row=6,column=4,sticky=Tk.W)
    label14=Tk.Label(root,text='   A\u2081   ')
    label14.grid(row=6,column=5,sticky=Tk.W)
    label15=Tk.Label(root,text='Rate')
    label15.grid(row=6,column=6,sticky=Tk.W)
    label16=Tk.Label(root,text='Slope')
    label16.grid(row=6,column=7,sticky=Tk.W)
    for i in range(a):
        xdata = time
        xdata = xdata[A:B]
        ydata = np.divide(quencher[i],buffer.mean())
        ydata = ydata[A:B]
        popt, pcov = curve_fit(func3, xdata, ydata, p0=[.001, .1, .7, .3],bounds=([0, .1, 0, 0], [100, 1, 1, 1]))
        residuals = ydata- func3(xdata, *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((ydata-np.mean(ydata))**2)
        label20=Tk.Label(root,text='%.0f' % quenchernum[i])
        label20.grid(row=7+i,column=0,sticky=Tk.W)
        r_squared = 1 - (ss_res / ss_tot)
        label21=Tk.Label(root,text='%.4f' % r_squared)
        label21.grid(row=7+i,column=1,sticky=Tk.W)
        t0=popt[0]
        label22=Tk.Label(root,text='%.4f' % t0)
        label22.grid(row=7+i,column=2,sticky=Tk.W)
        beta=popt[1]
        label23=Tk.Label(root,text='%.4f' % beta)
        label23.grid(row=7+i,column=3,sticky=Tk.W)
        Y0=popt[2]
        label24=Tk.Label(root,text='%.4f' % Y0)
        label24.grid(row=7+i,column=4,sticky=Tk.W)
        A1=popt[3]
        label25=Tk.Label(root,text='%.4f' % A1)
        label25.grid(row=7+i,column=5,sticky=Tk.W)
        rate=rate3(beta,t0)
        label26=Tk.Label(root,text='%.4f' % rate)
        label26.grid(row=7+i,column=6,sticky=Tk.W)
        slope=slope3(t0, beta, Y0, A1)
        label27=Tk.Label(root,text='%.4f' % slope)
        label27.grid(row=7+i,column=7,sticky=Tk.W)
    if c<a:
        c=a
    elif c>a:
        for i in range(c-a):
            label20=Tk.Label(root,text='- - - - -')
            label20.grid(row=7+a+i,column=0,sticky=Tk.W)
            label21=Tk.Label(root,text='- - - - -')
            label21.grid(row=7+a+i,column=1,sticky=Tk.W)
            label22=Tk.Label(root,text='- - - - -')
            label22.grid(row=7+a+i,column=2,sticky=Tk.W)
            label23=Tk.Label(root,text='- - - - -')
            label23.grid(row=7+a+i,column=3,sticky=Tk.W)
            label24=Tk.Label(root,text='- - - - -')
            label24.grid(row=7+a+i,column=4,sticky=Tk.W)
            label25=Tk.Label(root,text='- - - - -')
            label25.grid(row=7+a+i,column=5,sticky=Tk.W)
            label26=Tk.Label(root,text='- - - - -')
            label26.grid(row=7+a+i,column=6,sticky=Tk.W)
            label27=Tk.Label(root,text='- - - - -')
            label27.grid(row=7+a+i,column=7,sticky=Tk.W)
def Graph4():
    global bufferlines
    global Buffer_delete
    global time
    global a
    global quencherlines
    global Quencher_delete
    global b
    global buffer
    global quencher
    global c
    order = 1
    fs = np.fromstring(bufferlines[20], dtype=int, sep=",")
    cutoff = (fs/2)*.99
    Buffer_delete = entry1.get()
    Buffer_delete=[int(x) for x in Buffer_delete.split()]
    Buffer_delete=np.subtract(Buffer_delete,1)
    time=bufferlines[21]
    time=np.fromstring(time, dtype=float, sep=",")
    a=bufferlines[28]
    a=np.int(a)
    Quencher_delete = entry2.get()
    Quencher_delete=[int(x) for x in Quencher_delete.split()]
    Quencher_delete=np.subtract(Quencher_delete,1)
    b=quencherlines[28]
    b=np.int(b)
    f = Figure(figsize=(5, 4), dpi=100)
    aa = f.add_subplot(111)
    buffer=[]
    for i in range(fs[0]):
        if time[i]==0.002:
            A=i
        if time[i]<10000:
            B=i
    for i in range(a):
        if i in Buffer_delete:
            'nothing'
        else:
            n=i*2+33
            j=i+1
            buffer1=bufferlines[n]
            buffer1=np.fromstring(buffer1, dtype=float, sep=",")
            buffer2=np.max(buffer1)
            buffer3=np.subtract(buffer1,buffer2)
            buffer4=butter_lowpass_filter(buffer3, cutoff, fs, order)
            buffer1=np.add(buffer4,buffer2)
            buffer.append(buffer1)
            aa.plot(time,buffer1,label='%s B' % j)
    buffer=np.array(buffer)
    quencher=[]
    quenchernum=[]
    for i in range(b):
        if i in Quencher_delete:
            'nothing'
        else:
            n=i*2+33
            j=i+1
            quencher1=quencherlines[n]
            quencher1=np.fromstring(quencher1, dtype=float, sep=",")
            quencher2=np.max(quencher1)
            quencher3=np.subtract(quencher1,quencher2)
            quencher4=butter_lowpass_filter(quencher3, cutoff, fs, order)
            quencher1=np.add(quencher4,quencher2)
            quencher.append(quencher1)
            quenchernum.append(j)
            aa.plot(time,quencher1,label='%s Q' % j)
    quencher=np.array(quencher)
    a=[]
    a=quencher.shape
    a=a[0]
    aa.legend(loc='center left', bbox_to_anchor=(.9, 0.5))
    aa.set_xlabel('Time (s)')
    aa.set_ylabel('Relative Flourescence intensity (v)')
    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.draw()
    canvas.get_tk_widget().grid(row=5,column=0, columnspan=7)
    buffer=buffer.mean(0)
    label09=Tk.Label(root,text='Repeat')
    label09.grid(row=6,column=0,sticky=Tk.W)
    label10=Tk.Label(root,text='   r\u00b2   ')
    label10.grid(row=6,column=1,sticky=Tk.W)
    label11=Tk.Label(root,text='   t\u2080   ')
    label11.grid(row=6,column=2,sticky=Tk.W)
    label12=Tk.Label(root,text='Beta')
    label12.grid(row=6,column=3,sticky=Tk.W)
    label13=Tk.Label(root,text='   Y\u2080   ')
    label13.grid(row=6,column=4,sticky=Tk.W)
    label14=Tk.Label(root,text='   A\u2081   ')
    label14.grid(row=6,column=5,sticky=Tk.W)
    label15=Tk.Label(root,text='Rate')
    label15.grid(row=6,column=6,sticky=Tk.W)
    label16=Tk.Label(root,text='Slope')
    label16.grid(row=6,column=7,sticky=Tk.W)
    for i in range(a):
        xdata = time
        xdata = xdata[A:B]
        ydata = np.divide(quencher[i],buffer.mean())
        ydata = ydata[A:B]
        popt, pcov = curve_fit(func4, xdata, ydata, p0=[.001, .1, .7, .3],bounds=([0, .1, 0, 0], [100, 1, 1, 1]))
        residuals = ydata- func4(xdata, *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((ydata-np.mean(ydata))**2)
        label20=Tk.Label(root,text='%.0f' % quenchernum[i])
        label20.grid(row=7+i,column=0,sticky=Tk.W)
        r_squared = 1 - (ss_res / ss_tot)
        label21=Tk.Label(root,text='%.4f' % r_squared)
        label21.grid(row=7+i,column=1,sticky=Tk.W)
        t0=popt[0]
        label22=Tk.Label(root,text='%.4f' % t0)
        label22.grid(row=7+i,column=2,sticky=Tk.W)
        beta=popt[1]
        label23=Tk.Label(root,text='%.4f' % beta)
        label23.grid(row=7+i,column=3,sticky=Tk.W)
        Y0=popt[2]
        label24=Tk.Label(root,text='%.4f' % Y0)
        label24.grid(row=7+i,column=4,sticky=Tk.W)
        A1=popt[3]
        label25=Tk.Label(root,text='%.4f' % A1)
        label25.grid(row=7+i,column=5,sticky=Tk.W)
        rate=rate4(beta,t0)
        label26=Tk.Label(root,text='%.4f' % rate)
        label26.grid(row=7+i,column=6,sticky=Tk.W)
        slope=slope4(t0, beta, Y0, A1)
        label27=Tk.Label(root,text='%.4f' % slope)
        label27.grid(row=7+i,column=7,sticky=Tk.W)
    if c<a:
        c=a
    elif c>a:
        for i in range(c-a):
            label20=Tk.Label(root,text='- - - - -')
            label20.grid(row=7+a+i,column=0,sticky=Tk.W)
            label21=Tk.Label(root,text='- - - - -')
            label21.grid(row=7+a+i,column=1,sticky=Tk.W)
            label22=Tk.Label(root,text='- - - - -')
            label22.grid(row=7+a+i,column=2,sticky=Tk.W)
            label23=Tk.Label(root,text='- - - - -')
            label23.grid(row=7+a+i,column=3,sticky=Tk.W)
            label24=Tk.Label(root,text='- - - - -')
            label24.grid(row=7+a+i,column=4,sticky=Tk.W)
            label25=Tk.Label(root,text='- - - - -')
            label25.grid(row=7+a+i,column=5,sticky=Tk.W)
            label26=Tk.Label(root,text='- - - - -')
            label26.grid(row=7+a+i,column=6,sticky=Tk.W)
            label27=Tk.Label(root,text='- - - - -')
            label27.grid(row=7+a+i,column=7,sticky=Tk.W)
button1=Tk.Button(text='Select Directory',fg='black',command=Command1)
button2=Tk.Button(text='Select Buffer',fg='orange',command=Buffer1) 
button3=Tk.Button(text='Select Quencher',fg='yellow',command=Quencher1)
button4=Tk.Button(text='Stretched Exponentail fit(.2 sec at 2ms)',fg='black',command=Graph1)
button5=Tk.Button(text='Modified Exponentail fit(.2 sec at 2ms)',fg='black',command=Graph2)
button6=Tk.Button(text='KSV Stretched fit(at 2ms)',fg='black',command=Graph3)
button7=Tk.Button(text='KSV Modified fit(at 0ms)',fg='black',command=Graph4)
button8=Tk.Button(text='Reset',fg='black',command=Reset1)
label1=Tk.Label(root,text='Buffer Remove=')
label2=Tk.Label(root,text='Quencher Remove=')
entry1=Tk.Entry(root)
entry2=Tk.Entry(root)
button1.grid(row=0,column=0,sticky=Tk.W)
button2.grid(row=1,column=0,sticky=Tk.W)
button3.grid(row=2,column=0,sticky=Tk.W)
button4.grid(row=3,column=0,sticky=Tk.W,columnspan=2)
button5.grid(row=4,column=0,sticky=Tk.W,columnspan=2)
button6.grid(row=3,column=3,sticky=Tk.W,columnspan=2)
button7.grid(row=4,column=3,sticky=Tk.W,columnspan=2)
button8.grid(row=0,column=6,sticky=Tk.E)
label1.grid(row=1,column=3,sticky=Tk.W)
entry1.grid(row=1,column=4,sticky=Tk.W)
label2.grid(row=2,column=3,sticky=Tk.W)
entry2.grid(row=2,column=4,sticky=Tk.W)
root.mainloop()