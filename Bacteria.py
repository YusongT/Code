import math
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import glob

from scipy.optimize import leastsq



"""
Function
Snake_Splitter:
    Input a filename.
        open and read one file
        split the snakes
    return(xarray,yarray)
        xarray is an array consists of 2 or 3 data of snakes
        yarray is an array consists of 2 or 3 lists of 5-d lists:
        [data, which snake this is,f_or_s,pc,repeat]
"""

def Snake_Splitter(filename):
    print(filename)
    """analyze the information in the file name"""
    n=0
    index_of_letter=0
    #global re
    for letter in filename:
        if letter=='_':
            n=n+1
            if n==2:
                f_or_s=filename[index_of_letter+1]
                #print(f_or_s)
            if n==3:
                if filename[index_of_letter+2]=='p':
                    pc=filename[index_of_letter+1]
                elif filename[index_of_letter+3]=='p':
                    pc=filename[index_of_letter+1]+filename[index_of_letter+2]
                elif filename[index_of_letter+4]=='p':
                    pc=filename[index_of_letter+1]+filename[index_of_letter+2]+filename[index_of_letter+3]
                #print(pc)
            if n==4:
                
                re=filename[index_of_letter+1]
                if index_of_letter+7==len(filename):
                    #print('re two digits')
                    re=re+filename[index_of_letter+2]
        index_of_letter=index_of_letter+1
    
    
    #read the file
    filein=open(filename,"r")
    lines=filein.readlines()
    filein.close()
    
    #count the number of snakes
    index=0         #the total number of lines in the file
    snakenumber=0   #number of snakes
    startpoints=[]  #the index of the line that start to be the actual data
    
    for line in lines:
        index=index+1
        if line=='0\n':
            startpoints.append(index)
            snakenumber=snakenumber+1
    
    #end points of the snakes
    endpoints=[]
    for n in range(1,snakenumber):#the first startpoint-2 isn't end of any snakes, thus start from startpoint[1].
        endpoints.append(startpoints[n]-2)
    endpoints.append(index)
    
    
    xarray=[]
    yarray=[]
    #snakes
    for index_of_snake in range(0,len(startpoints)):
        n=0
        snakeline=[]
        for n in range(startpoints[index_of_snake],endpoints[index_of_snake]):
            snakeline.append(lines[n])
        
        splitted=[]
        for line in snakeline:
            splitted.append(line.split("\t"))
        
        xlist=[]
        ylist=[]
        for line in splitted:
            y=float(line[2])*0.06449
            y=(y**2)*math.pi
            ylist.append([y,index_of_snake+1,f_or_s,pc,re])
            #y is data in area: um^2
            #ylist is a list consist of lists: [data,index of snake,...]
            x=float(line[3])*3/60
            xlist.append(x)
        xarray.append(xlist)
        yarray.append(ylist)
    
    return(xarray,yarray)
    #xarray is an array consists of 2 or 3 data of snakes
    #yarray is an array consists of 2 or 3 lists of 5-d lists


"""
Function
Fitting: 
    Input xlist and ylist
        does least square fit, 
    plot it
    return k, b values.
"""

def Fitting(X,Y):
    def func(p,x):
        k,b=p
        return k*x+b
    
    def error(p,x,y):
        return func(p,x)-y
    
    p0=[1,20]
    
    Para=leastsq(error,p0,args=(X,Y))
    
    k,b=Para[0]
    #print("k=",k,"b=",b)
    #print("cost："+str(Para[1]))
    #print("The fitted line:")
    #print("y="+str(round(k,2))+"x+"+str(round(b,2)))
    """
    plt.figure(figsize=(8,6))#ratio of graph
    plt.scatter(X,Y,color="green",label="data",linewidth=2) 
    
    n=0
    while max(X)>n:
        #this n will be used to determine 
        #the range of the x axis of the plot later on.
        n=n+1
    
    #Plotting the fitted line.
    x=np.linspace(0,n,100) ##Draw 100 dots between the first two numbers.
    y=k*x+b ##The function
    
    plt.plot(x,y,color="red",label="fitted line",linewidth=2) 
    plt.legend(loc='lower right') 
    plt.xlabel('Time(h)')
    plt.ylabel('Area ($\mu$$m^2$)')
    plt.title('Least Square Fit')
    plt.show()
    """
    return(k,b)







"""main{}"""

"""
Make the data arrays
"""

list_of_blue=['#4dffff','#33bbff','#0040ff','#28288a']
list_of_orange=['#ffea00','#ffbf00','#e67300','#ff4000']
list_of_green=['#99e600','#3bb300','#006600','#003d10'] 
list_of_red=['#ffb3b3','#ff6666','#ff1f1f','#d60000','#8f0000']
# if only one red is needed, use #ff0000 (h=0,s=100,sl=50)
xarray=[]
yarray=[]

#filename=[]
filename=glob.glob(r'*_snakes_f_3pc_1.txt') 
#variable "filename" is the list of file names read in by "glob"
#wildcard character * and ? are used, e.g.: filename=glob.glob(r'*_snakes_?_4pc_?.txt')
for file_index in range(0,len(filename)):
    x,y=Snake_Splitter(filename[file_index])
    #loop around all files of which the file name are
    #in the list "filename".
    
    #then collect all x and ys of these files.
    
    for n in range(0,len(x)):
        #can't just xarray.append(x) 
        #as then the whole x is one element in xarray.
        xarray.append(x[n])
        yarray.append(y[n])
        #e.g. *_snakes_f_3pc_*.txt
        #xarray=[[1,2,3,...]
                #[1,2,3,...]
                #[1,2,3,...]]
        #yarray=[[[1,1,'f','3','4'],[2,1,f,3,4],[3,1,f,3,4],...]
                #[[1,2,f,3,4],[2,2,f,3,4],[3,2,f,3,4],...]
                #[[1,3,f,3,4],[.........],[.........],...]]
#e.g. print(yarray[0][0][0]) give:2.090898835624298


"""
log the data and append them in yarray (as when doing linear
fitting we need the logged data, the original data are of
course not linear):
"""

for n in range(0,len(yarray)):
    for m in range(0,len(yarray[n])):
        yarray[n][m].append(math.log(yarray[n][m][0]))
#now yarray[0][0][]:
#   0       1           2       3     4         5
#[data,which snake,'fast/slow','pc','no.',data after log]



"""
Linear fitting of the linear part of the first snake
"""

"""pick out the first snakes and the end points"""
X1=[]
Y1=[]
end_x1=[]
number_of_snake1=0
for n in range(0,len(xarray)):
    #print('snake1: ',yarray[n][0])
    index_of_snake=yarray[n][0][1]
    if index_of_snake==1:
        number_of_snake1=number_of_snake1+1
        X1.append(xarray[n])
        Y=[]
        for m in range(0,len(yarray[n])):
            Y.append(yarray[n][m])
        Y1.append(yarray[n])
        #print(len(X1[0]),len(Y1[0]))
    #X1 and Y1 are the data matrices for the first snakes
    elif index_of_snake==2:
        end_x1.append(xarray[n][0])


#print(end_x1)
end_points=[]
for n in range(0,len(X1)):
    end_point=0
    for m in range(0,len(X1[n])):
        if X1[n][m]<=end_x1[n]:
            end_point=end_point+1
    end_points.append(end_point)


K=[]
B=[]
xlambda=[]
#print('len(X1): ',len(X1))
for n in range(0,len(X1)):
    Ytemp=[]
    for m in range(0,len(Y1[n])):
        Ytemp.append(Y1[n][m][5])
    xlambda.append(Y1[n][0][3])
    """cut off from the beginning of the second snake"""
    xcut=[]
    ycut=[]
    for i in range(0,end_points[n]):
        xcut.append(X1[n][i])
        ycut.append(Ytemp[i])
    xcutarray=np.array(xcut)
    ycutarray=np.array(ycut)
    
    k,b=Fitting(xcutarray,ycutarray)
    K.append(k)
    B.append(b)
    #print('n: ',n)
#trusting that the order of the K,B and the order that corresponding file shows up in yarray is the same.
which_kb=0
for n in range(0,len(xarray)):
    if yarray[n][0][1]==1:
        #print(which_kb)
        #print('this is the first snake')
        for m in range(0,len(yarray[n])):
            yarray[n][m].append(K[which_kb])
            yarray[n][m].append(B[which_kb])
        which_kb=which_kb+1

#for n in range(0,len(yarray)):
    #if yarray[n][0][1]==1:
        #print('yarray[n][0]: ',yarray[n][0])

#now yarray[0][0][]:
#   0       1           2       3     4         5        6 7
#[data,which snake,'fast/slow','pc','no.',data after log,k,b]
# only the datas for the first snake 
# had been attatched the ks and bs.

#print('K: ',K)#K is the list of K values
#print('B: ',B)#B is the list of B values
#print(xlambda)#xlambda is the list of their x values


'''
"""Average and error bar of gradient and intercept"""
x_average=[]# the x list for the average
k_error=[]# ylist for error of k
b_error=[]# ylist for error of b
k_temp=[]# the list that stores the list of k correspond to same pc
b_temp=[]
k_average=[]#the list that store the sum of the k with same pc
b_average=[]#the list that store the sum of the b with same pc

for n in range(0,len(xlambda)):
    k=K[n]
    b=B[n]
    if n+1<len(xlambda):# when this isn't the last item in the list
        if xlambda[n+1]==xlambda[n]:# if the next would still be the same pc
            k_temp.append(k)
            b_temp.append(b)
        elif xlambda[n+1]!=xlambda[n]:      # if the next is a different pc
            k_temp.append(k)
            b_temp.append(b)
            k_average.append(np.mean(k_temp))# append the sum
            b_average.append(np.mean(b_temp))
            k_error.append(np.std(k_temp))
            b_error.append(np.std(b_temp))
            x_average.append(xlambda[n])
            k_temp=[]# reset the temporary lists to prepare for the next 
            b_temp=[]
    elif n+1>=len(xlambda):# for the last item in the x list
        k_temp.append(k)
        b_temp.append(b)
        k_average.append(np.mean(k_temp))
        b_average.append(np.mean(b_temp))
        k_error.append(np.std(k_temp))
        b_error.append(np.std(b_temp))
        x_average.append(xlambda[n])



"""
Now we have 2 list to be plotted on each graph
"""
"""K"""
plt.figure(figsize=(12,15))
plt.tick_params(labelsize=30)
plt.ylabel('Gradient',size=30)
plt.xlabel('Concentration (%)',size=30)
plt.title('Gradient vs Agar Concentration',size=30)
plt.scatter(xlambda,K,color='#1f77b4')

plt.scatter(x_average,k_average,color='#ff7f0e')
plt.errorbar(x_average,k_average,yerr=k_error,color='#ff7f0e',fmt='o',capsize=8)
#(add "elinewidth=5" change the width of the verticle line of the error bar thicker...)
#plt.savefig('gradient.png',dpi=300)
plt.show()

"""B"""
plt.figure(figsize=(12,15))
plt.tick_params(labelsize=30)
plt.ylabel('Intercept ($\mu$$m^2$)',size=30)
plt.xlabel('Concentration (%)',size=30)
plt.title('Intercept vs Agar Concentration',size=30)
plt.scatter(xlambda,B,color='#1f77b4')

plt.scatter(x_average,b_average,color='#ff7f0e')
plt.errorbar(x_average,b_average,yerr=b_error,color='#ff7f0e',fmt='o',capsize=8)
#plt.savefig('intercept.png',dpi=300)
plt.show()
'''




#Reminder: now yarray[0][0][]:
#   0       1           2       3     4         5        6 7
#[data,which snake,'fast/slow','pc','no.',data after log,k,b]



'''
#Code for plotting area versus time:
"""Plot different snakes with different colour in log scale"""

plt.figure(figsize=(15,12))
plt.tick_params(labelsize=30)
plt.ylabel('Area ($\mu$$m^2$)',size=30)
plt.xlabel('Time (hours)',size=30)
plt.title('Area of Colony vs Time',size=30)

fast_pc_list=['1.5','2','3','4']
slow_pc_list=['2','3.5','4']

for n in range(0,len(xarray)):
    if yarray[n][0][1]==1:
        if yarray[n][0][2]=='f':
            colour=list_of_blue[fast_pc_list.index(yarray[n][0][3])]
        else:
            colour=list_of_blue[slow_pc_list.index(yarray[n][0][3])]
    elif yarray[n][0][1]==2:
        if yarray[n][0][2]=='f':
            colour=list_of_orange[fast_pc_list.index(yarray[n][0][3])]
        else:
            colour=list_of_orange[slow_pc_list.index(yarray[n][0][3])]
    elif yarray[n][0][1]==3:
        if yarray[n][0][2]=='f':
            colour=list_of_green[fast_pc_list.index(yarray[n][0][3])]
        else:
            colour=list_of_green[slow_pc_list.index(yarray[n][0][3])]
    
    ylist=[]
    for m in range(0,len(yarray[n])):
        ylist.append(yarray[n][m][0])
    
    #print(colour)
    plt.semilogy(xarray[n],ylist,color=colour)

#plt.savefig('plot_s_final.png',dpi=300)
plt.show()
'''



'''
"""eliminate intercept for the first snake"""

plt.figure(figsize=(15,12))
plt.tick_params(labelsize=30)
plt.ylabel('log Area ($\mu$$m^2$)',size=30)
plt.xlabel('Time (hours)',size=30)
plt.title('Area Minus Intercept vs Time',size=30)

fast_pc_list=['1.5','2','3','4']
slow_pc_list=['2','3.5','4']

for n in range(0,len(xarray)):
    if yarray[n][0][1]==1:
        if yarray[n][0][2]=='f':
            colour=list_of_blue[fast_pc_list.index(yarray[n][0][3])]
        else:
            colour=list_of_blue[slow_pc_list.index(yarray[n][0][3])]
    elif yarray[n][0][1]==2:
        if yarray[n][0][2]=='f':
            colour=list_of_orange[fast_pc_list.index(yarray[n][0][3])]
        else:
            colour=list_of_orange[slow_pc_list.index(yarray[n][0][3])]
    elif yarray[n][0][1]==3:
        if yarray[n][0][2]=='f':
            colour=list_of_green[fast_pc_list.index(yarray[n][0][3])]
        else:
            colour=list_of_green[slow_pc_list.index(yarray[n][0][3])]
    
    ylist=[]
    if yarray[n][0][1]==1:
        #print(yarray[n][0])
        for m in range(0,len(yarray[n])):
            ylist.append(yarray[n][m][5]-yarray[n][m][7])
    else:
        for m in range(0,len(yarray[n])):
            ylist.append(yarray[n][m][5])
    plt.plot(xarray[n],ylist,color=colour)

#plt.savefig('plot_all_minus_intercept.png',dpi=300)
plt.show()# the 2nd and 3rd snakes are undisturbed
'''




"""Summing up three snakes for each file"""
#Reminder: now yarray[0][0][]:
#   0       1           2       3     4         5        6 7
#[data,which snake,'fast/slow','pc','no.',data after log,k,b]


xsumlist=[]#the xlist for the sum
ysumlist=[]

for n in range(0,len(xarray)):
    if yarray[n][0][1]==1:#when the current one is the 1st snake in the file
        #print('n: ',n)
        #print('len(yarray): ',len(yarray))
        if n+2<=len(yarray):#can't be the last list (maybe not needed)
            if n+2>=len(yarray) or yarray[n+2][0][1]!=3:# if there are 2 snakes
                print('2 snakes')
                #print(yarray[n+2][0][1])
                xend=min(xarray[n][len(xarray[n])-1],xarray[n+1][len(xarray[n+1])-1])
                xtemplist=[]
                ytemplist=[]
                for m in range(0,len(xarray[n])):
                    if xarray[n][m]<=xend:
                        #print(xarray[n][m])
                        #print(m)
                        xtemplist.append(xarray[n][m])
                        # sum up
                        if xarray[n][m]<xarray[n+1][0]:# before the second snake start
                            yitem=yarray[n][m]
                        elif xarray[n][m]>=xarray[n+1][0]:# after the second snake shows up
                            #print(m)
                            second=0
                            for i in range(0,len(xarray[n+1])):
                                if xarray[n+1][i]<xarray[n][m]:
                                    second=second+1# find the index for which xarray[n+1][]>=xarray[n][m]
                            #print(second)
                            y=yarray[n][m][0]+yarray[n+1][second][0]#sum up
                            
                            yitem=[y,'sum',yarray[n][m][2],yarray[n][m][3],yarray[n][m][4],math.log(y),yarray[n][m][6],yarray[n][m][7]]
                        ytemplist.append(yitem)
                #print('i: ',i)
                #print('xarray[n+1][i]: ',xarray[n+1][i])
                #print(xarray[n])
                #print(xarray[n+1])
                
                #print(n+1)
                #print(yarray[n+1][i][0])
                #print('yarray[n+1][0][0]: ',yarray[n+1][0][0])
                ysumlist.append(ytemplist)
                xsumlist.append(xtemplist)
            if n+2<len(yarray) and yarray[n+2][0][1]==3:# if there are 3 snakes
                #print('3 snakes')
                #print(yarray[n+2][0][1])
                xend=min(xarray[n][len(xarray[n])-1],xarray[n+1][len(xarray[n+1])-1],xarray[n+2][len(xarray[n+2])-1])
                xtemplist=[]
                ytemplist=[]
                for m in range(0,len(xarray[n])):
                    if xarray[n][m]<=xend:
                        #print(xarray[n][m])
                        #print(m)
                        xtemplist.append(xarray[n][m])
                        
                        if xarray[n][m]<xarray[n+1][0]:# before the second snake start
                            yitem=yarray[n][m]
                        elif xarray[n+2][0]>xarray[n][m]>=xarray[n+1][0]:# after the second, before the third
                            second=0
                            for i in range(0,len(xarray[n+1])):
                                if xarray[n+1][i]<xarray[n][m]:
                                    second=second+1# find the index for which xarray[n+1][]>=xarray[n][m]
                            y=yarray[n][m][0]+yarray[n+1][second][0]
                            yitem=[y,'sum',yarray[n][m][2],yarray[n][m][3],yarray[n][m][4],math.log(y),yarray[n][m][6],yarray[n][m][7]]
                        elif xarray[n][m]>=xarray[n+2][0]:
                            second=0
                            third=0
                            for i in range(0,len(xarray[n+1])):
                                if xarray[n+1][i]<xarray[n][m]:
                                    second=second+1# find the index for which xarray[n+1][]>=xarray[n][m]
                            for j in range(0,len(xarray[n+2])):
                                if xarray[n+2][j]<xarray[n][m]:
                                    third=third+1# find the index for which ..., for the 3rd snake
                            #print(third)
                            y=yarray[n][m][0]+yarray[n+1][second][0]+yarray[n+2][third][0]
                            yitem=[y,'sum',yarray[n][m][2],yarray[n][m][3],yarray[n][m][4],math.log(y),yarray[n][m][6],yarray[n][m][7]]
                        ytemplist.append(yitem)
                ysumlist.append(ytemplist)
                xsumlist.append(xtemplist)



plt.figure(figsize=(15,12))
plt.tick_params(labelsize=30)
plt.ylabel('Area ($\mu$$m^2$)',size=30)
plt.xlabel('Time (hours)',size=30)
plt.title('Area of Colony vs Time',size=30)

fast_pc_list=['1.5','2','3','4']
slow_pc_list=['2','3.5','4']
pc_list=['1.5','2','3','3.5','4']
for n in range(0,len(xarray)):
    if yarray[n][0][1]==1:
        if yarray[n][0][2]=='f':
            colour=list_of_blue[fast_pc_list.index(yarray[n][0][3])]
        else:
            colour=list_of_blue[slow_pc_list.index(yarray[n][0][3])]
    elif yarray[n][0][1]==2:
        if yarray[n][0][2]=='f':
            colour=list_of_orange[fast_pc_list.index(yarray[n][0][3])]
        else:
            colour=list_of_orange[slow_pc_list.index(yarray[n][0][3])]
    elif yarray[n][0][1]==3:
        if yarray[n][0][2]=='f':
            colour=list_of_green[fast_pc_list.index(yarray[n][0][3])]
        else:
            colour=list_of_green[slow_pc_list.index(yarray[n][0][3])]
    
    ylist=[]
    for m in range(0,len(yarray[n])):
        ylist.append(yarray[n][m][0])
    
    #print(colour)
    plt.semilogy(xarray[n],ylist,color=colour)

for n in range(0,len(xsumlist)):
    ylist=[]
    for m in range(0,len(ysumlist[n])):
        ylist.append(ysumlist[n][m][0])
    colour=list_of_red[pc_list.index(ysumlist[n][0][3])]
    plt.semilogy(xsumlist[n],ylist,color=colour)
plt.savefig('plot_07_snakes_f_3pc_1_sum.png',dpi=300)
plt.show()










'''
"""
standby 1: 
log axis
"""
plt.figure(figsize=(12,15))
plt.tick_params(labelsize=30)
plt.ylabel('Area (um^2)',size=30)
plt.xlabel('time (h)',size=30)
plt.title('Area versus time',size=30)

for n in range(0,len(xarray)):
    lines=yarray[n]
    ylist=[]
    
    for data in lines:
        ylist.append(data[0])
        
    index_of_snake=data[1]
    if index_of_snake==1:
        plt.semilogy(xarray[n],ylist,color='#1f77b4')
    elif index_of_snake==2:
        plt.semilogy(xarray[n],ylist,color='#ff7f0e')
    elif index_of_snake==3:
        plt.semilogy(xarray[n],ylist,color='#2ca02c')


plt.show()
'''



'''
"""
standby 2:
log the data instead of the axis
"""
plt.figure(figsize=(12,15))
plt.tick_params(labelsize=30)
plt.ylabel('log Area (um^2)',size=30)
plt.xlabel('time (h)',size=30)
plt.title('Area versus time',size=30)

datalist=[]

for n in range(0,len(xarray)):
    lines=yarray[n]
    ylist=[]
    
    for data in lines:
        a=data[0]
        a=math.log(a)
        ylist.append(a)
    
    datalist.append(ylist)
    index_of_snake=data[1]
    if index_of_snake==1:
        plt.plot(xarray[n],ylist,color='#1f77b4')
    elif index_of_snake==2:
        plt.plot(xarray[n],ylist,color='#ff7f0e')
    elif index_of_snake==3:
        plt.plot(xarray[n],ylist,color='#2ca02c')
        

plt.show()
'''