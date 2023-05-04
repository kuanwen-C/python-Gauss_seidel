import numpy as np
from copy import deepcopy

def check_diagonally(a):
    # Find diagonal coefficients
    diag=np.diag(np.abs(a))
    
    # Find row sum without diagonal coefficients
    off_diag=np.sum(np.abs(a),axis=1)-diag

    if np.all(diag > off_diag):
        print("diagonally dominant")
        return 1
    else:
        print("not diagonally dominant")
        return 0
#===========================================================Gauss Siedel============================================================     
def Gauss_Seidel(a, x, y, it, epsilon):
    converged = False                                 #initialize
    x_old = np.zeros(a.shape[0])  
    
    print("Iteration results")                             
    print(" k      ",end=" ")                          
    for i in range(1,len(a[0])+1):
        print("x",i,"      ",end=" ")
    print(" ")
    # You should to do...
    # Use the for loop to complete the iterative process
        # check if it is smaller than threshold 
        # assign the latest x value to the old value
    for k in range(1,it+1):   #跑iteration
        
        for i in range (0,len(a[0])):           #calculate xi
            temp=0                              #used to store the value of sigma
            for j in range (0,len(a[0])):       #calculate sigma
                if j!=i:                        
                    temp+=a[i][j]*x_old[j] 
            x[i]=(y[i]-temp)/a[i][i]
            
        dx=np.sqrt(np.dot(x-x_old,x-x_old))     #calculate error
        print(k,'      ',end=" ")             
        for i in range(0,len(a[0])-1):
            print('%.4f    '%x[i],end=" " )
        print('%.4f'%x[len(a[0])-1] )
        if (dx<epsilon):                        #determine if it converge
            converged=True
            print("Converged!")
            break
        x_old=deepcopy(x)
    if not converged:
        print("Not converged,increase the number of iteration")

    return x
#====================================Acquire user inputs================================================================================================ 
    num=int(input("A0*X0+A1*X1+A2*X2+........+An*Xn = B, How many varibales(n-1) are there in equation  :"))
coefficient=np.zeros((num,num),dtype=float)             #used a 2-D array to store num variables for num equations
B=np.zeros(num)                                         #used to store the ans for each equation 

#aquire the variables and ans for num equations
for num_of_equation in range(0,num):
    for num_of_variable in range (0,num):
        coefficient[num_of_equation][num_of_variable]=float(input("Please input A%d (For equation%d A0*X0+A1*X1+A2*X2+........+An*Xn = B)"%(num_of_variable,num_of_equation)))
    B[num_of_equation]=float(input("PLease input the ans B (For equation%d A0*X0+A1*X1+A2*X2+........+An*Xn = B):"%num_of_equation))

#show num equations acquired by the user 
for number_of_eq in range (0,num):
    print("Your input No.%d equation: "%number_of_eq,end=" ")
    for num_of_var in range (0,num-1):
        print("%.2f*X%d  + "%(coefficient[number_of_eq][num_of_var],num_of_var),end=" ")
    print("%.2f*X%d "%(coefficient[number_of_eq][num-1] ,num-1),end=" ")
    print("= %.2f "%B[number_of_eq])
#=====================================================================================================================================================
x = np.zeros(num)
check=check_diagonally(coefficient)

if check==1:
    x = Gauss_Seidel(coefficient, x, B, it=50, epsilon=0.0001)
    print('')
    print('Check')
    print('my solve:',x)
    x = np.linalg.solve(coefficient, B)
    print('np solve:',x)
else:
    x = Gauss_Seidel(coefficient, x, B, it=50, epsilon=0.0001)
    print('')
    print('Check')
    print('my solve:',x)
    x = np.linalg.solve(coefficient, B)
    print('np solve:',x)
    print("It's not diagonally dominant, so Gauss Seidel failed")
    
    #==End==============================================================================================================================================
    #by KuanWen-C
