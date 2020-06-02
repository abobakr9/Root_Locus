# importing the required module
import math
import matplotlib.pyplot as plt
import numpy as np
import sympy as sym
plt.xlim(-73,10)
z=[]
p=[0,-25,complex(-50,10),complex(-50,-10)]
s = sym.symbols('s')
chqEqn=s**4+125*s**3+5100*s**2+65000*s
#plot poles at the last of code

#Asymptotes(dotted line)
n=4
m=0
q=n-m
seta=180/q
segma= (p[0]+p[1]+p[2]+p[3])/q
segma=segma.real   #imag=0
endy = 80* math.sin(math.radians(seta))
endx = 80 * math.cos(math.radians(seta))
plt.plot([segma,segma+endx], [0,endy], color='black', linestyle='dashed',linewidth = 1) 
plt.plot([segma,segma-endx], [0,endy], color='black', linestyle='dashed',linewidth = 1)
plt.plot([segma,segma+endx], [0,-endy], color='black', linestyle='dashed', linewidth = 1) 
plt.plot([segma,segma-endx], [0,-endy], color='black', linestyle='dashed', linewidth = 1) 

#break out point
breakAwayPoint=[]
diffchqEqn=sym.diff(chqEqn,s)
a = sym.Poly(diffchqEqn, s)
sol = np.roots(a.all_coeffs())
for x in range(len(sol)):
    if(sol[x].imag ==0):
        breakAwayPoint.append(sol[x].real) #plot break away at the last of the code


#Find the intersections with the imaginary axis using Routh.
#from routh
k = sym.symbols('k')
aux=(4580*65000-125*k)/4580
sol = sym.solve(aux)
imagIntersect=math.sqrt(sol[0]/4580) #plot intersections at the last of the code

#Find Departure angles for complex poles
setadeparture=180-math.degrees(math.atan(10/(-25)))-math.degrees(math.atan(10/-50))-90
endy = 5 * math.sin(math.radians(180-setadeparture))
endx = 5 * math.cos(math.radians(180-setadeparture))
plt.plot([-50,-50+5], [10,10], color='black', linestyle='dashed',linewidth = 2) 
plt.plot([-50,-50-endx], [10,10+endy], color='black', linestyle='dashed',linewidth = 2)
plt.plot([-50,-50+5], [-10,-10], color='black', linestyle='dashed', linewidth = 2) 
plt.plot([-50,-50-endx], [-10,-10-endy], color='black', linestyle='dashed', linewidth = 2) 


#plotting root locus
for k in np.linspace(-10000,10000000,700):
  a = sym.Poly(chqEqn+k, s)
  sol = np.roots(a.all_coeffs())
  for x in range(len(sol)):
      plt.scatter([sol[x].real], [sol[x].imag], color= "black",  
            marker= "o", s=30)

plt.plot([0,-25], [0,0], color= "black", linewidth = 5)
degree=str(round(setadeparture,4))+'°'
plt.text(-49,12,degree,fontsize=10)
plt.text(-49,-13,degree,fontsize=10)

#plot poles
plt.scatter([0,-25,-50,-50], [0,0,10,-10], color= "red",  
            marker= "x", s=150)
plt.title('Root locus plot!')
#plot break away
plt.scatter(breakAwayPoint, [0,0,0], color= "red",  
            marker= "o", s=150)
#plot intersections with imag
plt.scatter([0,0], [imagIntersect,-imagIntersect], color= "yellow",  
            marker= "x", s=150)
plt.show()




