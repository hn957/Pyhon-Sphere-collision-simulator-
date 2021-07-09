# Copyright 2021 Hao Nguyen hhn@bu.edu
import sys
import numpy 
import math
from copy import deepcopy


velocity=[]
Position=[]
Radius=[]
Mass=[]
Name =[]
num_sphere=1

time=0
bounces=0
plog=[]
vlog=[]
collide_time=[]
event_type=[]
event_sphere=[]
elog=[]
mlog=[]

#the main loop, first call next Event to determine the earliest collision
def mainloop(velocity,Position,Radius,Mass,time,bounces,time_limit):
  

  while True:
    t,a,b,typ=nextEvent(velocity,Position,Radius,Mass)
    

    

    if (time+t)>time_limit:
      break

    if typ=="colliding":
      
      Forward(velocity,Position,t)
      Update(velocity,Position,Mass,a,b)
      
      plog.append(deepcopy(Position))  
      vlog.append(deepcopy(velocity))

      event_sphere.append([a,b])
      event_type.append(typ)
      time=time+t

      
      
      collide_time.append(time)
      bounces=bounces+1

      elog.append(Energy(velocity,Mass))
      mlog.append(deepcopy(Momentum(velocity,Mass)))








      
    elif typ=='reflecting':

      Forward(velocity,Position,t)
      Update_reflection(velocity,Position,Mass,a)


      
      plog.append(deepcopy(Position))
      vlog.append(deepcopy(velocity))

      event_sphere.append([a,b])
      event_type.append(typ)
      time=time+t
      collide_time.append(time)
      bounces=bounces+1

      elog.append(Energy(velocity,Mass))
      mlog.append(deepcopy(Momentum(velocity,Mass)))

      #break


    

  return bounces
  



def Update(velocity,Position,Mass,a,b):

  v1=numpy.array(velocity[a])
  v2=numpy.array(velocity[b])
  r1=numpy.array(Position[a])
  r2=numpy.array(Position[b])
  m1=Mass[a]
  m2=Mass[b]

  x=(2*m2/(m1+m2))*(numpy.dot(numpy.subtract(v1,v2),numpy.subtract(r1,r2))/numpy.dot(numpy.subtract(r1,r2),numpy.subtract(r1,r2)))
  vn1=numpy.subtract(v1,x*numpy.subtract(r1,r2))

  y=(2*m1/(m1+m2))*(numpy.dot(numpy.subtract(v2,v1),numpy.subtract(r2,r1))/numpy.dot(numpy.subtract(r2,r1),numpy.subtract(r2,r1)))
  vn2=numpy.subtract(v2,y*numpy.subtract(r2,r1))

  for i in range (0,3):
    velocity[a][i]=vn1[i]
    velocity[b][i]=vn2[i]

def Update_reflection(velocity,Position,Mass,a):
  b=0
  v1=numpy.array(velocity[a])
  v2=numpy.array(velocity[b])
  r1=numpy.array(Position[a])
  r2=numpy.array(Position[b])
  m1=Mass[a]
  m2=Mass[b]

  x=(2)*(numpy.dot(numpy.subtract(v1,v2),numpy.subtract(r1,r2))/numpy.dot(numpy.subtract(r1,r2),numpy.subtract(r1,r2)))
  vn1=numpy.subtract(v1,x*numpy.subtract(r1,r2))

  for i in range (0,3):
    velocity[a][i]=vn1[i]


def Forward(velocity,Position,t):
  for i in range(1,len(Position)):
    for j in range (0,3):
      Position[i][j]=Position[i][j]+t*velocity[i][j]

def equationroots( a, b, c):  
  
    # calculating discriminant using formula 
    dis = b * b - 4 * a * c  
    sqrt_val = math.sqrt(abs(dis))  
      
    # checking condition for discriminant 
    if dis > 0:  
        #print(" real and different roots ")  
        x=(-b + sqrt_val)/(2 * a)
        y=(-b - sqrt_val)/(2 * a) 

        return [min(x,y),max(x,y)]

        
      
    # when discriminant is less than 0 
    else: 
        #print("Complex Roots")  
        return [-1,-1]


def Reflection(velocity,Position,Radius,a):
  b=0
  pi=numpy.array(Position[a])
  pj=numpy.array(Position[b])
  vi=numpy.array(velocity[a])
  vj=numpy.array(velocity[b])

  x=numpy.subtract(pi,pj)
  y=numpy.subtract(vi,vj)
  z=numpy.dot(x,y)
  r=abs(Radius[a]-Radius[b])

  c=numpy.dot(x,x)-(r*r)
  b=2*numpy.dot(x,y)
  a=numpy.dot(y,y)

  t=equationroots(a,b,c)
  # print(ri,vi,rj,vj)
  # print(velocity[a],Position[a])
  return t


def collision_time(velocity,Position,Radius,a,b):
  pi=numpy.array(Position[a])
  pj=numpy.array(Position[b])
  vi=numpy.array(velocity[a])
  vj=numpy.array(velocity[b])

  x=numpy.subtract(pi,pj)
  y=numpy.subtract(vi,vj)
  z=numpy.dot(x,y)
  r=Radius[a]+Radius[b]

  c=numpy.dot(x,x)-(r*r)
  b=2*numpy.dot(x,y)
  a=numpy.dot(y,y)

  t=equationroots(a,b,c)
  # print(ri,vi,rj,vj)
  # print(velocity[a],Position[a])
  return t

def nextEvent(velocity,Position,Radius,Mass):
  t=math.inf
  s1=0
  s2=0
  typ=''

  
  for i in range(1,len(velocity)-1):

    for j in range (i+1,len(velocity)):


      t1,t2=collision_time(velocity,Position,Radius,i,j)

      if (t1<0.00000001) and (t1>0):
         t1=0

      if (t2<0.00000001) and (t2>0):
         t2=0

      if (t1<0) and (t1> -0.00000001):
        t1=0

      if (t2<0) and (t2>-0.00000001):
         t2=0

      

      if (t1<0) and (t2<0):
        continue

      #print(i,i+1,t,Collisioncheck(velocity,Position,t1,i,i+1))
      if (t1>=0) and Collisioncheck(velocity,Position,t1,i,j)==True:
          if (t1<t):
            t=t1
            s1=i
            s2=j
            typ="colliding"
          
          
      elif (t2>=0) and Collisioncheck(velocity,Position,t2,i,j)==True:
        if (t2<t):
          t=t2
          s1=i
          s2=j
          typ="colliding"

  
  for i in range (1,len(velocity)):
    #print(Reflectioncheck(velocity,Position,i),i)

    t1,t2=Reflection(velocity,Position,Radius,i)

    if (t1<0.00000001) and (t1>0):
       t1=0

    if (t2<0.00000001) and (t2>0):
       t2=0

    if (t1<0) and (t1> -0.00000001):
      t1=0

    if (t2<0) and (t2>-0.00000001):
       t2=0

    #print(t1,t2)

    if (t1<0) and (t2<0):
      continue


    if Reflectioncheck(velocity,Position,Radius,t1,i)==True and (t1>=0):
      if (t1<t):
        t=t1
        s1=i
        s2=0
        typ='reflecting'
    elif Reflectioncheck(velocity,Position,Radius,t2,i)==True and (t2>=0):
      if (t2<t):
        t=t2
        s1=i
        s2=0
        typ='reflecting'




  return (t,s1,s2,typ)





  #print(t,s1,s2)

def Collisioncheck(velocity,Position,t,a,b):
  pi=numpy.array(Position[a])
  pj=numpy.array(Position[b])

  vi=numpy.array(velocity[a])
  vj=numpy.array(velocity[b])

  for j in range (0,3):
    pi[j]=pi[j]+t*vi[j]
    pj[j]=pj[j]+t*vj[j]

  
  x=numpy.subtract(pi,pj)
  y=numpy.subtract(vi,vj)
  z=numpy.dot(x,y)
  
  




  if z<0:
    return True
  else:
    return False

def Reflectioncheck(velocity,Position,Radius,t,a):
  b=0
  pi=numpy.array(Position[a])
  pj=numpy.array(Position[b])

  vi=numpy.array(velocity[a])
  vj=numpy.array(velocity[b])

  
  for j in range (0,3):
    pi[j]=pi[j]+t*vi[j]
    pj[j]=pj[j]+t*vj[j]

  x=numpy.subtract(pi,pj)
  
  y=numpy.subtract(vi,vj)
  z=numpy.dot(x,y)

  #print(x,y,z)

  if z>0:

    return True
  else:
    return False

def Energy(velocity,Mass):
  e=0
  for i in range (1,len(velocity)):
    e=e+0.5*numpy.dot(velocity[i],velocity[i])*Mass[i]

  return e

def Momentum(velocity,Mass):
  m=numpy.array([0,0,0])
  for i in range (1,len(velocity)):
    x=numpy.array([Mass[i]*x for x in velocity[i]])
    m=numpy.add(m,x)

  return m.tolist()







try:
  Radius.append(float(sys.argv[1]))
  time_limit=float(sys.argv[2])
  Mass.append(numpy.inf)
  velocity.append([0,0,0])
  Position.append([0,0,0])
  Name.append("Universe")      #check for command line argument. return 1 if fail
except:
  print("command line failure")
  sys.exit(1)

print("Please enter the mass, radius, x/y/z position, x/y/z velocityand name of each sphere")
print("When complete, use EOF / Ctrl-D to stop entering")

#input block, check each stdin line, return 1 if fail
while True:
  input_str= sys.stdin.readline()
  if (input_str ==''):
    break

  try:
    data=input_str.split(' ')
    while '' in data:
      data.remove('')
    

    Mass.append(float(data[0]))
    Radius.append(float(data[1]))

    
    Position.append([float(data[2]),float(data[3]),float(data[4])])
    velocity.append([float(data[5]),float(data[6]),float(data[7])])

    Name.append(data[8].strip())
    num_sphere=num_sphere+1

  except:
    print("input fail")
    sys.exit(1) #stdin block#stding block#input block, check each stdin line, return 1 if fail 



#calculate initial energy and momentum
e1= Energy(velocity,Mass)
p1=Momentum(velocity,Mass)

#output block for initial conditions
print(" ")
print(f"Here are the initial conditions")
print(f"universe radius {Radius[0]:g}")
print(f"end simulation {time_limit:g}")

for i in range (1,num_sphere):
  print(f"{Name[i]} m={Mass[i]:g} R={Radius[i]:g} p=({Position[i][0]:g},{Position[i][1]:g},{Position[i][2]:g}) v=({velocity[i][0]:g},{velocity[i][1]:g},{velocity[i][2]:g})")

print(f"energy: {e1:g}")
print(f"momentum: ({p1[0]:g},{p1[1]:g},{p1[2]:g})") 
print(" ")
print("Here are the events")


bounces=mainloop(velocity,Position,Radius,Mass,time,bounces,time_limit)

#output block
for i in range (0,bounces):
  print(" ")
  print(f"time of event: {collide_time[i]:g} ")
  if (event_type[i]=="colliding"):
    print(f"{event_type[i]} {Name[event_sphere[i][0]]} {Name[event_sphere[i][1]]} ")
  else:
    print(f"{event_type[i]} {Name[event_sphere[i][0]]} ")
  for j in range (1,num_sphere):
    print(f"{Name[j]} m={Mass[j]:g} R={Radius[j]:g} p=({plog[i][j][0]:g},{plog[i][j][1]:g},{plog[i][j][2]:g}) v=({vlog[i][j][0]:g},{vlog[i][j][1]:g},{vlog[i][j][2]:g})")
  print(f"energy: {elog[i]:g}")
  print(f"momentum: ({mlog[i][0]:g},{mlog[i][1]:g},{mlog[i][2]:g})") 
  

  

#mainloop(velocity,Position,Radius,Mass,time,bounces,time_limit)
 


  

