import time,math,sys,json,pickle
import numpy as np
# DEFORM is a settlement calculation tool designed initially by URS Corp. in 1981 using FORTRAN program. This program translates the program to python for the convenience of further modifications and advancements.
print \
('\n\
|-------------------------DEFORM P1.0------------------------------|\n')
from time import localtime, strftime
current=strftime("%a, %d %b %Y %H:%M:%S", localtime())
print (current)
print \
('\n This Program calculates the settlement at any specified location \n\
on the surface of a multi-layered soil profile due to an array of \n\
rectangular loads on or below the surface. The program can also \n\
be used to calculate the stress distribution on any vertical profile.\n\
Copyright: URS Corporation\n\
March 2014\n\
Version description: converted from DEFORM1.FOR - Version 1.5\n')

#Functions
def corner(f,x,y,z,w,l):
    ww = w/2.
    ll = l/2.
    if (x-ww)<0:
      if (y<ll):
        a = ww+x
        b = ll+y
        i1 = f(a,b,z)
        a = ww-x
        i2 = f(a,b,z)
        b = ll-y
        i4 = f(a,b,z)
        a = ww+x
        i3 = f(a,b,z)
        return i1+i2+i3+i4
      elif (y>ll):
        a = ww+x
        b = ll+y
        i13 = f(a,b,z)
        a = ww-x
        i24 = f(a,b,z)
        b = y-ll
        i4 = f(a,b,z)
        a = ww+x
        i3 = f(a,b,z)
        return i13+i24-i3-i4
      else:
        a=ww+x
        i1=f(a,l,z)
        a=ww-x
        i2=f(a,l,z)
        return i1+i2
    elif (x-ww)==0:
      if (y<ll):
        b = y+ll
        i1 = f(w,b,z)
        b = ll-y
        i2 = f(w,b,z)
        return i1+i2
      elif (y>ll):
        b=ll+y
        i12=f(w,b,z)
        b=y-ll
        i2=f(w,b,z)
        return i12-i2
      else:
        return f(w,l,z)
    else:
      if (y>ll):
        a = ww+x
        b = ll+y
        i12 = f(a,b,z)
        b = ll-y
        i34 = f(a,b,z)
        a = x-ww
        i4 = f(a,b,z)
        b = ll+y 
        i2=f(a,b,z)
        return i12+i34-i2-i4
      elif (y<ll):
        a = ww+x
        b = ll+y
        i1234 = f(a,b,z)
        a = x-ww
        i24 = f(a,b,z)
        b = y-ll
        i4 = f(a,b,z)
        a = ww+x 
        i34 = f(a,b,z)
        return i1234-i24-i34+i4
      else:
        a = ww+x
        i12 = f(a,l,z)
        a = x-ww
        i2= f(a,l,z)
        return i12-i2

def fz(xz,yz,zz):
    if zz==0:
      return 0.25
    else:
      r1z = math.sqrt(xz*xz+zz*zz)
      r2z = math.sqrt(yz*yz+zz*zz)
      r3z = math.sqrt(xz*xz+yz*yz+zz*zz)
      az  = math.atan(xz*yz/zz/r3z)
      bz  = (xz*yz*zz/r3z)*(1./(r1z*r1z)+1./(r2z*r2z))
      return 0.5*(az+bz)/math.pi 

def fx(xx,yx,zx):
    if zx==0:
      return 0.25
    else:
      r1x = math.sqrt(xx*xx+zx*zx)
      r3x = math.sqrt(xx*xx+yx*yx+zx*zx)
      ax  = math.atan(xx*yx/zx/r3x)
      cx  = xx*yx*zx/r3x/r1x/r1x
      return 0.5*(ax-cx)/math.pi 

def fy(xy,yy,zy):
    if zy==0:
      return 0.25
    else:
      r2y = math.sqrt(yy*yy+zy*zy)
      r3y = math.sqrt(xy*xy+yy*yy+zy*zy)
      ay  = math.atan(xy*yy/zy/r3y)
      dy  = xy*yy*zy/r3y/r2y/r2y
      return 0.5*(ay-dy)/math.pi 

def localAxis(x,y,theta):
  xl = x*math.cos(theta*math.pi/180.) - y*math.sin(theta*math.pi/180.)
  yl = x*math.sin(theta*math.pi/180.) + y*math.cos(theta*math.pi/180.)
  return [xl,yl]

def localFactorC(theta):
  return math.cos(theta*math.pi/180)**2

def localFactorS(theta):
  return math.sin(theta*math.pi/180)**2



#Main Program
lines_list = sys.stdin.readlines()

calcType = int(lines_list[2].split()[0])
calcType -= 1
if calcType:
#TODO Calculate stress distribution
  print (calcType," Stress distribution code is under construction!")
else:
  iniVars=[]
  lCount = 5
  for i in range(3,lCount):
    for val in lines_list[i].split():
      iniVars.append(int(val))
  nSoilLayers = iniVars[0]
  nPrecon     = iniVars[1]
  soilArr     = [[0 for x in xrange(11)] for x in xrange(nSoilLayers)]
  lCount += nSoilLayers
  for nl in range(5,lCount):
    paraCount = 0
    for val in lines_list[nl].split():
      soilArr[nl-5][paraCount]=float(val)
      paraCount += 1

  for val in lines_list[lCount].split():
    yFinal = float(val)

  lCount += 1

  for val in lines_list[lCount].split():
    nLoadArea = int(val)

  areaArr     = [[0 for x in xrange(7)] for x in xrange(nLoadArea)]
  lCount     += 1
  for na in range(lCount, lCount+nLoadArea):
    paraCount = 0
    for val in lines_list[na].split():
      areaArr[na-lCount][paraCount]=float(val)
      paraCount += 1

  lCount     += nLoadArea
  for val in lines_list[lCount].split():
    nGrid     = int(val)

  varGrid = []
  matrixP = []
  lCount += 1
  if nGrid == 1:
    paraCount = 5
    for val in lines_list[lCount].split():
      varGrid.append(float(val))
    for x in np.arange(varGrid[0],varGrid[1]+varGrid[4],varGrid[4]):
      for y in np.arange(varGrid[2],varGrid[3]+varGrid[4],varGrid[4]):
        point=[]
        point.append(x)
        point.append(y)
        matrixP.append(point)

  else:
    paraCount = 2
    for np in range (lCount,len(lines_list)-1):
      point = []
      for val in lines_list[np].split():
        point.append(float(val))
      matrixP.append(point)

  hSL = []
  p0  = []
  pC  = []
  for nl in range(0,nSoilLayers):
    hSL.append(math.fabs(soilArr[nl][1]-soilArr[nl][0]))
    if nl==0:
      p0.append(hSL[nl]*0.5*soilArr[nl][2])
    else:
      p0.append(p0[nl-1]+hSL[nl-1]*0.5*soilArr[nl-1][2]+hSL[nl]*0.5*soilArr[nl][2])
    if (nPrecon==1):
      pC.append(soilArr[nl][5])
    else:
      pC.append(soilArr[nl][5]*p0[nl])

  resultPArr = []
  for point in matrixP:
   dSig = [[0 for x in xrange(3)] for x in xrange(nSoilLayers)]
   sb    = [0 for x in xrange(nSoilLayers)]
   pFin  = [0 for x in xrange(nSoilLayers)]
   cSet  = [0 for x in xrange(nSoilLayers)]
   eSet  = [0 for x in xrange(nSoilLayers)]
   sComp = [0 for x in xrange(nSoilLayers)]
   tCom  = [0 for x in xrange(nSoilLayers)]
   tSetC = 0
   tSetSC = 0
   tSetE = 0
   tSet  = 0
   resultP = []
   for na in range(0,nLoadArea):
     xl,yl=localAxis(point[0],point[1],areaArr[na][6])
     cxl,cyl=localAxis(areaArr[na][0],areaArr[na][1],areaArr[na][6])
     dx = math.fabs(xl-cxl)
     dy = math.fabs(yl-cyl)
     for nl in range(0,nSoilLayers):
       elevC = (soilArr[nl][0] + soilArr[nl][1])*0.5
       dz    = areaArr[na][2]-elevC
       if dz>=0:
         infX = corner(fx,dx,dy,dz,areaArr[na][3],areaArr[na][4])
         infY = corner(fy,dx,dy,dz,areaArr[na][3],areaArr[na][4])
         infZ = corner(fz,dx,dy,dz,areaArr[na][3],areaArr[na][4])
         dSig[nl][0]+=infX*areaArr[na][5]*localFactorC(areaArr[na][6])+infY*areaArr[na][5]*localFactorS(areaArr[na][6])
         dSig[nl][1]+=infY*areaArr[na][5]*localFactorC(areaArr[na][6])+infX*areaArr[na][5]*localFactorS(areaArr[na][6])
         dSig[nl][2]+=infZ*areaArr[na][5]

   for nl in range(0,nSoilLayers):
     if (dSig[nl][2]!=0):
       alpha = 0.5*(dSig[nl][0]+dSig[nl][1])/dSig[nl][2]
       if math.fabs(alpha)>1.:
         alpha = 1.
       sb[nl] = soilArr[nl][9]+alpha*(1.-soilArr[nl][9])
     pFin[nl] = p0[nl]+dSig[nl][2]
     if p0[nl] < pC[nl]:
       if pFin[nl] > pC[nl]:
         cSet[nl]=hSL[nl]*12.*(soilArr[nl][4]*math.log10(pC[nl]/p0[nl])+soilArr[nl][3]*math.log10(pFin[nl]/pC[nl]))
       else:
         cSet[nl]=hSL[nl]*12.*soilArr[nl][4]*math.log10(pFin[nl]/p0[nl])
     else:
       cSet[nl]=hSL[nl]*12.*soilArr[nl][3]*math.log10(pFin[nl]/p0[nl])
     eSet[nl] = hSL[nl]*12.*(dSig[nl][2]-soilArr[nl][7]*(dSig[nl][0]+dSig[nl][1]))/(soilArr[nl][6]*1000000.)
     cSet[nl] = cSet[nl]*sb[nl]
     sComp[nl] = hSL[nl]*12.*soilArr[nl][8]*math.log10(yFinal/soilArr[nl][10])
     tCom[nl] = cSet[nl]+eSet[nl]+sComp[nl]
     tSetC += cSet[nl]
     tSetSC += sComp[nl]
     tSetE += eSet[nl]
     tSet += tCom[nl]
     resultP.append(['{:^7.0f} {:^11.3f} {:^11.3f} {:^11.3f} {:^11.3f} {:^11.3f} {:^11.3f} {:^11.3f} {:^11.3f} {:^11.3f}'\
                    .format(nl+1, p0[nl],dSig[nl][2],dSig[nl][0],dSig[nl][1],sb[nl],cSet[nl],sComp[nl],eSet[nl],tCom[nl])])
     resultP.append('\n')
   resultP.append(['Total Compressive consolidation: ','{:>.3f}'.format(tSetC)])
   resultP.append(['Secondary Compression-20 Year: ','{:>.3f}'.format(tSetSC)])
   resultP.append(['Elastic Compression: ','{:>.3f}'.format(tSetE)])
   resultP.append(['Total Compression: ','{:>.3f}'.format(tSet)])
   resultPArr.append(resultP) 
  
  for x in range(len(resultPArr)):
    print ("Point: {:.3f}, {:.3f}\n".format(matrixP[x][0],matrixP[x][1]))
    print ("Soil    Initial   Vertical    Horizontal.X  Horizontal.Y    SBF      Consol.       Second.     Elastic      Total\n")
    print ("Layer   Stress  Stress delta  Stress delta  Stress delta              Comp.       Comp.20Yr      Comp.      Comp.\n")
    print ("         (psf)     (psf)        (psf)         (psf)                   (in)           (in)       (in)         (in)\n")
    for y in range(len(resultPArr[x])):
      print (' '.join(map(str,resultPArr[x][y])))
    print ('\n\n')

  print ('Stress distribution, Po, Pc, Pf')
  for nl in range(0,nSoilLayers):
    print ('No. {:3.0f}, P0: {:10.3f}, Pc: {:10.3f}, Pf: {:10.3f}\n'.format(nl+1,p0[nl],pC[nl],pFin[nl]))
  print ('\n\n')
