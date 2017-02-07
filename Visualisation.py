# -*- coding: utf-8 -*-
"""
Created on Sun Dec  4 14:00:08 2016

@author: reyansh
"""

import numpy as np
import pygame
import matplotlib.pyplot as plt
import time
import sys

print("Visualisation Program started.")


pygame.init() 

myfont = pygame.font.SysFont("Roboto", 20)

initial_time = time.time()

number = np.loadtxt('numberofparticles.txt', usecols=range(1))
number = int(number)

print('Number of particles in simulation = ' + str(number))



start_time = time.time()

X = []
for n in range(1, 2*number + 1):
    X.append(np.loadtxt('position.txt', usecols=range(n, n+1)))

X = np.array(X)

iters = np.size(X[0])

time_now = (time.time() - start_time)
print ('Position read in ' + str(time_now) + ' seconds.')


start_time = time.time()

Speed = []

for n in range(1, number+1):
    Speed.append(np.loadtxt('speed.txt', usecols=range(n, n+1)))

Speed = np.array(Speed)

time_now = (time.time() - start_time)
print ('Speed read in ' + str(time_now) + ' seconds.')


start_time = time.time()

Timeticks = []

for i in range(int(iters/50)):
    Timeticks.append(50*i)
    
Timeticks = np.array(Timeticks)

Timerange = np.linspace(0, 700, 1000)

for j in np.nditer(Timeticks):
    
    Hist = []
    
    for i in range(number):
        Hist.append(Speed[i][j])
    
    Hist = np.array(Hist)
    path = "/home/reyansh/Desktop/ExamDayPresentation/Regular/HistogramData/ "
    path+= str(int(j/50))    
    path += ".png"
    plt.figure()
    plt.hist(Hist, 40, normed = True)
    plt.xlim(0, 700)
    plt.ylim = (0, 0.010)
    ms_speed = sum(Hist*Hist/number)
    yrange = (2 * Timerange/ms_speed) * np.exp(-(Timerange ** 2)/ms_speed)
    plt.plot(Timerange, yrange, label = 'Maxwell-Boltzmann Distribution', color = 'red', lw = 3)
    plt.xlabel('Velocity (pixel s$^{-1}$)')
    plt.ylabel('Fraction of particles with velocity $v$' )    
    plt.legend(loc = "best")
    plt.savefig(path)
    plt.close()

time_now = (time.time() - start_time)
print ('Histograms created in ' + str(time_now) + ' seconds.')

Speed = Speed.reshape(np.size(Speed))


Speed = Speed*255/(max(Speed)) 


size = 5

start_time = time.time()

Pressure = []

Pressure.append((np.loadtxt('pressure.txt', usecols=range(1,2))))

Pressure = np.array(Pressure)
Pressure = Pressure.reshape(np.size(Pressure))

time_now = (time.time() - start_time)
print ('Pressure read in ' + str(time_now) + 'seconds.')


Even = []
Odd = []

for n in range(2*number):
    if (n%2 == 0):
        Even.append(n)
    else:
        Odd.append(n)

Even = np.array(Even)
Odd = np.array(Odd)


x_border_left = 5 
y_border_down = 5

width = 800 
height = 750 

Xcor = []
for n in np.nditer(Even):
    Xcor.append(X[n])

Ycor = []
for n in np.nditer(Odd):
    Ycor.append(X[n])

Xcor = np.array(Xcor)
Ycor = np.array(Ycor)


Xcor = Xcor.reshape(np.size(Xcor))
Ycor = Ycor.reshape(np.size(Ycor))

t = 0

clock = pygame.time.Clock()

pygame.display.set_caption('Hard disc simulation of gas')
screen=pygame.display.set_mode([1366, 768])

screen.fill([255,255,255])

for j in range(number):
    pygame.draw.circle(screen,[0,0,255],[int(Xcor[0 + j*iters]), int(Ycor[0 + j*iters])],size,size)

pygame.draw.rect(screen, [0, 0, 0], [x_border_left, y_border_down, width, height])
pygame.display.flip()

done = False

print("All done in total of " + str(time.time() - initial_time) + ' seconds.')

while not done:
    for event in pygame.event.get():
        if event.type==pygame.QUIT:
            done = True

    if (t >= iters):
        t = 0
    
    else:

        screen.fill([255,255,255])
        pygame.draw.rect(screen, [0,0,0], [x_border_left, y_border_down, width, height])
        img = pygame.image.load("/home/reyansh/Desktop/ExamDayPresentation/Regular/HistogramData/ " + str(int(t/50)) + ".png")
        img = pygame.transform.scale(img, (500, 400))        
        screen.blit(img, (900, 0))
        textsurface = myfont.render('Pressure = ' + str(Pressure[t]), True, (15,15,15))
        screen.blit(textsurface,(900, 400))
        textsurface = myfont.render('Volume = ' + str(width*height), True, (15,15,15))
        screen.blit(textsurface,(900, 450))
        textsurface = myfont.render('kT = ' + str(ms_speed/2), True, (15,15,15))
        screen.blit(textsurface,(900, 500))
        textsurface = myfont.render('Number of particles = ' + str(number), True, (15,15,15))
        screen.blit(textsurface,(900, 550))
        
        for j in range(number):
                                
                if (Speed[t + j*iters] <= int(0.25*255)):
                    pygame.draw.circle(screen,[255,0,0],[int(Xcor[t + j*iters]), int(Ycor[t + j*iters])],size,size)
                elif (Speed[t + j*iters] >int(0.25*255) and Speed[t + j*iters] <= int(0.5*255)):
                    pygame.draw.circle(screen,[255,246,2],[int(Xcor[t + j*iters]), int(Ycor[t + j*iters])],size,size)
                elif (Speed[t + j*iters] >int(0.5*255) and Speed[t + j*iters] <= int(0.75*255)):
                    pygame.draw.circle(screen,[72,162,234],[int(Xcor[t + j*iters]), int(Ycor[t + j*iters])],size,size)
                else:
                    pygame.draw.circle(screen,[255,255,255],[int(Xcor[t + j*iters]), int(Ycor[t + j*iters])],size,size)
                    
        

  
    pygame.display.flip()
    t += 1
    clock.tick(100)


pygame.quit()
sys.exit()
