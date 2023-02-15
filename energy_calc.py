#C:\\Users\\Jarek\\Documents\\Formula_e_enegry_usage\\energy_calc.py
import sys
import numpy as np 
import matplotlib.pyplot as plt 
import re
import functools
import operator
def Lagrange (x, distance , speed):
	#coefficients = [0, 0, 0]
	k = len(distance)
	def basis(j):
		p = [(x - distance[m])/(distance[j] - distance[m]) for m in range(k) if m != j]
		return functools.reduce(operator.mul, p)
		assert len(distance) != 0 and (len(speed) == len(distance))
	return sum(basis(j)*speed[j] for j in range(k))
def coef_from_lagrange_quadratic(distance,speed):
	if len(distance) != 3 :
		print ("not quadratic")
	coef = []
	coef.append(Lagrange(0, distance, speed))
	coef.append(((Lagrange(2, distance, speed) - Lagrange(0, distance, speed)) - 4*(Lagrange(1, distance, speed)- Lagrange(0, distance, speed)) )/(-2))
	coef.append(((Lagrange(2, distance, speed) - Lagrange(0, distance, speed)) - 2*(Lagrange(1, distance, speed) - Lagrange(0, distance, speed)) )/2)
	return coef
def roots_of_quadratic (a, b ,c):
	discriminant = b*b - 4*a*c
	if discriminant < 0:
		print ("lines do not intersect")
	elif discriminant == 0:
		return (-b + np.sqrt(discriminant))/(2*c)
	elif discriminant > 0:
		x1 = (-b + np.sqrt(discriminant))/(2*c)
		x2 = (-b - np.sqrt(discriminant))/(2*c)
		return [x1, x2]



car_data = []
i = 0 
file = open('C:\\Users\\Jarek\\Documents\\Formula_e_enegry_usage\\car.txt', 'r')
for line in file.readlines():
	input_data = re.findall(r"[-+]?(?:\d*\.*\d+)", line)
	car_data.append(eval(input_data[0])) 
	i += 1

mass   = car_data[0]
max_power  = car_data[1]
cd     = car_data[2]
niu    = car_data[3]
length = car_data[4]
v0     = car_data[5]
grip   = car_data[6]
dt     = 0.05
eff    = 0.95
# intitalizing arrays for plots
max_lap_time = 0
lap_time_cuttoff = 0.7
energy_graph_cutoff = 0
n1 = 4
n2 = 40
time_vs_energy = [[[] for _ in range(3)] for _ in range(n1)]
print (time_vs_energy)
        # trying to solve with Lagrange interpolation 
# mapping distance vs speed for breaking + calculating energy usage/recovery
distance_brk = [length]
speed_brk    = [v0]
energy_brk   = [0]
k = 0
while distance_brk[k] > 0 :
	speed_brk.append(speed_brk[k] + (grip*9.81*mass + cd*speed_brk[k]*speed_brk[k]/12.96)*dt*3.6/(mass))
	distance_brk.append(distance_brk[k] - (speed_brk[k + 1]+speed_brk[k])*dt/7.2)
	energy_brk.append(eff*min(max_power*3600/speed_brk[k] , grip*9.81*mass)*(speed_brk[k+1]+speed_brk[k])*dt/7.2 + energy_brk[k])
	k += 1

# mapping distance vs speed for coasting backwards from end 
distance_coast_end = [length]
speed_coast_end    = [v0 * 1.1]
l = 0
while distance_coast_end[l] > 0 :
	speed_coast_end.append(speed_coast_end[l] + (niu*9.81*mass + cd*speed_coast_end[l]*speed_coast_end[l]/12.96)*dt*3.6/(mass))
	distance_coast_end.append(distance_coast_end[l] - (speed_coast_end[l + 1]+speed_coast_end[l])*dt/7.2)
	l += 1
print (distance_coast_end[l-1], speed_coast_end[l-1])
for a in range(n1):
	power = max_power*(1 - 0.2*a)
	distance_acc = [0]
	speed_acc    = [v0]
	energy_acc   = [0]
	k = 0
	print(a)
# mapping acceleration
	while distance_acc[k] < length :
		speed_acc.append(speed_acc[k] + (min(power*3600/speed_acc[k] , grip*9.81*mass) - cd*speed_acc[k]*speed_acc[k]/12.96 - niu*mass*9.81)*dt*3.6/(mass))
		distance_acc.append(distance_acc[k] + (speed_acc[k+1]+speed_acc[k])*dt/7.2)
		energy_acc.append(min(power*3600/speed_acc[k] , grip*9.81*mass)*(speed_acc[k+1]+speed_acc[k])*dt/(7.2*eff) + energy_acc[k])
		k += 1
	#plt.plot(distance_acc,speed_acc)
	#plt.plot(distance_brk,speed_brk)
#plt.plot(distance_acc , speed_acc)
#plt.plot(distance_brk, speed_brk)
#plt.plot(distance_acc , energy_acc)
#plt.plot(distance_brk, energy_brk)
#plt.show()
# searching for the intersection of acceleration and breaking
	i = 1 
	k = 1 
	while True:
		while speed_brk[k] > speed_acc[i]:
			if distance_acc[i+1] + (length - distance_brk[k+1]) >= length:
				break
			i += 1
		while speed_brk[k] < speed_acc[i]:
			if distance_acc[i+1] + (length - distance_brk[k+1]) >= length:
				break
			k += 1
		if distance_acc[i+1] + (length - distance_brk[k+1]) >= length:
			break
	#print (speed_brk[k], speed_acc[i], distance_brk[k], distance_acc[i])
			# searching for the intersection of acceleration and coasting
	l = 1 
	i2 = 1 
	while distance_acc[i2] < distance_coast_end[l]:
		while speed_coast_end[l] > speed_acc[i2]:
			if distance_acc[i2] > distance_coast_end[l]:
				break
			i2 += 1
			#print (speed_coast_end[l] , speed_acc[i2],distance_coast_end[l] ,distance_acc[i2])
			#print (l , i2)
		while speed_coast_end[l] < speed_acc[i2]:
			if distance_acc[i2] > distance_coast_end[l]:
				break	
			l += 1
			#print (speed_coast_end[l] , speed_acc[i2], distance_coast_end[l] ,distance_acc[i2])
			#print (l , i2)
	i2 -= 1 
	                                                                 #### In my country there is a problem. And this problem is in the block

# calculating lift and coast 
	for b in range(n2):
		speed_coast   = [speed_acc[i - int((i-i2)*(b/n2))]]
		distance_coast = [distance_acc[i - int((i-i2)*(b/n2))]]
		#print (speed_coast, distance_coast ,'lift off point')
		n = 0
		while distance_coast[n] < length and speed_coast[n] > v0:
			#print (distance_coast[n], speed_coast[n])
			speed_coast.append(speed_coast[n] - (cd*speed_coast[n]*speed_coast[n]/12.96 + niu*mass*9.81)*dt*3.6/(mass))
			distance_coast.append(distance_coast[n] + (speed_coast[n + 1]+speed_coast[n])*dt/7.2)
			n += 1
		n = len(distance_coast)
		coast_coef = coef_from_lagrange_quadratic([distance_coast[0], distance_coast[n//2], distance_coast[n - 1]] , [speed_coast[0], speed_coast[n//2], speed_coast[n-1]])
		brk_coef = coef_from_lagrange_quadratic([distance_brk[k], distance_brk[k//2], distance_brk[0]] , [speed_brk[k], speed_brk[k//2], speed_brk[0]])
		x1_x2 = roots_of_quadratic (coast_coef[0] - brk_coef[0], coast_coef[1] - brk_coef[1], coast_coef[2] - brk_coef[2])
		x = max(x1_x2[0], x1_x2[1]) 
		xy = [x,coast_coef[0] + coast_coef[1]*x + coast_coef[2]*x*x]
		#print(xy , 'breaking point')
		#plt.plot(distance_coast,speed_coast)
		for x in range(len(distance_brk)) :
			if distance_brk[x + 1] < xy[0]: break
		energy_consump = energy_acc[i - int((i-i2)*(b/n2))] - energy_brk[x]
		for y in range(len(distance_coast)):
			if distance_coast[y + 1] > xy[0]: break
			#print(y)
		lap_time  = dt*(i - int((i-i2)*(b/n2)) + x + y)
		time_vs_energy[a][0].append(lap_time)
		time_vs_energy[a][1].append(energy_consump)
		time_vs_energy[a][2].append(y*dt)
		max_lap_time = max(max_lap_time, lap_time)
	#plt.show()
plt.figure(0)
for a in range(n1):
	plt.plot(time_vs_energy[a][0], time_vs_energy[a][1])
	for x in range(len(time_vs_energy[a][0])):
		if time_vs_energy[a][0][x] > max_lap_time*lap_time_cuttoff and time_vs_energy[a][1][x] > energy_graph_cutoff : 
			energy_graph_cutoff = time_vs_energy[a][1][x]
	#plt.plot(time_vs_energy[a][1], time_vs_energy[a][2], 'r-')
plt.xlabel('lap time, s')
plt.ylabel('Energy consumption , J', color='b')
plt.xlim(xmin = max_lap_time*lap_time_cuttoff)
plt.ylim(ymax = energy_graph_cutoff*1.1)
plt.figure(1)
for a in range(n1):
	plt.plot(time_vs_energy[a][0], time_vs_energy[a][2], 'b-')
	#plt.plot(time_vs_energy[a][1], time_vs_energy[a][2], 'r-')
plt.xlabel('Lap time, s')
plt.ylabel('Lift and coast time', color='r')
plt.show()
#print (lap_time, energy_consump)
