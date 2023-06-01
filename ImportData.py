from random import *
from math import *
import numpy as np

from CONSTANT import *

seed(random_seed)

def Import_data(H, data, n, Rs, Qmax, change):
	global n_default, Rs_default, Qmax_default

	

	file = ["bacgiang", "hanoi", "lamdong", "sonla", "thaibinh"]
	Dataset = [[[n+n_step*i,Rs, Qmax] for i in range(dataset_num)], [[n, Rs + Rs_step*i, Qmax] for i in range(dataset_num)], [[n, Rs, Qmax + Qmax_step*i] for i in range(dataset_num)]]

	file = file[data-1]

	with open(f"Data\\{file}.asc", "r") as f:
		f.readline()
		f.readline()
		xllcorner = float(f.readline()[9:-1])
		yllcorner = float(f.readline()[9:-1])
		cellsize = float(f.readline()[8:-1])
		NODATA_value = f.readline()
		data_asc = f.readlines()
		data_asc[0] = data_asc[0][13:]
		data_asc[0] = list(map(float, data_asc[0].split()))
		for i in range(1, len(data_asc)):
			data_asc[i] = list(map(float, data_asc[i].split()))
			data_asc[i-1].append(data_asc[i].pop(0))
		data_asc.pop()
		cell = int(H//25)
		data_asc = data_asc[-cell:]
		for i in range(len(data_asc)):
			data_asc[i] = data_asc[i][:cell]

	Dataset = Dataset[change-1]

	Targets, Qs = place_random(Dataset, data_asc)

	return Dataset, Targets, Qs, file, change - 1

def place_random(Dataset, data_asc):
	Targets = []
	Qs = []
	for j in range(len(Dataset)):
		n = Dataset[j][0]
		Qmax = Dataset[j][2]
		T = []
		for k in range(n):
			x, y = random()*H, random()*H
			z = data_asc[int(x//25)][int(y//25)]
			T.append([x,y,z])

		Qs.append([randint(1,Qmax) for _ in range(n)])
		Targets.append(T.copy())

	return Targets, Qs

def exportData(average_S, average_runtimeCov, Dataset, file, name, H, change):
	n = Dataset[0][0]
	Rs = Dataset[0][1]
	Q = Dataset[0][2]
	changes = ["n", "R", 'Qmax'] 

	with open(f"Result\\{name}\\{file}\\change {changes[change]} n{n} Rs{Rs} Q{Q} H{H} data.txt", "w") as f:
		if change == 0:
			f.write("changing n\n")
			for i in range(dataset_num):
				string = f'n = {Dataset[i][0]}, s1-1-{i+1}, {average_S["n"][i]}, {average_runtimeCov["n"][i]}\n'
				f.write(string)
		if change == 1:
			f.write("changing R\n")
			for i in range(dataset_num):
				string = f'R = {Dataset[i][1]}, s1-2-{i+1}, {average_S["R"][i]}, {average_runtimeCov["R"][i]}\n'
				f.write(string)
		if change == 2:
			f.write("changing Qmax\n")
			for i in range(dataset_num):
				string = f'Qmax = {Dataset[i][2]}, s1-3-{i+1}, {average_S["Q"][i]}, {average_runtimeCov["Q"][i]}\n'
				f.write(string)
