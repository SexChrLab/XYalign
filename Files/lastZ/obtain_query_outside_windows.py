import os
import sys

def obtain_query_outside_windows(window, target, query):
	output={}
	for each_window in window:
		each_window_query_in_target = []
		each_window_query_out_target = []
		each_window_query_in = []
		each_window_query_out = []
		for index in range(len(target)):
			target_coord = target[index]
			if target_coord[0] >= each_window[0] and target_coord[1] <= each_window[1]:
				query_coord = query[index]
				if (query_coord[0] >= each_window[0]) and (query_coord[1] <= each_window[1]):
					each_window_query_in.append(query_coord)
					each_window_query_in_target.append(target_coord)
				else:
					each_window_query_out.append(query_coord)
					each_window_query_out_target.append(target_coord)
		output[each_window]=[each_window_query_out]
	for keys in output.keys():
		for i in output[keys]:
			for j in i:
				if j[0] < j[1]:
					print "chrY",keys[0],keys[1],j[0],j[1]
				if j[0] > j[1]:
					print "chrY",keys[0],keys[1],j[1],j[0]


def main():
	window = []
	target = []
	query = []
	with open(sys.argv[1]) as f1:
		for line1 in f1:
     			line1 = line1.split("\t")
			window.append((int(line1[1]), int(line1[2])))
	
	with open(sys.argv[2]) as f2:
		for line2 in f2:
			line2 = line2.split("\t")
			target.append((int(line2[0]), int(line2[1])))
	
	with open(sys.argv[3]) as f3:
		for line3 in f3:
			line3 = line3.split("\t")
			query.append((int(line3[0]), int(line3[1])))

	obtain_query_outside_windows(window, target, query)

main()	
	
