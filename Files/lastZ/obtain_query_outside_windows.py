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
	with open(sys.argv[1]) as window_file:
		for each_line in window_file:
     			each_line = each_line.split("\t")
			window.append((int(each_line[1]), int(each_line[2])))
	
	with open(sys.argv[2]) as target_file:
		for each_line in target_file:
			each_line = each_line.split("\t")
			target.append((int(each_line[0]), int(each_line[1])))
	
	with open(sys.argv[3]) as query_file:
		for each_line in query_file:
			each_line = each_line.split("\t")
			query.append((int(each_line[0]), int(each_line[1])))

	obtain_query_outside_windows(window, target, query)

main()	
	
