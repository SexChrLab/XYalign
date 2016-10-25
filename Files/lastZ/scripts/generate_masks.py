import os
import sys

def generate_masks(window, target, query, xKb_buffer):
	output_target={}
	output_query={}
	for each_window in window:
		each_window_query_in_target = []
		each_window_query_out_target = []
		each_window_query_in = []
		each_window_query_out = []
		for index in range(len(target)):
			target_coord = target[index]
			if target_coord[0] >= each_window[0] and target_coord[1] <= each_window[1]:
				query_coord = query[index]
				window_range = (each_window[0] - xKb_buffer, each_window[1] + xKb_buffer)
				if (query_coord[0] >= window_range[0]) and (query_coord[1] <= window_range[1]):
					each_window_query_in.append(query_coord)
					each_window_query_in_target.append(target_coord)
				else:
					each_window_query_out.append(query_coord)
					each_window_query_out_target.append(target_coord)
		output_target[each_window]=each_window_query_out_target
		output_query[each_window]=each_window_query_out
	with open(sys.argv[5],'a') as target_outFile:
		target_outFile.write('chrY\n')
		for keys in output_target.keys():
			for i in output_target[keys]:
				target_outFile.write(str(i[0])+'\n')
				target_outFile.write(str(i[1])+'\n')
				target_outFile.write('NA\n')	
	with open(sys.argv[6],'a') as query_outFile:
		query_outFile.write('chrY\n')
		for keys in output_query.keys():
			for i in output_query[keys]:
				query_outFile.write(str(i[0])+'\n')
				query_outFile.write(str(i[1])+'\n')
				query_outFile.write('NA\n') 	

	with open(sys.argv[7],'a') as window_target_query_outFile:
		window_target_query_outFile.write(str("win_start"+'\t'+"win_end"+'\t'+"target_start"+'\t'+"target_end"+'\t'+"query_start"+'\t'+"query_end"+'\n'))
		for k in output_target:
			if len(output_target[k]) != 0:
				for index in range(0, len(output_target[k])):
					window_target_query_outFile.write(str(k[0])+'\t'+str(k[1])+'\t'+str(output_target[k][index][0])+'\t'+str(output_target[k][index][1])+'\t'+str(output_query[k][index][0])+'\t'+str(output_query[k][index][1])+'\n')
				
def main():

	""" From formatted lastZ output, return regions that are multimapped. 3 files are returned. File 1: coordinates for the target. File 2: coordinates for the query. File 1 and File 2 are in the rdotplot format to overlap the dotplot for visualization. File 3 contains the 10kb windows with the coordinates for the target and the query that are multimapped. Here, xKb_buffer is defined as the regions where if the query fell outside of that region, then it will be called multimapped. For example, if the window is (10000, 20000), and the target coordinates is (7000, 8000), and the query coordinates is (11000, 12000), then, this region will not be called multimapped if the buffer window is 50000. However, for the same window and target coordinates, if the query coordinates is (200000, 201000), then this region will be marked as "multimapped". This xKb buffer region can be specified by the user."""
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
	xKb_buffer = float(sys.argv[4])
	generate_masks(window, target, query, xKb_buffer)

main()	
	
