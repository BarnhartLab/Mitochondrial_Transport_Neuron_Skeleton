import numpy
import csv

def read_csv_file(filename, header=True):
	data = []
	with open(filename, newline='') as csvfile:
		csvreader = csv.reader(csvfile, delimiter=',')
		for row in csvreader:
			data.append(row)
	if header==True:
		out_header = data[0]
		out = data[1:]
		return out, out_header
	else:
		return out

def write_csv(data,header,filename):
	with open(filename, "w") as f:
		writer= csv.writer(f)
		writer.writerow(header)
		for row in data:
			writer.writerow(row)

def get_input_dict(row,header):
	input_dict = {}
	for r,h in zip(row,header):
		input_dict[h]=r
	return input_dict

def parse_track_file(track_file,input_dict):
	#xy calibration factor
	cal = float(input_dict['pixel_calibration'])
	#load file
	rows,header = read_csv_file(track_file)
	#get individual tracks, frames, and track numbers
	all_tracks = []
	track = []
	all_frames = []
	frames = []
	track_numbers = []
	n = 0
	while n<len(rows)-1:
		f = rows[n][int(input_dict['frame_index'])]
		x = rows[n][int(input_dict['x_index'])]
		y = rows[n][int(input_dict['y_index'])]
		track_num1 = rows[n][int(input_dict['track_index'])]
		track_num2 = rows[n+1][int(input_dict['track_index'])]
		if track_num1==track_num2:
			track.append([float(x),float(y)])
			frames.append(float(f))
		else:
			track_numbers.append(int(track_num1))
			track.append([float(x),float(y)])
			all_tracks.append(numpy.asarray(track)*cal)
			track = []
			frames.append(float(f))
			all_frames.append(numpy.asarray(frames))
			frames = []
		n = n+1
	all_tracks.append(numpy.asarray(track)*cal)
	all_frames.append(numpy.asarray(frames))
	track_numbers.append(int(track_num1))
	return all_tracks, all_frames, track_numbers

def load_mito_lengths(shape_file,input_dict):
	shape_measurements, header = read_csv_file(shape_file)
	SM = numpy.asarray(shape_measurements)
	shape_numbers = numpy.asarray(SM[:,int(input_dict['shape_index'])],dtype = 'int')
	mito_lengths = numpy.asarray(SM[:,int(input_dict['length_index'])],dtype = 'float')
	cal = float(input_dict['shape_calibration'])
	return mito_lengths*cal,shape_numbers
