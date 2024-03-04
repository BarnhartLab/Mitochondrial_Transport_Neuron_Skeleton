import numpy
import MitoTrackTools
import MitoTrackClass

label = 'primary'
parent_dir = '//'+label+'/'
output_dir = parent_dir+'measurements/'
input_csv = parent_dir+'inputs_'+label+'.csv'
rows,header = MitoTrackTools.read_csv_file(input_csv)

output_header = ['sample_name','mito_ID','time_interval','total_track_distance','total_track_time','average_speed','average_pause_free_speed','arrest_rate','number_of_stops','direction','mito_length']
OUT = []
output_header_byFly = ['sample_name','number_of_mitos','average_speed','average_pause_free_speed','average_number_of_stops','arrest_rate']
OUT_byFly = []

for row in rows[:]:
	fly = []
	input_dict = MitoTrackTools.get_input_dict(row,header)
	print('.....')
	print('analyzing sample '+input_dict['sample_name'])

	#get mito centroids for each mito track
	track_filename = parent_dir+input_dict['sample_name']+'/measurements/'+input_dict['track_file']+'.csv'
	print(track_filename)
	tracks,frames,track_numbers = MitoTrackTools.parse_track_file(track_filename,input_dict)
	print('found '+str(len(tracks))+' mito tracks')

	#get mito shape measurements, if they exist
	if input_dict['shapes_exist']=='TRUE':
		shape_filename = parent_dir+input_dict['sample_name']+'/measurements/'+input_dict['shape_file']+'.csv'
		lengths,shape_numbers = MitoTrackTools.load_mito_lengths(shape_filename,input_dict)
		print('found '+str(len(lengths))+' length measurements')
		if len(tracks)!=len(lengths):
			print('the number of tracks does not equal the number of length measurements!')
			shape = False
		else:
			shape = True
	else:
		shape = False
		print('there are no shape measurements for this sample')

	#generate mito track objects
	mitotracks = []
	if shape:
		for tr,f,tn,l,sn in zip(tracks,frames,track_numbers,lengths,shape_numbers):
			if tn!=sn:
				print('track and shape ID numbers do not match!!!')
				print('track ID number = '+str(tn))
				print('shape ID number ='+str(sn))
			mt = MitoTrackClass.MitoTrack(mito_number = tn, frames = f,centroids = tr,mito_length = l,t = float(input_dict['t']))
			mitotracks.append(mt)
	else:
		for tr,f,tn in zip(tracks,frames,track_numbers):
			mt = MitoTrackClass.MitoTrack(mito_number = tn, frames = f,centroids = tr,t = float(input_dict['t']))
			mitotracks.append(mt)

	for mt in mitotracks:
		min_pause = int(numpy.round(4/mt.t))
		mt_measurements = [input_dict['sample_name'],mt.mito_number,mt.t,mt.total_track_distance(),mt.total_track_time(),mt.average_speed(),mt.average_pause_free_speed(),mt.arrest_rate(min_pause = min_pause),numpy.sum(mt.stops(min_pause=min_pause)),mt.direction(),mt.mito_length]
		if not numpy.isnan(mt.average_pause_free_speed()):
			fly.append([mt.average_speed(),mt.average_pause_free_speed(),mt.arrest_rate(min_pause = min_pause)])
			OUT.append(mt_measurements)
	F = numpy.asarray(fly)
	FA = numpy.average(F,axis = 0)

	OUT_byFly.append([input_dict['sample_name'],len(F),FA[0],FA[1],FA[2]])


MitoTrackTools.write_csv(OUT,output_header,output_dir+'speed_and_shape_measurements-'+label+'.csv')
MitoTrackTools.write_csv(OUT_byFly,output_header_byFly,output_dir+'speed_measurements-'+label+'-byFly.csv')





	
