import numpy
import ResponseClassSimple
import ResponseTools
import utility
import os
from pylab import *

parent_dir = ''

input_csv = parent_dir+''
rows,header = ResponseTools.read_csv_file(input_csv)

for row in rows[:]:
	input_dict = ResponseTools.get_input_dict(row,header)
	print('analyzing sample ' + input_dict['sample_name'] + ' ' + input_dict['output_name'] + ' ' + input_dict['mask_name'] + ' ' + input_dict['stimulus_name'])

	sample_dir = parent_dir+input_dict['sample_name']
	image_dir = sample_dir+'/aligned_images/'
	mask_dir = sample_dir+'/masks/'
	stim_dir = sample_dir+'/stim_files/'
	plot_dir = sample_dir+'/plots/'
	output_dir = sample_dir+'/measurements/'

	if input_dict['aligned']=='TRUE':
		image_file = image_dir + input_dict['ch1_name']+'-aligned.tif'
	else:
		image_file = image_dir + input_dict['ch1_name']+'.tif'
	mask_file = mask_dir + input_dict['mask_name']+'.tif'
	stim_file = utility.get_file_names(stim_dir,file_type = 'csv',label = input_dict['stimulus_name'])[0]

	response_objects, stim_data, dataheader, labels = ResponseTools.extract_response_objects(image_file,mask_file,stim_file,input_dict)
	
	if input_dict['verbose']=='TRUE':
		parsed_stim_dir = stim_dir+'parsed/'
		if not os.path.exists(parsed_stim_dir):
			os.makedirs(parsed_stim_dir)
		utility.write_csv(stim_data,dataheader,parsed_stim_dir + 'parsed-'+input_dict['stimulus_name'] +'.csv')

	ResponseTools.smooth_responses(response_objects,input_dict)
	ResponseTools.save_raw_responses_csv(response_objects,output_dir+input_dict['output_name']+'-raw.csv')
	ResponseTools.plot_raw_responses(response_objects,plot_dir+input_dict['output_name'])
	print('extracted '+str(len(response_objects))+' response objects')
	
	ResponseTools.segment_individual_responses(response_objects,input_dict)
	ResponseTools.measure_average_dff(response_objects,input_dict)
	ResponseTools.save_individual_responses_csv(response_objects,output_dir+input_dict['output_name']+'-individual.csv')
	ResponseTools.save_average_responses_csv(response_objects,output_dir+input_dict['output_name']+'-average.csv')
	ResponseTools.plot_average_responses(response_objects,plot_dir+input_dict['output_name'])
	utility.save_tif(labels[0],mask_dir+input_dict['mask_name']+'-labels.tif')

	if input_dict['verbose']=='TRUE':
		parsed_stim_dir = stim_dir+'parsed/'
		if not os.path.exists(parsed_stim_dir):
			os.makedirs(parsed_stim_dir)
		utility.write_csv(stim_data,dataheader,parsed_stim_dir + 'parsed-'+input_dict['stimulus_name'] +'.csv')
