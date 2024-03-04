#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 10:00:35 2021

@author: Erin Barnhart
"""
import ResponseClassSimple
import numpy
from utility import *
import scipy.ndimage as ndimage
from pylab import *
from colour import Color
import scipy.stats as stats



def count_frames(filename,input_dict,threshold=1): 
	"""Reads in a stimulus output file and assigns an image frame number to each stimulus frame."""
	rows,header = read_csv_file(filename) #import stim file
	R = numpy.asarray(rows,dtype='float') # convert stim file list to an array
	
	"""replace intermediate values for AIN4 with 0 (if less than or equal to half the max value) 
	or the max value (if greater than half the max)"""
	max_v = numpy.max(R[:,-1])
	print(max_v)
	R[:,-1]=numpy.where(R[:,-1]<=max_v/2.,0,R[:,-1])
	R[:,-1]=numpy.where(R[:,-1]>max_v/2.,max_v,R[:,-1])

	#set up the output array
	output_array=numpy.zeros((R.shape[0],R.shape[1]+2))
	header.extend(['dt','frames']) 
	#calculate change in voltage signal
	vs = [0]
	vs.extend(R[1:,-1]-R[:-1,-1])
	#count image frames based on the change in voltage signal
	count_on = 0
	F_on = [0]
	count_off = 0
	F_off = [0]
	frame_labels = [0]
	n = 1
	while n<len(vs)-1:
		if vs[n]>vs[n-1] and vs[n]>vs[n+1] and vs[n] > threshold:
			count_on = count_on+1
			F_on.extend([count_on])
			F_off.extend([count_off])
		elif vs[n]<vs[n-1] and vs[n]<vs[n+1] and vs[n] < threshold*-1:
			count_off = count_off - 1
			F_off.extend([count_off])
			F_on.extend([F_on])
		else:
			F_on.extend([count_on])
			F_off.extend([count_off])
		frame_labels.extend([count_on*(count_on+count_off)])
		n=n+1
	frame_labels.extend([0])
	output_array[:,:R.shape[1]] = R
	output_array[:,-1] = frame_labels
	OAS = output_array[output_array[:,-1].argsort()]
	i1 = numpy.searchsorted(OAS[:,-1],1)
	OASc = OAS[i1:,:] #removes rows before image series start
	output_list = []
	n = 0
	frame=1
	same_frame = []
	while n<len(OASc):
		if int(OASc[n,-1])==frame:
			same_frame.append(list(OASc[n,:]))
		else:
			output_list.append(same_frame[0])
			same_frame = []
			same_frame.append(list(OASc[n,:]))
			frame=frame+1
		n=n+1
	output_list.append(same_frame[0])
	OL = numpy.asarray(output_list)
	gt = OL[:,int(input_dict['gt_index'])]
	dt = gt[1:]-gt[:-1]
	OL[1:,-2]=dt
	return OL,R,header

def find_dropped_frames(frames,time_interval,stim_data,stim_data_OG,gt_index):
	stim_frames = stim_data[-1,-1]
	print('number of image frames = ' +str(frames))
	print('number of stimulus frames = '+str(int(stim_frames)))
	if stim_frames != frames:
		print('uh oh!')
		target_T = frames * time_interval
		stim_T = numpy.sum(stim_data[:,-2])
		print('total time should be '+str(target_T)+' seconds')
		print('total time from stim file is '+str(stim_T)+' seconds')
		max_t_step = numpy.max(stim_data[:,-2])
		if numpy.round(max_t_step/time_interval)>=2:
			print('stimulus dropped at least one frame!')
			OUT = []
			num_df = 0
			for row in stim_data:
				if numpy.round(row[-2]/time_interval)>=2:
					num_df=num_df+1
					gt_dropped = row[gt_index]-time_interval
					stim_frame = numpy.searchsorted(stim_data_OG[:,gt_index],gt_dropped)
					print('frame dropped in row '+str(stim_frame)+' of original stim file (maybe)')
					OUT.append(stim_data_OG[stim_frame])
					OUT.append(row[:-2])
				else:
					OUT.append(row[:-2])
			print('found '+str(num_df)+' potential dropped frames')
		else:
			print('stim frames and image frames do not match, but no dropped frames found...double check the stim file')
	else:
		print('looks fine!')


def parse_stim_file(stim_info_array,frame_index = -1,rt_index = 1,st_index = None):
	"""Get frame numbers, global time, relative time per epoch, and stim_state (if it's in the stim_file)"""
	frames = stim_info_array[:,frame_index]
	rel_time = stim_info_array[:,rt_index]
	if st_index == 'None':
		stim_type = list(numpy.ones(len(frames),dtype = int))
	else:
		stim_type = stim_info_array[:,int(st_index)]
	return frames, rel_time, stim_type

def get_stim_position(stim_info_array,x_index = 3, y_index = 4):
	xpos = stim_info_array[:,x_index]
	ypos = stim_info_array[:,y_index]
	return xpos, ypos

def define_stim_state(rel_time,on_time,off_time):
	"""Define stimulus state (1 = ON; 0 = OFF) based on relative stimulus time."""
	stim_state = []
	for t in rel_time:
		if t>on_time and t<off_time:
			stim_state.extend([1])
		else:
			stim_state.extend([0])
	return stim_state

def segment_ROIs(mask_image):
	"""convert binary mask to labeled image"""
	labels = ndimage.measurements.label(mask_image)
	return labels[0]

def generate_ROI_mask(labels_image, ROI_int):
	return labels_image == ROI_int

def measure_ROI_fluorescence(image,mask):
	"""measure average fluorescence in an ROI"""
	masked_ROI = image * mask
	return numpy.sum(masked_ROI) / numpy.sum(mask)

def measure_ROI_ts(images,mask):
	out = []
	for image in images:
		out.append(measure_ROI_fluorescence(image,mask))
	return out

def measure_multiple_ROIs(images,mask_image):
	labels = segment_ROIs(mask_image)
	out = []
	num = []
	n = 1
	while n<=numpy.max(labels):
		mask = generate_ROI_mask(labels,n)
		out.append(measure_ROI_ts(images,mask))
		num.append(n)
		n=n+1
	return out,num,labels

def measure_one_ROI(images,mask_image):
	labels = segment_ROIs(mask_image)
	mask = labels >= 1
	out=[measure_ROI_ts(images,mask)]
	num=[1]
	return out,num,labels

def get_input_dict(row,header):
	input_dict = {}
	for r,h in zip(row,header):
		input_dict[h]=r
	return input_dict

def extract_response_objects(image_file,mask_file,stim_file,input_dict):
	"""inputs are file names for aligned images, binary mask, and unprocessed stimulus file
	outputs a list of response objects"""

	#read files
	I = read_tifs(image_file)
	mask = read_tifs(mask_file)

	#process stimulus file
	stim_data,stim_data_OG,header = count_frames(stim_file,input_dict)
	
	if (len(I))!=int(stim_data[-1][-1]):
		print("number of images does not match stimulus file")
		print('stimulus frames = ' + str(int(stim_data[-1][-1])))
		print('image frames = ' + str(len(I)))
		#stim_data = fix_dropped_frames(len(I),float(input_dict['time_interval']),stim_data,stim_data_OG,int(input_dict['gt_index']))

	print(I.shape)
	fmin = int(input_dict['fmin'])
	fmax = int(input_dict['fmax'])
	I = I[fmin:fmax]
	print(I.shape)
	stim_data_seg = stim_data[fmin:fmax]

	#get frames, relative time, stimulus type, and stimulus state from stim data
	fr,rt,st = parse_stim_file(stim_data_seg,rt_index = int(input_dict['rt_index']),st_index = input_dict['st_index'])
	ss = define_stim_state(rt,float(input_dict['on_time']),float(input_dict['off_time']))

	#get additional stim information for specific stim classes
	if input_dict['stim_class']=='mb':
		xpos,ypos = get_stim_position(stim_data_seg,x_index = int(input_dict['x_index']),y_index = int(input_dict['y_index']))

	#measure fluorscence intensities in each ROI
	if input_dict['one_ROI']=='TRUE':
		responses,num,labels = measure_one_ROI(I,mask)
	else:
		responses,num,labels = measure_multiple_ROIs(I,mask)

	print('number of ROIs = '+ str(numpy.max(num)))

	#load response objects
	response_objects = []
	for r,n in zip(responses,num):
		
		if input_dict['stim_class']=='mb':
			ro = ResponseClassSimple.movingBarResponse(F=r,stim_time = rt,stim_state = ss,ROI_num = n,stim_type = st,stim_x = xpos,stim_y = ypos)
		else:
			ro = ResponseClassSimple.Response(F=r,stim_time = rt,stim_state = ss,ROI_num = n,stim_type = st)

		ro.sample_name = input_dict['sample_name']
		ro.reporter_name = input_dict['reporter_name']
		ro.driver_name = input_dict['driver_name']
		ro.stimulus_name = input_dict['stimulus_name']
		ro.time_interval = float(input_dict['time_interval'])
		response_objects.append(ro)
	return response_objects, stim_data, header, labels

def extract_response_objects_from_arrays(F,stim_data,input_dict):
	"""inputs are numpy arrays containing raw fluorescence intensities (F) and stimulus data (stim_data)
        and a dictionary containing information about how to index the stim_data array (input_dict)"""
	#get frames, relative time, stimulus type, and stimulus state from stim data
	fr,rt,st = parse_stim_file(stim_data,rt_index = int(input_dict['rt_index']),st_index = input_dict['st_index'])
	ss = define_stim_state(rt,float(input_dict['on_time']),float(input_dict['off_time']))
	#load response objects
	response_objects = []
	for r,n in zip(F,numpy.arange(len(F))):
		ro = ResponseClassSimple.Response(F=r,stim_time = rt,stim_state = ss,ROI_num = n,stim_type = st)
		response_objects.append(ro)
	return response_objects

def smooth_responses(response_objects,input_dict):
	for ro in response_objects:
		ro.smooth(sigma = float(input_dict['sigma']))

def segment_individual_responses(response_objects,input_dict):
	for ro in response_objects:
		if input_dict['stim_class']=='mb':
			ro.segment_responses_and_bar_positions(int(input_dict['frames_before']),int(input_dict['frames_after']))
		else:
			ro.segment_responses(int(input_dict['frames_before']),int(input_dict['frames_after']))
		ro.measure_dff(int(input_dict['baseline_start']),int(input_dict['baseline_stop']))

def measure_average_dff(response_objects,input_dict):
	for ro in response_objects:
		ro.measure_average_dff(int(input_dict['epoch_length']))

def measure_average_stim_positions(response_objects,input_dict):
	for ro in response_objects:
		ro.average_stim_positions(int(input_dict['epoch_length']))

def map_RF_centers(response_objects):
	for ro in response_objects:
		ro.map_RF_center()

def align_mb_responses(response_objects,f_before = 20,f_after = 20):
	out = []
	for ro in response_objects:
		for a, st in zip(ro.average_dff,ro.stim_type_ind):
			max_index = numpy.argmax(a)
			if max_index>f_before and max_index+f_after<len(a):
				out.append([ro.ROI_num,st]+a[max_index-f_before:max_index+f_after])
	OUT = numpy.asarray(out)
	header = ['ROI','epoch_number']+list(numpy.arange(f_before+f_after))
	return OUT, header


def save_raw_responses_csv(response_objects,filename):
	OUT = []
	for ro in response_objects:
		out = [ro.sample_name,ro.driver_name,ro.reporter_name,ro.stimulus_name,ro.ROI_num]
		out.extend(ro.F)
		OUT.append(out)
	output_header = ['sample_name','driver_name','reporter_name','stimulus_name','ROI_num']
	output_header.extend(list(numpy.arange(0,ro.time_interval*len(ro.F),ro.time_interval)))
	write_csv(OUT,output_header,filename)

def save_individual_responses_csv(response_objects,filename):
	OUT = []
	for ro in response_objects:
		for st,dff in zip(ro.stim_type_ind,ro.dff):
			out = [ro.sample_name,ro.driver_name,ro.reporter_name,ro.stimulus_name,ro.ROI_num,st]
			out.extend(dff)
			OUT.append(out)
	output_header = ['sample_name','driver_name','reporter_name','stimulus_name','ROI_num','stim_type']
	output_header.extend(list(numpy.arange(0,ro.time_interval*len(dff),ro.time_interval)))
	write_csv(OUT,output_header,filename)

def save_average_responses_csv(response_objects,filename):
	OUT = []
	for ro in response_objects:
		for st,a in zip(ro.stim_type_ave,ro.average_dff):
			out = [ro.sample_name,ro.driver_name,ro.reporter_name,ro.stimulus_name,ro.ROI_num,st]
			out.extend(a)
			OUT.append(out)
	output_header = ['sample_name','driver_name','reporter_name','stimulus_name','ROI_num','stim_type']
	output_header.extend(list(numpy.arange(0,ro.time_interval*len(a),ro.time_interval)))
	write_csv(OUT,output_header,filename)

def save_RF_centers(response_objects,filename,threshold=0.1):
	OUT = []
	for ro in response_objects:
		if ro.RF_center[3]>=threshold:
			keep = 1
		else:
			keep = 0
		OUT.append([ro.sample_name,ro.driver_name,ro.reporter_name,ro.stimulus_name,ro.ROI_num]+ro.RF_center+[keep])
	output_header = ['sample_name','driver_name','reporter_name','stimulus_name','ROI_num','RF_center-x','RF_center-y','peak_response','min_peak']
	write_csv(OUT,output_header,filename)

def save_aligned_mb_responses(response_objects,filename):
	R, header = align_mb_responses(response_objects)
	write_csv(R,header,filename)

def plot_raw_responses(response_objects,filename):
	for ro in response_objects:
		plot(ro.F)
		savefig(filename+'-'+str(ro.ROI_num)+'-raw.png',dpi=300,bbox_inches = 'tight')
		clf()

def plot_average_responses(response_objects,filename):
	for ro in response_objects:
		for a in ro.average_dff:
			plot(a)
		savefig(filename+'-'+str(ro.ROI_num)+'-average.png',dpi=300,bbox_inches = 'tight')
		clf()

def plot_average_responses_by_stim_position(response_objects,filename):
	for ro in response_objects:
		for a,x,y in zip(ro.average_dff,ro.average_stim_x,ro.average_stim_y):
			if numpy.min(x)==numpy.max(x):
				plot(y,a)
			else:
				plot(x,a)
		savefig(filename+'-'+str(ro.ROI_num)+'-average_by_stim_position.png',dpi=300,bbox_inches = 'tight')
		clf()

def plot_RF_centers(response_objects,labels,filename,c1 = 'red',c2 = 'blue',screen_xlim=(-5,5),screen_ylim=(-4.5,8),threshold=0.1):
	color1 = Color(c1)
	color2 = Color(c2)
	color_range = list(color1.range_to(color2,len(response_objects)))
	
	D = []
	for ro in response_objects:
		D.append([ro.RF_center[0],ro.RF_center[1],ro.ROI_num,ro.RF_center[2],ro.RF_center[3]])
	D = numpy.asarray(D)
	D = D[numpy.argsort(D[:,1])] #sort by y pos
	#D = D[numpy.argsort(D[:,0])] #sort by x pos
	
	h,w = labels[0].shape
	im_rgb = np.zeros((h, w, 3), dtype=np.uint8)
	im_rgb[labels[0]==0,:] = [255,255,255]

	for d,c in zip(D,color_range):
		if d[4]>threshold:
			plot(d[1],d[0],'o',color = c.hex)
			im_rgb[labels[0]==d[2],:] = [c.red*255,c.green*255,c.blue*255]
		else: 
			im_rgb[labels[0]==d[2],:] = [200,200,200]

	ylim(screen_ylim)
	xlim(screen_xlim)
	ax = gca()
	ax.set_aspect(aspect = 1)
	savefig(filename+'.png',dpi=300,bbox_inches = 'tight')
	clf()

	save_tif(im_rgb,filename+'-labels.tif')

def plot_aligned_mb_responses(response_objects,filename):
	R, header = align_mb_responses(response_objects)
	for r in R:
		plot(r[2:],linewidth = 0.25,color = 'gray',alpha = 0.5)
	plot(numpy.average(R[:,2:],axis = 0),linewidth = 3)
	savefig(filename+'.png',dpi=300,bbox_inches = 'tight')
	clf()








