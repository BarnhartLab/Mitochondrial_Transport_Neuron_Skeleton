#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 10:09:32 2021

@author: Erin Barnhart
"""
import numpy as numpy
import scipy as scipy

class Response(object):
    """class attributes"""
    
    """instance attributes"""
    _instance_data = {'sample_name':None,
                       'ROI_num':None,
                       'reporter_name':None,
                       'driver_name':None,
                       'F':[],
                       'stimulus_name':None,
                       'stim_time':[],
                       'stim_state':[],
                       'stim_type':[],
                       'time_step':None,
                       'units':'seconds',
                       'smoothed':False}
    
    def __init__(self, **kws):
        
        """self: is the instance, __init__ takes the instace 'self' and populates it with variables"""
        
        for attr, value in self._instance_data.items():
            if attr in kws:
                value = kws[attr]
                setattr(self,attr,value)
            else:
                setattr(self,attr,value)
        
    """instance methods"""
    
    
    """find median fluorescence over time"""
    def med(self):
        return np.median(self.F)
        #self.median.append(median)

    """smooth raw fluoresence"""
    def smooth(self,sigma = 1):
        smoothed = scipy.ndimage.gaussian_filter1d(self.F,sigma)
        self.F = smoothed
        self.smoothed = True
        return smoothed


    """identify stim switch points"""
    def stimswitch(self):
        ON_indices = list(numpy.where(numpy.diff(self.stim_state)==1)[0]+1)
        OFF_indices = list(numpy.where(numpy.diff(self.stim_state)==-1)[0]+1)
        if ON_indices[0]<OFF_indices[0]:
            return ON_indices, OFF_indices
        else:
            return ON_indices,OFF_indices[1:]
    
    """break data into epochs"""
    def segment_responses(self,frames_before,frames_after):
        """points_before is the number of points before ON
        points_after is the number of points after OFF"""
        ON,OFF = self.stimswitch()
        r = []
        st_ind = []
        stim_type_ind = []
        for on, off in zip(ON,OFF):
            start = on - frames_before
            stop = off + frames_after
            if start>0 and stop<len(self.F)+1:
                r.append(self.F[start:stop])
                #print(len(self.F[start:stop]))
                st_ind.append(int(self.stim_type[on]))
        self.individual_responses = r
        self.stim_type_ind = st_ind
        return r, st_ind

    """get df/f"""
    def measure_dff(self,baseline_start,baseline_stop):
        #b = self.baseline_end
        ir = numpy.asarray(self.individual_responses)
        dff = []
        for i in ir:
            baseline = numpy.mean(i[baseline_start:baseline_stop])
            dff.append((list(i-baseline)/baseline))
        self.dff = dff
        return dff
    
    def measure_average_dff(self,epoch_length):
        #b = self.baseline_end
        A = []
        stim_type = 1
        while stim_type <= numpy.max(self.stim_type_ind):
            r = []
            for dff,st in zip(self.dff,self.stim_type_ind):
                if st == stim_type:
                    r.append(dff[:epoch_length])
            R = numpy.asarray(r)
            #print(R)
            A.append(list(numpy.average(R,axis = 0)))
            stim_type = stim_type +1
        self.average_dff = A
        #print(list(numpy.arange(1,numpy.max(self.stim_type_ind)+1)))
        self.stim_type_ave = list(numpy.arange(1,numpy.max(self.stim_type_ind)+1))
        return A

    def measure_stdev_dff(self,epoch_length):
        STDEV = []
        stim_type = 1
        while stim_type <= numpy.max(self.stim_type_ind):
            r = []
            for dff,st in zip(self.dff,self.stim_type_ind):
                if st == stim_type:
                    r.append(dff[:epoch_length])
            R = numpy.asarray(r)
            STDEV.append(list(numpy.std(R,axis = 0)))
            stim_type = stim_type +1
        self.stdev_dff = STDEV
        return STDEV

    def integrated_response(self,start_point,end_point):
        IR = []
        for aDFF in self.average_dff:
            IR.append(numpy.sum(aDFF[start_point:end_point]))
        return IR

class movingBarResponse(Response):

    _instance_data = {'sample_name':None,
                       'ROI_num':None,
                       'reporter_name':None,
                       'driver_name':None,
                       'F':[],
                       'stimulus_name':None,
                       'stim_time':[],
                       'stim_state':[],
                       'stim_type':[],
                       'time_step':None,
                       'units':'seconds',
                       'smoothed':False,
                       'stim_x':[],
                       'stim_y':[],
                       'RF_center':[]}

    def __init__(self, **kws):
        
        """self: is the instance, __init__ takes the instace 'self' and populates it with variables"""
        super().__init__(**kws)
        
        for attr, value in self._instance_data.items():
            if attr in kws:
                value = kws[attr]
                setattr(self,attr,value)
            else:
                setattr(self,attr,value)

    """break data into epochs, include bar positions"""
    def segment_responses_and_bar_positions(self,frames_before,frames_after):
        """points_before is the number of points before ON
        points_after is the number of points after OFF"""
        ON,OFF = self.stimswitch()
        r = []
        x = []
        y = []
        st_ind = []
        stim_type_ind = []
        for on, off in zip(ON,OFF):
            start = on - frames_before
            stop = off + frames_after
            if start>0 and stop<len(self.F)+1:
                r.append(self.F[start:stop])
                x.append(self.stim_x[start:stop])
                y.append(self.stim_y[start:stop])
                #print(len(self.F[start:stop]))
                st_ind.append(int(self.stim_type[on]))
        self.individual_responses = r
        self.stim_x_ind = x
        self.stim_y_ind = y
        self.stim_type_ind = st_ind
        return r, x, y, st_ind

    """get average bar positions"""
    def average_stim_positions(self,epoch_length):
        AX = []
        AY = []
        stim_type = 1
        while stim_type <= numpy.max(self.stim_type_ind):
            all_x = []
            all_y = []
            for x,y,st in zip(self.stim_x_ind,self.stim_y_ind,self.stim_type_ind):
                if st == stim_type:
                    all_x.append(x[:epoch_length])
                    all_y.append(y[:epoch_length])
            X = numpy.asarray(all_x)
            Y = numpy.asarray(all_y)
            AX.append(list(numpy.average(X,axis = 0)))
            AY.append(list(numpy.average(Y,axis = 0)))
            stim_type = stim_type +1
        self.average_stim_x = AX
        self.average_stim_y = AY
        return AX, AY

    """find the receptive field center (x and y coordinates) based on response amplitudes"""
    def map_RF_center(self):
        y_max = []
        x_max = []
        max_amp = []
        for a,x,y in zip(self.average_dff,self.average_stim_x,self.average_stim_y):
            max_amp.append([numpy.max(a)])
            max_index = numpy.argmax(a)
            if numpy.min(x)==numpy.max(x):
                #x coodinate is not changing, so bar is moving in y
                y_max.append(y[max_index])
            else:
                x_max.append(x[max_index])
        center_x = numpy.average(x_max)
        center_y = numpy.average(y_max)
        peak = numpy.average(max_amp)
        min_peak = numpy.min(max_amp)
        self.RF_center = [center_x,center_y,peak,min_peak]
        return center_x, center_y, peak,min_peak
