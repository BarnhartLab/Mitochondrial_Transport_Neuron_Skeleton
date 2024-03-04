import numpy as numpy
from scipy import signal

class MitoTrack(object):

    """instance attributes"""
    """centroids is an Nx2 numpy array containing x and y coordinates of the mito centroid over time.
    't' is the length of time that elapsed between measurements of the mito centroid position.
    'normal' is a 1x2 numpy array that defines direction of positive displacement."""

    _instance_data = {'sample_name':None,
                       'mito_number':None,
                       'branch_type':None,
                       'branch_radius':None,
                       'xy_units':None,
                       't_units':None,
                       'frames':[],
                       'centroids':[],
                       'mito_length':None,
                       't':1,
                       'normal':[0,-1],
                       'stimulus':[]}

    
    def __init__(self, **kws):
        
        """self: is the instance, __init__ takes the instace 'self' and populates it with variables"""
        
        for attr, value in self._instance_data.items():
            if attr in kws:
                value = kws[attr]
                setattr(self,attr,value)
            else:
                setattr(self,attr,value)
        
    """instance methods"""
    def total_track_distance(self):
        displacement_vectors = self.centroids[1:] - self.centroids[:-1]
        distances = numpy.sqrt((displacement_vectors ** 2).sum(axis=1))
        return numpy.sum(distances)

    def total_track_time(self):
        return len(self.centroids)*self.t

    def velocity(self):
        """Calculate the velocity of the mito at every time point."""
        displacement_vectors = self.centroids[1:] - self.centroids[:-1]
        distances = numpy.sqrt((displacement_vectors ** 2).sum(axis=1))
        if self.normal:
            dot_products = (self.normal * displacement_vectors).sum(axis=1)
            return numpy.sign(dot_products) * distances / self.t
        else:
            print('normal vector required to calculate velocity')

    def average_velocity(self,first=0,last=-1):
        """Calculate the average velocity for a mito track."""
        return numpy.average(self.velocity()[first:last])

    def direction(self):
        """Get mito direction based on average velocity."""
        if self.average_velocity()>0:
            return 'anterograde'
        if self.average_velocity()<0:
            return 'retrograde'
        else:
            return 'no net displacement'

    def speed(self):
        """Calculate the speed of the mito at every time point."""
        displacement_vectors = self.centroids[1:] - self.centroids[:-1]
        distances = numpy.sqrt((displacement_vectors ** 2).sum(axis=1))
        #print(distances)
        return distances/self.t

    def average_speed(self, first = 0, last = -1):
        """Calculate the average speed for a mito track."""
        return numpy.average(self.speed()[first:last])

    def median_speed(self, first = 0, last = -1):
        """Calculte the median speed for a mito track."""
        return numpy.median(self.speed()[first:last])

    def pause_free_speed(self,threshold = 0.1):
        """get instant speeds above threshold (ie pause-free speeds)"""
        out = []
        speed = self.speed()
        mito_state = self.mito_state(threshold=threshold)
        for s,ms in zip(speed,mito_state):
            if ms==1:
                out.append(s)
        return(out)

    def average_pause_free_speed(self,threshold=0.1,first = 0, last = -1):
        return numpy.average(self.pause_free_speed(threshold=threshold)[first:last])


    def mito_state(self, threshold=0.1):
        """Classify the mito as moving [1] or stationary [0] at each time point"""
        out = []
        for x in self.speed():
            if x >= threshold:
                out.extend([1])
            else:
                out.extend([0])
        return out

    def stops(self,threshold = 0.1, min_pause = 1):
        """Identify time points where mito arrests motility."""
        mito_state = self.mito_state(threshold=threshold)
        out = []
        n = 1
        while n<len(mito_state)-min_pause:
            if mito_state[n-1]==1 and numpy.sum(mito_state[n:n+min_pause])==0:
                out.extend([1])
            else:
                out.extend([0])
            n = n+1
        return out

    def arrest_rate(self,threshold = 0.1,min_pause = 1):
        """Calculate the mitochondrila arrest rate (stops per unit time)."""
        total_time = numpy.sum(self.mito_state(threshold=threshold))*self.t
        total_stops = numpy.sum(self.stops(threshold=threshold,min_pause=min_pause))
        return total_stops/total_time

    



