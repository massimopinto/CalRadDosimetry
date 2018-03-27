import numpy as np
import matplotlib.pyplot as pyplt
import matplotlib
import math
from scipy import stats
import sys
from scipy.signal import butter, filtfilt, argrelextrema


dirname = "."
dir_img = "./img/"
dir_data = "./sample_data/"
PREFIX_LENGTH = 14 # To strip off the "dir_data" part of the filename
POSTFIX_LENGTH = -10

script = sys.argv[0]
filename = sys.argv[1]

Core_data = np.loadtxt(filename, skiprows=1)
#data = Core_data[-5000:-10][:]
data = Core_data[-10000:][:]
#data = Core_data[10000:][:]
# filtration of noise inspired by https://stackoverflow.com/questions/28536191/how-to-filter-smooth-with-scipy-numpy
# See also the documentation on filtfilt here
# https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.signal.filtfilt.html

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filtfilt(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    sm_y = filtfilt(b, a, data)
    return sm_y

# Parameters used for the filtration of noisy data
cutoff = 1500
fs = 50000

time = data[:,0]
resistance = data[:,1]

#Smoothing of resistance data using filtfilt.
y_smooth = butter_lowpass_filtfilt(resistance, cutoff, fs=60000)

# By finite differences, the first order derivative of the smoothed y axes for each mean value of x over your array is given by:
dy=np.diff(y_smooth,1) # note that the first derivative is calculated over y_smooth, *not* resistance
dx=np.diff(time,1)

dy_raw = np.diff(resistance,1)
dy_smooth_gradient=np.gradient(y_smooth,4)
dx_gradient=np.gradient(time,4)

#dx=np.gradient(time,1)
yfirst=dy/dx
yfirst_raw=dy_raw/dx
yfirst_gradient = dy_smooth_gradient/dx_gradient
yfirst_gradient = butter_lowpass_filtfilt(yfirst_gradient, cutoff=1500, fs=60000)

# [...] the corresponding values of x are 
xfirst=0.5*(time[:-1]+time[1:])

# Should one need to smooth the first derivative as well
yfirst_smooth = butter_lowpass_filtfilt(yfirst, cutoff=1500, fs=10000)

# For the second order, do the same process again
dyfirst=np.diff(yfirst,1)
dxfirst=np.diff(xfirst,1)
ysecond=dyfirst/dxfirst

dyfirst_gradient = np.gradient(yfirst_gradient,3)
dxfirst_gradient = np.gradient(time,3)
ysecond_gradient = dyfirst_gradient/dxfirst_gradient
ysecond_gradient = butter_lowpass_filtfilt(ysecond_gradient, cutoff=1500, fs=60000)

maxima = argrelextrema(ysecond_gradient, np.greater)
minima = argrelextrema(ysecond_gradient, np.less)

# the corresponding values of x are: 
xsecond=0.5*(xfirst[:-1]+xfirst[1:])

pyplt.figure(figsize=(10,5))
pyplt.plot(time, resistance)
#pyplt.show()
pyplt.savefig(dir_img+"glance_"+filename[PREFIX_LENGTH:POSTFIX_LENGTH]+".png")

# create a grid for the plots and arrangements of plots
f = pyplt.figure(figsize=(20,10))
gs = matplotlib.gridspec.GridSpec(6, 6, hspace=0.8, wspace=0.4)
f.suptitle("Analysis of "+filename, fontsize=14)

ax1 = f.add_subplot(gs[:2,:3])
ax2 = f.add_subplot(gs[2:4,:3])
ax3 = f.add_subplot(gs[4:6,:3])
ax4 = f.add_subplot(gs[:2,3:6])
ax5 = f.add_subplot(gs[2:4,3:6])
ax6 = f.add_subplot(gs[4:6,3:6])

#Populating the plots
#ax1: original data plus smoothing
ax1.plot(time,resistance)
ax1.set_ylabel("resistance [Ohm]")
ax1.set_title("Raw Data and filt.filt'ed")
ax1.plot(time,y_smooth, color='y', linewidth=3.0)

ax2.plot(xfirst,yfirst)
ax2.set_title("First Derivative on filtered data with np.diff")

ax3.plot(xsecond,ysecond)
ax3.set_xlabel("Time [s]")
ax3.set_title("Second Derivative on filtered data with np.diff(array,1)")

ax4.set_title("First Derivative on raw data with np.diff(array,1)")
ax4.plot(xfirst,yfirst_raw)

ax5.plot(time, yfirst_gradient)
ax5.set_title("(filtered) First Derivative on filtered data with np.gradient(array,2)")

ax6.set_xlabel("Time [s]")
ax6.plot(time, ysecond_gradient)
ax6.scatter(time[maxima],ysecond_gradient[maxima], marker='o', linestyle='None', s=40, c='red')
ax6.scatter(time[minima],ysecond_gradient[minima], marker='o', linestyle='None', s=40, c='black')

ax6.set_title("(filtered) Second Derivative on filtered data with np.gradient(array,2) and filt.filt'ered")

pyplt.savefig(dir_img+"/"+"derivate_filters_"+filename[PREFIX_LENGTH:POSTFIX_LENGTH]+".png")
# pyplt.show()

# ## Selection of "correct" maxima and minima in the second derivative
# At this point we have a series of maxima and minima, but only the **larger** maxima and minima, in module, are the ones poiting at the end and the beginning of the heating runs, respectively. The smaller maxima and minima in between are artefact of the filtering that was applied to the data. What you really need then is to make a selection of those maxima and minima that will point you to the heating runs.
# Let's first try with a selection based on the value of the maxima and the minima.

ymins = ysecond_gradient[minima] # subselection of ysecond_gradient
tmins = time[minima] # subselection of time
ymaxs = ysecond_gradient[maxima]
tmaxs = time[maxima]

# The threshold value might depend on the noise and this procedure may need generalization.
min_index = ysecond_gradient[minima] < -2.0e-5
max_index = ysecond_gradient[maxima] > 2.0e-5

t_beginnings = tmins[min_index]
t_ends = tmaxs[max_index]

pyplt.figure(figsize=(20,10))
pyplt.plot(time, ysecond_gradient)
pyplt.title("Selection of second derivative's maxima and minima in "+filename, fontsize=14)
pyplt.ylabel("Second derivative of resistance with respec to time")
pyplt.xlabel("time [seconds]")
pyplt.scatter(tmaxs[max_index],ymaxs[max_index], marker='o', s=100, c='red')
pyplt.scatter(tmins[min_index],ymins[min_index], marker='o', s=100, c='black')
pyplt.savefig(dir_img+"/"+"selection_of_maxima_"+filename[PREFIX_LENGTH:POSTFIX_LENGTH]+".png")
#pyplt.plot

# ## Identification of the beginning and the end of each run on the resistance plots
# The _minima_ of the second derivative, after appropriate polishing for the artefacts introduced by differentation of noisy data, are useful to identify the _beginning_ of a haeting run, while the _maxima_ correspond to the _end_ of the heating run. However maxima and minima of the second derivative are much smaller arrays than the original `time` and `resistance` arrays, so first one needs to locate the index values of these minima and maxima on the original arrays. 


time_mins = np.isin(time,t_beginnings) # finds the position of the real minima back in the orginal 'time' array
time_maxs = np.isin(time,t_ends)

pyplt.figure(figsize=(10,5))
pyplt.title("Localization of the heating runs in "+filename)
pyplt.ylabel("Resistance [Ohm]")
pyplt.xlabel("Time [seconds]")
pyplt.plot(time,resistance)
pyplt.plot(time,y_smooth, color='y', linewidth=2.0)
pyplt.scatter(time[time_mins],resistance[time_mins], marker ='o', s=100, c='black')
pyplt.scatter(time[time_maxs],resistance[time_maxs], marker='o', s=100, c='red')
pyplt.savefig(dir_img+"/"+"localizing_starts_ends_"+filename[PREFIX_LENGTH:POSTFIX_LENGTH]+".png")

# ## Linear regressions and extrapolations to the mid-run
# As common in graphite calorimetry for radiation dosimetry, the _relative resistance change_ during the heating run is determined via a _forward extrapolation_ of the pre-heating trend, and a _backwards extrapolation_ of the post-heating trend, both at _mid-run_. These trends are determined via linear fitting of the corresponding trend data over a time span of 120". The linear fitting data, the extrapolations, and further computations will be saved in a proper array (or list?)

assert len(t_beginnings) == len(t_ends), 'Number of minima and maxima is not equal. Cannot define heating runs.'
heating_times = t_ends - t_beginnings
pos = heating_times > 0 
assert pos.all(), 'Heating time returned negative in one heating cycle. Cannot calculate deltaR/R'

def find_nearest(array,value):
    '''Returns the index of the element of an array that is neares to a value'''
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        #return array[idx-1]
        return idx-1
    else:
        #return array[idx]
        return idx


def generate_fit_regions(array, fit_margin, span=120, direction='backwards'):
    '''Returns a new two-dimensional array with indices of where the linear fits will take place.
    
    Parameters
    ----------
    array : ndarray
        The array with all the times, intended in seconds
    fit_ends : ndarray
        The array containing all the indices where the fit ends (or starts)
    span : int
        The interval over which to carry out the linear regressions. Defaults to 2' (120").
    direction : string
        Specifies whether moving forward or backwards (default) on the array
    
    Returns
    -------
    out: ndarray
        Array of dtype=int with the indices defining the regions for (later) linear regressions
    '''
    assert len(fit_margin) > 0, 'Cannot define any fitting region'
    fit_array = np.empty([2,len(fit_margin)], dtype=int)
    if direction == 'forward':
        for i in range(0, len(fit_margin)):
            fit_array[0,i] = np.searchsorted(array,fit_margin[i])
            fit_array[1,i] = find_nearest(array, array[fit_array[0,i]]+span)
    elif direction == 'backwards':
        for i in range(0,len(fit_margin)):
            fit_array[1,i] = np.searchsorted(array,fit_margin[i])
            fit_array[0,i] = find_nearest(array, array[fit_array[1,i]]-span)
    return fit_array



fit_pre_heat = generate_fit_regions(time, t_beginnings, 120, 'backwards')
fit_post_heat = generate_fit_regions(time, t_ends, 120, 'forward')

pyplt.figure(figsize=(10,5))

pyplt.title("Localization of the fitting regions in "+filename)
pyplt.ylabel("Resistance [Ohm]")
pyplt.xlabel("Time [seconds]")
pyplt.plot(time,resistance, linewidth=1.0)
pyplt.plot(time,y_smooth, color='y', linewidth=3.0)
pyplt.scatter(time[time_mins],resistance[time_mins], marker ='o', s=100, c='black')
pyplt.scatter(time[fit_pre_heat[0,]],resistance[fit_pre_heat[0,]], marker ='o', s=100, c='black')
pyplt.scatter(time[time_maxs],resistance[time_maxs], marker='o', s=100, c='red')
pyplt.scatter(time[fit_post_heat[1,]],resistance[fit_post_heat[1,]], marker ='o', s=100, c='red')

pyplt.savefig(dir_img+"/"+"localize_fitting_regions_"+filename[PREFIX_LENGTH:POSTFIX_LENGTH]+".png")


# ### Organization of output data in columns 
# 
# | Heating starts (s) | Heating ends (s) | Midpoint (s) | delta_T (s) | pre-slope | pre-intercept | Resistance-pre | post-slope | post-intercept | Resistance-post | delta_R | R_mid_average | delta_R_over_R | slopes_control |
# | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---|
# | 1600 | 1720 | 1660 | 120 | 1.04e-4 | 5e-5| 10004.5 | bla | bla | 10004.22 | 0.255 | 10004.3 | 2.5 e-5 | 0.30367314 | 
# | 1900 | 2020 | 1960 | 120 | 1.05e-4 | 5.2e-5| 10004.2 | bla | bla | 10002.12 |0.245 | 10004.24 | 2.45 e-5 | 0.2772314 |



# Preparation of data for the linear fits and for their extrapolations.
x1 = fit_pre_heat[0,] # indices of the beginning of the linear fit, pre-heating drifts
x2 = np.searchsorted(time,time[time_mins]) # indices of the end of the linear fit, pre-heating drifts
x3 = np.searchsorted(time,time[time_maxs]) # indices of the beginning of the linear fit, post-heating drifts
x4 = fit_post_heat[1,] # indices of the end of the linear fit, post-heating drifts
delta_T = time[fit_post_heat[0,]] - time[fit_pre_heat[1,]] # heating run durations (seconds)
half_runs = time[fit_pre_heat[1,]] + 0.5*delta_T # midpoints of heating run (seconds) 
xmid = np.empty([len(half_runs)], dtype=int) # indices of the midpoints
for i in range(0,len(xmid)):
    xmid[i] = find_nearest(time, half_runs[i])

assert ((len(x1) == len(x2)) & (x1 < x2).all()), '''Number of linear fit extremes in pre-heating drifts 
is not equal or left limit(s) > right limit(s)'''
assert ((len(x3) == len(x4)) & (x3 < x4).all()), '''Number of linear fit extremes in post-heating drifts 
is not equal or left limit(s) > right limit(s)'''

fit_results_pre = np.empty([5,len(x1)], dtype=float)
fit_results_post= np.empty([5,len(x1)], dtype=float)
for i in range(0, len(x1)):
    fit_results_pre[:,i]  = stats.linregress(time[x1[i]:x2[i]],resistance[x1[i]:x2[i]])
    fit_results_post[:,i] = stats.linregress(time[x3[i]:x4[i]],resistance[x3[i]:x4[i]])

R_mid_pre  = fit_results_pre[0,:]*time[xmid]+fit_results_pre[1,:]
R_mid_post = fit_results_post[0,:]*time[xmid]+fit_results_post[1,:]

#Organize results in a numpy array    
results = np.vstack((t_beginnings, t_ends, half_runs, delta_T, fit_results_pre[0,:],fit_results_pre[1,:],R_mid_pre))    
results = np.vstack((results, fit_results_post[0,:], fit_results_post[1,:], R_mid_post))

delta_R = abs(results[6,]-results[9,])
R_mid_average = 0.5*abs(results[6,]+results[9,])
delta_R_over_R = delta_R / R_mid_average
slopes_control = abs(abs(results[4,]-abs(results[7,]))) / (0.5*(abs(results[4,]+abs(results[7,])))) # relative difference of the pre-drift and post-drift slopes. An indicator of thermal equilirbium

results = np.vstack((results, delta_R, R_mid_average, delta_R_over_R, slopes_control))

# create a list of numpy arrays of variable lenghts, each corresponding to a best fit line
line_pre = [] # empty
line_post = [] 
for i in range(0, len(x2)):
    line_pre.append(fit_results_pre[0,i]*time[x1[i]:xmid[i]] + fit_results_pre[1,i])
    line_post.append(fit_results_post[0,i]*time[xmid[i]:x4[i]] + fit_results_post[1,i])

results_fig = pyplt.figure(figsize=(20,10))
results_grid = matplotlib.gridspec.GridSpec(2, 1, hspace=0.2, wspace=0.1)
ax1 = results_fig.add_subplot(results_grid[0,0])
ax2 = results_fig.add_subplot(results_grid[1,0])

ax1.plot(time,resistance, linewidth=1.0)
for i in range(0, len(x2)):
    ax1.plot(time[x1[i]:xmid[i]], line_pre[i], color='black', linewidth=3.0)
    ax1.plot(time[xmid[i]:x4[i]], line_post[i], color='red', linewidth=3.0)
#ax1.plot(time,resistance)
ax1.set_ylabel("resistance [Ohm]")
ax1.set_title("Automated analysis of file "+filename, fontsize=16)

ax2.set_ylim(2e-5,2.8e-5)
ax2.scatter(time[time_mins],delta_R_over_R, s=200, c='black')
ax2.set_ylabel("deltaR / R [Ohm / Ohm]")
ax2.set_xlabel("time [seconds]")

pyplt.savefig(dir_img+"/"+"final_results_"+filename[PREFIX_LENGTH:POSTFIX_LENGTH]+".png")
