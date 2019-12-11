from __future__ import print_function
def butter_bandpass(fs, lowcut=None, highcut=None, order=10):
    """
    Parameters
    ----------
    fs : float
        Sample Rate
    lowcut : float
        Low pass frequency cutoff
    highcut : float
        High pass frequency cutoff
    order : int
        Butterworth filter order [Default = 10, i.e., 10th order Butterworth filter]

    Reference
    ---------
    http://scipy-cookbook.readthedocs.io/items/ButterworthBandpass.html
    """
    from scipy.signal import butter
    nyq = 0.5 * fs
    if lowcut and highcut:
        high = highcut / nyq
        low = lowcut / nyq
        sos = butter(order, [low, high], analog=False, btype='band', output='sos')
    elif highcut:
        high = highcut / nyq
        sos = butter(order, high, btype='low', output='sos')
    elif lowcut:
        low = lowcut / nyq
        sos = butter(order, low, btype='high', output='sos')
    else:
        print('Error -- must supply lowcut, highcut or both')
        sos = None

    return sos

def butter_bandpass_filter(data, fs=None, lowcut=None, highcut=None, order=10):
    """
    Buterworth high, low or band pass filter.
    
    Uses scipy.signal.sosfiltfilt -- A forward-backward digital filter using cascaded 
    second-order sections.

    Parameters
    ----------
    data : array
        Data array
    fs : float
        Sample Rate
    lowcut : float
        Low pass frequency cutoff
    highcut : float
        High pass frequency cutoff
    order : int
        Butterworth filter order [Default = 10, i.e., 10th order Butterworth filter]

    Reference
    ---------
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.sosfiltfilt.html
    http://scipy-cookbook.readthedocs.io/items/ButterworthBandpass.html
    """
    from scipy.signal import sosfiltfilt
    sos = butter_bandpass(fs, lowcut=lowcut, highcut=highcut, order=order)
    return sosfiltfilt(sos, data)



