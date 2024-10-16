import os

from collections import OrderedDict

from astropy.io import fits

log_headers = {
    'date-beg': 'DATE-BEG',
    'effective_exposure_time': 'EXPTIMEE',
    'object_name': 'OBJNAME',
    'object_type': 'OBJTYPE',
    'ra': 'RA',
    'dec': 'DEC',
    'filter_yj': 'FILTER1',
    'filter_hk': 'FILTER2',
    'filter_aux': 'FILTER3',
    'aperture': 'FILTER4',
}



def append_log(log_name, ):
    if not os.path.exists(log_name):
        with open(log_name, 'w') as f:
            f.write(','.join(log_headers.keys()))

