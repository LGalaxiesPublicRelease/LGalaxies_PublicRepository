# numpy dtype for SFHbins
import numpy as np

sfh_struct_dtype = np.dtype([
('Snapnum',np.int32,1),
('Bin',np.int32,1),
('LookbackTime',np.float64,1),
('dt',np.float64,1),
('Nbin',np.int32,1)
])

