# pyGITRM: A Python interface for GITRM

### simTransSource

Given a native paraolid file, the source code and makefile required to convert to a simmetrix (smd) file are in this folder. Please read     
https://github.com/jguterl/pyGITRM/wiki/Simmetrix--to-Omega%E2%80%90h-mesh-conversion#basics

### Omega-h mesh conversion

The next step is converting the simmetrix mesh to an omega-h mesh. Please read    
https://github.com/jguterl/pyGITRM/wiki/Simmetrix--to-Omega%E2%80%90h-mesh-conversion#convert-simmetrix-sms-mesh-to-omega-h-sms

### Python script: write_fields.py   

Given the OEDGE data in a .npy format where the field information and the mesh information (the x and y co-ordinates of the four vertices of the quadrilateral) are in oder:    
a) It creates a netCDF4 file for the fields which GITRm understands    
b) It creates a lu.nc file which is essentially the OEDGE mesh information GITRm understands (This is required as fields are transferred directly from OEDGE to GITRm. without need to convert to a rectilinear grid). 


