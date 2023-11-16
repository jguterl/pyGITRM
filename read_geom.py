import io, libconf
print("Hello Reading Geometry file")
with io.open("gitrGeometryFromGitrm.cfg") as f:
	config = libconf.load(f)
x1=config['geom']['x1']
print(x1)
