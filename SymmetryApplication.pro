TEMPLATE = subdirs

CONFIG += \
	warn on

SUBDIRS += Util \
		SimpleMesh \
		Harris3D \
		RotationalSymmetry \
		SymArch

Harris3D.depends	+= Util SimpleMesh
RotationalSymmetry	+= Util SimpleMesh Harris3D
SymArch			+= RotationalSymmetry
