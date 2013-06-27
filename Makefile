### Write a proper make file

### Build eye. Then cython eye. Then install python egg?

all: cython

cython: SchematicEye/eye SchematicEye/eye.so
	cd SchematicEye/ && make python && \
	echo "leaving SchematicEye" &&\
	cd ../
