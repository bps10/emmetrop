### Write a proper make file

### Build eye. Then cython eye. Then install python egg?

all: cython

cython: SchematicEye/eye SchematicEye/eye.so
	cd SchematicEye/ && git pull origin master && \
	make python && \
	echo "leaving SchematicEye" &&\
	cd ../
