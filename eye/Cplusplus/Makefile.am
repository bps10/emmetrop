# Copyright (C) 2010-2011 Free Software Foundation, Inc
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

noinst_PROGRAMS = SchematicEye

SchematicEye_SOURCES = SchematicEye.cc
SchematicEye_LDADD = $C:/Cplusplus/goptical-1.0/goptical_core/src/libgoptical.la
#SchematicEye_CPPFLAGS = -I$C:/Cplusplus/goptical-1.0/goptical_core/src


AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

bin_PROGRAMS = SchematicEye
hello_SOURCES = src/hello.h src/hello.cc src/main.cc \
                src/foo/foo.h src/foo/foo.cc

dist_noinst_SCRIPTS = autogen.sh