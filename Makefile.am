## Makefile.am -- Process this file with automake to produce Makefile.in
## Copyright (C) 2000 Gary V. Vaughan
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

MAINTAINERCLEANFILES 	= Makefile.in aclocal.m4 configure

INCLUDES		= -I$(MEGAPOLY_VECMATH_INCDIR) \
				-I$(MEGAPOLY_UTIL_INCDIR) \
				-I$(top_srcdir)

bin_PROGRAMS	= mesh2slim

mesh2slim_SOURCES	= mesh2slim.cxx \
			PIEM.cxx PIEM.hxx \
			PIEMHandler.hxx \
			SlimSimplify.cxx SlimSimplify.hxx \
			SlimBall.hxx \
			SlimTree.hxx \
			SlimTreeIO.hxx \
			SlimTreeApp.hxx \
			MeshR.hxx MeshRApp.hxx \
			MeshRDG.hxx \
			SMFRIO.hxx RIO.hxx \
			SVD.cpp SVD.h \
			jacobi.h nrutil.h
mesh2slim_LDADD	= -lm
