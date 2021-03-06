#
# Teem: Tools to process and visualize scientific data and images
# Copyright (C) 2009--2019  University of Chicago
# Copyright (C) 2008, 2007, 2006, 2005  Gordon Kindlmann
# Copyright (C) 2004, 2003, 2002, 2001, 2000, 1999, 1998  University of Utah
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# (LGPL) as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# The terms of redistributing and/or modifying this software also
# include exceptions to the LGPL that facilitate static linking.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#

#-----------------------------------------------------------------------------
#
# TeemConfig.cmake - Teem CMake configuration file for external projects.
#
# This file is configured by Teem and used by the TeemUse.cmake module
# to load Teem's settings for an external project.

# The directory of TeemConfig.cmake is, by definition, Teem_DIR.
# (this_dir == Teem_DIR)
#
get_filename_component(this_dir "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component(Teem_ROOT_DIR "${this_dir}/../../.." ABSOLUTE)

# CMake files required to build client applications that use Teem.
set(Teem_BUILD_SETTINGS_FILE "/scratch/shannon/c/aether/Projects/BOS/image-generation/analysis/src/teem-build/TeemBuildSettings.cmake")
set(Teem_USE_FILE "")

# The Teem directories.
set(Teem_EXECUTABLE_DIRS "${Teem_ROOT_DIR}/bin")
set(Teem_LIBRARY_DIRS "${Teem_ROOT_DIR}/lib")
set(Teem_INCLUDE_DIRS "${Teem_ROOT_DIR}/include")

# The Teem libraries.
set(Teem_LIBRARIES "teem")

# The C flags added by Teem to the cmake-configured flags.
set(Teem_REQUIRED_C_FLAGS "")

# The Teem version number
set(Teem_VERSION_MAJOR "1")
set(Teem_VERSION_MINOR "12")
set(Teem_VERSION_PATCH "0")

# Is Teem using shared libraries?
set(Teem_BUILD_SHARED_LIBS "ON")

# The list of tools in teem
set(Teem_TOOLS "")

# The Teem library dependencies.
if(NOT TARGET teem)
  include("${Teem_ROOT_DIR}/lib/cmake/teem/TeemLibraryDepends.cmake")
endif()
