#
# Makefile for Cubic Spline interpolating polynomial
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# written by Conrad Shyu (shyu4751@yahoo.com)
#
# Author's comments:
# ------------------
# Department of Physics
# University of Idaho, Moscow, ID 83844
#
# revised on September 3, 2008
#
all: spline

spline:
	g++ -O3 -I. driver.cpp spline.cpp -o spline -lm

clean:
	rm -f spline
