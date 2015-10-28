#!/usr/bin/env Rscript
#; -*- mode: R;-*-
# =========================================================
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with iRAP.  If not, see <http://www.gnu.org/licenses/>.
#
#
# =========================================================

# Install all packages required
repo<-"http://www.stats.bris.ac.uk/R/"
install.packages("optparse",repos=repo)
install.packages("parallel",repos=repo)

q(status=0)
