# =========================================================
#
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

function run_wrapper {
    local timit="/usr/bin/time -a "
    local tags=$1
    local logfile=$2
    shift 2
    if [ ! -e $logfile ]; then
	echo "Label Time Seconds Memory"|tr " " "\t" > $logfile
    fi
    $timit -o $logfile --format "$tags\t%E\t%e\t%M" bash -c "$*"
}
