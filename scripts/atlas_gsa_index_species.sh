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

###########################################################
# Expression Atlas specific code
DIRS="/nfs/production3/ma/home/atlas3-production/analysis/differential/microarray/experiments /nfs/production3/ma/home/atlas3-production/analysis/differential/rna-seq/experiments"
contrast_details="/nfs/public/ro/fg/atlas/experiments/contrastdetails.tsv"

###########################################################
species=$1

if [ "$species-" == "-" ]; then
    echo "Usage: gsa_index_species.sh species/organism"
    exit 1
fi

if [ ! -e $contrast_details  ]; then
    echo "$contrast_details not found"
    exit 1
fi
 
species_no_gaps=`echo $species|tr " " "_"`

rm -f $species_no_gaps.tsv.files.txt

echo "INFO: gathering files from $species_no_gaps"
touch $species_no_gaps.tsv.files.txt
#touch $species_no_gaps.xml.files.txt
grep "$species" /nfs/public/ro/fg/atlas/experiments/contrastdetails.tsv|cut -f 1|uniq| ( while  read n; do
											     for d in $DIRS; do
												 find $d/$n -maxdepth 2 -name "*analytics.tsv" -print  >> $species_no_gaps.tsv.files.txt 2>/dev/null;
											     done;
											 done)

echo "INFO: gathering files from $species_no_gaps...done"
echo "Files: $species_no_gaps.tsv.files.txt"

#
#mkdir -p db
SCRIPT=$(readlink -f $0)
`dirname $SCRIPT`/gsa_prepare_data.R -c 4 -i $species_no_gaps.tsv.files.txt  -x $species_no_gaps.tsv.files.txt -o $species_no_gaps.po

exit
