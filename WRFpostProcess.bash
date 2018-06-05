#!/bin/bash
set -e

usage="Usage: bash postProcess.bash inputfile outputfile"

: ${NCCOPY:=/home/mikarant/github/cozoc/netcdf/bin/nccopy}

input_file=$1
output_file=$2

#input_file=~/wrfout_d01_0001-01-01_00\:00\:00_INTRP.nc
#output_file=~/temp1.nc

# Extract the wanted variables from the input file
IGNORE_ATT_COORDINATES=1 cdo exprf,myexpr $input_file $output_file


# Extract the time axis from the input file (rather crudely, but didn't figure out better way)
ncks -A -v XTIME $input_file xtime.tmp.nc
ncks -A xtime.tmp.nc $output_file
rm -rf xtime.tmp*

# Remove unnecessary attributes from the variables
ncatted -O -a ,,d,, $output_file

# Add units and long names to the variables
arr=(
LEV  "Pa"       "Pressure levels"
U    "m s**-1"  "U velocity"
V    "m s**-1"  "V velocity"
T    "K"        "Temperature"
Z    "m"        "Geopotential height"
SP   "hPa"      "Surface pressure"
Q    "K s**-1"  "Diabatic heating"
FU   "m s**-2"  "Frictional U-tendency"
FV   "m s**-2"  "Frictional U-tendency"
)

cmd="ncatted -O"
i=0
while (( i < ${#arr[@]} )); do
    cmd="$cmd -a units,${arr[$i]},c,c,\"${arr[$((i+1))]}\" -a long_name,${arr[$i]},c,c,\"${arr[$((i+2))]}\""
    (( i += 3 ))
done
cmd="$cmd $output_file"

echo $cmd

eval $cmd


# Remove the long history-attribute
ncatted -O -a history,global,d,, $output_file
ncatted -O -a history_of_appended_files,global,d,, $output_file

# Rename dimensions
ncrename -d vlevs,lev -d west_east,x -d south_north,y ${output_file}

# Convert the output file to netcdf-4 format
$NCCOPY -k nc4 $output_file ${output_file}4
rm -rf $output_file

