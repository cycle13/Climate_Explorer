
# Shell script to run the idl code makeERApentads_FEB2016.pro with a read in value
# run using ./run_makeERApentads_FEB2016.pro

# Work out what the name of the script is
echo "Command: $0 $*"
battorun=$1
echo "bat file to run: $battorun"

progtorun=$2
echo "Program to run: $progtorun"

param=$3
echo "Param to run: $param"

# run the IDL code 
/opt/ukmo/idl/ukmo/bin/tidl << EOF
waveoff
@$battorun 
$progtorun,'$param'
exit
EOF
