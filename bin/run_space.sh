#!/bin/bash

SPACE_HOME="/apps/dibig_tools/tenX/"
STEPS="1234"
LANE="*"
MASK="Y28N*,I10,N10,Y91N*"

#source /apps/dibig_tools/1.0/lib/sh/utils.sh

usage() {
    echo "run_space.sh - run spaceranger pipeline"
    echo
    echo "Usage: run_space.sh [options] params master"
    echo
    echo "Options should be specified before the required arguments. Options:"
    echo
    echo " -s S | Perform only the steps specified in the string S. Default: $STEPS"
    echo " -h   | Print this help message."
    echo
    echo "Steps:"
    echo 
    echo " 1 | Generate fastq files (spaceranger 'mkfastq' tool). Creates a directory"
    echo "     called Fastq and saves fastq files in it."
    echo " 2 | Generate gene counts (spaceranger 'count' tool). Creates one directory"
    echo "     for each sample."
    echo " 3 | Aggregate samples    (spaceranger 'aggr' tool). Writes aggregated data"
    echo "     to a folder named 'Aggregated'."
    echo " 4 | Generate final report. Creates the directory specified with the NAME variable"
    echo "     and saves HTML and .cloupe files to it, plus an index.html file linking"
    echo "     all output files in a table. Finally, zips this folder."
    echo
    echo "Params file:"
    echo
    echo "This file should contain variable assignments in bash format (var=value, no spaces),"
    echo "one per line. The following variables should be defined:"
    echo
    echo "  NAME (project name, used to name final report directory)"
    echo "  GENOME (directory containing 10X reference for this genome)"
    echo "  RUNDIR (Illumina run directory, containing the RunInfo.xml file)"
    echo "  MASK (base mask for bcl2fastq, optional, default: $MASK)"
    echo "  LANE (lane containing samples, optional, default: $LANE)"
    echo
    echo "Master file:"
    echo
    echo "This should be a comma-delimited file with one line for each sample, and the"
    echo "following columns:"
    echo
    echo "  Sample name"
    echo "  10X barcode (e.g. SI-TT-A1)"
    echo "  Path to slide image"
    echo "  Slide ID (e.g. V11F01-278)"
    echo "  Slide area"
    echo
}

# If submit is available, use it, otherwise fall back to sbatch
# Submit is here: 
command -v submit >/dev/null 2>&1
if [[ $? == 0 ]];
then
  SUBMIT=submit
else
  SUBMIT=sbatch
fi

# Process command-line options
while getopts "s:h" opt; do
    case $opt in
	s)
	    STEPS="$OPTARG"
	    ;;
	h)
	    usage
	    exit 0
	    ;;
    esac
done
shift $((OPTIND-1))
PARAMS=$1
MASTER=$2

if [[ "$PARAMS" == "" ]];
then
  usage
  exit 1
fi

write_samplesheet() {
  TMPIFS=$IFS
  IFS=,
  echo "Lane,Sample,Index"
  while read SMP IDX IMG SLIDE AREA;
  do
    echo "${LANE},${SMP},${IDX}"
  done < $MASTER
  IFS=$TMPIFS
}

run_mkfastq() {
    # Check that we have a valid run directory
  if [[ ! -f $RUNDIR/RunInfo.xml ]];
  then
    echo "Error: file RunInfo.xml not found - please specify path to a valid Illumina run directory with the -d option."
    exit 2
  fi

  # Call mkfastq
  echo "Starting step mkfastq."
  write_samplesheet > samplesheet.csv
  $SUBMIT -W ${SPACE_HOME}/scripts/space_mkfastq.qsub Fastq samplesheet.csv $RUNDIR $MASK
  echo "Step mkfastq terminated."
}

run_count() {
  echo "Starting step count."
  IFS=,
  nj=0
  while read SMP IDX IMG SLIDE AREA;
  do
    $SUBMIT -W ${SPACE_HOME}/scripts/space_count.qsub $SMP $IMG $GENOME $SLIDE $AREA Fastq --reorient-images &
    nj=$((nj+1))
  done < $MASTER
  wait
  echo "Step count terminated."
}

write_aggr() {
  TMPIFS=$IFS
  IFS=,
  echo "library_id,molecule_h5,cloupe_file"
  while read SMP IDX IMG SLIDE AREA;
  do
    echo "${SMP},${SMP}/outs/molecule_info.h5,${SMP}/outs/cloupe.cloupe"
  done < $MASTER
  IFS=$TMPIFS
}

run_aggr() {
  echo "Starting step aggr."
  write_aggr > aggr.csv
  $SUBMIT -W ${SPACE_HOME}/space_aggr.qsub Aggregated aggr.csv &
  wait
  echo "Step aggr terminated."
}

make_report() {
  mkdir -p $NAME
  HTML=${NAME}/index.html
  cat > ${HTML} <<EOF
<!DOCTYPE html>
<HTML>
<HEAD>
<TITLE>10X Visium analyis</TITLE>
<STYLE>
TABLE {
  width:60%;
  border: 2px solid black;
  border-collapse: collapse;
}
TR {
  border-top: 1px solid grey;
}
</STYLE>
</HEAD>
<BODY>
<CENTER>
<H1>10X Visium analyis</H1>
<TABLE>
<THEAD>
<TR><TH>Sample</TH><TH>Summary</TH><TH>Loupe<SUP>*</SUP></TH></TR>
</THEAD>
<TBODY>
EOF

  TMPIFS=$IFS
  IFS=,
  while read SMP IDX IMG SLIDE AREA;
  do
    cp -v ${SMP}/outs/cloupe.cloupe ${NAME}/${SMP}.cloupe
    cp -v ${SMP}/outs/web_summary.html ${NAME}/${SMP}.summary.html
    echo "<TR><TH>${SMP}</TH><TD><A href='${SMP}.summary.html'>${SMP}.summary.html</A></TD><TD><A href='${SMP}.cloupe'>${SMP}.cloupe</A></TD></TR>" >> ${HTML}
  done < $MASTER
  cp -v Aggregated/outs/cloupe.cloupe ${NAME}/Aggregated.cloupe
  cp -v Aggregated/outs/web_summary.html ${NAME}/Aggregated.summary.html
  echo "<TR><TH>Aggregated</TH><TD><A href='Aggregated.summary.html'>Aggregated.summary.html</A></TD><TD><A href='Aggregated.cloupe'>Aggregated.cloupe</A></TD></TR>" >> ${HTML}

  cat >> ${NAME}/index.html <<EOF
</TBODY>
</TABLE>
<BR>
<B>*</B> Loupe files should be opened with the <A href='https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest'>10X Loupe Browser</A>.
</CENTER>
</BODY>
</HTML>
EOF

  zip -r ${NAME}.zip ${NAME}
}

source $PARAMS

if [[ $STEPS == *1* ]]; then run_mkfastq; fi
if [[ $STEPS == *2* ]]; then run_count; fi
if [[ $STEPS == *3* ]]; then run_aggr; fi
if [[ $STEPS == *4* ]]; then make_report; fi

