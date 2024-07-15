#!/bin/bash

set -e

SPACE_HOME="/apps/dibig_tools/tenX/"
STEPS="12345"
LANE="*"
MASK="Y28N*,I10,N10,Y91N*"
MODE=visium
#source /apps/dibig_tools/1.0/lib/sh/utils.sh

usage() {
cat <<EOF

run_space.sh - run spaceranger pipeline

Usage: run_space.sh [options] params master

Options should be specified before the required arguments. Options:

 -s S | Perform only the steps specified in the string S. Default: $STEPS
 -h   | Print this help message.

Steps:

 1 | Generate fastq files (spaceranger 'mkfastq' tool). Creates a directory
     called Fastq and saves fastq files in it.
 2 | Generate gene counts (spaceranger 'count' tool). Creates one directory
     for each sample.
 3 | Aggregate samples    (spaceranger 'aggr' tool). Writes aggregated data
     to a folder named 'Aggregated'.
 4 | Correct aggregate using Harmony.
 5 | Generate final report. Creates the directory specified with the NAME variable
     and saves HTML and .cloupe files to it, plus an index.html file linking
     all output files in a table. Finally, zips this folder.

Params file:

This file should contain variable assignments in bash format (var=value, no spaces),
one per line. The following variables should be defined:

  NAME (project name, used to name final report directory)
  GENOME (directory containing 10X reference for this genome)
  PANEL (file containing panel probles, optional)
  RUNDIR (Illumina run directory, containing the RunInfo.xml file)
  MASK (base mask for bcl2fastq, optional, default: $MASK)
  LANE (lane containing samples, optional, default: $LANE)
  MODE (visium, cytassist, or hd, default: $MODE)

Master file:

This should be a comma-delimited file with one line for each sample. For visium
mode, it should have the following columns:

  Sample name
  10X barcode (e.g. SI-TT-A1)
  Path to slide image
  Slide ID (e.g. V11F01-278)
  Slide area

For cytassist or hd mode:

  Sample name
  Path to HE image file
  Path to CytAssist image file
  Slide ID
  Slide area

EOF
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
    $SUBMIT -W ${SPACE_HOME}/scripts/space_count.qsub $PARAMS $SMP $IMG $SLIDE $AREA Fastq --reorient-images true &
    nj=$((nj+1))
  done < $MASTER
  wait
  echo "Step count terminated."
}

run_count_cyta() {
    echo "Starting step count (CytAssist)."
    IFS=,
    nj=0
#  while read SMP IDX IMG CYTIMG SLIDE AREA;
#  do
#    $SUBMIT -W ${SPACE_HOME}/scripts/space_count_cyta.qsub $PARAMS $SMP $IMG $CYTIMG $SLIDE $AREA --reorient-images true &
#    nj=$((nj+1))
#  done < $MASTER
#  wait

  # In HD mode cloupe files are links
  if [[ "$MODE" == "hd" ]];
  then
      while read SMP IDX IMG CYTIMG SLIDE AREA; do
	  echo "[${SMP}]"
	  pushd ${SMP}/outs/
	  cp $(readlink *.cloupe) cloupe.cloupe
	  popd
      done < $MASTER
  fi
  
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
  $SUBMIT -W ${SPACE_HOME}/scripts/space_aggr.qsub Aggregated aggr.csv &
  wait
  echo "Step aggr terminated."
}

run_harmony() {
    echo "Starting step harmony."
    mkdir -p harmony
    cut -f 1 -d , $MASTER | xargs $SUBMIT -W ${SPACE_HOME}/scripts/Rbatch.qsub ${SPACE_HOME}/R/harmony.R
    echo "Step harmony terminated."
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
  width:90%;
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
<H1>${NAME} - 10X Visium analyis</H1>
<TABLE>
<THEAD>
<TR>
  <TH>Sample</TH>
  <TH>Summary</TH>
  <TH>Loupe<SUP>*</SUP></TH>
  <TH>Spots<BR>under tissue</TH>
  <TH>Mean reads<BR>per spot</TH>
  <TH>Median genes<BR>per spot</TH>
  <TH>Fraction of spots<BR>under tissue</TH>
  <TH>Fraction reads in<BR>spots under tissue</TH>
  <TH>Number of reads</TH>
</TR>
</THEAD>
<TBODY>
EOF

  TMPIFS=$IFS
  IFS=,
  ALLMET=""
  case $MODE in
      visium|cytassist)
	  repflag=-v
	  ;;
      hd)
	  repflag=-d
	  ;;
  esac
  while read SMP IDX IMG SLIDE AREA;
  do
    met=${SMP}/outs/metrics_summary.csv
    ALLMET="${ALLMET} ${met}"
    metrics=$(${SPACE_HOME}/scripts/getMetrics.py $repflag $met)
    cp -v ${SMP}/outs/cloupe.cloupe ${NAME}/${SMP}.cloupe
    cp -v ${SMP}/outs/web_summary.html ${NAME}/${SMP}.summary.html
    echo "<TR><TH>${SMP}</TH><TD><A target='_blank' href='${SMP}.summary.html'>${SMP}.summary.html</A></TD><TD><A href='${SMP}.cloupe'>${SMP}.cloupe</A></TD>" >> ${HTML}
    echo $metrics >> ${HTML}
    echo "</TR>" >> ${HTML}
  done < $MASTER
  echo ${ALLMET} | xargs ${SPACE_HOME}/scripts/getMetrics.py -s $repflag > ${NAME}/metrics_summary.csv
  echo "<TR>" >> ${HTML}
  if [[ $STEPS == 3 ]]; then
      cp -v Aggregated/outs/cloupe.cloupe ${NAME}/Aggregated.cloupe
      cp -v Aggregated/outs/web_summary.html ${NAME}/Aggregated.summary.html
      echo "<TH>Aggregated</TH><TD><A href='Aggregated.summary.html' download='Y'>Aggregated.summary.html</A></TD><TD><A href='Aggregated.cloupe'>Aggregated.cloupe</A></TD>" >> ${HTML}
  else
      echo "<TD colspan=3>(aggregation not performed)</TD>" >> ${HTML}
  fi
  echo "<TD align='center' colspan='6'>Download <A href='${NAME}/metrics_summary.csv'>metrics summary file</A> (CSV format)</TD></TR>" >> ${HTML}

  cat >> ${HTML} <<EOF
</TBODY>
</TABLE>
<BR>
<B>*</B> Loupe files should be opened with the <A href='https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest'>10X Loupe Browser</A>.
EOF

  if [[ $STEPS == *4* ]]; then
      mkdir -p ${NAME}/harmony/
      cp harmony/* ${NAME}/harmony/
      cat >> ${HTML} <<EOF
<H2>Batch effect correction (Harmony)</H2>
<TABLE>
  <TBODY>
    <TR><TD colspan='2'><P>Batch effect correction was performed with <A href=''>Harmony</A>. The following plots show the distribution of samples before and after batch correction:</P></TD></TR>
    <TR>
      <TD align='center'><A href='harmony/pre-correction.png'><IMG width='400' src='harmony/pre-correction.png' /></A></TD>
      <TD align='center'><A href='harmony/corrected-ident.png'><IMG width='400' src='harmony/corrected-ident.png' /></A></TD>
    <TR><TD colspan='2'><P>Use the following links to download the corrected UMAP and cluster files:
<UL>
<LI><A href='harmony/corrected_umap.csv' download>corrected_umap.csv</A></LI>
<LI><A href='harmony/corrected_clusters.csv' download>corrected_clusters.csv</A></LI>
</UL>
These files should be imported into the Loupe Browser using the <A href='https://www.10xgenomics.com/resources/analysis-guides/correcting-batch-effects-in-visium-data#import-csv-in-the-loupe-browser' target='_blank'>procedure described here</A>.</P></TD>
    </TR>
  </TBODY>
</TABLE>
EOF
  fi

cat >> ${NAME}/index.html <<EOF
</CENTER>
</BODY>
</HTML>
EOF

  zip -r ${NAME}.zip ${NAME}
}

source $PARAMS

if [[ $STEPS == *1* ]]; then run_mkfastq; fi
if [[ $STEPS == *2* ]]; then
    case $MODE in
	cytassist|hd)
	    run_count_cyta
	    ;;
	*)
	    run_count
	    ;;
    esac
fi
if [[ $STEPS == *3* ]]; then run_aggr; fi
if [[ $STEPS == *4* ]]; then run_harmony; fi
if [[ $STEPS == *5* ]]; then make_report; fi

