#!/bin/bash

# Change this if scripts are moved somewhere else!
TENX_HOME=/apps/dibig_tools/tenX

# Specify default versions of cellranger and R
# Note: these should start with "/".
export R_VERSION="/4.0"

STEPS=1235
TENX_RUN_DIR=.
TENX_RUN_ID=Fastq
REPORT=Report
AGGRDIR=Aggregated
# Organism ID
TENX_ORG=GRCh38
# Base mask
BM="Y28N*,I10,I10,Y91N*"
# Chemistry
CHEM=""
# Differential analysis
DIFF=""

# If submit is available, use it, otherwise fall back to sbatch
# Submit is here: 
command -v submit >/dev/null 2>&1
if [[ $? == 0 ]];
then
  SUBMIT=submit
else
  SUBMIT=sbatch
fi

function usage() {
cat <<EOF
    Usage: run_atac_10x.sh [options] sampleSheet

    Options should be specified before the sampleSheeet. Options:

     -s S | Perform only the steps specified in the string S. Default: $STEPS
            See below for a description of steps.
     -o O | Specify organism name. Default: $TENX_ORG. Possible values:
            GRCh38, hg19, hg19_and_mm10, mm10, ercc92, Rnor96.
     -d D | Illumina run directory. Should be the path to the directory that
            contains the RunInfo.xml file. Default: '${TENX_RUN_DIR}'.
     -b B | Base mask for bcl2fastq (default: Y28N*,I10,I10,Y91N*).
     -f F | Store fastq files in folder F. Default: $TENX_RUN_ID.
     -a A | Store aggregated dataset in folder A. Default: $AGGRDIR.
     -c C | Use chemistry C. Default: auto-detect.
     -x X | Perform differential analysis on the pairs of samples listed in
            file X (a tab-delimited file with a pair of sample names in each row).
     -r R | Write final report to folder R. Default: $REPORT.
          | The report folder will be zipped to a zip file with the same name.
     -h   | Print this help message.

    Steps:

     1 | Generate fastq files (cellranger-atac 'mkfastq' tool). Creates the directory
         specified with the -f option and saves fastq files in it.
     2 | Generate gene counts (cellranger-atac 'count' tool). Creates one directory
         for each sample.
     3 | Aggregate samples    (cellranger-atac 'aggr' tool). Writes aggregated data
         to a folder named 'Aggregated'.
     4 | Reanalyze samples    (cellranger-atac 'reanalize' tool).
     5 | Generate final report. Creates the directory specified with the -r option,
         and saves HTML and .cloupe files to it, plus an index.html file linking
         all output files in a table. Finally, zips this folder.

    Sample sheet format:

    The sample sheet should be a comma-delimited text file, with each line having the
    format 'lane,samplename,barcode'. Example:

     Lane,Sample,Index
     *,PBMC_1,SI-GA-A1
     *,PBMC_2,SI-GA-B1
     *,PBMC_3,SI-GA-C1

    Return values:

     0 | Pipeline terminated successfully.
     1 | Help message printed.
     2 | Error - missing sample sheet, or RunInfo.xml not found.

EOF
}

while getopts "s:r:o:f:d:b:c:x:a:h" opt; do
    case $opt in
	s)
	    STEPS="$OPTARG"
	    ;;
	r)
            REPORT="$OPTARG"
            ;;
	o)
	    TENX_ORG="$OPTARG"
	    ;;
	f)
	    TENX_RUN_ID="$OPTARG"
	    ;;
	d)
	    TENX_RUN_DIR="$OPTARG"
	    ;;
	b)
	    BM="$OPTARG"
	    ;;
	c)
	    CHEM="$OPTARG"
	    ;;
	x)
	    DIFF="$OPTARG"
	    ;;
	a)
	    AGGRDIR="$OPTARG"
	    ;;
	h)
	    usage
	    exit 1
	    ;;
    esac
done
shift $((OPTIND-1))

# Name of sample sheet file. Example:
# *,PBMC_1,SI-GA-A1
# *,PBMC_2,SI-GA-B1
# *,PBMC_3,SI-GA-C1

# Sample sheet for mkfastq
TENX_SAMPLE_SHEET=$1

if [[ -z "$TENX_SAMPLE_SHEET" ]];
then
  usage
  exit 1
fi

# Check that we have what we need...

if [[ ! -f "$TENX_SAMPLE_SHEET" ]];
then
    echo "Error: sample sheet $TENX_SAMPLE_SHEET not found!"
    exit 2
fi

# Let's get started...
#source /apps/dibig_tools/1.0/lib/sh/utils.sh

SAMPLES=$(grep -i -v ^lane ${TENX_SAMPLE_SHEET} | cut -d , -f 2)
NSAMPLES=$(grep -i -c -v ^lane ${TENX_SAMPLE_SHEET})

echo $NSAMPLES samples: $SAMPLES
if [[ -f options.txt ]];
then
  echo Using options from: options.txt
fi



function step1() {

  # Check that we have a valid run directory
  if [ ! -f $TENX_RUN_DIR/RunInfo.xml ];
  then
    echo "Error: file RunInfo.xml not found - please specify path to a valid Illumina run directory with the -d option."
    exit 2
  fi

  # Call mkfastq
  echo "Starting step mkfastq."
  rm -fr $TENX_RUN_ID
  $SUBMIT -W $TENX_HOME/scripts/mkfastq.qsub $TENX_RUN_ID $TENX_SAMPLE_SHEET $TENX_RUN_DIR $BM
  echo "Step mkfastq terminated."
}

function step2() {
  # Call count

  OPTS=""
  if [[ -f options.txt ]];
  then
    OPTS=$(grep ^count options.txt | cut -f 2)
    SLURMOPTS=$(grep ^slurm options.txt | cut -f 2)
  fi
  echo "Starting step count."
  for smp in $SAMPLES;
  do
    rm -rf $smp
    $SUBMIT -W $SLURMOPTS $TENX_HOME/scripts/atac_count.qsub $smp $TENX_ORG $TENX_RUN_ID $CHEM $OPTS &
  done
  wait
  echo "Step count terminated."
}

function make_aggr_table() {
    IN=$1
    OUT=$2

    # IN = Sample sheet used in mkfastq step
    # OUT = Sample sheet needed by aggr step

    rm -f $OUT
    HERE=$(pwd)
    echo "library_id,fragments,cells" > $OUT
    for s in $( grep -i -v ^lane $IN | cut -f 2 -d , );
    do
	FRAG="${HERE}/$s/outs/fragments.tsv.gz"
	CELL="${HERE}/$s/outs/singlecell.csv"
	if [[ -f $FRAG ]];
	then
	    echo "${s},${FRAG},${CELL}" >> $OUT
	else
	    echo "Warning: fragments file $FRAG does not exist!"
	    exit 1
	fi
    done
}

function step3() {
    # Call aggr

    echo "Starting step aggr."
    rm -rf $AGGRDIR
    TENX_AGGR_SHEET=${TENX_SAMPLE_SHEET%.*}.aggr.csv
    make_aggr_table $TENX_SAMPLE_SHEET $TENX_AGGR_SHEET
    $SUBMIT -W $TENX_HOME/scripts/do_atac_aggr.qsub $AGGRDIR $TENX_AGGR_SHEET $TENX_ORG
    echo "Step aggr terminated."
}

function write_index() {
    REP=${REPORT}/index.html
    cat > ${REP} <<EOF
<!DOCTYPE html>
<HTML>
<HEAD>
  <TITLE>${REPORT} - 10X Analysis Report</TITLE>
  <STYLE>
TABLE {
  border-collapse: collapse;
  border: 2px solid black;
}
TD {
  padding: 8px;
  border: 1px solid black;
}
TH {
  padding: 8px;
}
</STYLE>
</HEAD>
<BODY>
<CENTER>
<H1>${REPORT} - 10X Analysis Report</H1>
<BR>
<H2><I>Quality Control</I></H2>
<TABLE>
<TR><TH>Sample</TH><TH>Report</TH><TH>Cloupe</TH><TH>Number of cells</TH><TH>Mean fragments / cell</TH><TH>Frac fragments over targets</TH><TH>Frac fragments in peaks</TH><TH>Number of fragments</TH></TR>
EOF

    ALLMET=""
    for smp in $SAMPLES;
    do
      MET=${smp}/outs/summary.csv
      ALLMET="$ALLMET $MET"
      metrics=$(${TENX_HOME}/scripts/getMetrics.py -a $MET)
      echo "<TR><TD>$smp</TD><TD><A href='${smp}_summary.html' target='${smp}_report'>${smp}_summary.html</A></TD><TD><A href='${smp}.cloupe'>${smp}.cloupe</A></td>" >> ${REP}
      echo $metrics >> ${REP}
      echo "</TR>" >> ${REP}
    done
    ${TENX_HOME}/scripts/getMetrics.py -a -s $ALLMET > ${REPORT}/summary.csv
    echo "<TR><TD><B>Aggregated</B></TD><TD><A href='Aggregated_summary.html' target='${smp}_report'>Aggregated_summary.html</A></TD><TD><A href='Aggregated.cloupe'>Aggregated.cloupe</A></TD><TD colspan=5 align='center'><A href='metrics_summary.csv'><B>Metrics Summary</B></A></TD></TR>" >> ${REP}
    cat >> ${REP} <<EOF
</TABLE>
<BR><BR>
</CENTER>
</BODY>
</HTML>
EOF
}

function step5() {
    echo "Create output directory."
    rm -fr ${REPORT}
    mkdir ${REPORT}
    cp -v Aggregated/outs/cloupe.cloupe ${REPORT}/Aggregated.cloupe
    cp -v Aggregated/outs/web_summary.html ${REPORT}/Aggregated_summary.html
    for smp in $SAMPLES; 
    do
	cp -v $smp/outs/cloupe.cloupe ${REPORT}/${smp}.cloupe
	cp -v $smp/outs/web_summary.html ${REPORT}/${smp}_summary.html
	# cp -v ${smp}.seurat.html ${smp}.rds ${REPORT}
	# ${TENX_HOME}/scripts/markers.py ${smp}.markers.csv ${REPORT}/${smp}.markers.xlsx
    done
    write_index
    zip -r ${REPORT}.zip ${REPORT}/
    echo "Report created."
}

# Main

if [[ $STEPS == *1* ]]; then step1; fi
if [[ $STEPS == *2* ]]; then step2; fi
if [[ $STEPS == *3* ]]; then step3; fi
if [[ $STEPS == *5* ]]; then step5; fi

exit 0
