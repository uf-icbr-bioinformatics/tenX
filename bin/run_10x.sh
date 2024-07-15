#!/bin/bash

# Change this if scripts are moved somewhere else!
TENX_HOME=/apps/dibig_tools/tenX

# Specify default versions of cellranger and R
# Note: these should start with "/".
#export CELLRANGER_VERSION="/7.0.1"
export CELLRANGER_VERSION="/8.0.0"
export R_VERSION="/4.3"

STEPS=123456
TENX_RUN_DIR=.
TENX_RUN_ID=Fastq
REPORT=Report
AGGRDIR=Aggregated
OPTFILE="options.txt"
# Organism ID
TENX_ORG=GRCh38
# Base mask
BM="Y28N*,I10,I10,Y91N*"
# Chemistry
CHEM=""
# Differential analysis
DIFF=""
# ParseBiosciences support
PARSE="FALSE"

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
  * Usage: run_10x.sh [options] sampleSheet

    Options should be specified before the sampleSheeet. Options:

     -s S | Perform only the steps specified in the string S. Default: $STEPS
            See below for a description of steps.
     -v V | Read options from file V (default: options.txt).
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
     -p   | Dataset is from ParseBiosciences pipeline (only run steps 4-6).
     -h   | Print this help message.

  * Steps:

     1 | Generate fastq files (cellranger 'mkfastq' tool). Creates the directory
         specified with the -f option and saves fastq files in it.
     2 | Generate gene counts (cellranger 'count' tool). Creates one directory
         for each sample.
     3 | Aggregate samples    (cellranger 'aggr' tool). Writes aggregated data
         to a folder named 'Aggregated'.
     4 | Run Seurat on all samples (including the aggregate).
     5 | Perform differential analysis within each cluster.
     6 | Generate final report. Creates the directory specified with the -r option,
         and saves HTML and .cloupe files to it, plus an index.html file linking
         all output files in a table. Finally, zips this folder.

  * Sample sheet format:

    The sample sheet should be a comma-delimited text file, with each line having the
    format 'lane,samplename,barcode'. Example:

     Lane,Sample,Index
     *,PBMC_1,SI-GA-A1
     *,PBMC_2,SI-GA-B1
     *,PBMC_3,SI-GA-C1

    Optionally, a fourth column called Condition can be added, to group replicates of
    the same condition.

    Return values:

     0 | Pipeline terminated successfully.
     1 | Help message printed.
     2 | Error - missing sample sheet, or RunInfo.xml not found.


  * Step-specific options:

    Additional options can be specified in the options.txt file (or the file specified
    with the -v option). This should be a tab-delimited file in which the first column
    contains a key and the second column contains the options to be passed to the appropriate
    process. These keys are currently defined:

    slurm  - used for SLURM job submission
    count  - passed to cellranger count (step 2)
    seurat - passed to seurat (steps 4 and 5)
    diff   - passed to seurat for differential analysis (step 5)

EOF
}

while getopts "s:r:o:f:d:b:c:x:a:v:ph" opt; do
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
	v)
	    OPTFILE="$OPTARG"
	    ;;
	p)
	    PARSE="TRUE"
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
if [[ -f $OPTFILE ]];
then
  echo Using options from: $OPTFILE
fi

if [[ "$PARSE" == "TRUE" ]];
then
  PLATFORM="ParseBio"
else
  PLATFORM="10X"
fi
echo "Running in ${PLATFORM} mode."

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
  if [[ -f $OPTFILE ]];
  then
    OPTS=$(grep ^count $OPTFILE | cut -f 2)
    SLURMOPTS=$(grep ^slurm $OPTFILE | cut -f 2)
  fi
  echo "Starting step count."
  for smp in $SAMPLES;
  do
    rm -rf $smp
    $SUBMIT -W $SLURMOPTS $TENX_HOME/scripts/do_count.qsub $smp $TENX_ORG $TENX_RUN_ID $CHEM $OPTS &
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
    echo "sample_id,molecule_h5" > $OUT
    for s in $( grep -i -v ^lane $IN | cut -f 2 -d , );
    do
	H5="$s/outs/molecule_info.h5"
	if [[ -f $H5 ]];
	then
	    echo "$s,$H5" >> $OUT
	else
	    echo "Warning: h5 file $H5 does not exist!"
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
    $SUBMIT -W $TENX_HOME/scripts/do_aggr.qsub $AGGRDIR $TENX_AGGR_SHEET
    echo "Step aggr terminated."
}

function step4() {
    echo "Starting step Seurat."

    mingenes=500
    maxgenes=5000
    maxmt=5
    mtpatt="^mt-"
    nfeatures=2000
    ntop=12
    dimensions=10
    resolution=0.5
    minpct=0.25
    logfc=0.25
    celldex="FALSE"
    parse=${PARSE}

    if [[ -f $OPTFILE ]];
    then
      vars=$(grep ^seurat $OPTFILE | cut -f 2)
      eval $vars
    fi

    mkdir -p R markers rds
    for smp in $SAMPLES;
    do
	rmd=${smp}.seurat.Rmd
	cat > ${rmd} <<EOF
---
output: html_document
params:
  funcs: "${TENX_HOME}/R/seuratfuncs.R"
  data_dir: "${smp}/outs/filtered_feature_bc_matrix/"
  sample: "${smp}"
  mingenes: ${mingenes}
  maxgenes: ${maxgenes}
  maxmt: ${maxmt}
  mtpatt: "${mtpatt}"
  method: "LogNormalize"
  nfeatures: ${nfeatures}
  ntop: ${ntop}
  dimensions: ${dimensions}
  resolution: ${resolution}
  minpct: ${minpct}
  logfc: ${logfc}
  parse: ${parse}
  celldex: "${celldex}"
---

EOF
	cat ${TENX_HOME}/R/seurat.Rmd.template >> $rmd
	$SUBMIT -W ${TENX_HOME}/scripts/Rbatch.qsub $TENX_HOME/R/Render.R $rmd &
    done
    wait
    echo "Step seurat terminated."
}

function step5() {
    if [[ ! -f "$DIFF" ]];
    then
	echo "Error: file ${DIFF} does not exist, cannot perform differential analysis."
	return
    fi

    echo "Starting differential analysis step."
    mingenes=500
    maxgenes=5000
    maxmt=5
    logfc=0.25
    pvalue=0.05
    resolution=0.5
    mtpatt="^mt-"
    mincells=50
    ntop=12
    parse=${PARSE}
    celldex="FALSE"
    
    if [[ -f $OPTFILE ]];
    then
      vars=$(grep ^seurat $OPTFILE | cut -f 2)
      eval $vars
      vars=$(grep ^diff $OPTFILE | cut -f 2)
      eval $vars
    fi

    mkdir -p _diff
    while read cond1 cond2; do
	smp1=$($TENX_HOME/scripts/cond_samples.py $TENX_SAMPLE_SHEET $cond1)
	smp2=$($TENX_HOME/scripts/cond_samples.py $TENX_SAMPLE_SHEET $cond2)
	rmd=${cond1}.vs.${cond2}.integ.Rmd
	cat > ${rmd} <<EOF
---
output: html_document
params:
  funcs: ${TENX_HOME}/R/seuratfuncs.R
  sample1: ${smp1}
  sample2: ${smp2}
  condition1: ${cond1}
  condition2: ${cond2}
  mingenes: ${mingenes}
  maxgenes: ${maxgenes}
  maxmt: ${maxmt}
  mtpatt: ${mtpatt}
  resolution: ${resolution}
  logfc: ${logfc}
  pvalue: ${pvalue}
  mincells: ${mincells}
  ntop: ${ntop}
  parse: ${parse}
  celldex: "${celldex}"
---
EOF
	cat ${TENX_HOME}/R/integ.Rmd.template >> $rmd
	echo $cond1 vs $cond2
	$SUBMIT -W -o --mem=50G ${TENX_HOME}/scripts/Rbatch.qsub ${TENX_HOME}/R/Render.R $rmd &
    done < ${DIFF}
    wait
}

function check_seurat() {
    for smp in $SAMPLES; do
	if [[ ! -f ${smp}.seurat.html ]]; then return 1; fi
    done
    return 0
}

function write_index() {
    REP=${REPORT}/index.html
    cat > ${REP} <<EOF
<!DOCTYPE html>
<HTML>
<HEAD>
  <TITLE>${REPORT} - ${PLATFORM} Analysis Report</TITLE>
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
<H1>${REPORT} - ${PLATFORM} Analysis Report</H1>
<BR>
<H2><I>Quality Control</I></H2>
<TABLE>
<TR><TH>Sample</TH><TH>Report</TH><TH>Cloupe</TH><TH>Number of cells</TH><TH>Mean reads / cell</TH><TH>Median genes / cell</TH><TH>Fraction reads in cells</TH><TH>Number of reads</TH></TR>
EOF

    ALLMET=""
    for smp in $SAMPLES;
    do
      MET=${smp}/outs/metrics_summary.csv
      ALLMET="$ALLMET $MET"
      metrics=$(${TENX_HOME}/scripts/getMetrics.py $MET)
      echo "<TR><TD>$smp</TD><TD><A href='${smp}_summary.html' target='${smp}_report'>${smp}_summary.html</A></TD><TD><A href='${smp}.cloupe'>${smp}.cloupe</A></td>" >> ${REP}
      echo $metrics >> ${REP}
      echo "</TR>" >> ${REP}
    done
    ${TENX_HOME}/scripts/getMetrics.py -s $ALLMET > ${REPORT}/metrics_summary.csv
    if [[ "$PLATFORM" == "10X" ]];
    then
      echo "<TR><TD><B>Aggregated</B></TD><TD><A href='Aggregated_summary.html' target='${smp}_report'>Aggregated_summary.html</A></TD><TD><A href='Aggregated.cloupe'>Aggregated.cloupe</A></TD><TD colspan=5 align='center'><A href='metrics_summary.csv'><B>Metrics Summary</B></A></TD></TR>" >> ${REP}
    fi
    cat >> ${REP} <<EOF
</TABLE>
<BR><BR>
EOF

    check_seurat
    if [[ "$?" == "0" ]]; then
	cat >> ${REP} <<EOF
<H2><I>Analysis</I></H2>
<TABLE>
<TR><TH>Sample</TH><TH>Seurat Report</TH><TH>Seurat RDS</TH><TH>Marker genes</TH></TR>
EOF

	for smp in $SAMPLES;
	do
	    echo "<TR><TD>$smp</TD><TD><A href='${smp}.seurat.html' target='_blank'>${smp}.seurat</A></TD><TD><A href='${smp}.rds' target='_blank'>${smp}.rds</A></TD><TD><A href='${smp}.markers.xlsx' target='_blank'>${smp}.markers.xlsx</A></TD></TR>" >> ${REP}
	done

	echo "</TABLE>" >> ${REP}
    fi
    
    if [[ $DIFF ]];
    then
	cat >>${REP} <<EOF
<BR><BR>
<H2><I>Differential Analysis</I></H2>
<TABLE>
<TR><TH>Test</TH><TH>Ctrl</TH><TH>Report</TH><TH>Clusters</TH><TH>Markers</TH><TH>DE genes</TH>
EOF
	while read test ctrl; do
	    report="${test}.vs.${ctrl}.integ.html"
	    markers="${test}_${ctrl}.markers.csv"
	    markerx="${test}_${ctrl}.markers.xlsx"
	    zip="${test}.vs.${ctrl}.zip"
	    nclust=$(unzip -Z -1 $zip | wc -l)
	    cp -v $report $zip ${REPORT}
	    ${TENX_HOME}/scripts/markers.py _diff/$markers ${REPORT}/$markerx
	    echo "<TR><TD align='center'>$test</TD><TD align='center'>$ctrl</TD><TD><A href='${report}' target='_blank'>${report}</A></TD><TD align='right'>${nclust}</TD><TD><A href='${markerx}' target='_blank'>${markerx}</A></TD><TD><A href='${zip}' target='_blank'>${zip}</A></TD></TR>" >> ${REP}
	done < $DIFF
	echo "</TABLE>" >> ${REP}
    fi
cat >>${REP} <<EOF
<BR><BR>
</CENTER>
</BODY>
</HTML>
EOF
}

function step6() {
    echo "Create output directory."
    rm -fr ${REPORT}
    mkdir ${REPORT}
    cp -v Aggregated/outs/count/cloupe.cloupe ${REPORT}/Aggregated.cloupe
    cp -v Aggregated/outs/web_summary.html ${REPORT}/Aggregated_summary.html
    for smp in $SAMPLES; 
    do
	cp -v $smp/outs/cloupe.cloupe ${REPORT}/${smp}.cloupe
	cp -v $smp/outs/web_summary.html ${REPORT}/${smp}_summary.html
	if [[ -f ${smp}.seurat.html ]]; then
	    cp -v ${smp}.seurat.html ${smp}.rds ${REPORT}
	    ${TENX_HOME}/scripts/markers.py ${smp}.markers.csv ${REPORT}/${smp}.markers.xlsx
	fi
    done
    write_index
    zip -r ${REPORT}.zip ${REPORT}/
    echo "Report created."
}

# Clean environment
#loaded_modules="$(module -t --redirect list)"
#module purge
#module load dibig_tools

module load python/3

rm -f *.done

if [[ $STEPS == *1* ]]; then step1; fi
if [[ $STEPS == *2* ]]; then step2; fi
if [[ $STEPS == *3* ]]; then step3; fi
if [[ $STEPS == *4* ]]; then step4; fi
if [[ $STEPS == *5* ]]; then step5; fi
if [[ $STEPS == *6* ]]; then step6; fi

#module load "$loaded_modules"

exit 0
