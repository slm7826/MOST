#!/bin/bash
set -e # exit on errors
set -u # exit on undefined vars


# store full command line for future reference
cmdline="$0 $@"

# ---- default values
verbosity=0
     tool="$HOME/work/notes/lm4p2/sfcx-experiments/crashes/plot-watch-output.py"
   outdir=.

print_help () {
cat <<EOF
Usage:
`basename $0` [-v] -o output-dir crashlog
   -v = increase verbosity
   -o = use specified output directory
EOF
}

# print error message and exit with specified error code
die () {
   echo "${0##*/} :: ERROR :: $1" >&2
   exit ${2:-1}
}

# parse command-line arguments
while getopts ":hvo:" Option
do
  case "$Option" in
    h) print_help ; exit 1 ;;
    v) verbosity=$(($verbosity+1)) ;;
    o) outdir="$OPTARG" ;;
  [?]) print_help >& 2 ; die "Illegal command-line option \"-$OPTARG\"" ;;
  esac
done
shift $(( $OPTIND - 1 ))

if ((verbosity > 0)); then
  tool="$tool --print"
fi

crashlog="$1"
[[ -e $crashlog ]] || die "file \"$crashlog\" does not exist"
[[ -f $crashlog ]] || die "file \"$crashlog\" is not a regular file"
tool="$tool -i $crashlog"

[[ -e $outdir ]] || die "output directory \"$outdir\" does not exist"
[[ -d $outdir ]] || die "\"$outdir\" is not a directory"

$tool --e cosz --save="$outdir/cosz.pdf"
$tool --ylim=250:370 --e atmos_T --e cana_T+delta_Tc --e grnd_T+delta_Tg --e vegn_T+delta_Tv --save="$outdir/T.pdf"
$tool --e Ha0  --e land_sens --save="$outdir/sens.pdf"
$tool --e fog+delta_fog --save="$outdir/fog.pdf"
$tool --e atmos_wind --e ustar --e u_sfc --e ustar_sfc --save="$outdir/winds.pdf"
$tool --e vegn_fsw --e fswg --e flwg --save="$outdir/radiation.pdf"
$tool --e con_g_h --save="$outdir/con_g_h.pdf"
$tool --e drag_q --save="$outdir/drag_q.pdf"
$tool --e drag_q  --e atmos_wind --save="$outdir/drag-vs-wind.pdf"
