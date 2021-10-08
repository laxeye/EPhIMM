#!/bin/bash

# EPhIMM - Express Phylogenetic Inference based on Multiple Markers
# Written by Aleksei Korzhenkov, 2019
# https://github.com/laxeye/EPhIMM

#Get a path to script
EPHIMMPATH=$(dirname $(readlink -f $0))

#Defaults
MRKDIR="$EPHIMMPATH/markers/"
OUTPUTFOLDER="output-$(date +%F)"  #Default name includes current date
EVALUE=1e-6   # e-value for HMMER or BLAST
AUTO=0        # 1 to perform autoedit of alignment
MAXGAPFRACTION=50  # Set maximum gap fraction in column to leave it in the alignment.
PREDICT=0     # 1 to predict protein coding genes. Input will be genome assemblies *.fna.
USEBLAST=0    # 1 to use BLASTp homology search instead of HMMER (default).
MARKEREXT="hmm"  # Extension of markers. Default is "hmm" for HMMER.
MAXMULTICOPYSHARE=0.5  # Max number of genomes, having more than one copy of marker gene.
PROTEOMEEXT="faa"
GENOMEEXT="fna"
QUIET=0
THREADS=1
FORCEBH=0
PREDICTEDEXT="faa"
USENUCL=0

#Parse args
while (( "$#" )); do
  case "$1" in
    -m|--markers)
      if [[ -d $2 ]]; then
        MRKDIR=$2
      fi
      shift 2
      ;;
    -s|--sequences)
      if [[ -d $2 ]]; then
        MRKDIR=$2
      fi
      MARKEREXT="faa"
      USEBLAST=1
      shift 2
      ;;
    -g|--max-gap-fraction)
      if [[ $2 -gt 0 && $2 -lt 100 ]]; then
        MAXGAPFRACTION=$2
        AUTO=1
      else
        echo "Error: Max gap fraction should be greater than 0 and less than 100"
        exit 1
      fi
      shift 2
      ;;
    -e|--evalue)
      EVALUE=$2
      shift 2
      ;;
    -o|--output)
      OUTPUTFOLDER=$2
      shift 2
      ;;
    -p|--predict-proteins)
      PREDICT=1
      shift
      ;;
    -b|--besthit)
      FORCEBH=1
      shift
      ;;
    -x|--extension)
      if [[ PREDICT -eq 1 ]]; then
        GENOMEEXT=$2
      else
        PROTEOMEEXT=$2
      fi
      shift 2
      ;;
    -d|--nucleotide)
      PREDICTEDEXT="fna"
      USENUCL=1
      shift
      ;;
    -a|--auto)
      AUTO=1
      shift
      ;;
    -t|--threads)
      THREADS=$2
      shift 2
      ;;
    -q|--quiet)
      QUIET=1
      shift
      ;;
    -h|--help)
      echo -e "\nEPhIMM v 0.3 - Express Phylogenetic Inference based on Multiple Markers"
      echo -e "Written by Aleksei Korzhenkov, 2019-2021"
      echo -e "For new versions please check https://github.com/laxeye/EPhIMM\n"
      echo -ne "Usage: $0 [-p|--predict-proteins] [-e|--evalue <N>] "
      echo -ne "[-g|--max-gap-fraction <0-100>] [-o|--output <folder>] "
      echo -ne "{[-m|--markers <folder>] | [-s|--sequences <folder>] " 
      echo -e "[-d|--nucleotide] [-a|--auto] [-t|--threads <N>]"
      echo "or"
      echo -e "$0 [-h|--help]\n"
      exit 0
      ;;
    --) # end argument parsing
      shift
      break
      ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      echo "Possible flags: -m, -e, -g, -o, -p, -a, -t, -x " >&2
      exit 1
      ;;
    *) # preserve positional arguments
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac
done

function log {
  if [[ $QUIET = 1 ]]; then
    echo -e "#$(date +"%T")\t$1" >> $logfile
  else
    echo -e "#$(date +"%T")\t$1" | tee -a $logfile
  fi
}

maxthreads=`grep -c "^processor" /proc/cpuinfo`
if [[ $THREADS -gt $maxthreads ]]; then
    THREADS=$maxthreads
fi

if [[ -d $OUTPUTFOLDER ]]; then
    rm -fr $OUTPUTFOLDER
fi
mkdir -p $OUTPUTFOLDER
logfile="$OUTPUTFOLDER/fetching.log"

log "#Analysys started"
# Predict CDS
if [[ $PREDICT -eq 1 ]]; then
    genomecount=`ls *$GENOMEEXT | wc -l`
    if [[ $genomecount -lt 3 ]]; then
      log "Error: too low number of genomes - $genomecount. Please check input files and extension."
      exit 1
    fi
    log "Protein prediction and file indexing started"
    for file in *fna;
      do aafile=cds.${file/$GENOMEEXT/faa};
      nucfile=cds.${file/$GENOMEEXT/fna}
      if ! [[ -f $aafile ]]; then
        prodigal -o /dev/null -i $file -a $aafile -d $nucfile -q;
        esl-sfetch --index $aafile > /dev/null
        esl-sfetch --index $nucfile > /dev/null
      fi;
    done
else
  for aafile in *faa; do
    esl-sfetch --index $aafile > /dev/null;
  done
fi

log "Processing files"
genomecount=`ls *faa | wc -l`

MARKERLIST=$OUTPUTFOLDER/markers.list
ls $MRKDIR/ | grep $MARKEREXT$ > $MARKERLIST
markercount=`wc -l < $MARKERLIST`

if [[ $markercount -eq 0 ]]; then
    log "Error! No markers found in $MRKDIR"
    exit 1
else
    log "#$markercount markers found in $MRKDIR"
fi

log "#Searching and extracting marker genes"

for marker in $(cat $MARKERLIST); do
  mrkpath=$MRKDIR$marker
  marker=${marker/.$MARKEREXT/}

  if [[ -d $marker ]]; then
      rm -r -f $marker
  fi
  mkdir $marker;

  for file in *faa; do
    result=$marker/${file%.*}.out

    if [[ $USEBLAST -eq 0 ]]; then
      hmmsearch --cpu $THREADS -E $EVALUE -o /dev/null --tblout /dev/stdout $mrkpath $file | grep -v "^#" > $result;
    else
      blastp -evalue $EVALUE -outfmt "6 sseqid qseqid slen qlen length pident ppos evalue bitscore" -max_target_seqs 5 -query $mrkpath -subject $file -out $result;
    fi

    lines=`wc -l < $result`

    case $lines in
        0)
            log "Marker gene \"$marker\" not found in $file."
            echo -e ">$marker.${file/.faa/}\nXXXXXX" > $marker/${file/.faa/.fasta};
            echo "$marker.${file/.faa/}" >> $marker/missing.txt;
            ;;
        1)
            esl-sfetch -o $marker/${file/.faa/.fasta} -n $marker.${file/.faa/} ${file/.faa/.$PREDICTEDEXT} $(awk '{print $1}' $result) > /dev/null;
            ;;
        *)
            log "Marker gene \"$marker\" has many copies in $file."
            if [[ $FORCEBH -eq 1 ]]; then
              esl-sfetch -o $marker/$file -n $marker.${file/.faa/} $file $(awk '{print $1}' $result | head -1) > /dev/null;
            else
              echo "$marker.${file/.faa/}" >> $marker/multicopy.txt;
              let n=1
              if [[ $USENUCL -eq 1 ]]; then
                for line in $(awk '{print $1}' $result)
                  do esl-sfetch -n $marker.${file/.faa/}.$n ${file/.faa/.$GENOMEEXT} $line >> $marker/${file/.faa/.fasta};
                  ((n++));
                done;
              else
                for line in $(awk '{print $1}' $result)
                  do esl-sfetch -n $marker.${file/.faa/}.$n $file $line >> $marker/${file/.faa/.fasta};
                  ((n++));
                done;
              fi
            fi
            ;;
    esac

  done;

done;

log "Working with multicopy markers"
for marker in $(cat $MARKERLIST); do
  marker=${marker/.$MARKEREXT/}
  if [[ -f $marker/multicopy.txt ]] ; then 
    cat $marker/*fasta > $marker/$marker.1.fasta ;
    if [[ -f $marker/missing.txt ]] ; then
      java -jar $EPHIMMPATH/BioKotlin.jar RemoveByName $marker/$marker.1.fasta $marker/missing.txt > $marker/$marker.fasta ;
      rm -f $marker/$marker.1.fasta
    else
      mv $marker/$marker.1.fasta $marker/$marker.fasta ;
    fi
    mafft --quiet --auto $marker/$marker.fasta > $marker/aln.fasta ;
    java -jar $EPHIMMPATH/BioKotlin.jar DistMeanProtJC $marker/aln.fasta | grep -f $marker/multicopy.txt > $marker/multicopy.dist ;
    rm -f $marker/aln.fasta
    for mcopy in $(cat $marker/multicopy.txt); do
      grep $mcopy $marker/multicopy.dist | sort -k 2 | cut -f 1 | head -1 > $marker/keep.txt
      genome=${mcopy/$marker./}
      java -jar $EPHIMMPATH/BioKotlin.jar ExtractByName $marker/$genome.fasta $marker/keep.txt > $marker/$genome.1.fasta
      mv $marker/$genome.1.fasta $marker/$genome.fasta
    done
  fi
done

log "Collecting markers"

echo -e "#Genome\tMissing markers\tMultiple copies" > $OUTPUTFOLDER/genomes.stats.txt
for file in *faa;
    do missing=`grep $file $logfile | grep -c "not found"`;
    multi=`grep $file $logfile | grep -c "many copies"`;
    echo -e "$file\t$missing\t$multi" >> $OUTPUTFOLDER/genomes.stats.txt;
done;

echo "#Genomes with more than 50% missing marker genes:" > $OUTPUTFOLDER/bad.genomes.txt
grep -v "^#" $OUTPUTFOLDER/genomes.stats.txt | awk '$2 > "'"$markercount"'" * 0.5 ' >> $OUTPUTFOLDER/bad.genomes.txt
echo "#Genomes with share of multi-copy marker genes more than $MAXMULTICOPYSHARE:" >> $OUTPUTFOLDER/bad.genomes.txt
grep -v "^#" $OUTPUTFOLDER/genomes.stats.txt | awk '$3 > "'"$markercount"'" * "'"$MAXMULTICOPYSHARE"'" ' >> $OUTPUTFOLDER/bad.genomes.txt

echo -e "#Marker\tMissing markers\tMultiple copies" > $OUTPUTFOLDER/marker.stats.txt
for marker in $(cat $MARKERLIST); do
  marker=${marker/.$MARKEREXT/}
  missing=`grep "$marker" $logfile | grep -c "not found"`
  multi=`grep "$marker" $logfile | grep -c "many"`
  echo -e "$marker\t$missing\t$multi" >> $OUTPUTFOLDER/marker.stats.txt;
done

echo "#Markers missing in more than 50% genomes:" > $OUTPUTFOLDER/bad.markers.txt
grep -v "^#" $OUTPUTFOLDER/marker.stats.txt | awk '$2 > "'"$genomecount"'" * 0.5 ' >> $OUTPUTFOLDER/bad.markers.txt
echo "#Markers with multiple copies in more than $MAXMULTICOPYSHARE genomes:" >> $OUTPUTFOLDER/bad.markers.txt
grep -v "^#" $OUTPUTFOLDER/marker.stats.txt | awk '$3 > "'"$genomecount"'" * "'"$MAXMULTICOPYSHARE"'" ' >> $OUTPUTFOLDER/bad.markers.txt

for file in *faa; do
    if ! grep -q "$file" $OUTPUTFOLDER/bad.genomes.txt ; then
        echo ">${file/.faa/}" > $OUTPUTFOLDER/$file;
        for marker in $(cat $MARKERLIST); do
          marker=${marker%.*}
          cat $marker/${file/.faa/.fasta} | sed 's/\*//' | grep -v ">" | grep -v "^X\+$"  >> $OUTPUTFOLDER/$file;
        done
    else
        echo "Genome $file will not be used because of incompleteness or contamination."
    fi
done

cat $OUTPUTFOLDER/*faa > $OUTPUTFOLDER/all.markers.fasta

log "Concatenated markers stored in \"$OUTPUTFOLDER/all.markers.fasta\""


log "Performing MSA"

mafft --thread -1 --auto  $OUTPUTFOLDER/all.markers.fasta | sed 's/X/-/g' > $OUTPUTFOLDER/all.aligned.fasta

alnfilename="$OUTPUTFOLDER/all.aligned.e.fasta"

if [[ $AUTO -eq 1 ]]; then

  java -jar $EPHIMMPATH/BioKotlin.jar AlignmentClearGaps $OUTPUTFOLDER/all.aligned.fasta $MAXGAPFRACTION > $alnfilename
  log "Performing autoedit of MSA"

else

  echo "It's strictly recommended to check alignment file \"$OUTPUTFOLDER/all.aligned.fasta\""
  echo "You may automaticly delete all columns with more than $MAXGAPFRACTION% gaps"

  read -p "Press \"N\" if You like to edit file manually.
  Press \"Y\" if You like to clean gaps automaticly: " response

  if [[ $response -eq "Y" ]]; then
      java -jar $EPHIMMPATH/BioKotlin.jar AlignmentClearGaps $OUTPUTFOLDER/all.aligned.faa $MAXGAPFRACTION > $alnfilename
      log "Performing autoedit of MSA"
  else
      read -p "Provide a filename of edited alignment:" alnfilename
  fi

fi

log "Performing phylogenetic inference"

if [[ $USENUCL -eq 1 ]]; then
  fasttree -log $OUTPUTFOLDER/fastree.log -gamma -gtr $alnfilename > $OUTPUTFOLDER/all.markers.nwk
else
  fasttree -log $OUTPUTFOLDER/fastree.log -gamma -lg $alnfilename > $OUTPUTFOLDER/all.markers.nwk
fi

log "Workflow finished"
