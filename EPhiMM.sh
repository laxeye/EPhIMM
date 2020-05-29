#!/bin/bash

# EPhIMM - Express Phylogenetic Inference based on Multiple Markers
# Written by Aleksei Korzhenkov, 2019
# https://github.com/laxeye/EPhIMM

#Get a path to script
EPHIMMPATH=$(dirname $(readlink -f $0))
#Defaults
MRKDIR="$EPHIMMPATH/markers/"
HMMEVALUE=1e-6
AUTO=0
#Set maximum gap fraction in column to leave it in alignment
MAXGAPFRACTION=50
OUTPUTFOLDER="output-$(date +%F)"
PREDICT=0

#Parse args
while (( "$#" )); do
  case "$1" in
    -m|--markers)
      if [[ -d $2 ]]; then
        MRKDIR=$2
      fi
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
      HMMEVALUE=$2
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
    -a|--auto)
      AUTO=1
      shift
      ;;
    -t|--threads)
      THREADS=$2
      shift 2
      ;;
    -h|--help)
      echo -e "\nEPhIMM v 0.1 - Express Phylogenetic Inference based on Multiple Markers"
      echo -e "Written by Aleksei Korzhenkov, 2019-2020"
      echo -e "For new versions please check https://github.com/laxeye/EPhIMM\n"
      echo -n "Usage: $0 [-p|--predict-proteins] [-e|--evalue <N>] [-g|--max-gap-fraction <0-100>] [-o|--output <directory>] "
      echo -e "[-m|--markers <directory>] [-a|--auto] [-t|--threads <N>]"
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
      echo "Possible flags: -m, -e, -g, -o, -p, -a " >&2
      exit 1
      ;;
    *) # preserve positional arguments
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac
done


#Using half of available CPU threads for hmmsearch
if [[ -z $THREADS ]]; then
  THREADS=1
fi
maxthreads=`grep -c "^processor" /proc/cpuinfo`
if [[ $THREADS -gt $maxthreads ]]; then
    THREADS=$maxthreads
fi

if [[ -d $OUTPUTFOLDER ]]; then
    rm -r $OUTPUTFOLDER
fi
mkdir -p $OUTPUTFOLDER
logfile="$OUTPUTFOLDER/fetching.log"

echo "#Analysys started at $(date)" | tee $logfile

# Predict CDS
if [[ $PREDICT -eq 1 ]]; then
    echo -e "#$(date +"%T")\tProtein prediction started" | tee -a $logfile
    for file in *fna;
      do aafile=${file/fna/faa};
      if ! [[ -f $aafile ]]; then
        prodigal -o /dev/null -i $file -a $aafile -q;
      fi;
    done
fi

# Indexing fasta files
echo -e "#$(date +"%T")\tIndexing files" | tee -a $logfile
for file in *faa; do esl-sfetch --index $file > /dev/null; done

echo "#Directory with HMM profiles of marker genes: $MRKDIR" | tee -a $logfile

markercount=`ls $MRKDIR/*hmm | wc -l`
genomecount=`ls *faa | wc -l`

if [[ $markercount -eq 0 ]]; then
    echo "No HMM profiles found in $MRKDIR" | tee -a $logfile
else
    echo "#$markercount HMM profiles found" | tee -a $logfile
fi

echo -e "#$(date +"%T")\tSearching and extracting marker genes" | tee -a $logfile

for marker in $(ls $MRKDIR/);
    do if [[ -d $marker ]]; then
        rm -r -f $marker
    fi
    mkdir $marker;
    mrkpath=$MRKDIR$marker

    for file in *faa;
    do hmmsearch --cpu $THREADS -E $HMMEVALUE -o /dev/null --tblout /dev/stdout $mrkpath $file | grep -v "^#" > $marker/${file/faa/out};

    lines=`wc -l < $marker/${file/faa/out}`

    case $lines in
        0)
            echo -e "Marker gene \"$marker\" not found in $file!" | tee -a $logfile;
    		    echo -e ">$marker.${file/.faa/}\nXXXXXX" > $marker/$file;
            echo "$marker.${file/.faa/}" >> $marker/missing.txt;
        		;;
        1)
            esl-sfetch -o $marker/$file -n $marker.${file/.faa/} $file $(awk '{print $1}' $marker/${file/faa/out}) > /dev/null;
            ;;
        *)
            echo -e "Marker gene \"$marker\" has many copies in $file! You should inspect sequences manually!" | tee -a $logfile;
            echo "$marker.${file/.faa/}" >> $marker/multicopy.txt;
            let n=1
            for line in $(awk '{print $1}' $marker/${file/faa/out})
                do esl-sfetch -n $marker.${file/.faa/}.$n $file $line >> $marker/$file;
                ((n++));
            done;
            ;;
    esac

    done;
done;


for marker in $(ls $MRKDIR/); do
  if [[ -f $marker/multicopy.txt ]] ; then 
    cat $marker/*faa > $marker/$marker.1.faa ;
    if [[ -f $marker/missing.txt ]] ; then
      java -jar $EPHIMMPATH/BioKotlin.jar RemoveByName $marker/$marker.1.faa $marker/missing.txt > $marker/$marker.faa ;
      rm -f $marker/$marker.1.faa
    else
      mv $marker/$marker.1.faa $marker/$marker.faa ;
    fi
    mafft --quiet --auto $marker/$marker.faa > $marker/aln.faa ;
    java -jar $EPHIMMPATH/BioKotlin.jar DistMeanProtJC $marker/aln.faa | grep -f $marker/multicopy.txt > $marker/multicopy.dist ;
    rm -f $marker/aln.faa
    for mcopy in $(cat $marker/multicopy.txt); do
      grep $mcopy $marker/multicopy.dist | sort -nk 2 | cut -f 1 | head -1 > $marker/keep.txt
      genome=${mcopy/$marker./}
      java -jar $EPHIMMPATH/BioKotlin.jar ExtractByName $marker/$genome.faa $marker/keep.txt > $marker/$genome.1.faa
      mv $marker/$genome.1.faa $marker/$genome.faa
    done
  fi
done

echo -e "#$(date +"%T")\tCollecting statistics" | tee -a $logfile

echo -e "#Genome\tMissing markers\tMultiple copies" > $OUTPUTFOLDER/genomes.stats.txt
for file in *faa;
    do missing=`grep $file $logfile | grep -c "not found"`;
    multi=`grep $file $logfile | grep -c "many copies"`;
    echo -e "$file\t$missing\t$multi" >> $OUTPUTFOLDER/genomes.stats.txt;
done;

echo "#Genomes with more than 50% missing marker genes:" > $OUTPUTFOLDER/bad.genomes.txt
grep -v "^#" $OUTPUTFOLDER/genomes.stats.txt | awk '$2 > "'"$markercount"'"/2' >> $OUTPUTFOLDER/bad.genomes.txt
echo "#Genomes with more than 20% multi-copy marker genes:" >> $OUTPUTFOLDER/bad.genomes.txt
grep -v "^#" $OUTPUTFOLDER/genomes.stats.txt | awk '$3 > "'"$markercount"'"/5' >> $OUTPUTFOLDER/bad.genomes.txt

echo -e "#Marker\tMissing markers\tMultiple copies" > $OUTPUTFOLDER/hmm.stats.txt
for marker in $(ls $MRKDIR/);
    do missing=`grep "$marker" $logfile | grep -c "not found"`
    multi=`grep "$marker" $logfile | grep -c "many"`
    echo -e "$marker\t$missing\t$multi" >> $OUTPUTFOLDER/hmm.stats.txt;
done

echo "#Markers missing in more than 50% genomes:" > $OUTPUTFOLDER/bad.hmms.txt
grep -v "^#" $OUTPUTFOLDER/hmm.stats.txt | awk '$2 > "'"$genomecount"'"/2' >> $OUTPUTFOLDER/bad.hmms.txt
echo "#Markers with multiple copies in more than 20% genomes:" >> $OUTPUTFOLDER/bad.hmms.txt
grep -v "^#" $OUTPUTFOLDER/hmm.stats.txt | awk '$3 > "'"$genomecount"'"/5' >> $OUTPUTFOLDER/bad.hmms.txt

for file in *faa;
do
    if ! grep -q "$file" $OUTPUTFOLDER/bad.genomes.txt ; then
        echo ">${file/.faa/}" > $OUTPUTFOLDER/$file;
        cat *hmm/$file | sed 's/\*//' | grep -v ">" | grep -v "^X\+$"  >> $OUTPUTFOLDER/$file;
    else
        echo "Genome $file will not be used because of incompleteness or contamination."
    fi
done

cat $OUTPUTFOLDER/*faa > $OUTPUTFOLDER/all.markers.faa

echo "Concatenated markers stored in \"$OUTPUTFOLDER/all.markers.faa\""

echo -e "#$(date +"%T")\tPerforming MSA" | tee -a $logfile

mafft --thread -1 --auto  $OUTPUTFOLDER/all.markers.faa | sed 's/X/-/g' > $OUTPUTFOLDER/all.aligned.faa

alnfilename="$OUTPUTFOLDER/all.aligned.e.faa"

if [[ $AUTO -eq 1 ]]; then

  java -jar $EPHIMMPATH/BioKotlin.jar AlignmentClearGaps $OUTPUTFOLDER/all.aligned.faa $MAXGAPFRACTION > $alnfilename
  echo -e "#$(date +"%T")\tPerforming autoedit of MSA" | tee -a $logfile

else

  echo "It's strictly recommended to check alignment file \"$OUTPUTFOLDER/all.aligned.faa\""
  echo "You may automaticly delete all columns with more than $MAXGAPFRACTION% gaps"

  read -p "Press \"N\" if You like to edit file manually.
  Press \"Y\" if You like to clean gaps automaticly: " response

  if [[ $response -eq "Y" ]]; then
      java -jar $EPHIMMPATH/BioKotlin.jar AlignmentClearGaps $OUTPUTFOLDER/all.aligned.faa $MAXGAPFRACTION > $alnfilename
      echo -e "#$(date +"%T")\tPerforming autoedit of MSA" | tee -a $logfile
  else
      read -p "Provide a filename of edited alignment:" alnfilename
  fi

fi

echo -e "#$(date +"%T")\tPerforming phylogenetic inference" | tee -a $logfile

fasttree -gamma -lg $alnfilename > $OUTPUTFOLDER/all.markers.nwk

echo -e "#$(date +"%T")\tWorkflow finished" | tee -a $logfile

