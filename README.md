# Rekonstruktion der Verwandschaftsbeziehung von Chordatieren mit Hilfe genomischer Daten

Zum Erstellen des Stammbaums der Chordatiere in dieser Vitrine benutzen wir Sequenzdaten aus der [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/)-Datenbank des [National Center for Biotechnology Information (NCBI)](https://www.ncbi.nlm.nih.gov/), welche eine sehr umfangreiche Sammlung an kuratierten Genomdaten unterschiedlichster Organismen enthält. Diese Datenbank, die aufgrund des technischen Fortschritts auf dem Gebiet der Genomsequenzierung [rapide wächst](https://www.ncbi.nlm.nih.gov/refseq/statistics/) und Stand 2022 weit über 100,000 komplette Genome enthält. erlaubt es Sequenzdaten sowohl der genomischen DNA, als auch der Aminosäuresequenz von transkribierten Genen zu beziehen. In diesem Fall werden DNA-codons, das sind Dreiergruppen von aufeinander folgenden Nukleotiden die die Bausteinen von Proteinen codieren, bioinformatisch in die entsprechenden Aminosäuren übersetzt. Die Sequenz der Aminosäuren eignet sich gegenüber der Nukleotidsequenz von genomischer DNA besser, um weit entfernte Verwandtschaften zu rekonstruieren. In unserem Fall fokussieren wir uns auf das Genom von Mitochondrien, welche Zellorganellen darstellen, die bakteriellen Ursprungs sind und in allen eukaryotischen Lebewesen die Energieversorgung der Zellen gewährleisten.

![Mito](https://upload.wikimedia.org/wikipedia/commons/6/64/Cell_structure_%2813080952404%29.jpg)

Zur Rekonstruktion der Verwandschaftsbäume werden leistungsstarke Computer mit ausreichend Arbeitsspeicher (>100GB) und starken Prozessen mit vielen Prozessorkernen (>100 cores) zur parallelisierten Analyse benötigt. Außerdem sollte das benutzte Betriebssystem auf UNIX basieren, damit die einzelnen Befehle über die Commandline an den Computer übermittelt werden können. Wer sich mit UNIX und dem Grundlagen bioinformatischen Arbeitens näher beschäftigen möchte findet im Internet viele nützliche Tutorials und Ressourcen, wie zum Beispiel [hier](http://www.ee.surrey.ac.uk/Teaching/Unix/index.html) und [hier](https://practicalcomputing.org/).

Im Folgenden stellen wir die fundamentalen Analyseschritte zur Erstellung einer Phylogenie der Chordatiere (wie in Vitrine dargestellt). Da die Zahl der Organismen in der RefSeq Datenbank ständig zunimmt, können die Ergebnisse vom bereits generierten Baum in der Vitrine abweichen.

## (1) Download der Sequenzdaten

Zunächst wird ein Datensatz mit sämtlichen verfügbaren mitochondrialen Daten heruntergeladen. Dazu fokussieren wir uns, wie bereits eingangs erwähnt auf Aminosäuresequenzen sämtlicher kodierender Gene in den mitochondrialen Genomen. Wir laden zusätzlich die DNA Sequenzen herunter. Basierend auf der Metainformation der einzelnen Sequenzen erstellen wir eine Datenbank mit den Namen aller vorhandenen Organismen auf verschiedenen taxonomischen levels.

Benötigte zusätzliche Programme:

-   [NCBI-edirect](https://www.ebi.ac.uk/Tools/msa/muscle/Tools/NCBIedirect)

```bash
## make directory
mkdir /media/inter/mkapun/projects/EukMitGenomeTree/data

## go to directory
cd /media/inter/mkapun/projects/EukMitGenomeTree/data

## download aminoacid sequence dataset
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.protein.faa.gz

## download DNA sequence dataset
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.1.genomic.fna.gz

## simplify the header of the FASTA to only retain the name of the Species
gunzip -c mitochondrion.1.1.genomic.fna.gz \
  | awk '{if ($1~"^>") {print $1"_"$2"_"$3} else {print}}' \
  > mitochondrion.1.1.genomic_fixed.fasta

gunzip -c mitochondrion.1.1.genomic.fna.gz \
  | awk '{if ($1~"^>") {print substr($1,2)}}' \
  > mitochondrion.1.1.genomic_fixed.list

### now get the Taxonomy Table for each RefSeq
module load Tools/NCBIedirect
while read -r line
do
  #echo $line
  ID=`esearch -db nucleotide -query ${line} < /dev/null |esummary | xtract -pattern TaxId  -element TaxId `
  ID2=`efetch -db taxonomy -id ${ID} -format xml | xtract -pattern Taxon -tab "," -first TaxId ScientificName \
    -group Taxon -KING "(-)" -PHYL "(-)" -CLSS "(-)" -ORDR "(-)" -FMLY "(-)" -GNUS "(-)" \
    -block "*/Taxon" -match "Rank:kingdom" -KING ScientificName \
    -block "*/Taxon" -match "Rank:phylum" -PHYL ScientificName \
    -block "*/Taxon" -match "Rank:class" -CLSS ScientificName \
    -block "*/Taxon" -match "Rank:order" -ORDR ScientificName \
    -block "*/Taxon" -match "Rank:family" -FMLY ScientificName \
    -block "*/Taxon" -match "Rank:genus" -GNUS ScientificName \
    -group Taxon -tab "," -element "&KING" "&PHYL" "&CLSS" "&ORDR" "&FMLY" "&GNUS"`
  echo ${line}","$ID2
done < mitochondrion.1.1.genomic_fixed.list > mitochondrion.1.1.genomic_fixed_taxon.list

## convert the table to rename the tips of the Trees with the Species names
awk -F "," '{split($1,a,"."); print a[1]".:"$3}' mitochondrion.1.1.genomic_fixed_taxon.list \
  > mitochondrion.1.1.genomic_fixed_taxon.txt

## Isolate the class taxonomic level to color code the tree
awk -F "," '$1!~/""/{split($3,a," "); print a[1]"_"a[2]"\t"$6}' mitochondrion.1.1.genomic_fixed_taxon.list \
  | grep -v "^_"  > mitochondrion.1.1.genomic_fixed_taxon.colors
```

## (2) Filtern des Datensatzes

Der Datensatz `mitochondrion.1.protein.faa.gz` enthält nun die Aminosäuresequenz mitochondrialer Gene von allen verfügbaren Organismen. Mit Hilfe eines eigens für diesen Zweck geschriebenen Skripts ([proteins2genome.py](scripts/proteins2genome.py)) in der Programmiersprache _Python_ isolieren wir nun Sequenzendaten der Chrordatieren und verwerfen alle anderen Daten. Außerdem reduzieren wir den Datensatz auf Gene, die in mindestens 90% aller Taxa vorkommen. Außerdem verbinden wir für jedes Taxon die Aminosäurensequenzen der einzelnen Gene in der immer gleichen Reihenfolge zu einer langen Einzelsequenz.

```bash
## here, we reduce the FASTA file to contain only genes that are present in 90% of all taxa that belong to the Chordates
mkdir -p /media/inter/mkapun/projects/EukMitGenomeTree/results_AA/Chordata


cd /media/inter/mkapun/projects/EukMitGenomeTree

python /media/inter/mkapun/projects/EukMitGenomeTree/scripts/proteins2genome.py \
  --TaxList data/mitochondrion.1.1.genomic_fixed_taxon.list \
  --Tax ${Taxon} \
  --FreqTH 0.90 \
  --input data/mitochondrion.1.protein.faa.gz  \
  > results_AA/${Taxon}/mitochondrion.1.protein_${Taxon}.fasta
```

## (3) Alignment der Sequenzdaten

Nun sind die Aminosäuresequenzen der einzelnen Taxa entsprechend der Gene geordnet und zu einer langen Kette von Aminosäuren verbunden. Diese Kette berücksichtigt aber keine "strukturellen" Mutationen, wie Insertionen oder Deletionen, die zur Verschiebung der relativen Positionen der gleichen Aminosäuren in verschiedenen Taxa führen kann. Durch ein sogenanntes Alignment werden zusammengehörige Aminosäuren in Spalten übereinander geordnet, wobei jede Zeile die Sequenzdaten eines Taxons enthält. Positionen mit "Löchern" aufgrund von strukturellen Mutationen werden mit einem `-` Symbol aufgefüllt. Dieser kritische Schritt ist notwendig um gleiche Merkmale, in unserem Fall, Aminosäuren an der gleichen Sequenzposition, über Taxa hinweg miteinander zu vergleichen und so genetische Unterschiede zu erkennen. Die zugrunde liegenden Algorithmen sind sehr rechenintensiv, weshalb wir die Rechenlast auf 200 Prozessorkerne aufteilen.  

MAFFT verändert die Metainformation der ursprünglichen Daten, was zu einem Problem führt, wenn die ursprünglichen Taxon-Namen für die Benennung der Astenden des Stammbaums benutzt werden sollen. Ein eigenes Skript ([fixIDAfterMafft.py](scripts/fixIDAfterMafft.py)) wird deshalb benutzt um die ursprüngliche Metainformation wieder herzustellen. In einem weiteren Schritt entfernen wir mit einem weiteren Schritt ([reduceAln2FASTA.py](scripts/reduceAln2FASTA.py)) sämtliche Positionen an den mehr als 50% aller Individuen ein `-` Symbol enthalten. Dies soll Verhindern, dass eine hohe Anzahl an nicht informativen Positionen die Rekonstruktion von Verwandtschaftsverhältnisse im nächsten Schritt verzerrt.  

Benötigte zusätzliche Programme:

-   [MAFFT](https://mafft.cbrc.jp/alignment/software/)

```bash
## carry out the alignment with MAFFT
mafft \
  --thread 200 \
  --auto \
  results_AA_stringent/Chordata/mitochondrion.1.protein_Chordata.fasta \
  > results_AA_stringent/Chordata/mitochondrion.1.protein_Chordata_aln_full.fasta

## Fix ID's after MAFFT :-(((((
python  /media/inter/mkapun/projects/EukMitGenomeTree/scripts/fixIDAfterMafft.py \
  --Alignment /media/inter/mkapun/projects/EukMitGenomeTree/results_AA_stringent/Chordata/mitochondrion.1.protein_Chordata_aln_full.fasta \
  --input /media/inter/mkapun/projects/EukMitGenomeTree/results_AA_stringent/Chordata/mitochondrion.1.protein_Chordata.fasta \
  > /media/inter/mkapun/projects/EukMitGenomeTree/results_AA_stringent/Chordata/mitochondrion.1.protein_Chordata_aln_full_fixed.fasta

## Only retain Position where less than 50% of all taxa have gaps
python /media/inter/mkapun/projects/EukMitGenomeTree/scripts/reduceAln2FASTA.py \
  --input results_AA_stringent/Chordata/mitochondrion.1.protein_Chordata_aln_full_fixed.fasta  \
  --threshold 0.5 \
  > results_AA_stringent/Chordata/mitochondrion.1.protein_Chordata_aln.fasta

## replace ambiguous AA with gaps
sed -i '/^>/! s/[BJZX]/\-/g' /media/inter/mkapun/projects/EukMitGenomeTree/results_AA_stringent/Chordata/mitochondrion.1.protein_Chordata_aln.fasta
```

Das Hintergrundbild in der Vitrine und das untere Bild sind Beispiele für ein solches Alignment mit DNA-Sequenzen

![Alignment](data/Alignment_small.jpg)

## (4) Stammbaum-Rekonstruktion mit Hilfe des Maximum-Likelihood Verfahrens

Der aufbereitete Datensatz kann nun benutzt werden um die Verwandtschaftsbeziehung der einzelnen Taxa anhand von Sequenzunterschiede abzuschätzen. Dazu wird ein maximum-likelihood Verfahren angewandt, bei der die Wahrscheinlichkeit verschiedener Baum-topologien, also Verwandtschaftsbeziehungen zwischen den einzelnen Taxa, verglichen werden. Ziel ist es den Stammbau zu finden, der am besten zu der Datenmatrix passt wenn man spezifische Annahmen zum Ablauf der Evolution macht, z.B.: wie häufig eine Aminosäure auf Grund von Mutationen durch eine andere ersetzt wird. Dieser Rechenschritt ist äußerst rechenintensiv, da sehr viele verschiedenen Baumtopologien miteinander verglichen werden.

Benötigte zusätzliche Programme:

-   [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/)

```bash

```
