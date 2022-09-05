# Rekonstruktion der Verwandschaftsbeziehung von Chordatieren mit Hilfe genomischer Daten

In diesem bioinformatischen Projekt benutzen wir Sequenzdaten aus der [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/)-Datenbank des [National Center for Biotechnology Information (NCBI)](https://www.ncbi.nlm.nih.gov/), welche eine sehr umfangreiche Sammlung an kuratierten Genomdaten unterschiedlichster Organismen enthält. Diese Datenbank, die aufgrund des technischen Fortschritts auf dem Gebiet der Genomsequenzierung rapide wächst und Stand 2022 weit über 100,000 komplette Genome enthält. erlaubt es Sequenzdaten sowohl der genomischen DNA, als auch der Aminosäuresequenz von transkribierten Genen zu beziehen. In diesem Fall werden DNA-codons, das sind Dreiergruppen von aufeinanderfolgenden Nukleotiden die die Bausteinen von Proteinen codieren, bioinformatisch in die entsprechenden Aminosäuren übersetzt. Die Sequenz der Aminosäuren eignet sich gegenüber der Nukleotid-sequenz von genomischer DNA besser, um weit entfernte Verwandschaften zu rekonstruieren. In unserem Fall fokussieren wir uns auf das Genom von Mitochondrien, welche Zellorganellen darstellen, die bakteriellen Ursprungs sind und in allen eukaryotischen Lebewesen die Energieversorgung der Zellen gewährleisten.

![Mito](https://upload.wikimedia.org/wikipedia/commons/6/64/Cell_structure_%2813080952404%29.jpg)

Zur Rekonstruktion der Verwandschaftsbäume werden leistungsstrake Computer mit ausreichend Arbeitsspeicher (>100GB) und starken Prozessen mit vielen Prozessorkernen (>24 cores) zur parallelisierten Analyse benötigt. Ausserdem sollte das benutzte Betriebssystem auf UNIX basieren, damit die einzelnen Befehle über die Commandline an den Computer übermittelt werden können. Wer sich mit UNIX und dem Grundlagen bioinformtischen Arbeitens näher beschäftigen möchte findet im Internet viele nützliche Tutorials, wie zum Beispiel [hier](http://www.ee.surrey.ac.uk/Teaching/Unix/index.html) und [hier](<>).

Im folgenden stellen wir die fundamentalen Analyseschritte zur Erstellung einer Phylogenie der Chordatiere (wie in Vitrine dargestellt). Da die Zahl der Organismen in der RefSeq Datenbank ständig zunimmt, können die Ergebnisse vom gezeigten Baum abweichen.

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

Der Datensatz `mitochondrion.1.protein.faa.gz` enthält nun die Aminosäuresequenz von mitochondrialer Gene von allen verfügbaren Organismen. Mit Hilfe eines eignes für diesen Zweck geschriebenen Pythonskripts isolieren wir nun Sequenzendaten der Chrordatieren und verwerfen alle anderen Daten. Ausserdem
