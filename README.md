# Rekonstruktion der Verwandtschaftsbeziehungen von Chordatieren mit Hilfe genomischer Daten

Zum Erstellen eines "Stammbaums" (genauer eines Verwandtschaftsbaums) der Chordatiere in dieser Vitrine benutzen wir Sequenzdaten aus der öffentlich zugänglichen [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/)-Datenbank des [National Center for Biotechnology Information (NCBI)](https://www.ncbi.nlm.nih.gov/), welche eine sehr umfangreiche Sammlung an kuratierten Genomdaten unterschiedlichster Organismen enthält. Diese Datenbank, die aufgrund des technischen Fortschritts auf dem Gebiet der Genomsequenzierung [rapide wächst](https://www.ncbi.nlm.nih.gov/refseq/statistics/) und Stand 2022 weit über 100.000 Genome enthält, erlaubt es, sowohl genomische DNA-Sequenzen als auch Aminosäuresequenzen von transkribierten Genen zu beziehen. In letzterem Fall werden DNA-Codons, das sind Dreiergruppen von aufeinanderfolgenden Nukleotiden, die die Bausteine von Proteinen kodieren, bioinformatisch in die entsprechenden Aminosäuren übersetzt. Die Sequenz der Aminosäuren eignet sich gegenüber der Nukleotidsequenz von genomischer DNA besser, um weit entfernte Verwandtschaften zu rekonstruieren. In unserem Fall fokussieren wir auf das Genom von Mitochondrien, das sind Zellorganellen in eukaryotischen Lebewesen, welche die Energieversorgung der Zellen gewährleisten.

![Mito](https://upload.wikimedia.org/wikipedia/commons/6/64/Cell_structure_%2813080952404%29.jpg)

Zur Rekonstruktion der Verwandtschaftsbäume werden leistungsstarke Computer mit ausreichendem Arbeitsspeicher (>100GB) und starken Computerchips mit vielen Prozessorkernen (>100 Cores) zur parallelisierten Analyse benötigt. Außerdem sollte das benutzte Betriebssystem auf UNIX basieren, damit die einzelnen Befehle über die Commandline an den Computer übermittelt werden können. Wer sich mit UNIX und dem Grundlagen bioinformatischen Arbeitens näher beschäftigen möchte, findet im Internet viele nützliche Tutorials und Ressourcen, wie zum Beispiel  [hier](http://www.ee.surrey.ac.uk/Teaching/Unix/index.html) und [hier](https://practicalcomputing.org/).

Im Folgenden stellen wir die fundamentalen Analyseschritte zur Erstellung einer Phylogenie der Chordatiere (wie in der Vitrine dargestellt) vor. Da die Zahl der Organismen in der RefSeq Datenbank ständig zunimmt, können die Ergebnisse vom bereits generierten Baum in der Vitrine abweichen.

## (1) Download der Sequenzdaten

Zunächst wird ein Datensatz mit sämtlichen verfügbaren mitochondrialen Daten heruntergeladen. Dazu fokussieren wir, wie bereits eingangs erwähnt, auf Aminosäuresequenzen sämtlicher Protein-codierender Gene in den mitochondrialen Genomen. Wir laden zusätzlich die DNA-Sequenzen herunter. Basierend auf der Metainformation der einzelnen Sequenzen erstellen wir eine Datenbank mit den Namen aller im Datensatz vorhandenen Organismen (wir sprechen von "Taxon", Mehrzahl "Taxa"), wobei die verschiedenen taxonomischen Ebenen inkludiert sind (wie, zum Beispiel Art, Gattung, Familie, Ordnung, etc.).

Benötigte zusätzliche Programme:

-   [NCBI-edirect](https://www.ebi.ac.uk/Tools/msa/muscle/Tools/NCBIedirect)

```bash
## make directory
mkdir ~/PhylogenyDIY/data

## go to directory
cd ~/PhylogenyDIY/data

## download aminoacid sequence dataset
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.protein.faa.gz

## download DNA sequence dataset
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.1.genomic.fna.gz

## simplify the header of the DNA FASTA file and only retain the name of the Species
gunzip -c mitochondrion.1.1.genomic.fna.gz \
  | awk '{if ($1~"^>") {print $1"_"$2"_"$3} else {print}}' \
  > mitochondrion.1.1.genomic_fixed.fasta

gunzip -c mitochondrion.1.1.genomic.fna.gz \
  | awk '{if ($1~"^>") {print substr($1,2)}}' \
  > mitochondrion.1.1.genomic_fixed.list

### now get the taxonomy table for each RefSeq entry
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

## convert the table to rename the tips of the trees with species names
awk -F "," '{split($1,a,"."); print a[1]".:"$3}' mitochondrion.1.1.genomic_fixed_taxon.list \
  > mitochondrion.1.1.genomic_fixed_taxon.txt

## isolate the class taxonomic level to color code the tree
awk -F "," '$1!~/""/{split($3,a," "); print a[1]"_"a[2]"\t"$6}' mitochondrion.1.1.genomic_fixed_taxon.list \
  | grep -v "^_"  > mitochondrion.1.1.genomic_fixed_taxon.colors
```

## (2) Filtern des Datensatzes

Der Datensatz `mitochondrion.1.protein.faa.gz` enthält nun die Aminosäuresequenzen mitochondrialer Gene aller derzeit verfügbaren Taxa. Mit Hilfe eines eigens für diesen Zweck geschriebenen Skripts  ([proteins2genome.py](scripts/proteins2genome.py)) in der Programmiersprache _Python_ isolieren wir nun Sequenzdaten der Chordatiere und verwerfen alle anderen Daten. Außerdem reduzieren wir den Datensatz auf Gene, die in mindestens 90% aller Taxa vorkommen, und verbinden die Aminosäuresequenzen der einzelnen Gene in der gleichen Reihenfolge zu einer langen Einzelsequenz für jedes Taxon in dem Datensatz.

```bash
## make new directory
mkdir -p ~/PhylogenyDIY/results

cd ~/PhylogenyDIY

## here, we reduce the FASTA file to contain only genes that are present in 90% of all taxa that belong to the Chordates and concatenate the sequence of all retained genes per taxon
python ~/PhylogenyDIY/scripts/proteins2genome.py \
  --TaxList data/mitochondrion.1.1.genomic_fixed_taxon.list \
  --Tax Chordata \
  --FreqTH 0.90 \
  --input data/mitochondrion.1.protein.faa.gz  \
  > results/mitochondrion.1.protein_Chordata.fasta
```

## (3) Alignment der Sequenzen

Nun sind die Aminosäuresequenzen der einzelnen Taxa entsprechend der Gene geordnet und zu einer langen Kette von Aminosäuren verbunden. Der Sequenzvergleich wird durch ein sogenanntes Alignment bewerkstelligt, in dem die Sequenzen der einzelnen Taxa untereinandergeschrieben werden. Änderungen in der Aminosäuresequenz (= ein Austausch einzelner Aminosäuren) werden so sichtbar. Nun müssen noch Mutationen, die zu unterschiedlichen Längen führen, berücksichtigt werden. Insertionen oder Deletionen von Aminosäuren würden im Sequenzvergleich zur Verschiebung der relativen Positionen von Aminosäuren gleichen Ursprungs führen. Positionen mit "Löchern" aufgrund von Deletionen werden mit einem `-` Symbol (man nennt dies "gap") aufgefüllt. Aber auch Mutationen, die zu zusätzlichen Aminosäuren führen (Insertionen) müssen im Alignment durch die Einführung von "Gaps" kompensiert werden. Dieser kritische Schritt ist notwendig, um Merkmale gleichen Ursprungs, in unserem Fall Aminosäuren an der gleichen Sequenzposition, über Taxa hinweg miteinander zu vergleichen und so genetische Unterschiede erkennen zu können. Die zugrunde liegenden Algorithmen sind sehr rechenintensiv, weshalb wir die Rechenlast auf 200 Prozessorkerne aufteilen.

MAFFT verändert die Metainformation der ursprünglichen Daten, was zu einem Problem führt, wenn die ursprünglichen Taxon-Namen für die Benennung der Astenden des Stammbaums benutzt werden sollen. Ein eigenes Skript ([fixIDAfterMafft.py](scripts/fixIDAfterMafft.py)) wird deshalb benutzt, um die ursprüngliche Metainformation wieder herzustellen. In einem weiteren Schritt ([reduceAln2FASTA.py](scripts/reduceAln2FASTA.py)) entfernen wir Positionen, an denen mehr als 50% aller Individuen ein `-` Symbol enthalten. Dies soll verhindern, dass eine hohe Anzahl an nicht-informativen Positionen die Rekonstruktion der Verwandtschaftsverhältnisse verzerrt.

Benötigte zusätzliche Programme:

-   [MAFFT](https://mafft.cbrc.jp/alignment/software/)

```bash
conda activate mafft-7.487

cd ~/PhylogenyDIY

## carry out the alignment with MAFFT
mafft \
  --thread 200 \
  --auto \
  results/mitochondrion.1.protein_Chordata.fasta \
  > results/mitochondrion.1.protein_Chordata_aln_full.fasta

## fix ID's after MAFFT alignment
python  ~/PhylogenyDIY/scripts/fixIDAfterMafft.py \
  --Alignment ~/PhylogenyDIY/results/mitochondrion.1.protein_Chordata_aln_full.fasta \
  --input ~/PhylogenyDIY/results/mitochondrion.1.protein_Chordata.fasta \
  > ~/PhylogenyDIY/results/mitochondrion.1.protein_Chordata_aln_full_fixed.fasta

## only retain Position where less than 50% of all taxa have gaps
python ~/PhylogenyDIY/scripts/reduceAln2FASTA.py \
  --input results/mitochondrion.1.protein_Chordata_aln_full_fixed.fasta  \
  --threshold 0.5 \
  > results/mitochondrion.1.protein_Chordata_aln.fasta

## replace ambiguous AA with gaps
sed -i '/^>/! s/[BJZX]/\-/g' ~/PhylogenyDIY/results/mitochondrion.1.protein_Chordata_aln.fasta
```

Das Hintergrundbild in der Vitrine und das untere Bild sind Beispiele für ein solches Alignment basierend auf DNA-Sequenzen.

![Alignment](data/Alignment_small.jpg)

## (4) Stammbaum-Rekonstruktion mit Hilfe des Maximum-Likelihood-Verfahrens

Der aufbereitete Datensatz kann nun benutzt werden, um die Verwandtschaftsbeziehungen der einzelnen Taxa anhand von Sequenz-Unterschieden abzuschätzen. Dazu wird ein Maximum-Likelihood-Verfahren angewandt, bei der die Wahrscheinlichkeiten verschiedener Baum-Topologien, also möglicher Verwandtschaftsbeziehungen zwischen den einzelnen Taxa, verglichen werden. Ziel ist es, den Stammbaum zu finden, welcher am besten zu der Sequenzdatenmatrix passt, wenn man spezifische Annahmen zum evolutionären Ablauf macht, z.B. wie häufig eine Aminosäure auf Grund von Mutationen durch eine andere ersetzt wird. Dieser Rechenschritt ist äußerst rechenintensiv, da sehr viele verschiedene Baum-Topologien miteinander verglichen werden. Wir verteilen die Rechenlast deshalb wieder auf 200 Prozessorkerne.

Benötigte zusätzliche Programme:

-   [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/)

```bash
module load Phylogeny/RAxML-2.8.10

## make new directory
mkdir ~/PhylogenyDIY/results/raxml
cd ~/PhylogenyDIY/results/raxml

## run ML tree reconstruction
raxmlHPC-PTHREADS-SSE3 \
  -m PROTGAMMAWAG  \
  -N 20 \
  -p 772374015 \
  -n Chordata_const \
  -s ../mitochondrion.1.protein_Chordata_aln.fasta \
  -# 3 \
  -T 200
```

Im Anschluss müssen noch die Taxon-IDs aus den ursprünglichen Metadaten mit Hilfe eines zusätzlichen Skripts ([RenameTreeLeaves_new.py](scripts/RenameTreeLeaves_new.py)) durch korrekte Artnamen ersetzt werden. Außerdem bestimmen wir mit Hilfe eines weiteren Skripts ([MakeOutgroup.py](scripts/MakeOutgroup.py)), welche Taxa zu den _Hyperoartia_, _Ascidiacea_ und _Leptocardii_ gehören, die wir als Außengruppen für die visuelle Darstellung des Baums definieren.

```bash
python ~/PhylogenyDIY/scripts/RenameTreeLeaves_new.py \
  --input ~/PhylogenyDIY/results/raxml/RAxML_bestTree.Chordata_const \
  > ~/PhylogenyDIY/results/raxml/RAxML_bestTree_renamed.Chordata_const

outgroup=`python ~/PhylogenyDIY/scripts/MakeOutgroup.py --tree ~/PhylogenyDIY/results/iqtree_const/mitochondrion.1.protein_Chordata_aln.fasta_renamed.parstree --taxa ~/PhylogenyDIY/data/mitochondrion.1.1.genomic_fixed_taxon.list --list Hyperoartia,Ascidiacea,Leptocardii`
```

Schlussendlich erstellen wir mit Hilfe der Programmiersprache _R_ eine visuelle Darstellung des besten Stammbaum-Modells, welches mit RAxML identifiziert wurde

```R
# load necessary R libraries
library('ggtree')
library('gridExtra')
library('ggrepel')
library('ape')
library('ggplot2')
library('phangorn')
library('ggimage')
library('dplyr')
library('plotly')

## load tree file and root with outgroup taxa
tree<-read.tree('~/PhylogenyDIY/results/raxml/RAxML_bestTree_renamed.Chordata_const')
tree<-root(tree,outgroup=c($outgroup))

## load color information for highlighting different orders in different colors
Col=read.table('~/PhylogenyDIY/data/mitochondrion.1.1.genomic_fixed_taxon.colors',
header=F)
colnames(Col)<-c('tip','cat')

## plot tree
phylo.tree <- ggtree(tree,
layout='roundrect',
lwd=.1,
branch.length='none')+
theme_tree2()+
theme_bw()+
xlab('av. subst./site') +
theme(axis.title.y=element_blank(),
  axis.text.y=element_blank(),
axis.ticks.y=element_blank())+
theme(legend.position='bottom') +
scale_colour_discrete('Orders')+
theme(legend.title = element_text(size=10))+
theme(legend.text = element_text(size=8))+
guides(color = guide_legend(override.aes = list(size = 3)))

phylo.tree <- phylo.tree  %<+% Col+
geom_tiplab(aes(color=cat), size = 0.2)

## export tree
ggsave(filename='~/PhylogenyDIY/results/raxml/tree_rect.pdf',
  phylo.tree,
  width=10,
  height=30,limitsize=F)
```

Der Baum, der unterhalb sowie in der Geschichtsvitrine gezeigt wird, wurde mit Hilfe von Graphikprogrammen nachbearbeitet. So wurde, z.B. die ursprüngliche Legende entfernt und durch eine vertikale Beschriftung der Klassen mit deutschen Namen ersetzt. Außerdem wurden die Äste mancher Taxon-Gruppen um die X-Achse rotiert, was den Baum nicht verändert, aber die Lesbarkeit verbessert.

![Tree](data/Tree_rect.jpg)
