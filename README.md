#### What is EffectorP?

Fungal plant pathogens secrete effector proteins that modulate the host cell to facilitate infection. 
Computational effector candidate identification and subsequent functional characterization delivers valuable insights 
into plant-pathogen interactions. However, effector prediction in fungi has been challenging due to a lack of unifying
sequence features such as conserved N-terminal sequence motifs. Fungal effectors are commonly predicted from secretomes 
based on criteria such as small size and cysteine-rich, which suffers from poor accuracy.

EffectorP is a machine learning method for fungal effector prediction in secretomes and has been trained to distinguish secreted proteins 
from secreted effectors in plant-pathogenic fungi.
EffectorP improves fungal effector prediction from secretomes based on a robust signal of sequence-derived properties.
EffectorP 2.0 achieves an accuracy of 89%, compared with 82% for EffectorP 1.0 and 59.8% for a small size classifier.

#### What is EffectorP not?

EffectorP is not a tool for secretome prediction. 

EffectorP has been trained to find fungal effectors in secretomes, 
so please run it on a FASTA file of secreted fungal proteins to test if they are predicted effectors. It is recommended 
to use tools such as SignalP or Phobius	to predict first if a protein is likely to be secreted.
Alternatively, experimentally determined secretomes instead of computationally predicted secretomes can be submitted to EffectorP. 

#### Running EffectorP

You can submit secreted fungal proteins to the webserver at http://effectorp.csiro.au/.

Alternatively, you can install EffectorP on your machine to run it locally. 
 
All training and evaluation data can be found [here](http://effectorp.csiro.au/data.html).

#### Installing EffectorP 

EffectorP has been written in Python and uses pepstats from the EMBOSS software and the WEKA 3.8.1 software. To get EffectorP to work on your local machine, you need to install the EMBOSS and WEKA softwares from source. Both are already provided in the EffectorP distribution to ensure that compatible versions are used. **Effector from version 2.0.1 inclusive uses Python 3.** 

0. Download the latest release from this github repo (or alternatively you can clone the github repo and skip step 1).

1. Make sure EffectorP has the permission to execute. Then unpack LOCALIZER in your desired location
```
tar xvf EffectorP-2.0.1.tar.gz
chmod -R 755 EffectorP-2.0.1/
cd EffectorP-2.0.1
```

2. For the EMBOSS installation, you need to switch to the Scripts directory and unpack, configure and make. Alternatively, if you are on a computer cluster and EMBOSS is already installed, you can change the variable PEPSTATS_PATH in the EffectorP.py script to the EMBOSS directory that contains pepstats on the machine you are using.
```
cd Scripts
tar xvf emboss-latest.tar.gz
cd EMBOSS-6.5.7/
./configure
make
cd ../ 
```

3. For WEKA, you need to simply unzip the file weka-3-8-1.zip
```
unzip weka-3-8-1.zip
```
If you are having troube installing EMBOSS, please see [here](http://emboss.sourceforge.net/download/) for help.
If you are having troube installing WEKA, please see [here](https://www.cs.waikato.ac.nz/~ml/weka/index.html) for help. 

4. Test if EffectorP is working
```
python EffectorP.py -i Effector_Testing.fasta
```

#### EffectorP output format
Run this to get a feel for the output format:
```
python EffectorP.py -i Effector_Testing.fasta
```

#### Citation for EffectorP 
 
Sperschneider J, Dodds PN, Gardiner DM, Singh KB, Taylor JM (2018) Improved prediction of fungal effector proteins from secretomes with EffectorP 2.0. Molecular Plant Pathology. [Link to paper](https://bsppjournals.onlinelibrary.wiley.com/doi/full/10.1111/mpp.12682)
