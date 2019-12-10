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

#### Citation for EffectorP 2.0:
 
 Sperschneider J, Dodds PN, Gardiner DM, Singh KB, Taylor JM (2018) Improved prediction of fungal effector proteins from secretomes with EffectorP 2.0. Molecular Plant Pathology. Abstract

#### Running EffectorP

You can submit secreted fungal proteins to the webserver at http://effectorp.csiro.au/.

Alternatively, you can install EffectorP on your machine to run it locally. For detailed installation instructions see here: http://effectorp.csiro.au/software.html
