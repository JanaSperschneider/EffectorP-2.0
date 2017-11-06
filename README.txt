-----------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------
EffectorP: predicting fungal effector proteins from secretomes using machine learning
Copyright (C) 2015-2016 Jana Sperschneider	
Contact: jana.sperschneider@csiro.au
-----------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------
Installation instructions for EffectorP 1.0
-----------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------
EffectorP relies on two tools, the EMBOSS software and the WEKA 3.6 software. These have been shipped 
with the EffectorP 1.0.tar.gz file, but they need to be installed by the user. 

1) Extract the EffectorP 1.0 archive:

-----------------------------------------
tar xvf EffectorP_1.0.tar.gz
cd EffectorP_1.0 
-----------------------------------------

2) Install EMBOSS

-----------------------------------------
cd Scripts
tar xvf emboss-latest.tar.gz 
cd EMBOSS-6.5.7/
./configure
make
cd ../
-----------------------------------------

3) Install WEKA: simply unzip the file weka-3-6-12.zip

-----------------------------------------
unzip weka-3-6-12.zip
-----------------------------------------

3) Run EffectorP

To test that EffectorP is working, type the following command in the working directory EffectorP_1.0/Scripts

-----------------------------------------
python EffectorP.py -i Effector_Testing.fasta
-----------------------------------------

Note that EffectorP runs under Python 2.x, not under Python 3.x.
If you are having troube installing EMBOSS, please see here for help: http://emboss.sourceforge.net/download/
If you are having troube installing WEKA, please see here for help: http://www.cs.waikato.ac.nz/~ml/weka/index.html

