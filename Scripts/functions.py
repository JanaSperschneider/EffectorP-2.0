#!/usr/bin/env python
"""
    Improved prediction of fungal effector proteins from secretomes with EffectorP 2.0
    Copyright (C) 2017-2018 Jana Sperschneider	

    This program is free software; you can redistribute it and/or modify  
    it under the terms of the GNU General Public License as published by  
    the Free Software Foundation; either version 3 of the License, or     
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Please also see the CSIRO Disclaimer provided with EffectorP (LICENCE.txt).

    Contact: jana.sperschneider@csiro.au
"""
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
import os
import sys
import re
import io
import getopt
import regex as re
import random
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
def usage():
    """ Function: usage()

        Purpose:  Print helpful information for the user.        
        
        Input:    None.
    
        Return:   Print options for running EffectorP 2.0.       
    """
    print '''
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# EffectorP :: predicting fungal effector proteins from secretomes using machine learning
# EffectorP 2.0 (November 2017); http://effectorp.csiro.au/
# Copyright (C) 2017-2018 Jana Sperschneider, CSIRO.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    '''
    print "Usage for EffectorP 2.0: ", 
    print "python EffectorP.py [-options] -i <input_file>"
    print 
    print "where basic options are:"
    print "-h : show brief help on version and usage" 
    print 
    print "options for output format:"
    print "-s : short output format that provides predictions for all proteins as one tab-delimited table [default long format]"
    print
    print "options directing output:"
    print "-o <f> : direct output to file <f>, not stdout"
    print "-E <f> : save predicted effectors to FASTA file <f>"        
    print "-N <f> : save predicted non-effectors to FASTA file <f>"        
    print
    print "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
    print
    sys.exit()    

    return
# -----------------------------------------------------------------------------------------------------------
def scan_arguments(commandline):
    """ Function: scan_arguments()

        Purpose:  Scan the input options given to the EffectorP program.        
        
        Input:    Input options given by the user.
    
        Return:   Parsed options.
    """
    try:
        opts, args = getopt.getopt(commandline, "hso:E:N:i:", ["help"])        
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    FASTA_FILE = None
    short_format = False
    output_file = None
    effector_output = None
    noneffector_output = None

    i_count, o_count, E_count, N_count, P_count = 0, 0, 0, 0, 0
   
    for opt, arg in opts:
        if opt in ("-o"):
            output_file = arg
            o_count += 1
        elif opt in ("-s"):
            short_format = True
        elif opt in ("-i"):
            FASTA_FILE = arg
            i_count += 1
        elif opt in ("-E"):
            effector_output = arg
            E_count += 1
        elif opt in ("-N"):
            noneffector_output = arg
            N_count += 1
        elif opt in ("-h", "--help"):
            usage()
        else:
            assert False, "unhandled option"

    if i_count > 1 or o_count > 1 or E_count > 1 or N_count > 1:
       usage()

    return FASTA_FILE, short_format, output_file, effector_output, noneffector_output
# -----------------------------------------------------------------------------------------------------------
models_bayes = ['/MODEL_FILES/PathogenSecretomes_Bayes/trainingdata_samegenomes_iteration25_ratio3_bayes.model',
                '/MODEL_FILES/PathogenSecretomes_Bayes/trainingdata_samegenomes_iteration69_ratio3_bayes.model',
                '/MODEL_FILES/PathogenSecretomes_Bayes/trainingdata_samegenomes_iteration79_ratio3_bayes.model',
                '/MODEL_FILES/PathogenSecretomes_Bayes/trainingdata_samegenomes_iteration36_ratio3_bayes.model',
                '/MODEL_FILES/PathogenSecretomes_Bayes/trainingdata_samegenomes_iteration44_ratio3_bayes.model',
                '/MODEL_FILES/PathogenSecretomes_Bayes/trainingdata_samegenomes_iteration49_ratio3_bayes.model',
                '/MODEL_FILES/PathogenSecretomes_Bayes/trainingdata_samegenomes_iteration10_ratio3_bayes.model',
                '/MODEL_FILES/PathogenSecretomes_Bayes/trainingdata_samegenomes_iteration32_ratio3_bayes.model',
                '/MODEL_FILES/PathogenSecretomes_Bayes/trainingdata_samegenomes_iteration47_ratio3_bayes.model',
                '/MODEL_FILES/PathogenSecretomes_Bayes/trainingdata_samegenomes_iteration46_ratio3_bayes.model',

                '/MODEL_FILES/AnimalPathogenSecretomes_Bayes/trainingdata_samegenomes_iteration47_ratio3_bayes.model',
                '/MODEL_FILES/AnimalPathogenSecretomes_Bayes/trainingdata_samegenomes_iteration24_ratio3_bayes.model',
                '/MODEL_FILES/AnimalPathogenSecretomes_Bayes/trainingdata_samegenomes_iteration97_ratio3_bayes.model',
                '/MODEL_FILES/AnimalPathogenSecretomes_Bayes/trainingdata_samegenomes_iteration12_ratio3_bayes.model',
                '/MODEL_FILES/AnimalPathogenSecretomes_Bayes/trainingdata_samegenomes_iteration51_ratio3_bayes.model',

                '/MODEL_FILES/SaprophyteSecretomes_Bayes/trainingdata_samegenomes_iteration85_ratio3_bayes.model',
                '/MODEL_FILES/SaprophyteSecretomes_Bayes/trainingdata_samegenomes_iteration73_ratio3_bayes.model',
                '/MODEL_FILES/SaprophyteSecretomes_Bayes/trainingdata_samegenomes_iteration46_ratio3_bayes.model',
                '/MODEL_FILES/SaprophyteSecretomes_Bayes/trainingdata_samegenomes_iteration10_ratio3_bayes.model',
                '/MODEL_FILES/SaprophyteSecretomes_Bayes/trainingdata_samegenomes_iteration62_ratio3_bayes.model',
                '/MODEL_FILES/SaprophyteSecretomes_Bayes/trainingdata_samegenomes_iteration2_ratio3_bayes.model',
                '/MODEL_FILES/SaprophyteSecretomes_Bayes/trainingdata_samegenomes_iteration70_ratio3_bayes.model',
                '/MODEL_FILES/SaprophyteSecretomes_Bayes/trainingdata_samegenomes_iteration60_ratio3_bayes.model',
                '/MODEL_FILES/SaprophyteSecretomes_Bayes/trainingdata_samegenomes_iteration87_ratio3_bayes.model',
                '/MODEL_FILES/SaprophyteSecretomes_Bayes/trainingdata_samegenomes_iteration21_ratio3_bayes.model']

models_J48 = ['/MODEL_FILES/PathogenSecretomes_Trees/trainingdata_samegenomes_iteration66_ratio3_J48.model',
                '/MODEL_FILES/PathogenSecretomes_Trees/trainingdata_samegenomes_iteration95_ratio3_J48.model',
                '/MODEL_FILES/PathogenSecretomes_Trees/trainingdata_samegenomes_iteration15_ratio3_J48.model',
                '/MODEL_FILES/PathogenSecretomes_Trees/trainingdata_samegenomes_iteration54_ratio3_J48.model',
                '/MODEL_FILES/PathogenSecretomes_Trees/trainingdata_samegenomes_iteration50_ratio3_J48.model',
                '/MODEL_FILES/PathogenSecretomes_Trees/trainingdata_samegenomes_iteration81_ratio3_J48.model',
                '/MODEL_FILES/PathogenSecretomes_Trees/trainingdata_samegenomes_iteration91_ratio3_J48.model',
                '/MODEL_FILES/PathogenSecretomes_Trees/trainingdata_samegenomes_iteration3_ratio3_J48.model',
                '/MODEL_FILES/PathogenSecretomes_Trees/trainingdata_samegenomes_iteration89_ratio3_J48.model',
                '/MODEL_FILES/PathogenSecretomes_Trees/trainingdata_samegenomes_iteration2_ratio3_J48.model',

                '/MODEL_FILES/AnimalPathogenSecretomes_Trees/trainingdata_samegenomes_iteration69_ratio3_J48.model',
                '/MODEL_FILES/AnimalPathogenSecretomes_Trees/trainingdata_samegenomes_iteration9_ratio3_J48.model',
                '/MODEL_FILES/AnimalPathogenSecretomes_Trees/trainingdata_samegenomes_iteration86_ratio3_J48.model',
                '/MODEL_FILES/AnimalPathogenSecretomes_Trees/trainingdata_samegenomes_iteration81_ratio3_J48.model',
                '/MODEL_FILES/AnimalPathogenSecretomes_Trees/trainingdata_samegenomes_iteration29_ratio3_J48.model',     

                '/MODEL_FILES/SaprophyteSecretomes_Trees/trainingdata_samegenomes_iteration35_ratio3_J48.model',
                '/MODEL_FILES/SaprophyteSecretomes_Trees/trainingdata_samegenomes_iteration24_ratio3_J48.model',
                '/MODEL_FILES/SaprophyteSecretomes_Trees/trainingdata_samegenomes_iteration87_ratio3_J48.model',
                '/MODEL_FILES/SaprophyteSecretomes_Trees/trainingdata_samegenomes_iteration13_ratio3_J48.model',
                '/MODEL_FILES/SaprophyteSecretomes_Trees/trainingdata_samegenomes_iteration22_ratio3_J48.model',
                '/MODEL_FILES/SaprophyteSecretomes_Trees/trainingdata_samegenomes_iteration45_ratio3_J48.model',
                '/MODEL_FILES/SaprophyteSecretomes_Trees/trainingdata_samegenomes_iteration97_ratio3_J48.model',
                '/MODEL_FILES/SaprophyteSecretomes_Trees/trainingdata_samegenomes_iteration27_ratio3_J48.model',
                '/MODEL_FILES/SaprophyteSecretomes_Trees/trainingdata_samegenomes_iteration46_ratio3_J48.model',
                '/MODEL_FILES/SaprophyteSecretomes_Trees/trainingdata_samegenomes_iteration28_ratio3_J48.model']
# -----------------------------------------------------------------------------------------------------------
ARFF_HEADER = '''@RELATION effectors
@ATTRIBUTE Tiny NUMERIC
@ATTRIBUTE Small NUMERIC
@ATTRIBUTE Aliphatic NUMERIC
@ATTRIBUTE Aromatic NUMERIC
@ATTRIBUTE Nonpolar NUMERIC
@ATTRIBUTE Polar NUMERIC
@ATTRIBUTE Charged NUMERIC
@ATTRIBUTE Basic NUMERIC
@ATTRIBUTE Acidic NUMERIC
@ATTRIBUTE A NUMERIC
@ATTRIBUTE C NUMERIC
@ATTRIBUTE D NUMERIC
@ATTRIBUTE E NUMERIC
@ATTRIBUTE F NUMERIC
@ATTRIBUTE G NUMERIC
@ATTRIBUTE H NUMERIC
@ATTRIBUTE I NUMERIC
@ATTRIBUTE K NUMERIC
@ATTRIBUTE L NUMERIC
@ATTRIBUTE M NUMERIC
@ATTRIBUTE N NUMERIC
@ATTRIBUTE P NUMERIC
@ATTRIBUTE Q NUMERIC
@ATTRIBUTE R NUMERIC
@ATTRIBUTE S NUMERIC
@ATTRIBUTE T NUMERIC
@ATTRIBUTE V NUMERIC
@ATTRIBUTE W NUMERIC
@ATTRIBUTE Y NUMERIC
@ATTRIBUTE MolecularWeight NUMERIC
@ATTRIBUTE Charge NUMERIC
@ATTRIBUTE GRAVY NUMERIC
@ATTRIBUTE Exposed NUMERIC
@ATTRIBUTE Disorder NUMERIC
@ATTRIBUTE Hydro NUMERIC
@ATTRIBUTE Bulky NUMERIC
@ATTRIBUTE Interface NUMERIC
@ATTRIBUTE class {effector,non-effector}
@DATA
'''
# -----------------------------------------------------------------------------------------------------------
# Taken from http://web.expasy.org/protscale/pscale/Hphob.Doolittle.html
GRAVY_DIC = {
'A': 1.8, 
'R': -4.5,
'N': -3.5,
'D': -3.5,  
'C':  2.5,  
'Q': -3.5,  
'E': -3.5,  
'G': -0.4,  
'H': -3.2,  
'I':  4.5,  
'L':  3.8,  
'K': -3.9,  
'M':  1.9,  
'F':  2.8,  
'P': -1.6,  
'S': -0.8,  
'T': -0.7,  
'W': -0.9,  
'Y': -1.3,  
'V':  4.2}
# -----------------------------------------------------------------------------------------------------------
# Hydrophobicity (Fauchere and Pliska, 1983)
HYDRO_DIC = {
'R': -1.01,
'K': -0.99,  
'D': -0.77,  
'E': -0.64,  
'N': -0.6,
'Q': -0.22,  
'S': -0.04,  
'G': -0.0,  
'H': 0.13,  
'T': 0.26,  
'A': 0.31, 
'P': 0.72,  
'Y': 0.96,  
'V': 1.22,
'C': 1.54,  
'L': 1.7,  
'F': 1.79,  
'I': 1.8,  
'M': 1.23 ,  
'W': 2.25,  
}
# -----------------------------------------------------------------------------------------------------------
# Taken from http://www.cprofiler.org/help.html
# Surface exposure (Janin, 1979), these are free energy values
EXPOSED_DIC = {
'A': 0.3, 
'R': -1.4,
'N': -0.5,
'D': -0.6,  
'C': 0.9,  
'Q': -0.7,  
'E': -0.7,  
'G': 0.3,  
'H': -0.1,  
'I': 0.7,  
'L': 0.5,  
'K': -1.8,  
'M': 0.4,  
'F': 0.5,  
'P': -0.3,  
'S': -0.1,  
'T': -0.2,  
'W': 0.3,  
'Y': -0.4,  
'V': 0.6
}
# Disorder propensity (Dunker et al., 2001)
DISORDER_DIC = {
'A': 1.0,
'R': 1.0,
'S': 1.0,
'Q': 1.0,
'E': 1.0,
'G': 1.0,
'K': 1.0,
'P': 1.0,
'D': 0.0,
'H': 0.0,
'M': 0.0,
'T': 0.0,
'N': -1.0,
'C': -1.0,
'I': -1.0,
'L': -1.0,
'F': -1.0,
'W': -1.0,
'Y': -1.0,
'V': -1.0,
}
# Bulkiness (Zimmerman et al., 1968)
BULKY_DIC = {
'G' : 3.4,
'S' : 9.47,	
'A' : 11.5, 	
'D' : 11.68,	
'N' : 12.82,	
'C' : 13.46,	
'E' : 13.57,	
'H' : 13.69,	
'R' : 14.28,
'Q' : 14.45,	
'K' : 15.71,	
'T' : 15.77,
'M' : 16.25,	
'P' : 17.43,
'Y' : 18.03,	
'F' : 19.8,	
'I' : 21.4,
'L' : 21.4, 	
'V' : 21.57,	
'W' : 21.67 			 			 			 			 		
}
# Interface propensity (Jones and Thornton, 1997)
INTERFACE_DIC = {
'A': -0.17, 
'R': 0.27,
'N': 0.12,
'D': -0.38,  
'C': 0.43,  
'Q': -0.11,  
'E': -0.13,  
'G': -0.07,  
'H': 0.41,  
'I': 0.44,  
'L': 0.4,  
'K': -0.36,  
'M': 0.66,  
'F': 0.82,  
'P': -0.25,  
'S': -0.33,  
'T': -0.18,  
'W': 0.83,  
'Y': 0.66,  
'V': 0.27}
# -----------------------------------------------------------------------------------------------------------
def GRAVY(sequence):
    # The GRAVY value for a peptide or protein is calculated as the sum of hydropathy values of all the amino acids, 
    # divided by the number of residues in the sequence. 

    gravy = 0.0
    for aa in sequence:
	    if aa.upper() in GRAVY_DIC:
	        gravy += GRAVY_DIC[aa.upper()]

    gravy = gravy/len(sequence)

    return gravy
# -----------------------------------------------------------------------------------------------------------
def HYDRO(sequence):

    hydro = 0.0
    for aa in sequence:
	    if aa.upper() in HYDRO_DIC:
	        hydro += HYDRO_DIC[aa.upper()]

    hydro = hydro/len(sequence)

    return hydro
# -----------------------------------------------------------------------------------------------------------
def DISORDER(sequence):
    
    disorder = 0.0
    for aa in sequence:
	    if aa.upper() in DISORDER_DIC:
	        disorder += DISORDER_DIC[aa.upper()]

    disorder = disorder/len(sequence)

    return disorder
# -----------------------------------------------------------------------------------------------------------
def EXPOSED(sequence):
    
    exposed = 0.0
    for aa in sequence:
	    if aa.upper() in EXPOSED_DIC:
	        exposed += EXPOSED_DIC[aa.upper()]

    exposed = exposed/len(sequence)

    return exposed
# -----------------------------------------------------------------------------------------------------------
def BULKY(sequence):
    
    bulky = 0.0
    for aa in sequence:
	    if aa.upper() in BULKY_DIC:
	        bulky += BULKY_DIC[aa.upper()]

    bulky = bulky/len(sequence)

    return bulky
# -----------------------------------------------------------------------------------------------------------
def INTERFACE(sequence):
    
    interface = 0.0
    for aa in sequence:
	    if aa.upper() in INTERFACE_DIC:
	        interface += INTERFACE_DIC[aa.upper()]

    interface = interface/len(sequence)

    return interface
# -----------------------------------------------------------------------------------------------------------
def get_seqs_ids_fasta(FASTA_FILE):
    """ Function: get_seqs_ids_fasta()

        Purpose:  Given a FASTA format file, this function extracts
                  the list of identifiers and the list of sequences 
                  in the order in which they appear in the FASTA file.
              
        Input:    Path to FASTA format file.
    
        Return:   List of identifiers and list of sequences in the order 
                  in which they appear in the FASTA file.
    """ 
    identifiers = []
    sequences = []

    with open(FASTA_FILE) as f: 
        content = f.readlines()

        for position, line in enumerate(content):
            if '>' in line:
                identifiers.append(line)
                seq = []
                following_lines = content[position + 1:]
                for next_line in following_lines:
                    if '>' not in next_line:
                        seq.append(next_line.strip())
                    else:
                        break
                sequence = "".join(seq)
                sequence = sequence.replace('*', '')
                sequences.append(sequence)

    return identifiers, sequences
# -----------------------------------------------------------------------------------------------------------
def write_FASTA_short_ids(f_output, ORIGINAL_IDENTIFIERS, ORIGINAL_SEQUENCES):
    """ Function: write_FASTA_short_ids()

        Purpose:  Given a list of identifiers and the corresponding list 
                  of sequence, write these to a FASTA file using short
                  identifiers such as protein1, protein2, .... This is 
                  done because some programs like pepstats do not like 
                  long identifier names as input.
              
        Input:    Path to desired FASTA format output file, list of 
                  identifiers and list of corresponding sequences.
    
        Return:   List of short identifiers.
    """

    with open(f_output, 'w') as f:
        SHORT_IDENTIFIERS = []
        # Change identifiers to protein1, protein2, ...
        # and write to temporary file
        SET = zip(ORIGINAL_IDENTIFIERS, ORIGINAL_SEQUENCES)
    
        for index,  (identifier, sequence) in enumerate(SET):
            short_id = '>protein' + str(index)
            SHORT_IDENTIFIERS.append(short_id)
            f.writelines(short_id + '\n')
            f.writelines(sequence + '\n')

    return SHORT_IDENTIFIERS
# -----------------------------------------------------------------------------------------------------------
def pepstats(SHORT_IDENTIFIERS, SEQUENCES, pepstats_file):
    """ Function: pepstats()

        Purpose:  Given a set that contains the list of identifiers and 
                  the corresponding list of sequences, scan the given 
                  pepstats result file to extract protein properties.
              
        Input:    Set that contains the list of identifiers and 
                  the corresponding list of sequences and peptstats 
                  result file. 
    
        Return:   Dictionary of protein properties for each protein in the 
                  set.
    """

    pepstats_dic = {}

    with open(pepstats_file) as f: 
        content = f.readlines()
        for start, line in enumerate(content):
            if 'PEPSTATS of ' in line:
                TARGET_ID = line.split('PEPSTATS of ')[1]
                TARGET_ID = TARGET_ID.split('from 1 to')[0]
                TARGET_ID = TARGET_ID.strip()
                sequence = None

                for (identifier, seq) in zip(SHORT_IDENTIFIERS, SEQUENCES):
                    if identifier.replace('>', '') == TARGET_ID:
                        sequence = seq.strip()

                if sequence:
                    length = float(len(sequence))
                    # Amino acid frequencies in the sequence
                    amino_acid_frequencies = []
                    amino_acid_frequencies.append(100.0*sequence.count('A')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('C')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('D')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('E')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('F')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('G')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('H')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('I')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('K')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('L')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('M')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('N')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('P')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('Q')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('R')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('S')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('T')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('V')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('W')/length)
                    amino_acid_frequencies.append(100.0*sequence.count('Y')/length)
                    # Extract molecular weight
                    mwline = content[start + 2:start + 3]
                    molecular_weight = float(re.findall("\d+.\d+", str(mwline))[0])
                    # Extract charge
                    charge_line = content[start + 3:start + 4]
                    charge = float(re.findall("[-+]?\d+.\d+", str(charge_line))[1])
                    # Extract amino acid class frequencies

		    # Watch out, in the pepstats software, if isoelectric point == None, an 
                    # extra line will be introduced
                    start_aas = content[start:].index('Property\tResidues\t\tNumber\t\tMole%\n')
                    perline = content[start + start_aas + 1:start + start_aas + 10]

                    tiny = float(re.findall("\d+.\d+", str(perline[0]))[-1])
                    small = float(re.findall("\d+.\d+", str(perline[1]))[-1])
                    aliphatic = float(re.findall("\d+.\d+", str(perline[2]))[-1])
                    aromatic = float(re.findall("\d+.\d+", str(perline[3]))[-1])
                    non_polar = float(re.findall("\d+.\d+", str(perline[4]))[-1])
                    polar = float(re.findall("\d+.\d+", str(perline[5]))[-1])
                    charged = float(re.findall("\d+.\d+", str(perline[6]))[-1])
                    basic = float(re.findall("\d+.\d+", str(perline[7]))[-1])
                    acidic = float(re.findall("\d+.\d+", str(perline[8]))[-1])
                    amino_acid_classes = []
                    amino_acid_classes.append(tiny)
                    amino_acid_classes.append(small)
                    amino_acid_classes.append(aliphatic)
                    amino_acid_classes.append(aromatic)
                    amino_acid_classes.append(non_polar)
                    amino_acid_classes.append(polar)
                    amino_acid_classes.append(charged)
                    amino_acid_classes.append(basic)
                    amino_acid_classes.append(acidic)
                    # Store values in dictionary
                    pepstats_dic[TARGET_ID] = molecular_weight, charge, amino_acid_classes, amino_acid_frequencies, length, GRAVY(sequence), EXPOSED(sequence), DISORDER(sequence), HYDRO(sequence), BULKY(sequence), INTERFACE(sequence)

                else:
                    print 'There was an error scanning the pepstats file.'
                    print 'Could not find corresponding sequence for identifier', TARGET_ID
                    sys.exit()

    return pepstats_dic
# -----------------------------------------------------------------------------------------------------------
def write_weka_input(weka_input, SHORT_IDENTIFIERS, pepstats_dic):
    """ Function: write_weka_input()

        Purpose:  Given the query identifiers and pepstats-calculated
                  protein features, write the input arff file for WEKA. 
              
        Input:    WEKA arff file name, query identifiers and pepstats dictionary.                  
    
        Return:   None. 
    """   
    with open(weka_input, 'w') as f:
        # Create a list of features for each protein
        X = [[] for __ in xrange(len(SHORT_IDENTIFIERS))]

        for protein_position, TARGET_ID in enumerate(SHORT_IDENTIFIERS):
            TARGET_ID = TARGET_ID.replace('>', '')
            TARGET_ID = TARGET_ID.strip()

            molecular_weight, charge, amino_acid_classes, amino_acid_frequencies, length, GRAVY, Exposed, Disorder, Hydro, Bulky, Interface = pepstats_dic[TARGET_ID]

            X[protein_position] = amino_acid_classes + amino_acid_frequencies + [molecular_weight, charge, GRAVY, Exposed, Disorder, Hydro, Bulky, Interface]

        # Write protein feature data to WEKA arff file
        f.writelines(ARFF_HEADER)
        for index, vector in enumerate(X):
            for feature in vector:
                f.writelines(str(feature) + ',')
            f.writelines('?\n')

    return
# -----------------------------------------------------------------------------------------------------------
def parse_weka_output(file_input, ORIGINAL_IDENTIFIERS, SEQUENCES):
    """ Function: parse_weka_output()

        Purpose:  Given the WEKA output file and the query identifiers and sequences, 
                  parse the predicted class for each protein from the WEKA output. 
                  Write the predicted effectors to a FASTA file.
              
        Input:    WEKA output file and the query identifiers and sequences.                  
    
        Return:   The set of predicted effectors only as well as all predictions. 
    """    
    predicted_effectors, predicted_noneffectors, predictions = [], [], []

    with open(file_input) as f:

        content = f.readlines()
        #print content

        content_start = content.index('    inst#     actual  predicted error prediction (Tiny,Small,Aliphatic,Aromatic,Nonpolar,Polar,Charged,Basic,Acidic,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,MolecularWeight,Charge,GRAVY,Exposed,Disorder,Hydro,Bulky,Interface)\n') 
        
        content = content[content_start + 1:]

        for line in content:
            if line.strip():
                position = line.split()[0]
                prediction = line.split()[2]
                prob = float(line.split()[3])        
 
                # WEKA output counts from position 1, our identifiers are counted from zero
                identifier = ORIGINAL_IDENTIFIERS[int(position) - 1]
                sequence = SEQUENCES[int(position) - 1]

                if 'non-eff' in prediction:                               
                    noneffector = identifier.strip()
                    noneffector = noneffector.replace('>', '')  
                    predictions.append((noneffector, 'Non-effector', prob, sequence))
                    predicted_noneffectors.append((noneffector, prob, sequence))
                else:                    
                    effector = identifier.strip()
                    effector = effector.replace('>', '')                                               
                    predictions.append((effector, 'Effector', prob, sequence))
	                # Append predicted effector to list of predicted effectors
                    predicted_effectors.append((effector, prob, sequence))

    return predicted_effectors, predicted_noneffectors, predictions
# -----------------------------------------------------------------------------------------------------------
def short_output(predictions):
    """ Function: short_output()

        Purpose:  Given the WEKA predictions for each protein, write  
                  string that contains the short output format.
              
        Input:    WEKA predictions for each protein.                  
    
        Return:   String that contains predictions for all proteins as tab-delimited table.
    """
    # Output predictions for all proteins as tab-delimited table
    short_output_string = '# Identifier \t Prediction \t Probability \n'
    for protein, pred, prob, sequence in predictions:    
        short_output_string += protein + '\t' + pred + '\t' + str(prob) + '\n'            

    return short_output_string
# -----------------------------------------------------------------------------------------------------------
def long_output(ORIGINAL_IDENTIFIERS, predicted_effectors):
    """ Function: long_output()

        Purpose:  Given the predicted effectors and identifiers for the test set,  
                  write string that contains the long output format.
              
        Input:    Predicted effectors and identifiers of test set.                  
    
        Return:   String that contains list of predicted effectors with posterior probabilites
                  and a short statistic on the percentage of predicted effectors in the test set.
    """
    # Output predicted effectors for long format
    long_output_string = '-----------------\n'
    long_output_string += 'Predicted effectors:\n\n'
    for effector, prob, sequence in predicted_effectors:
        long_output_string += effector + '| Effector probability:' + str(prob) + '\n'

    long_output_string += '-----------------\n\n'
    long_output_string += 'Number of proteins that were tested: ' + str(len(ORIGINAL_IDENTIFIERS)) + '\n' 
    long_output_string += 'Number of predicted effectors: ' + str(len(predicted_effectors)) + '\n' 
    long_output_string += '\n' + '-----------------' + '\n' 
    long_output_string += str(round(100.0*len(predicted_effectors)/len(ORIGINAL_IDENTIFIERS), 1)) + ' percent are predicted to be effectors.'  
    long_output_string += '\n' + '-----------------' + '\n'

    return long_output_string
# -----------------------------------------------------------------------------------------------------------           
