#! /usr/bin/python
"""
    EffectorP 2.0: predicting fungal effector proteins from secretomes using machine learning
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
import functions
import subprocess
import errno
import uuid
import shutil
import tempfile
import numpy as np
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# Main Program starts here
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
if __name__ == '__main__': 
    # -----------------------------------------------------------------------------------------------------------
    SCRIPT_PATH = sys.path[0]
    # Change the path to WEKA to the appropriate location on your computer
    WEKA_PATH = SCRIPT_PATH + '/weka-3-8-1/weka.jar'
    PEPSTATS_PATH = SCRIPT_PATH + '/EMBOSS-6.5.7/emboss/'
    # -----------------------------------------------------------------------------------------------------------
    # Check that the path to the WEKA software exists
    path_exists = os.access(WEKA_PATH, os.F_OK)
    if not path_exists:
        print
        print "Path to WEKA software does not exist!"
        print "Check the installation and the given path to the WEKA software %s in EffectorP.py (line 47)."%WEKA_PATH
        print
        sys.exit()
    # -----------------------------------------------------------------------------------------------------------
    # Check that the path to the EMBOSS software exists for pepstats
    path_exists = os.access(PEPSTATS_PATH, os.F_OK)
    if not path_exists:
        print
        print "Path to EMBOSS software does not exist!"
        print "Check the installation and the given path to the EMBOSS software %s in EffectorP.py (line 48)."%PEPSTATS_PATH
        print
        sys.exit()
    # -----------------------------------------------------------------------------------------------------------
    commandline = sys.argv[1:]
    # -----------------------------------------------------------------------------------------------------------
    if commandline:
        FASTA_FILE, short_format, output_file, effector_output, noneffector_output = functions.scan_arguments(commandline)
	# If no FASTA file was provided with the -i option
        if not FASTA_FILE:
            print
            print 'Please specify a FASTA input file using the -i option!'
            functions.usage()
    else:
        functions.usage()
    # -----------------------------------------------------------------------------------------------------------
    # Temporary folder name identifier that will be used to store results
    RESULTS_PATH = tempfile.mkdtemp() + '/'
    # -----------------------------------------------------------------------------------------------------------
    # Check if FASTA file exists
    try:
        open(FASTA_FILE, 'r') 
    except IOError as (errno, strerror):
        print "Unable to open FASTA file:", FASTA_FILE  #Does not exist OR no read permissions
        print "I/O error({0}): {1}".format(errno, strerror)
        sys.exit()
    # -----------------------------------------------------------------------------------------------------------
    # Try to create folder where results will be stored
    try:
        os.mkdir(RESULTS_PATH)
    except OSError as exception:        
        if exception.errno != errno.EEXIST:
            raise
    # -----------------------------------------------------------------------------------------------------------
    # Extract the identifiers and sequences from input FASTA file
    ORIGINAL_IDENTIFIERS, SEQUENCES = functions.get_seqs_ids_fasta(FASTA_FILE)
    SEQUENCES = [seq.upper() for seq in SEQUENCES]
    # -----------------------------------------------------------------------------------------------------------
    print '-----------------'
    print 
    print "EffectorP 2.0 is running for", len(ORIGINAL_IDENTIFIERS), "proteins given in FASTA file", FASTA_FILE
    print
    # -----------------------------------------------------------------------------------------------------------
    # Write new FASTA file with short identifiers because pepstats can't handle long names
    f_output = RESULTS_PATH + 'short_ids.fasta'
    SHORT_IDENTIFIERS = functions.write_FASTA_short_ids(f_output, ORIGINAL_IDENTIFIERS, SEQUENCES)
    # -----------------------------------------------------------------------------------------------------------
    # Call pepstats
    print 'Call pepstats...'
    ProcessExe = PEPSTATS_PATH + 'pepstats'
    ParamList = [ProcessExe, '-sequence', RESULTS_PATH + 'short_ids.fasta', 
              '-outfile', RESULTS_PATH +  'pepstats.out']

    try:
        Process = subprocess.Popen(ParamList, shell=False)
        sts = Process.wait()
        cstdout, cstderr = Process.communicate()

        if Process.returncode:
            raise Exception("Calling pepstats returned %s"%Process.returncode)
        if cstdout:
            pass
        elif cstderr:
            sys.exit()
    except:
        e = sys.exc_info()[1]
        print "Error calling pepstats: %s" % e
        sys.exit()
    print 'Done.'
    print
    # -----------------------------------------------------------------------------------------------------------
    # Parse pepstats file
    print 'Scan pepstats file'
    pepstats_dic = functions.pepstats(SHORT_IDENTIFIERS, SEQUENCES, RESULTS_PATH +  'pepstats.out')
    print 'Done.'
    print
    # -----------------------------------------------------------------------------------------------------------
    # Write the WEKA arff file for classification of the input FASTA file
    weka_input = RESULTS_PATH + 'weka.arff'    
    functions.write_weka_input(weka_input, SHORT_IDENTIFIERS, pepstats_dic)
    # -----------------------------------------------------------------------------------------------------------
    # Ensembl averaging approach, use seq0,seq1,seq2... as keys in case there are duplicate FASTA identifiers
    ensembl_votes = {}
    # -----------------------------------------------------------------------------------------------------------
    # Call WEKA models for classification of input FASTA file
    # -----------------------------------------------------------------------------------------------------------
    models = functions.models_bayes + functions.models_J48
    # -----------------------------------------------------------------------------------------------------------
    print 'Ensemble classification'
    for model in functions.models_bayes:
        #--------------------------------------------------------------
        ParamList = ['java', '-cp', WEKA_PATH, 'weka.classifiers.bayes.NaiveBayes', '-l', SCRIPT_PATH + model, '-T', RESULTS_PATH + 'weka.arff', '-p', 'first-last']

        with open(RESULTS_PATH + 'Predictions.txt', 'wb') as out:
            try:
                Process = subprocess.Popen(ParamList, shell=False, stdout=out)
                sts = Process.wait()
                cstdout, cstderr = Process.communicate()

                if Process.returncode:
                    raise Exception("Calling WEKA returned %s"%Process.returncode)
                if cstdout:
                    pass
                elif cstderr:
                    sys.exit()
            except:
                e = sys.exc_info()[1]
                print "Error calling WEKA: %s" % e
                sys.exit()
        #-------------------------------------------------------------- 
        # Parse the WEKA output file
        file_input = RESULTS_PATH + 'Predictions.txt'
        predicted_effectors, predicted_noneffectors, predictions = functions.parse_weka_output(file_input, ORIGINAL_IDENTIFIERS, SEQUENCES)
        
        for index, (ident, prediction, prob, seq) in enumerate(predictions):

            short_ident = 'seq' + str(index)

            if short_ident in ensembl_votes:
                previous_predictions = ensembl_votes[short_ident] 
                ensembl_votes[short_ident] = previous_predictions + [(prediction, prob)]
            else:
                ensembl_votes[short_ident] = [(prediction, prob)]
        #-------------------------------------------------------------- 
    for model in functions.models_J48:
        #--------------------------------------------------------------
        ParamList = ['java', '-cp', WEKA_PATH, 'weka.classifiers.trees.J48', '-l', SCRIPT_PATH + model, '-T', RESULTS_PATH + 'weka.arff', '-p', 'first-last']

        with open(RESULTS_PATH + 'Predictions.txt', 'wb') as out:
            try:
                Process = subprocess.Popen(ParamList, shell=False, stdout=out)
                sts = Process.wait()
                cstdout, cstderr = Process.communicate()

                if Process.returncode:
                    raise Exception("Calling WEKA returned %s"%Process.returncode)
                if cstdout:
                    pass
                elif cstderr:
                    sys.exit()
            except:
                e = sys.exc_info()[1]
                print "Error calling WEKA: %s" % e
                sys.exit()
        #-------------------------------------------------------------- 
        # Parse the WEKA output file
        file_input = RESULTS_PATH + 'Predictions.txt'
        predicted_effectors, predicted_noneffectors, predictions = functions.parse_weka_output(file_input, ORIGINAL_IDENTIFIERS, SEQUENCES)
        
        for index, (ident, prediction, prob, seq) in enumerate(predictions):

            short_ident = 'seq' + str(index)

            if short_ident in ensembl_votes:
                previous_predictions = ensembl_votes[short_ident] 
                ensembl_votes[short_ident] = previous_predictions + [(prediction, prob)]
            else:
                ensembl_votes[short_ident] = [(prediction, prob)]
    print 'Done.'
    print
    #--------------------------------------------------------------
    # Soft voting
    #--------------------------------------------------------------
    ensemble_predictions, predicted_effectors, predicted_noneffectors, predicted_weakeffectors = [], [], [], []
    #--------------------------------------------------------------
    for index, (ident, prediction, prob, seq) in enumerate(predictions):
 
        short_ident = 'seq' + str(index)
        yes_prob, no_prob = [], []

        for vote, prob in ensembl_votes[short_ident]:

            if vote == 'Non-effector':
                no_prob.append(prob)
                yes_prob.append(1.0 - prob)                        

            if vote == 'Effector':
                yes_prob.append(prob)
                no_prob.append(1.0 - prob)        
             
        yes_prob = sum(yes_prob)/float(len(models))
        no_prob = sum(no_prob)/float(len(models))
        
        yes_prob, no_prob = round(yes_prob,3), round(no_prob,3)

        if yes_prob > no_prob:
            if yes_prob > 0.55:
                prediction = 'Effector'
                prob = round(yes_prob,3)
                predicted_effectors.append((ident, prob, seq))
            else:
                prediction = 'Unlikely effector'
                prob = round(yes_prob,3)
                predicted_weakeffectors.append((ident, prob, seq))
        else:
            prediction = 'Non-effector'
            prob = round(no_prob,3)
            predicted_noneffectors.append((ident, prob, seq))

        ensemble_predictions.append((ident, prediction, prob, seq))
    #--------------------------------------------------------------
    # If user wants the stdout output directed to a specified file
    if output_file:

        with open(output_file, 'wb') as out:
            # Short format: output predictions for all proteins as tab-delimited table
            if short_format:
                out.writelines(functions.short_output(ensemble_predictions))
            # If the user wants to see the long format, output additional information and stats
            else:
                out.writelines(functions.short_output(ensemble_predictions))
                out.writelines(functions.long_output(ORIGINAL_IDENTIFIERS, predicted_effectors))                
        print 'EffectorP results were saved to output file:', output_file  

    else:
        # Short format: output predictions for all proteins as tab-delimited table to stdout
        if short_format:
            print functions.short_output(ensemble_predictions)
        # If the user wants to see the long format, output additional information and stats
        else:
            print functions.short_output(ensemble_predictions)
            print functions.long_output(ORIGINAL_IDENTIFIERS, predicted_effectors)
    # -----------------------------------------------------------------------------------------------------------
    # If the user additionally wants to save the predicted effectors in a provided FASTA file
    if effector_output:
        with open(effector_output, 'w') as f_output:
            for effector, prob, sequence in predicted_effectors:
                f_output.writelines('>' + effector + ' | Effector probability: ' + str(prob) + '\n')
                f_output.writelines(sequence + '\n')  
    if noneffector_output:
        with open(noneffector_output, 'w') as f_output:
            for effector, prob, sequence in predicted_noneffectors:
                f_output.writelines('>' + effector + ' | Non-effector probability: ' + str(prob) + '\n')
                f_output.writelines(sequence + '\n')  
            for effector, prob, sequence in predicted_weakeffectors:
                f_output.writelines('>' + effector + ' | Unlikely effector probability: ' + str(prob) + '\n')
                f_output.writelines(sequence + '\n')  
    # -----------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------
    # Clean up and delete temporary folder that was created
    shutil.rmtree(RESULTS_PATH)
    # -----------------------------------------------------------------------------------------------------------
    try:
        sys.stdout.close()
    except:
        pass
    try:
        sys.stderr.close()
    except:
        pass
    # -----------------------------------------------------------------------------------------------------------

