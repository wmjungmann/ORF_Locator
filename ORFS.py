"""

Written by: Matthew Jungmann
07/03/2024

Program that takes a fasta file of an E-Coli pal gene sequence as an input. It then generates the complimentary strand and reverse complimentary 
strand of the primary sequence. After then converting these sequences to amino acid sequences, it will find and display the frame number, the 
nucleotide and amino acid start and stop positions, the length of the ORFs in nucleotide/amino positions, and display the DNA/AA sequence for 
all potential realistic ORFs.

When locating ORFs, the minimum length can be set inside of 

"""

# Matthew Jungmann
# This program will determine and display the position of all the Start amino Acids
# all the stop amino acids for the fist, second, and third reading frames of the primary strand.
# the start is 'M' and the stop is '_'

import sys

codon_table = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

def ConvertComp(dna_seq):
    complimentary_strand = ''
    
    # Convert nucleotides to their counterpart (T<->A | G<->C) and add to new complimentary strand
    for nucleotide in dna_seq:
        if nucleotide == 'T':
            complimentary_strand += 'A'
        if nucleotide == 'A':
            complimentary_strand += 'T'
        if nucleotide == 'G':
            complimentary_strand += 'C'
        if nucleotide == 'C':
            complimentary_strand += 'G'
    
    return complimentary_strand


def CreateAminos(read_frames):
    
    amino_list = []

    # Convert each frame to a list of codons
    for frame in read_frames:
        amino_seq = ''
        
        # Cycle through the frame in sets of 3 nucleotides simulating codons
        for n in range(0, len(frame), 3):
           codon = frame[n:n+3]
           
           # Check if the codon substring is in the codon table dictionary
           if codon in codon_table:
               amino_acid = codon_table[codon]
               amino_seq += amino_acid
        
        amino_list.append(amino_seq)
        
    
    return amino_list

def GenerateFrames(primary, complimentary):
    orf = []
    
    # Reverse the complimentary sequence so it can be converted from the start of the string instead of the end
    complimentary_rev = complimentary[::-1]
    
    # Get frames +1, +2, +3
    for start_nt in range(0,3):
       frame = primary[start_nt:len(primary)]
       orf.append(frame)
       
    print("\n")
       
    # Get frames -1, -2, -3
    for start_nt in range(0,3):
       frame = complimentary_rev[start_nt:len(complimentary_rev)]
       orf.append(frame)
    
    return orf
        

def FindOrfs(rf, rf_num, frame_num, min_length):
    
    true_orfs = []
    start = 0
    stop = 0
    potential = ""
    nuc_start_pos = []
    nuc_stop_pos = []
    aa_start_pos = []
    aa_stop_pos = []
    
    
    # Check if the reading frame is from the primary or complimentary strand by checking reading fram enumber
    if rf_num > 0:
        
        # Set the frame offset and length of the reading frame in nucleotides
        frame_offset = rf_num - 1
        
        # Cycle through amino acids (AA)
        while start < len(rf):
            j = 0
            potential = ''
            
            # Check if current AA is equal to start codon 
            if rf[start] == "M":
                
                # Add start codon to potential ORF
                potential += rf[start]
                stop = start + 1

                # Starting from the codon after M, check each AA
                while stop < len(rf):
                    
                    # Check if current AA is equal to a stop AA
                    if rf[stop] == '*':
                        
                        # Check that potential ORF is longer than the minimum requirement
                        if len(potential) >= min_length:
                            # Find amino acid start and stop position
                            aa_start_position = str(start + 1)
                            aa_stop_position = str(stop + 1)
                            
                            # Find nucleotide start and stop position using the offset of the current reading frame
                            start_position = str(((((start + 1) * 3) - 2) + frame_offset))
                            stop_position = str((((stop + 1) * 3) + frame_offset))
                            
                            # Add ORF and start and stop positions to lists
                            true_orfs.append(potential)
                            nuc_start_pos.append(start_position)
                            nuc_stop_pos.append(stop_position)
                            aa_start_pos.append(aa_start_position)
                            aa_stop_pos.append(aa_stop_position)
                            
                            # Set start position to the last stop codon position for next search
                            start = stop
                            
                            break
                        
                        # If potential ORF is not of sufficient length, stop searching for stop codon and restart search for start codon
                        else:
                            break
                        
                    # If current codon does not equal stop codon
                    else:
                        # Check if the current AA is the last AA in the reading frame
                        if stop == (len(rf) - 1):
                            
                            # Check that the potential ORF is longer than the minimum requirement
                            if len(potential) >= min_length:
                                
                                aa_start_position = str(start + 1)
                                aa_stop_position = ">" + str(stop + 1)
                                
                                # Find nucleotide start and stop position using the offset of the current reading frame, 
                                # add symbol to show stop position is out of view
                                start_position = str(((((start + 1) * 3) - 2) + frame_offset))
                                stop_position = ">" + str((((stop + 1) * 3) + frame_offset))
                                
                                # Add ORF and start and stop positions to lists
                                true_orfs.append(potential)
                                nuc_start_pos.append(start_position)
                                nuc_stop_pos.append(stop_position)
                                aa_start_pos.append(aa_start_position)
                                aa_stop_pos.append(aa_stop_position)
                                
                                # Set start position to the last stop codon position for next search
                                start = stop
                                break
                            
                            # If potential ORF is not of sufficient length, stop searching for stop codon and restart search for start codon
                            else:
                                break
                            
                        # If current codon is not the final codon in the sequence, move to the next codon
                        else:
                            potential += rf[stop]
                            stop += 1
            
            start += 1
    else:
        # Set start position to the end of the sequence
        start = len(rf) - 1 
        frame_offset = rf_num + 1
        
        # Cycle through amino acids (AA)
        while start >= 0:
            stop = 0
            potential = ''
            
            # Check if current AA is equal to start codon 
            if rf[start] == "M":
                
                # Set stop to the next AA in the sequence
                potential += rf[start]
                stop = start - 1

                # Cycle through AA until the end of the sequence
                while stop >= 0:
                    
                    # Check if current AA is equal to stop codon 
                    if rf[stop] == '*':
                        
                        # Check if the potential ORF reaches the minimum length requirement of 20
                        if len(potential) >= min_length:
                            
                            # Register amino acid start and stop positions
                            aa_start_position = str(start + 1)
                            aa_stop_position = str(stop + 1)
                            
                            # Find nucleotide start and stop position using the offset of the current reading frame, 
                            # add symbol to show stop position is out of view
                            start_position = str(((((start + 1) * 3) + 2) + frame_offset))
                            stop_position = str((((stop + 1) * 3) + frame_offset))
                            
                            # Reverse the potential string to display the ORF the correct direction from right to left
                            potential = potential[::-1]
                            
                            # Add ORF and start and stop positions to lists
                            true_orfs.append(potential)
                            nuc_start_pos.append(start_position)
                            nuc_stop_pos.append(stop_position)
                            aa_start_pos.append(aa_start_position)
                            aa_stop_pos.append(aa_stop_position)
                            
                            start = stop
                            break
                        
                        # If potential ORF is not of sufficient length, stop searching for stop codon and restart search for start codon
                        else:
                            break
                    else:
                        # Check if stop position has reached the end of the sequence
                        if stop == 0:
                            if len(potential) >= min_length:
                                # Register amino acid start and stop positions
                                aa_start_position = str(start + 1)
                                aa_stop_position = ">" + str(stop + 1)

                                # Find nucleotide start and stop position using the offset of the current reading frame, 
                                # add symbol to show stop position is out of view
                                frame_bounds = rf_num + 4
                                start_position = str(((((start + 1) * 3) + 2) + frame_offset))
                                stop_position = ">" + str(frame_bounds)
                                
                                # Add ORF and start and stop positions to lists
                                true_orfs.append(potential)
                                nuc_start_pos.append(start_position)
                                nuc_stop_pos.append(stop_position)
                                aa_start_pos.append(aa_start_position)
                                aa_stop_pos.append(aa_stop_position)
                                
                                start = stop
                                break
                            
                            # If potential ORF is not of sufficient length, stop searching for stop codon and restart search for start codon
                            else:
                                break
                        
                        # If current codon is not the final codon in the sequence, move to the next codon
                        else:
                            potential += rf[stop]
                            stop -= 1
                
            start -= 1
        
    
    count = 0
    if len(true_orfs) > 0:  
    
        
        
        for orf in true_orfs:
            # Find length of ORF in amino acid sequence 
            AA_length = len(orf)
            
            # Convert to length in nucleotide sequence
            nucleotide_len = (AA_length * 3) + 3    
            

            # Display the ORF information
            print("\n{:<5} | {:<8} | {:<8} | {:<8} | {:<8} | {:<8} | {:<15}".format('ORF', 'Frame', 'AA Start', 'AA Stop', 'NT Start', 'NT Stop', 'Length (nt|aa)'))
            print("-----------------------------------------------------------------------------")
            print("{:<5} | {:<8} | {:<8} | {:<8} | {:<8} | {:<8} | {:<3}|{:<3}".format(str(frame_num), str(rf_num), str(aa_start_pos[count]), str(aa_stop_pos[count]) , str(nuc_start_pos[count]), str(nuc_stop_pos[count]), str(nucleotide_len), str(AA_length)))
            print("-----------------------------------------------------------------------------")
            
            # Display the ORF
            print("The AA sequence of the ORF: ")
            print(orf + "\n")
            
            
            count += 1
            frame_num += 1
    
    # Return frame num to keep track of ORFS across differnt reading frames
    return frame_num


                


#******************************* main function ********************************************************

def main():

    # Prompt the user for the file and open the file for reading
    FileName = input("Enter the name of the file: ")
    
    # Prompt user to enter the minimum length of ORF 
    min_orf_length = input("Enter the minimum length of ORF (Amino Acid): ")
    
    # Check that input length is valid, prompt again for input if not
    while min_orf_length.isnumeric() == False:
        print("\nInvalid input, please enter a number!")
        min_orf_length = input("Enter the minimum length of ORF (Amino Acid): ")
    
    # Change input to integer
    min_orf_length = int(min_orf_length)
    
    # Open the file in reading text mode, error check to see if file exists
    try:
        Fp1 = open(FileName,'r')

        # Read in Descriptor Line
        descriptor_line = Fp1.readlines(1)

        # Read and print the DNA contents: assume it starts on 2nd line
        dna_lines = Fp1.read()
     
    # Code to exit "gracefully" if file does not exist
    except IOError:
        print("error unable to read file or file does not exist!!!")
        print("Exiting the program")
        stop = input("press enter to exit....")
        Fp1.close()
        sys.exit(1)


    # Use split (\n) /join ('') to form a contiguous sequence
    list_seq = dna_lines.split('\n')
    
    primary_strand = ('').join(list_seq)
    
    # Generate complimentary strand
    complimentary_strand = ConvertComp(primary_strand)
    
    # Frames are returned in order +1, +2, +3, -1, -2, -3
    dna_frames = GenerateFrames(primary_strand, complimentary_strand)
    
    # Convert nucleotide reading frames to amino acids
    amino_frames = CreateAminos(dna_frames)
    
    # Divide amino acid frames into 5'-3' and 3'-5' frames
    amino_frames_prim = amino_frames[0:3]
    amino_frames_rev = amino_frames[3:7]
    
    
    # Initialise numbers to track frames when finding orfs
    rf_num = 0
    nuc_num = 0
    frame_num = 1
    
    print("~ORFs in reading frames 1, 2, 3 are read left to right~")
    print("~ORFs in reading frames -1, -2, -3 are read right to left~\n")
    

    # Cycle thorugh frames 1 to 3, locating ORFs and relevant information
    for frame in amino_frames_prim:
        rf_num += 1
        frame_num = FindOrfs(frame, rf_num, frame_num, min_orf_length)
        nuc_num += 1
    
    
    # Cycle through frames -1 to -3, locating ORFs and relevant information
    rf_num = 0
    nuc_num = 0
    for frame in amino_frames_rev:
        rf_num -= 1
        frame = frame[::-1]
        frame_num = FindOrfs(frame, rf_num, frame_num, min_orf_length)
    
    
    stop = input("\n press return to finish....")


"""****************** test plan ********************************
  Run the program and compare output to output found on the NCBI ORF finder. 
    - input incorrect file name to ensure program correctly identifies this and ends gracefully
    - input invalid input for minimum ORF length to ensure program deals with this accordingly
    - Comapare lengths, start and stop positions, and ORF amino acids with those found on the NCBI ORF finder
   

"""

#**************** execute program **************************

main()
