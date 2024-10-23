# Methodology

The method to find all “potential” realistic open reading frames (ORF) in the RNA sequence of SARS-CoV-2 can be split into four distinct steps: 
 - Creating the complimentary strand
 - Define the reading frames
 - Convert each reading frame from nucleotides to amino acids 
 - Identify realistic ORFs withing each reading frame  

To create the complimentary strand, we take the primary strand of RNA which is read 5’-3’ and convert each nucleotide to its known counterpart – if a particular nucleotide is guanine in the primary strand, it will be cytosine in the complimentary strand and vice versa. Likewise, if a particular nucleotide is adenine in the primary strand, it will be thymine in the complimentary strand, and vice versa. This complimentary strand is 3’-5’, meaning it is read right to left, as opposed to left to right as is done with the primary strand. After this is completed, we can check that this conversion has taken place correctly by verifying that the first and last three nucleotides from the primary strand are compliments of the first and last three nucleotides from the complimentary strand.

![image](https://github.com/user-attachments/assets/35885726-d37b-4084-885e-3418d3d63a5a)

Now the complimentary strand has been created, we now need to define the reading frames for both the primary and complimentary strands. As we can not be sure where the relevant genetic sequence begins, we must account for every possible outcome. We know that each codon is three nucleotides in length, therefore, three reading frames for each strand must be defined to account for every combination of nucleotides. For the primary strand, which is read left to right, we can do this by start reading frame 1 from nucleotide 1, reading frame 2 from nucleotide 2, and reading frame 3 from nucleotide 3. 

![image](https://github.com/user-attachments/assets/5832da6c-b3bd-4a7e-a461-bdb1b8abc34b)

For the complimentary strand, as we are reading from right to left, we begin from the end of the sequence, hence, reading frame -1 starts from the final nucleotide in the sequence, reading frame -2 from the penultimate nucleotide, and reading frame -3 from the third last nucleotide in the sequence.

![image](https://github.com/user-attachments/assets/b11460e9-8417-4b0e-beb5-6aa2cd324005)

With all the reading frames now defined, we must convert them to amino acids to facilitate the locating of start and stop codons of ORFs. To do so, we must check each group of three nucleotides from the beginning of each reading frame and translate them to their amino acid equivalent.  We do so by comparing each sequence to a known dictionary of amino acid nucleotide counterparts. This amino acid is then added to the list representing the amino acid sequence, and the next three nucleotides are then converted. Once all nucleotides in the frame are accounted for, the amino acid sequence is complete. This is method is then carried out for all the primary strand reading frames and complimentary strand reading frames.  Finally, we must now locate the potential ORFs from these generated amino acid sequences. We know that each ORF begins with a start codon and ends with a stop codon, therefore, we must try to locate these positions in the amino acid sequences. Additionally, we know that normally ORFs are only considered realistic if they are a certain number of amino acids in length. For this method, I will be using 20 amino acids as the minimum number for an ORF to be considered realistic. 

To locate ORFs for the primary strand reading frames, we analyse each amino acid in the sequence from start to finish. Start codons are represented in amino acid sequences by the letter “M”, while stop codons for this method will be represented by the symbol “\*”. To find the first start codon, we analyse each amino acid in the sequence until we encounter an “M” codon. We note the position of this codon, and then starting from the discovered start codon, we again sequentially look at each codon looking for a “\*” codon. If a “\*” codon is found, we must then check if the length of the potential ORF is greater than the length we determined before our search. If it is, we note down the position of the stop codon. The amino acids from the start codon to just before the stop codon are considered a “potential” realistic ORF. This method must now be repeated for all primary strand reading frames.

![image](https://github.com/user-attachments/assets/0d7cbe9e-26ce-46e7-9bfe-9e52670be7a4)

To locate to position of the start and stop codon for the complimentary strand reading frames, we must work backwards as the frames are 3’-5’ meaning they are read right to left. The same method applies, sequentially checking if an amino acid equals a start or stop codon. As well as checking whether the length of the ORF meets the minimum requirement to be considered a “potential” realistic ORF.

![image](https://github.com/user-attachments/assets/7b498511-99dd-4846-a375-6ff30ce5c783)

# Evaluation

To evaluate my program using SARS-CoV-2 genome, I set the minimum length of a potential ORF to 200 to limit the number of ORFs identified. For the evaluation to be thorough, there must be data that I will compare my results to. This data was collected from the NCBI ORF finder, the limit was set to 600nt.

![image](https://github.com/user-attachments/assets/252198a4-98b3-41ae-8aaf-4428642e6f1d)

Using this genome, my program successfully identified all “potential” realistic ORFs. This is a very good result considering the significant difference between SARS-CoV-2 genome and the e-coli genomes it was originally designed for. While this result was sufficient, the program was not fully tested as there were no complimentary strand ORFs detected. I decided to change the minimum length of the detected frames to 30nt to really test the program to its fullest extent.  Even with this significant change in the number of ORFs, the program success rate remained at 100%. This can be seen by comparing an ORF detected by my program to the same ORF from the NCBI ORF finder.

## NCBI ORF Finder
![image](https://github.com/user-attachments/assets/4cfec871-38f8-4d86-9305-6d4f0e641253)

![image](https://github.com/user-attachments/assets/b7a60a03-b9c3-4aa5-933d-0cb4cb68f383)

## ORF_Locator
![image](https://github.com/user-attachments/assets/04e76a76-7f5d-44a7-a5d3-0042762c03b0)


There is a significant limitation to my program that must be mentioned. For previously mentioned comparisons with the NCBI finder, I was only checking for “ATG” codons as start codons. If we look for alternative start codons, the accuracy of my program drops significantly – my program continues to detect 399 ORFs whereas the finder identifies 758 ORFs. This is something that would have to be addressed for the program to be up to the standard of other online software.

Overall, my created program can successfully detect “potential” realistic ORFs from a SARS-CoV-2 genome, across a range of different minimum lengths. It does have a glaring limitation, as it only detects ORFs that begin with an “ATG” codon. ORFs that begin with alternative start codons are not detected.





