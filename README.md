## Pipeline used in:

# Phylogenetic Analysis of Allotetraploid Species Using Polarized Genomic Sequences
J.L. Leal, Pascal Milesi, Eva Hodková, Jarkko Salojärvi, and Martin Lascoux

(forthcoming)

____


If you made it to this page, you are likely looking for the script used to phase allotetraploid sequences using the genomic polarization protocol.

The script is called "01s_polarizeTETRA.py" and is written in python3. You can find it in the folder "03_Polarize_Allotetraploid" above.

Run this script using the following command:

```bash
python3 01s_polarizeTETRA.py [input.fasta] [No. species] [refSequence ID] [polyploid ID]
```

Input arguments:
- [input.fasta]: MSA in fasta format. Sequences are assumed to be aligned and coded using IUPAC nomenclature. 
  If your MSA includes gaps (-), recode them as masked sites (N) before running this script.                   
- [No. species]: Number of species in MSA.
- [refSequence ID]: Reference sequence ID (polarizing sequence), as shown in MSA.
- [polyploid ID]: Allopolyploid sequence ID, as shown in MSA.

NOTE: When running this script over a very large number of loci, you may occasionally get the following error message "Error: check number of sequences in fasta file". This happens if some of your MSAs are missing some sequences. This could have happened for a variety of reasons. For instance, trimAl, a program used to clean MSAs, will remove any sequence that contains only masked sites (NNNNN). If only a few of your MSAs are missing sequences, you can ignore the error messages you get when running the polarization script.

If you use this script, please cite the paper listed above.
