## Pipeline used in:

# Phylogenetic Analysis of Allotetraploid Species Using Polarized Genomic Sequences
J.L. Leal, Pascal Milesi, Eva Hodková, Jarkko Salojärvi, and Martin Lascoux

(forthcoming)

____


If you made it to this page, you are likely looking for the script used to phase allotetraploid sequences using our polarization protocol.

The script is called "01s_polarizeTETRA.py" and is written in python3. You can find it in the folder "03_Polarize_Allotetraploid" above.

Run this script using the following command:

```bash
python3 01s_polarizeTETRA.py [input.fasta] [No. species] [refSequence ID] [polyploid ID]
```

Input arguments:
- <input.fasta>: MSA in fasta format. Sequences are assumed to be aligned and coded using IUPAC nomenclature. 
  If your MSA includes gaps (-), recode them as masked sites (N) before running this script.                   
- <No. species>: Number of species in MSA.
- [refSequence ID>]: Reference sequence ID (polarizing sequence), as shown in MSA.
- [polyploid ID]: Allopolyploid sequence ID, as shown in MSA.

If you use this script, please cite the paper listed above.
