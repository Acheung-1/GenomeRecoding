# GenomeRecoding

Users input a segments of a virus genome. This application optimizes the locations of CpGs using dynamic programming and works as follows:

Identifies the current locations of all CpGs (Cs followed by Gs)
Identifies all locations that can become a CpG (by mutating the nucleotide sequence such that it encodes for the same exact amino acid sequence)
Outputs a DNA sequence that will encode the same amino acid as the user input

This outputted sequence will optimize the CpG spacing to be 13 nucleotides away from each other. There is evidence in the scientific literature that the virus will be weakened.
