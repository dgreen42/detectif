A tool that takes a genetic sequence and identifies motifs

squences will be no less than 3 nucleotides. The cap is somewhat arbitrary, however there should be a default to prevent the program from going forever

 - default cap will be 20 nucleotides (chars)

The patterns will be derived from the sequence. This means that there will be initial pattern to start from

Algorithm:

1. Start (size = 3) window at the beginning of the sequence
2. Save new patterns in a hashmap with a count starting at 1. If a pattern already exists add to the count.
3. Once the whole sequence has been searched save results to a file: 
 - If the size <= 20 (or specifies size)  go to step 1 and add 1 to the window
 - if the size > 20 (or specified size) finish
