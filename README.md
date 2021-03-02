# suitename
Suitename - a tool for classifying backbone "suite" conformations (linkages between ribose sugars) in an RNA molecule.

Suitename takes an input file as its first argument, or if none provided, reads the standard input. It writes its results to standard output. It classifies the conformation of each suite into one of several dozen predefined clusters that have been determined by years of study. 

Two forms of input are supported:
  1. A list of the 6 dihedral angles $/alpha$, $/beta$, $/gamma$, $/delta$, $/epsilon$, $/zeta$) for each residue of the RNA molecule. Suitename will re-parse the dihedral angles in the residues to obtain the 7 dihedral angles in each suite ($/delta/-1, $/epsilon/-1, $/zeta/-1, $/alpha$, $/beta$, $/gamma$, $/delta$), and then operate on those. This is the default input.

      Each line of this format describes one residue, with fields separated by colons. The first several (default 6) are ID information, the remainder are angles. Sample:
1: A:   7: : :  U:-75.533:-154.742:48.162:80.895:-148.423:-159.688

  2. A kinemage file providing a list of 7 or 9 dihedral angles in each suite. Mark this by using the --suitein command line flag. 9 angles include $/chi$-1 and $/chi$.

Three forms of output are supported:
  1. A text report, showing the classification of each suite into a cluster or outlier, and the neatness of fit ("suiteness") of the suite into that cluster. Suiteness is the cosine of the normalized distance of a suite datapoint from the power 3 hyperellipsoid boundary of the cluster in toward its center. A statistical summary is included at the end. This is the default output format.
  2. A kinemage file, which will display a 7D data point for each suite in the data. Points are color-coded according to the clusters to which they have been assigned. Each cluster is displayed as a colored ring surrounding its defined center. Specify this by using the --kinemage command line flag.
  3. A brief string showing only the cluster assignments. It consists of three characters per suite - base identity (uc) and 2-character number-letter name of the suite cluster (e.g., C1aG1gU1aA1aA1cG). Specify this by using the --string command line flag.

Many other command line options are available; type suitename --help to display them.
