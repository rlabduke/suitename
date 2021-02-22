# suitename
Suitename - a tool for classifying "suites" (linkages between bases) in an RNA molecule.

Suitename takes an input file as its first argument, or if none provided, reads the standard input. It writes its results to standard output. It classifies the conformation of each suite into one of several dozen predefined clusters that have been determined by years of study. 

Two forms of input are supported:
1. A list of residues in "dangle" format. Suitename will analyze the angles in the residues to form suites, and then operate on those. This is the default.
2. A kinemage file providing a list of suites. Mark this by using the --suitein the command line flag.

Three forms of output are supported:
1. A report, showing the classification of each suite into a cluster, and the neatness of fit ("suiteness") of the suite into that cluster. By default, a statistical summary is included at the end. This is the default.
2. A kinemage file, which will display a data point for each suite in the data. Points are color-coded according to the clusters to which they have been assigned. Each cluster is also displayed as a colored ring surrounding its defined center. Specify this by using the --kinemage command line flag.
3. A brief string showing only the cluster assignments. It consists of three characters per suite. Specify this by using the --string command line flag.

Many other command line options are available; type suitename --help to display them.