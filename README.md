# deuterocol

A suite of membrane protein structure comparison tools

## Core

### Deuterocol1 ###

This script collects available indices and coordinates files for any number of TCDB superfamilies/families/subfamilies/systems.
It can also attempt to generate a negative control of comparable size to an existing set of positive control structures.

`--fams`
	A list of TC-IDs corresponding to the superfamilies, families, subfamilies, or systems of interest

`-o OUTDIR`
	Where to store the resulting data (default: deuterocol1)

`--invert-inclusive`
	Creates a negative control from excluded subfamilies belonging to known superfamilies

`--invert-inclusive-size FACTOR`
	Multiplies the effective size of the positive control by this value to generate a negative control of adjusted size

`--min-tms`
	Minimum number of TMSs for negative control proteins

`--max-tms`
	Maximum number of TMSs for negative control proteins

`tmdatadir`
	Where TMS assignment data is stored

### Deuterocol2 ###

This script prepares a directory structure for bulk structure comparisons

`--fams1, --fams2`
	Lists of TC-IDs corresponding to the superfamilies, families, subfamilies, or systems to compare against each other

`--d1dir`
	Path to Deuterocol1 directory

`--tmdatadir`
	Path to TMdata directory

`--outdir`
	Path to write to

`--allow-internal`
	Allow self-vs-self comparisons, e.g. 1.A.1 vs 1.A.1

### superpose ###

<!--CITEME-->
This wrapper for SSM Superpose runs superpose on all scheduled alignments for a given Deuterocol2 output directory.
Can be safely run multiple times on the same directory without clobbering

`d2dir`
	Path to Deuterocol2 directory

### tmalign ###

This wrapper for TMalign runs TMalign on all scheduled alignments for a given Deuterocol2 output directory.
Can be safely run multiple times on the same directory without clobbering

`d2dir`
	Path to Deuterocol2 directory

### winnower ###

This script selects the best superposition(s) for each pair of polypeptides after filtering by coverage, length, or quality.

`infile`
	A single superposition TSV or a Deuterocol2 output path

`outfile`
	Where to write the surviving alignments

`-d`
	Distance threshold for coverage calculations (default: 4ang)

`--min-cov`
	Minimum coverage (0...1.)

`-s`
	How many extra alignments to keep for each pair of structures (default:0)

`--tmdatadir`
	Path to TMdata directory

## TMdatabase management

These scripts retrieve various TMS predictions from various databases.
Downloading pregenerated copies of the TMData files is recommended over running the scripts individually.

### opm\_dbtool ###

*(Deprecated)* This script turns OPM database dumps into per-chain TMS assignments

### pdbtm\_dbtool ###

This script fetches and unpacks PDBTM TMS assignments

### superfamily\_dbtool ###

This script parses the text on the TCDB [http://tcdb.org/superfamily.php](list of superfamilies) in order to allow for automated negative control generation


