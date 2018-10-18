# deuterocol

A suite of membrane protein structure comparison tools

## Core

### Deuterocol1 ###

This script collects available indices and coordinates files for any number of TCDB superfamilies/families/subfamilies/systems.
It can also attempt to generate a negative control of comparable size to an existing set of positive control structures.

### Deuterocol2 ###

This script prepares a directory structure for bulk structure comparisons


### superpose ###

<!--CITEME-->
This wrapper for SSM Superpose runs superpose on all scheduled alignments for a given Deuterocol2 output directory.

### tmalign ###

This wrapper for TMalign runs TMalign on all scheduled alignments for a given Deuterocol2 output directory.

### winnower ###

This script selects the best superposition(s) for each pair of polypeptides after filtering by coverage, length, or quality.

## TMdatabase management

These scripts retrieve various TMS predictions from various databases.
Downloading pregenerated copies of the TMData files is recommended over running the scripts individually.

### opm\_dbtool ###

*(Deprecated)* This script turns OPM database dumps into per-chain TMS assignments

### pdbtm\_dbtool ###

This script fetches and unpacks PDBTM TMS assignments

### superfamily\_dbtool ###

This script parses the text on the TCDB [http://tcdb.org/superfamily.php](list of superfamilies) in order to allow for automated negative control generation


