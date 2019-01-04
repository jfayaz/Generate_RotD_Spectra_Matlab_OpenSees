###################################################################################################################################################################################
## Author : Jawad Fayaz  (jfayaz@uci.edu)
################################################################################################################################################################################

proc ReadGMFile {inFilename outFilename dt NumPts} {
	upvar $dt DT;
	upvar $NumPts NPTS;
	if [catch {open $inFilename r} inFileID] {
		puts stderr "Cannot open $inFilename for reading";
	} else {
		set outFileID	[open $outFilename w];
		set flag	0;
		foreach line [split [read $inFileID] \n] {
			if {[llength $line] == 0} {
				continue;
			} elseif {$flag == 2} {
				puts $outFileID $line;
			} else {
				foreach word [split $line] {
					if {$flag == 1 && $word != "DT="} {
						set NPTS1	$word;
						set NPTS2	[split $NPTS1 ,];
						set NPTS	[lindex $NPTS2 0];
					};
					if {$flag == 2} {
						set DT		$word;
						break;
					};
					if {[string match $word "NPTS="] == 1} {
						set flag	1;
					};
					if {[string match $word "DT="] == 1} {
						set flag	2;
					};
				};
			};
		};
		close $outFileID;
		close $inFileID;
	};
};
