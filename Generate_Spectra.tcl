###################################################################################################################################################################################
## Author : Jawad Fayaz  (jfayaz@uci.edu)
##################################################################################################################################################################################
## RUN THE MATLAB FILE (THIS IS A FUNCTION CALLED IN THE MATLAB FILE)
###################################################################################################################################################################################
## This code develops the RotD50 Sa and RotD100 Sa Spectra of the Bi-Directional 
## Ground Motion records provided in the 'GM' folder which must be in current folder. 
## The two directions of the ground motion record must be named as 'GM1i' and 'GM2i',
## where 'i' is the ground motion number which goes from 1 to 'n', 'n' being the total
## number of ground motions for which the Spectra needs to be generated. The extension
## of the files must be '.AT2'
## 
## For example: If the Spectra of two ground motion records are required, 4 files with
## the following names must be provided in the given 'GM' folder:
##     'GM11.AT2' - Ground Motion 1 in direction 1 (direction 1 can be either one of the bi-directional GM as we are rotating the ground motions it does not matter) 
##     'GM21.AT2' - Ground Motion 1 in direction 2 (direction 2 is the other direction of the bi-directional GM)
##     'GM12.AT2' - Ground Motion 2 in direction 1 (direction 1 can be either one of the bi-directional GM as we are rotating the ground motions it does not matter)  
##     'GM22.AT2' - Ground Motion 2 in direction 2 (direction 2 is the other direction of the bi-directional GM)
## The Ground Motion file must be a vector file with 4 header lines.The first 3 lines can have
## any content, however, the 4th header line must be written exactly as per the following example:
##     'NPTS=  15864, DT= 0.0050' of 0.25 sec}
###################################################################################################################################################################################
set No_of_GMs 2
for {set iM 1} {$iM <= $No_of_GMs} {incr iM 1} { 
puts "Running Ground Motion $iM"
file mkdir Results/RESULTS$iM
for {set k 1} {$k <= 72} {incr k 1} { 
set j 0
for {set i 2} {$i <= 100} {incr i 2} {
set j [expr $j +1]
set T_step($j) $i
}
for {set i 110} {$i <= 200} {incr i 10} {
set j [expr $j +1]
set T_step($j) $i
}
for {set i 225} {$i <= 500} {incr i 25} {
set j [expr $j +1]
set T_step($j) $i
}
set Tp [expr $T_step($k)]
model BasicBuilder -ndm 3 -ndf 6;		
set pi 		[expr 3.141]
set dir 	[pwd]
set g 		386.2
set GMinter 0
# Geometry
set T		[expr double($Tp)/100] 
set L 		[expr 1]; 		
set d   	[expr 2];
set r   	[expr $d/2]
set A 		[expr $pi*($r**2)];
set E  		[expr 1.0];
set Iy 		[expr $pi*($r**4)/4];
set G  		[expr 1.0];
set J 		[expr $pi*($r**4)/2];
set Iz  	[expr $pi*($r**4)/4];
set K   	[expr 3*$E*$Iz/($L**3)];
set M   	[expr ($K*($T**2)/(4*($pi**2)))];
set Tn  	[expr 2*$pi/pow($K/$M,0.5)];
# Nodal coordinates:
node 1 0 0 0;			
node 2 0 0 $L;
	
# Boundary Conditions
fix 1 1 1 1 1 1 1; 			
# Nodal masses:
mass 2 $M $M 0. 0. 0. 0.;	
# ##### ELEMENTS 
# Geometric transformation: performs a linear geometric transformation of stiffness and resisting force from the basic system to the global-coordinate system
set TransfTag 1; 			
set TransfType Linear ;		
geomTransf $TransfType $TransfTag 0 -1 0 ; 	
 	
# Element connectivity:
element elasticBeamColumn 12 1 2 $A $E $G $J $Iy $Iz $TransfTag;			
# Periods
set wb				  [eigen  1];
set wwb			      [lindex $wb 0];
set Tb			      [expr 2*$pi/sqrt($wwb)];
puts [format "Running Period = %.2f for Ground Motion %.0f" $Tb $iM]
# Damping
set xDamp 		0.05;					
set MpropSwitch 1.0;
set KcurrSwitch 0.0;
set KcommSwitch 1.0;
set KinitSwitch 0.0;
set lambdaN 	[eigen 1]; 	   																
set lambdaI 	[lindex $lambdaN [expr 0]];
set omegaI 		[expr pow($lambdaI,0.5)];
set alphaM 		[expr $MpropSwitch*$xDamp*(2*$omegaI)/($omegaI)];
set betaKcurr 	[expr $KcurrSwitch*2.*$xDamp/($omegaI)];
set betaKcomm 	[expr $KcommSwitch*2.*$xDamp/($omegaI)];
set betaKinit 	[expr $KinitSwitch*2.*$xDamp/($omegaI)];
rayleigh 		$alphaM $betaKcurr $betaKinit $betaKcomm;
# ##### RECORDERS 
set st _
recorder Node -file $dir/Results/RESULTS$iM/Node$iM$st$T.out -time -node 2 -dof 1 2 3 4 5 6 disp;		
# ##### ANALYSIS 
set dtfact 					1;
set Tol						1.0e-6;
set maxNumIter				100;
set printFlag				0;
set TestType				EnergyIncr;
set NewmarkGamma			0.50;
set NewmarkBeta				0.25;
set algorithmType			NewtonLineSearch;
constraints Penalty 		1.e18 1.e18;
numberer 					RCM;
system SparseGeneral 		-piv;
test 						$TestType $Tol $maxNumIter $printFlag;
algorithm 					$algorithmType;
integrator 					TRBDF2;
analysis 					Transient;
# ## GENERATING G3 FILES
source 						ReadGMFile.tcl
set iEQ 					[expr $iM];
set iGMinput 				"GM1$iEQ GM2$iEQ"
foreach GMinput 			$iGMinput {
	set inFile 				"$dir/GMs/$GMinput.AT2";
	set outFile 			"$dir/GMs/$GMinput.g3";
	ReadGMFile				$inFile $outFile dt NumPts;
};
set DtAnalysis				[expr $dt/$dtfact];
set TmaxAnalysis			[expr $dt*($NumPts-1)];
set Nsteps					[expr int($TmaxAnalysis/$DtAnalysis)];
set iGMfile		            "$dir/GMs/GM1$iEQ.g3 $dir/GMs/GM2$iEQ.g3 $dir/GMs/GM1$iEQ.g3 $dir/GMs/GM2$iEQ.g3";
set iGMdirection			"1 1 2 2";
set iGMfact1				[expr cos($GMinter*$pi/180)];
set iGMfact2				[expr -sin($GMinter*$pi/180)]; 
set iGMfact3				[expr sin($GMinter*$pi/180)];
set iGMfact4				[expr cos($GMinter*$pi/180)];
set iGMfact 				"$iGMfact1 $iGMfact2 $iGMfact3 $iGMfact4"
set IDloadTag 				4
set iloop 					"1 2 3 4"
foreach Sloop $iloop GMdirection $iGMdirection GMfile $iGMfile GMfact $iGMfact {
	incr IDloadTag 
	set AccelSeries 			"Series -dt $dt -filePath $GMfile -factor [expr $GMfact*$g]"
	pattern UniformExcitation 	$IDloadTag $GMdirection -accel $AccelSeries;
};
for {set ik 1} {$ik <= $Nsteps} {incr ik 1} {
	set ok	  [analyze 1 $DtAnalysis];
	
	if {$ok != 0} {
        
		set DtAnalysis 0.0001
		puts "Trying dt = 0.0001..."
		
		for {set ikDummy 1} {$ikDummy <= [expr int($dt/$DtAnalysis)]} {incr ikDummy 1} {
			set ok [analyze 1 $DtAnalysis]
	if {$ok != 0} {
		puts "Trying Bisection ...";
		algorithm NewtonLineSearch <-type Bisection>;
		test $TestType 1.0e-5 $maxNumIter $printFlag;
		set ok [analyze 1 $DtAnalysis];
		test $TestType $Tol $maxNumIter $printFlag;
		algorithm $algorithmType;
	};
	if {$ok != 0} {
		puts "Trying Secant ...";
		algorithm NewtonLineSearch <-type Secant>;
		test $TestType 1.0e-5 $maxNumIter $printFlag;
		set ok [analyze 1 $DtAnalysis];
		test $TestType $Tol $maxNumIter $printFlag;
		algorithm $algorithmType;
	};
	if {$ok != 0} {
		puts "Trying Newton ...";
		algorithm Newton;
		test $TestType 1.0e-5 $maxNumIter $printFlag;
		set ok [analyze 1 $DtAnalysis];
		test $TestType $Tol $maxNumIter $printFlag;
		algorithm $algorithmType;
	};
	if {$ok != 0} {
		puts "Trying Modified Newton ...";
		algorithm ModifiedNewton -initial;
		test $TestType 1.0e-5 $maxNumIter $printFlag;
		set ok [analyze 1 $DtAnalysis];
		test $TestType $Tol $maxNumIter $printFlag;
		algorithm $algorithmType;
	};
	if {$ok != 0} {
		puts "Trying Secant Newton...";
		algorithm SecantNewton;
		test $TestType 1.0e-5 $maxNumIter $printFlag;
		set ok [analyze 1 $DtAnalysis];
		test $TestType $Tol $maxNumIter $printFlag;
		algorithm $algorithmType;
	};
	if {$ok != 0} {
		puts "Trying more iterations...";
		test $TestType 1.0e-5 1000 $printFlag;
		set ok [analyze 1 $DtAnalysis];
		test $TestType $Tol $maxNumIter $printFlag;
	};
	if {$ok != 0} {
		puts "Trying tolerance 1.0e-3 ...";
		test $TestType 1.0e-3 $maxNumIter 0;
		set ok [analyze 1 $DtAnalysis];
		test $TestType $Tol $maxNumIter $printFlag;
	};
	if {$ok != 0} {
		puts "Trying Bisection ...";
		algorithm NewtonLineSearch <-type Bisection>;
		test $TestType 1.0e-3 $maxNumIter $printFlag;
		set ok [analyze 1 $DtAnalysis];
		test $TestType $Tol $maxNumIter $printFlag;
		algorithm $algorithmType;
	};
	if {$ok != 0} {
		puts "Trying Secant ...";
		algorithm NewtonLineSearch <-type Secant>;
		test $TestType 1.0e-3 $maxNumIter $printFlag;
		set ok [analyze 1 $DtAnalysis];
		test $TestType $Tol $maxNumIter $printFlag;
		algorithm $algorithmType;
	};
	if {$ok != 0} {
		puts "Trying Newton ...";
		algorithm Newton;
		test $TestType 1.0e-3 $maxNumIter $printFlag;
		set ok [analyze 1 $DtAnalysis];
		test $TestType $Tol $maxNumIter $printFlag;
		algorithm $algorithmType;
	};
	if {$ok != 0} {
		puts "Trying Modified Newton ...";
		algorithm ModifiedNewton -initial;
		test $TestType 1.0e-3 $maxNumIter $printFlag;
		set ok [analyze 1 $DtAnalysis];
		test $TestType $Tol $maxNumIter $printFlag;
		algorithm $algorithmType;
	};
	if {$ok != 0} {
		puts "Trying Secant Newton...";
		algorithm SecantNewton;
		test $TestType 1.0e-3 $maxNumIter $printFlag;
		set ok [analyze 1 $DtAnalysis];
		test $TestType $Tol $maxNumIter $printFlag;
		algorithm $algorithmType;
	};
	if {$ok != 0} {
		puts "Trying more iterations...";
		test $TestType 1.0e-3 1000 $printFlag;
		set ok [analyze 1 $DtAnalysis];
		test $TestType $Tol $maxNumIter $printFlag;
	};
				if {$ok != 0} {
					set Nstepsmax [expr $ik-1];
					break;
				};	
		};                
			if { $ok != 0 } {
				puts "Analysis did not converge; terminate simulation ..."
				break;
			};
           		set DtAnalysis [expr $dt/$dtfact];    			                  
	};
		
};
wipe all 
};
wipe all 
};

