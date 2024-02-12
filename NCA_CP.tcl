spinsys {
	channels 15N 13C 
	nuclei 15N 13C
	dipole 1 2 988.16 0 0 0
}

par {
	proton_frequency 	600e6
	start_operator   	I1x
	detect_operator  	I2p
	np	           		24
	crystal_file    	rep168
	spin_rate        	14e3
	gamma_angles     	3
	sw		   			spin_rate/2
	verbose	   			1101
	method	   			direct
	variable rfN     	spin_rate*2.5
	variable rfC     	spin_rate*1.5
	variable tsw		1.0e6/sw
}

proc pulseq {} {
    global par

	# 15N to 13C CP simulation; shaped pulses used; powers defined in 'main'
	acq_block { pulse_shaped $par(tsw) $par(rect) $par(tacn80) }
	
}

proc main {} {
	global par

	# Open file of different B1 scaling values & split into array "lines"
	set fp1 [open B1_values.txt]
	set lines [split [read $fp1] "\n"]
	close $fp1

	# Prepare shaped pulse for 13C channel, based on Bruker 'tacn80'
	set par(tacn80) [load_shape tacn80]
	set par(tacn80) [shape_dup $par(tacn80) 0 [expr $par(rfC)/100]]
	
	# Create output file for transfer efficiencies
	set fp2 [open B1_dependence.out w+]
	
	# Loop for each line in the file read at the beginning
	foreach line $lines {
		
		# Set a rectangular shaped pulse for 15N, scale by B1 from "line"
		set par(rect) [shape_create 1000 -ampl [expr $par(rfN)*$line] -phase 0]
		
		# Run SIMPSON simulation with the adjusted 15N power level and save
		set f [fsimpson]
		fsave $f dB1_$line.fid
		
		# Find maximum in the CP buildup curve, i.e. transfer efficiency
		set max [findex $f 1]
		for {set j 2} {$j <= [fget $f -np]} {incr j} {
			if {[findex $f $j -re] > $max} {
				set max [findex $f $j -re]
			}
		}
		puts $fp2 "$line $max"
		
		# Free memory
		funload $f
		free_shape $par(rect)
	}
	
	# All simulations ended--close the output file created
	close $fp2
}
