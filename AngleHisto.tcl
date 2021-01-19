package require struct
# reference orientation needs to be in structure 0
# trajectory needs to be in structure 1
set firstFrame 0
set lastFrame 199
set firstres 334
set lastres 430
set filename Rotmatrix.dat


# open file for data output
set outfile [open $filename w]

# loop over all frames
for {set f $firstFrame} {$f<=$lastFrame} {incr f} {

    #obtain transformation matrix
    set sel0 [atomselect 0 "backbone and resid $firstres to $lastres"]
    set sel1 [atomselect 1 "backbone and resid $firstres to $lastres" frame $f]
    set M [measure fit $sel0 $sel1]
    # $sel0 move $M

    #write result to file
    puts $outfile $M
}

# close file
close $outfile
