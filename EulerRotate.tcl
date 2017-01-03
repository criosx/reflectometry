set filename "protein"
set betastart  0
set betaend    180
set betainc    5
set gammastart 0
set gammaend   360
set gammainc   5

for {set beta $betastart} {$beta <= $betaend} {set beta [expr {$beta + $betainc}]} {
    for {set gamma $gammastart} {$gamma < $gammaend} {set gamma [expr {$gamma + $gammainc}]} {
        set fname $filename
        append fname ".pdb"
        animate read pdb $fname waitfor all
        set nf [molinfo top get numframes]
        for {set fr 0} {$fr < $nf} {incr fr} {
            set repres [atomselect top all frame $fr]
            $repres move [transaxis z $gamma]
            $repres move [transaxis x $beta]
            #rotate z by $gamma
            #rotate x by $beta
        }
        set fname $filename
        append fname "_tilt" $beta "_orie" $gamma ".pdb"
        animate write pdb $fname waitfor all
        animate delete all
    }
}
