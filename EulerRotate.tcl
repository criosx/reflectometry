set filename "1OGS"
set betastart -90
set betaend    90
set betainc    5
set gammastart 0
set gammaend   360
set gammainc   5

for {set beta $betastart} {$beta <= $betaend} {set beta [expr {$beta + $betainc}]} {
    for {set gamma $gammastart} {$gamma < $gammaend} {set gamma [expr {$gamma + $gammainc}]} {
        set fname $filename
        append fname ".pdb"
        mol new $fname
        set repres [atomselect top all]
        $repres move [transaxis z $gamma]
        $repres move [transaxis x $beta]
        #rotate z by $gamma
        #rotate x by $beta
        set fname $filename
        append fname "_tilt" $beta "_orie" $gamma ".pdb"
        $repres writepdb $fname
        mol delete all
    }
}
