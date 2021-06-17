#!bin/bash
solve-field --scale-units arcsecperpix --scale-high 0.4 --scale-low 0.2 --no-plots --new-fit ~/Documents/NovaAstrometry/solved_$1 ~/Documents/NovaAstrometry/$1 --continue

n="$1"
shortn=$(echo $n | sed -r 's/.{5}$//')

rm $shortn-indx.xyls
rm $shortn.wcs
rm $shortn.solved
rm $shortn.rdls
rm $shortn.match
rm $shortn.corr
rm $shortn.axy
rm $1
mv ~/Documents/NovaAstrometry/solved_$1 ~/Documents/NovaAstrometry/$1
