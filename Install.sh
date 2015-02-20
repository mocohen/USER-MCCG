# Install/unInstall package files in LAMMPS

echo ''
echo 'Add MCCG '

if (test $1 = 1) then
    
    cp -p fix_mccg.cpp ../
    cp -p fix_mccg.h   ../
    

elif (test $1 = 0) then

    rm fix_mccg.cpp
    rm fix_mccg.h



fi