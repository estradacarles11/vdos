#!/bin/bash

echo "///////////////////////////////////////////////////////////////////////////////"
echo "//                                                                           //"
echo "//                                                                           //"
echo "//                                                                           //"
echo "//                              WELCOME TO VDOS                              //"
echo "//                                    by                                     //"
echo "//                           Carles Estrada Girbau                           //"
echo "//                                                                           //"
echo "//                                                                           //"
echo "//                                                                           //"
echo "///////////////////////////////////////////////////////////////////////////////"

if test -e "vdos.param";then
    echo "Ready to read vdos.param"
else
    echo
    echo "You need to create a file named \"vdos.param\" as:"
    echo
    echo "system_label=si54"
    echo
    echo "histogram_steps=70"
    echo
    echo "reciprocal_space_divisions=25-25-25"
    echo
    echo "output=si54.vdos"
    exit
fi

which grep > /dev/null
if [ $? -eq 0 ]; then
    echo 
else
    echo "Please install grep to proceed using:"
    echo
    echo "sudo apt-get install grep"
    exit
fi

which sed > /dev/null
if [ $? -eq 0 ]; then
    echo 
else
    echo "Please install grep to proceed using:"
    echo
    echo "sudo apt-get install grep"
    exit
fi

if grep -q "system_label=" vdos.param; then
    systemlabel=$(grep "system_label=" vdos.param | sed 's/system_label=// ')
    echo "System label: $systemlabel"
else
    echo "No system_label"
    exit
fi

if grep -q "histogram_steps=" vdos.param; then
    hsteps=$(grep "histogram_steps=" vdos.param | sed 's/histogram_steps=// ')
    echo "Histogram steps: $hsteps"
else
    echo "No histogram steps, using default."
    hsteps=100
fi

if grep -q "reciprocal_space_divisions=" vdos.param; then
    line=$(grep "reciprocal_space_divisions=" vdos.param)
    if [[ $line == "reciprocal_space_divisions="*"-"*"-"* ]]; then
        n=$(echo $line | sed 's/reciprocal_space_divisions=// ')
        IFS=- read h k l <<< $n
    else
        h=$(grep "reciprocal_space_divisions=" vdos.param | sed 's/reciprocal_space_divisions=// ')
        k=h
        l=h
    fi
    echo "Reciprocal space divisions: $h $k $l"
else
    echo "No reciprocal_space_divisions, using defaults."
    h=100
    k=100
    l=100
fi

if grep -q "output=" vdos.param; then
    output=$(grep "output=" vdos.param | sed 's/output=// ')
    echo "Output: $output"
else
    echo "No output, using default."
    output="$systemlabel.vdos"
fi

which python > /dev/null
if [ $? -eq 0 ]; then
    echo "Executing the script..."

    python vdos.py -s $hsteps -nh $h -nk $k -nl $l -o $output $systemlabel > $systemlabel.vdosout

    echo "Done!"
else
    echo "Please install python2.7 to proceed using:"
    echo
    echo "sudo apt-get install python"
    exit
fi

which gnuplot > /dev/null
if [ $? -eq 0 ]; then
    echo 
else
    echo "Please install gnuplot to proceed using:"
    echo
    echo "sudo apt-get install gnuplot"
    exit
fi

if test -e "$output";then
    echo "Generating plot..."
else
    echo "Please, make sure you have installed this python modules:"
    echo
    echo "argparse"
    echo "numpy"
    echo "time"
    exit
fi

if test -e "vdos.pdf";then
    rm vdos.pdf
fi

if test -e "$systemlabel.pdf";then
    rm $systemlabel.pdf
fi

gnuplot -e "filename='$systemlabel.vdos'" vdos.gp
mv vdos.pdf $systemlabel-$hsteps-$h-$k-$l.pdf

echo "Done!"

exit
