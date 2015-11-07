make clean
mv output output_temp
mkdir output
tar cvzf lgal.tgz --exclude '.*' --exclude '*.jpg' --exclude '*.pdf'\
    Idl/ \
    Makefile \
    Makefile_compilers \
    Makefile_options \
    My_Makefile_options \
    README.txt \
    awk/ \
    code/ \
    devel/ \
    docs/ \
    input/*.* \
    output \
    run/
rmdir output
mv output_temp output

# Add these back in for a standalone backup
#    CoolFunctions/ \
#    PhotTables/ \
#    PhotTables_2.0/ \
#    SpecPhotTables/ \

