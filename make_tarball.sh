make clean
mv output output_temp
mkdir output
tar cvzf lgal.tgz --exclude '.*' --exclude '*.jpg' --exclude '*.pdf'\
    Auxcode/ \
    H2frac \
    Makefile \
    Makefile_compilers \
    Makefile_options \
    Makefile_options_MCMC \
    My_Makefile_options \
    README.txt \
    code/ \
    docs/ \
    input/ \
    output/
rmdir output
mv output_temp output

# Add these back in for a standalone backup
#    CoolFunctions/ \
#    MergerTrees \
#    PhotTables/ \
#    PhotTables_2.0/ \
#    SpecPhotTables/ \
#    YieldTables/ \

