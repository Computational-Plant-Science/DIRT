python3 /opt/DIRT/main.py "$INPUT" $INDEX $MASK_THRESHOLD $EXCISED_ROOTS 1 1 25.1 0 0 0 $WORKDIR /opt/DIRT/traits.csv
cat $INDEX/output.csv >> results.csv
find . -type f -name "*.png" -exec cp {} $WORKDIR \;
