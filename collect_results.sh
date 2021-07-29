cat $INDEX/output.csv >> results.csv
find . -type f -name "*.png" -exec cp {} $WORKDIR \;
