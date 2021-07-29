cat $INDEX/output.csv >> results.csv
find . -type f | grep -i .png | xargs -i cp {} $WORKDIR
