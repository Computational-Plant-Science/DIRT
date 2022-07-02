python3 /opt/DIRT/main.py "$INPUT" $INDEX $MASK_THRESHOLD $EXCISED_ROOTS 1 1 $MARKER_DIAMETER 0 0 0 $WORKDIR /opt/DIRT/traits.csv

# copy individual CSVs to working dir and concatenate CSV output into a single file
results="$WORKDIR/output.$INDEX.csv"
cp "$INDEX/output.csv" $results
all_results="$WORKDIR/results.csv"
if [ -s $all_results ]
then
        tail -qn +2 $results >> $all_results
else
        cat $results > $all_results
fi

# copy all pngs from nested output dirs to working dir
find $INDEX -type f -name "*.png" -exec cp {} $WORKDIR \;
