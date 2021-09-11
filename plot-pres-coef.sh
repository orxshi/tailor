# Merge all files into 'upper1' if the second column (y) is positive or zero.
awk '$2>=0' pres_coef-* > upper1
# Sort upper for the first column (x)
sort -k 1,1 upper1 > upper
# Remove intermediate file
rm upper1
# Merge all files into 'lower1' if the second column (y) is negative.
awk '$2<0' pres_coef-* > lower1
# Sort upper for the first column (x)
sort -k 1,1 lower1 > lower
# Remove intermediate file
rm lower1
# Plot
gnuplot -persist -e "plot \"upper\" using 1:(\$4*-1) with lines notitle, \"lower\" using 1:(\$4*-1) with lines notitle"
