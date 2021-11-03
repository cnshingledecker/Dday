#
# Remove old output files and directories
./clean.sh

mkdir ab
mkdir csv

# Run the program
./monacoGeneralization
# Copy the output files to the appropriate folders, if needed
if [ -f "H2.ab" ]; then
  mv *.ab ab/ 
fi
if [ -f "H2.csv" ]; then
  mv *.csv csv/
fi
echo "Done!"


