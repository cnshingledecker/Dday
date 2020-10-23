# Remove old output files and directories
if [ -d "ab" ]; then
   rm -rf ab 
else
  echo "ab directory doesn't exist"
fi

if [ -d "csv" ]; then
   rm -rf csv 
else
  echo "csv directory doesn't exist"
fi

if [ -f "H2.ab" ]; then
  rm *.ab 
else
  echo "no old *.ab files to remove"
fi

if [ -f "H2.csv" ]; then
  rm *.csv 
else
  echo "no old *.csv files to remove"
fi

mkdir ab
mkdir csv

# Run the program
./monaco

# Copy the output files to the appropriate folders, if needed
if [ -f "H2.ab" ]; then
  mv *.ab ab/ 
fi

if [ -f "H2.csv" ]; then
  mv *.csv csv/
fi

echo "Done!"


