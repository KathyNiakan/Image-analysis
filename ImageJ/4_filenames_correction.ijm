// This script accompanies the manuscript Simon et. al., (2025)
// Repository available on https://ctr.uniofcam.dev/ctr-bioinformatics/niakan-lab/simon-et-al-2025

// This macro script will rename .tif files split to Ch and Z single files 
// and change the number of leading zeros for use with CellProfiler

// Input:	.tif single Ch and Z files output from 3_split_files
// Output:	.tif files renamed with matching number of leading zeros

dir = getDirectory("Please select a directory");
print(dir)
filenames = getFileList(dir);

WrongPattern1 = ".*Z[0-9]{1,2}\\.tif$";
WrongPattern2 = ".*Z[0-9]{1,2}_C.*\\.tif$";

for(i = 0; i < filenames.length; i++) {
	newFileName = filenames[i];
	while (matches(newFileName, WrongPattern1) | matches(newFileName, WrongPattern2)) {
		print("Wrong Pattern Matched");
		newFileName = replace(newFileName, "Z", "Z0");
	} 
	File.rename(dir+filenames[i], dir+newFileName);
}

print("Done")
