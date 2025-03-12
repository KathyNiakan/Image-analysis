// This script accompanies the manuscript Simon et. al., (2025)
// Repository available on https://ctr.uniofcam.dev/ctr-bioinformatics/niakan-lab/simon-et-al-2025

// This macro script will take .tif and stardist.tif files and split to Ch and Z single files
// Each image sequence will be saved in a seperate folder for each embryo
// Files can then be renamed for use in CellProfiler pipeline

// Inputs:  .tif multichannel Z stack output from 1_LIF_series_exporter
//			.tif stardist segmentation output from 2_batch_stardist
// Output:	.tif single Ch and Z files in separate folders for each embryo

setBatchMode(true);
print("\\Clear");
dir = getDirectory("Choose file dir");
print("dir: "+dir);
list = getFileList(dir);
print("Number o files: "+list.length);
fs = File.separator;

save_dir=dir+fs+"split_files";
File.makeDirectory(save_dir);
print("Saving files in: "+save_dir);


WrongPattern1 = ".*Z[0-9]{1,2}\\.tif$";
WrongPattern2 = ".*Z[0-9]{1,2}_C.*\\.tif$";


for(i=0;i<list.length;i++) {
	if(endsWith(list[i], ".tif") == true && endsWith(list[i], "_stardist.tif") == false) {
		print("File nr: "+i+1);
		run("Bio-Formats Importer", "open=["+dir+fs+list[i]+"] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
		title = File.nameWithoutExtension();
		print("Image title: "+title);
		run("Bio-Formats Exporter", "save=["+save_dir+fs+title+fs+title+".tif] write_each_z_section write_each_channel compression=Uncompressed");
		run("Close All");
		print("done");
	};
	if(endsWith(list[i], "_stardist.tif") == true) {
		print("File nr: "+i+1);
		run("Bio-Formats Importer", "open=["+dir+fs+list[i]+"] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
		title = File.nameWithoutExtension();
		titleshort = replace(title, "_stardist", "");
		print("Image title: "+title);
		run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
		run("Bio-Formats Exporter", "save=["+save_dir+fs+titleshort+fs+title+".tif] write_each_z_section compression=Uncompressed");
		run("Close All");
		print("done");
	};
	showProgress(-i/list.length);
};
print("All done");
setBatchMode(false);