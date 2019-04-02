echo "Starting AssemblyTools."

echo "Removing old files."
rm -rf ncbi-genomes-2018-12-20

echo "Decompressing new files."
tar zxf genome_assemblies.tar

echo "Renaming files."
python /Users/viniWS/Bio/masters/scripts/assembly_tools/renamer.py -b $1 

echo "Moving genomic files to folder."
rm -rf genomic_fna/
mkdir genomic_fna/
mv $1/*/*genomic.fna genomic_fna/

echo "Starting Prodigal."
python /Users/viniWS/Bio/masters/scripts/assembly_tools/prodigal.py -i genomic_fna/
