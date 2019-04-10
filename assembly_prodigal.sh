echo "Starting AssemblyTools."

echo "Removing old files."
rm -rf /Users/viniWS/Bio/masters/test_data/ncbi/ncbi-genomes-2018-12-20

echo "Decompressing new files."
tar -C /Users/viniWS/Bio/masters/test_data/ncbi/ -zxf /Users/viniWS/Bio/masters/test_data/ncbi/genome_assemblies.tar

echo "Renaming files."
python /Users/viniWS/Bio/masters/scripts/assembly_tools/renamer.py -i /Users/viniWS/Bio/masters/test_data/ncbi/ncbi-genomes-2018-12-20

echo "Moving genomic files to folder."
rm -rf /Users/viniWS/Bio/masters/test_data/ncbi/genomic_fna/
mkdir /Users/viniWS/Bio/masters/test_data/ncbi/ncbi-genomes-2018-12-20/genomic_fna/
mv /Users/viniWS/Bio/masters/test_data/ncbi/ncbi-genomes-2018-12-20/*/*genomic.fna /Users/viniWS/Bio/masters/test_data/ncbi/ncbi-genomes-2018-12-20/genomic_fna/

#echo "Starting Prodigal."
#python /Users/viniWS/Bio/masters/scripts/assembly_tools/prodigal.py -i /Users/viniWS/Bio/masters/test_data/ncbi/ncbi-genomes-2018-12-20/genomic_fna/
