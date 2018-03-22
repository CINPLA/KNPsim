file="neuron_input_1.h5"
if [ ! -f "$file" ]
then
    echo "$0: File '${file}' not found. Fetching from figshare.";
    wget https://ndownloader.figshare.com/files/10622146
    mv 10622146 neuron_input_1.rar
    unrar e neuron_input_1.rar
fi
python run_3d_column.py
python run_3d_column_nofield.py
