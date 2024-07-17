exec=$TCHEM_INSTALL_PATH/examples/TChem_RateofProgress.x
input=$TCHEM_INSTALL_PATH/examples/runs/atmospheric_chemistry/uci_col/uci_v2_test3.yaml
inputFile=$TCHEM_INSTALL_PATH/examples/runs/atmospheric_chemistry/uci_col/input_conditions_multi_col.yaml
run_this="$exec --chemfile=$input \
          --inputfile=$inputFile \
          --outputfile=rate_of_progress.dat" 

echo $run_this
eval $run_this
