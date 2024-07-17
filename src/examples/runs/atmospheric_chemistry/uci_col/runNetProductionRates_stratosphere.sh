exec=$TCHEM_INSTALL_PATH/examples/TChem_NetProductionRates.x
input=$TCHEM_INSTALL_PATH/examples/runs/atmospheric_chemistry/uci_col/uci_explicit_mech.yaml
inputFile=$TCHEM_INSTALL_PATH/examples/runs/atmospheric_chemistry/uci_col/input_conditions_explicit_part_multi_col.yaml
run_this="$exec --chemfile=$input \
          --inputfile=$inputFile \
          --outputfile=net_production_rates_stratosphere.dat" 

echo $run_this
eval $run_this
