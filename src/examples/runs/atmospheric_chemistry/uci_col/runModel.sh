exec=$TCHEM_INSTALL_PATH/examples/TChem_AtmosphericChemistryE3SM.implicit_euler.x
input=$TCHEM_INSTALL_PATH/examples/runs/atmospheric_chemistry/uci_col/uci_v2_test3.yaml
inputFile=$TCHEM_INSTALL_PATH/examples/runs/atmospheric_chemistry/uci_col/input_conditions_multi_col.yaml
run_this="$exec --chemfile=$input \
          --inputfile=$inputFile \
          --outputfile=full_gas.dat \
          --time-iterations-per-interval=100 \
          --tol-time=1e-6 \
          --dtmin=1800 \
          --dtmax=1800 \
          --tend=1800 \
          --atol-newton=1e-18 \
          --rtol-newton=1e-8 \
          --max-newton-iterations=20 \
          --max-time-iterations=20000"

echo $run_this
eval $run_this
