exec=$TCHEM_INSTALL_PATH/examples/TChem_AtmosphericChemistryE3SM.explicit_euler.x
input=$TCHEM_INSTALL_PATH/examples/runs/atmospheric_chemistry/uci_col/uci_explicit_mech.yaml
inputFile=$TCHEM_INSTALL_PATH/examples/runs/atmospheric_chemistry/uci_col/input_conditions_explicit_part_multi_col.yaml
run_this="$exec --chemfile=$input \
          --inputfile=$inputFile \
          --outputfile=full_gas_stratosphere.dat \
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
