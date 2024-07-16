exec=$TCHEM_INSTALL_PATH/examples/TChem_AtmosphericChemistryE3SM.x
input=$TCHEM_INSTALL_PATH/examples/runs/atmospheric_chemistry/e3sm_v2/e3sm_v2.yaml
run_this="$exec --chemfile=$input \
          --outputfile=full_gas.dat \
          --time-iterations-per-interval=100 \
          --tol-time=1e-6 \
          --dtmin=1e-20 \
          --dtmax=10 \
          --tend=1.08e4 \
          --atol-newton=1e-18 \
          --rtol-newton=1e-8 \
          --max-newton-iterations=20 \
          --max-time-iterations=20000"

echo $run_this
eval $run_this
