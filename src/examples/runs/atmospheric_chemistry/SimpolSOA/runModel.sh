exec=$TCHEM_INSTALL_PATH/examples/TChem_AerosolChemistry_CVODE.x

run_this="$exec --chemfile=test.yaml \
	  --aerofile=mechanism_aero.yaml \
          --inputfile_particles=scenario_conditions_particle.yaml \
	  --outputfile=full_gas.dat \
          --use-cvode=true \
	  --time-iterations-per-interval=10 \
          --max-time-iterations=100\
	  --tol-time=1e-3 \
          --atol-time=1e-12 \
	  --dtmin=1e-20 \
          --dtmax=10\
          --tend=1000 \
          --atol-newton=1e-12 \
          --rtol-newton=1e-8 \
          --max-newton-iterations=20 \
          --max-time-iterations=200000"

echo $run_this
eval $run_this
