if [ -z "${CAMP_INSTALL_PATH}"]; then
   echo "Set CAMP_INSTALL_PATH to run CAMP test"
   exit 1
fi
mkdir out
exec=$CAMP_INSTALL_PATH/test_chemistry_cb05cl_ae5
echo "$exec"
eval "$exec"
