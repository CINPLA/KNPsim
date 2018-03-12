mkdir make_the_figures/data
cd examples/step_concentration
sh run_application_1.sh
mv *.h5 ../../make_the_figures/data
cd ../point_source
sh run_application_2.sh
mv *.h5 ../../make_the_figures/data
cd ../hay_model
sh run_application_3_fast.sh
mv *.h5 ../../make_the_figures/data
